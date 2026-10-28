#include "amici/backwardproblem.h"

#include "amici/edata.h"
#include "amici/exception.h"
#include "amici/forwardproblem.h"
#include "amici/model.h"
#include "amici/solver.h"
#include "amici/steadystateproblem.h"

#include <ranges>

namespace amici {
constexpr realtype conv_thresh = 1.0;

BackwardProblem::BackwardProblem(ForwardProblem& fwd)
    : model_(fwd.model)
    , solver_(fwd.solver)
    , edata_(fwd.edata)
    , t_(fwd.get_final_time())
    , discs_main_(fwd.get_discontinuities())
    , main_final_state_(fwd.get_final_simulation_state())
    , x_main_initial_(fwd.get_initial_simulation_state().sol.x)
    , dJydx_(fwd.get_adjoint_updates(*model_, *edata_))
    , dJzdx_(fwd.get_dJzdx())
    , preeq_problem_(fwd.get_preeq_problem())
    , posteq_problem_(fwd.get_posteq_problem())
    , presim_result(fwd.get_presimulation_result())
    , ws_(model_, solver_)
    , simulator_(model_, solver_, &ws_) {}

void BackwardProblem::workBackwardProblem() {

    if (model_->nx_solver <= 0
        || solver_->get_sensitivity_order() < SensitivityOrder::first
        || solver_->get_sensitivity_method() != SensitivityMethod::adjoint
        || model_->nplist() == 0) {
        return;
    }

    handle_postequilibration();

    // handle main simulation

    // If we have posteq, infinity timepoints were already treated
    int it = model_->nt() - 1;
    while (it >= 0 && std::isinf(model_->get_timepoint(it))) {
        --it;
    }

    model_->restore_trajectory_state(main_final_state_);
    ws_.discs_ = discs_main_;
    ws_.nroots_ = compute_nroots(discs_main_, model_->ne);
    simulator_.run(
        t_, model_->t0(), it, model_->get_timepoints(), &dJydx_, &dJzdx_
    );

    bool const asa_preeq = preeq_problem_
                           && preeq_problem_->get_solver()
                                      ->get_sensitivity_method_pre_equilibration()
                                  == SensitivityMethod::adjoint;
    if (asa_preeq) {
        // Reinitialization at the start of the main simulation. With forward
        // pre-equilibration, this is part of the state sensitivities instead.
        // ASA presimulation with reinitialization is not supported.
        model_->reinitialize_b(
            model_->t0(), x_main_initial_, ws_.xB_, ws_.xQB_
        );
    }

    // handle presimulation
    if (edata_ && edata_->t_presim > 0) {
        ConditionContext cc(
            model_, edata_, FixedParameterContext::presimulation
        );
        model_->restore_trajectory_state(presim_result.final_state_);
        ws_.discs_ = presim_result.discs;
        ws_.nroots_ = compute_nroots(ws_.discs_, model_->ne);
        simulator_.run(
            model_->t0(), model_->t0() - edata_->t_presim, -1, {}, &dJydx_,
            &dJzdx_
        );
    }

    // store pre-pre-equilibration adjoint state and quadrature still needed
    // for computing sllh in ReturnData
    xB_pre_preeq_ = ws_.xB_;
    xQB_pre_preeq_ = ws_.xQB_;

    // handle pre-equilibration
    if (asa_preeq) {
        auto preeq_solver = preeq_problem_->get_solver();

        ConditionContext cc2(
            model_, edata_, FixedParameterContext::preequilibration
        );

        auto const t0 = std::isnan(model_->t0_preeq()) ? model_->t0()
                                                       : model_->t0_preeq();

        auto const& preeq_result = preeq_problem_->get_result();
        model_->restore_trajectory_state(preeq_result.final_state_);

        // If there were no discontinuities or no simulation was performed,
        // we can use the steady-state shortcuts.
        // If not we need to do the regular backward integration.
        if (preeq_problem_->get_steady_state_status()[1]
                == SteadyStateStatus::not_run
            || preeq_result.discs.empty()) {
            preeq_solver->update_and_reinit_states_and_sensitivities(model_);

            auto preeq_final_fwd_state = preeq_result.final_state_;
            preeq_problem_bwd_.emplace(
                *solver_, *model_, preeq_final_fwd_state.sol, &ws_
            );
            preeq_problem_bwd_->run(t0);
        } else if (preeq_problem_->get_steady_state_status()[1]
                   != SteadyStateStatus::not_run) {
            // backward integration of the pre-equilibration problem
            // (only if preequilibration was done via simulation)
            ws_.discs_ = preeq_result.discs;
            ws_.nroots_ = compute_nroots(ws_.discs_, model_->ne);
            EventHandlingBwdSimulator preeq_simulator(
                model_, preeq_solver, &ws_
            );
            preeq_simulator.run(
                preeq_result.final_state_.sol.t,
                preeq_result.initial_state_.sol.t, -1, {}, &dJydx_, &dJzdx_
            );
        } else {
            Expects(
                false
                && "Unhandled preequilibration case in "
                   "BackwardProblem::workBackwardProblem()"
            );
        }
    }

    // Leave the model in its state at the end of the main simulation, which
    // is used for computing the outputs in ReturnData.
    model_->restore_trajectory_state(main_final_state_);
}

void BackwardProblem::handle_postequilibration() {
    if (!posteq_problem_) {
        return;
    }

    // initialize xB - only process the post-equilibration timepoints
    for (int it = 0; it < model_->nt(); it++) {
        if (std::isinf(model_->get_timepoint(it))) {
            for (int ix = 0; ix < model_->nxtrue_solver; ix++)
                ws_.xB_[ix] += dJydx_[ix + it * model_->nx_solver];
        }
    }

    auto const& posteq_result = posteq_problem_->get_result();
    model_->restore_trajectory_state(posteq_result.final_state_);

    // If there were no discontinuities or no simulation was performed,
    // we can use the steady-state shortcuts.
    // If not we need to do the regular backward integration.
    if (posteq_problem_->get_steady_state_status()[1]
            == SteadyStateStatus::not_run
        || posteq_result.discs.empty()) {

        auto final_state = posteq_problem_->get_final_simulation_state();
        posteq_problem_bwd_.emplace(*solver_, *model_, final_state.sol, &ws_);
        posteq_problem_bwd_->run(model_->t0());

        // re-initialize state vectors for main simulation
        model_->initialize_b(ws_.xB_, ws_.dxB_, ws_.xQB_, true);
    } else if (posteq_problem_->get_steady_state_status()[1]
               != SteadyStateStatus::not_run) {
        // backward integration of the post-equilibration problem
        // (only if post-equilibration was done via simulation)
        ws_.discs_ = posteq_result.discs;
        ws_.nroots_ = compute_nroots(ws_.discs_, model_->ne);
        EventHandlingBwdSimulator posteq_simulator(model_, solver_, &ws_);
        posteq_simulator.run(
            posteq_result.final_state_.sol.t,
            posteq_result.initial_state_.sol.t, -1, {}, &dJydx_, &dJzdx_
        );
    } else {
        Expects(
            false
            && "Unhandled post-equilibration case in "
               "BackwardProblem::handlePostequilibration()"
        );
    }
}

void EventHandlingBwdSimulator::handle_event_b(
    Discontinuity const& disc, std::vector<realtype> const* dJzdx
) {
    // Undo events in the reverse of the order they were applied going
    // forward (which is not necessarily their index order).
    // All events of this discontinuity share the event time, so the terms
    // depending on the event-time sensitivity are accumulated and applied
    // at the trigger state below. The same holds for the dependence on the
    // state at the trigger time for events using values from the trigger
    // time.
    std::vector<realtype> stauB(model_->nJ, 0.0);
    // This also applies to the secondary discontinuities triggered by this
    // one, which were handled right before.
    if (!stauB_secondary_.empty()) {
        Expects(!disc.event_applications.empty());
        stauB = std::move(stauB_secondary_);
        stauB_secondary_.clear();
    }
    // Jump in the right-hand side due to a Heaviside update not followed by
    // any event application (adjoint of the corresponding `sx` update).
    if (disc.xdot_jump.size() > 0) {
        Expects(disc.event_applications.empty());
        for (int iJ = 0; iJ < model_->nJ; ++iJ)
            for (int ix = 0; ix < model_->nxtrue_solver; ++ix)
                stauB.at(iJ) += ws_->xB_[ix + iJ * model_->nxtrue_solver]
                                * disc.xdot_jump[ix];
    }
    AmiVector xB_at_trigger(ws_->xB_.size(), solver_->get_sun_context());
    xB_at_trigger.zero();
    // Events of later discontinuities may use the state at which this one
    // was reached.
    if (!xB_at_trigger_earlier_.empty()) {
        xB_at_trigger.copy(xB_at_trigger_earlier_.front());
        xB_at_trigger_earlier_.erase(xB_at_trigger_earlier_.begin());
    }
    for (auto const& application :
         std::ranges::reverse_view(disc.event_applications)) {
        auto const offset = application.trigger_disc_offset;
        Expects(offset >= 0);
        while (gsl::narrow<int>(xB_at_trigger_earlier_.size()) < offset) {
            xB_at_trigger_earlier_.emplace_back(
                ws_->xB_.size(), solver_->get_sun_context()
            );
            xB_at_trigger_earlier_.back().zero();
        }
        model_->add_adjoint_event_update(
            ws_->xB_, ws_->xQB_,
            offset == 0 ? xB_at_trigger : xB_at_trigger_earlier_.at(offset - 1),
            stauB, application.ie, t_, application.x_post,
            application.xdot_post, application.xdot_eval, application.xdot_pre,
            application.x_eval, application.h_assign, application.dx_post,
            application.h_eval
        );
    }
    for (int i = 0; i < ws_->xB_.size(); ++i)
        ws_->xB_[i] += xB_at_trigger[i];
    // There is no earlier discontinuity at this time
    Expects(disc.secondary_event || xB_at_trigger_earlier_.empty());

    // Event outputs were evaluated at the state at which the discontinuity
    // was reached, also for events without state updates and for the events
    // filled in at the last timepoint.
    for (int ie = 0; ie < model_->ne; ie++) {
        if (disc.root_info[ie] != 1)
            continue;

        auto const iroot = --ws_->nroots_[ie];
        Expects(iroot >= 0);
        if (model_->nz == 0 || iroot >= model_->n_max_event())
            continue;

        auto const offset = (ie * model_->n_max_event() + iroot)
                            * model_->nx_solver * model_->nJ;
        for (int ix = 0; ix < model_->nxtrue_solver; ++ix) {
            for (int iJ = 0; iJ < model_->nJ; ++iJ) {
                ws_->xB_[ix + iJ * model_->nxtrue_solver]
                    += (*dJzdx)[offset + iJ + ix * model_->nJ];
            }
        }
    }

    // restore the state at which the discontinuity was reached
    auto state = model_->get_model_state();
    state.h = disc.h_at_trigger;
    state.total_cl = disc.total_cl_at_trigger;
    model_->set_model_state(state);

    // The event time doesn't depend on the parameters for events at the
    // initial time. Otherwise, it's defined by the triggering root function
    // (if several were triggered simultaneously, the last one, as for
    // forward sensitivities). Secondary discontinuities share the event time
    // of the discontinuity that triggered them.
    if (disc.secondary_event) {
        stauB_secondary_ = std::move(stauB);
    } else if (!disc.initial_event && !disc.event_applications.empty()) {
        int ie_tau = -1;
        for (int ie = 0; ie < model_->ne; ie++) {
            if (disc.root_info[ie] == 1)
                ie_tau = ie;
        }
        if (ie_tau >= 0) {
            model_->add_adjoint_event_time_update(
                ws_->xB_, ws_->xQB_, stauB, ie_tau, t_, disc.x_at_trigger,
                disc.dx_at_trigger
            );
        }
    }
}

void EventHandlingBwdSimulator::handle_datapoint_b(
    int const it, std::vector<realtype> const* dJydx
) {
    // The solver wasn't reset yet, as xB_ is necessary for solver setup.
    // For initial time point (we are integrating backwards!), diagnosis needs
    // to be stored outside this function.
    if (it < model_->nt() - 1)
        solver_->store_diagnosis_b(ws_->which);

    Expects(dJydx != nullptr);
    for (int ix = 0; ix < model_->nxtrue_solver; ix++) {
        for (int iJ = 0; iJ < model_->nJ; iJ++)
            // we only need the 1:nxtrue_solver (not the nx_true) slice here!
            ws_->xB_[ix + iJ * model_->nxtrue_solver]
                += (*dJydx)[iJ + (ix + it * model_->nx_solver) * model_->nJ];
    }
}

realtype EventHandlingBwdSimulator::get_next_t(int const it) {
    if (it < 0 && ws_->discs_.empty()) {
        throw AmiException(
            "No more timepoints (it=%d, ie=%d) available at %f. This should "
            "not happen, please report a bug including this stacktrace at "
            "https://github.com/AMICI-dev/AMICI/issues/new/choose",
            it, ws_->discs_.size(), this->t_
        );
    }

    if (!ws_->discs_.empty()
        && (it < 0 || ws_->discs_.back().time > model_->get_timepoint(it))) {
        double const tdisc = ws_->discs_.back().time;
        return tdisc;
    }

    return model_->get_timepoint(it);
}

BwdSimWorkspace::BwdSimWorkspace(
    gsl::not_null<Model*> model, gsl::not_null<Solver const*> solver
)
    : model_(model)
    , xB_(model_->nx_solver, solver->get_sun_context())
    , dxB_(model_->nx_solver, solver->get_sun_context())
    , xQB_(model_->nJ * model_->nplist(), solver->get_sun_context()) {}

void EventHandlingBwdSimulator::run(
    realtype const t_start, realtype const t_end, realtype it,
    std::vector<realtype> const& timepoints, std::vector<realtype> const* dJydx,
    std::vector<realtype> const* dJzdx
) {
    Expects(model_->nz == 0 || dJzdx != nullptr);
    Expects(t_start >= t_end);
    Expects(it < 0 || t_start >= timepoints[it]);
    Expects(it < 0 || t_end <= timepoints.front());

    t_ = t_start;

    // If there is nothing to integrate, only the updates at `t_start` are
    // applied. The backward solver is not set up in this case, as IDAS
    // requires that the forward problem was integrated before.
    bool const integrate = t_start > t_end;

    // datapoint at t_start?
    if (it >= 0 && timepoints[it] == t_start) {
        handle_datapoint_b(it, dJydx);
        if (integrate)
            solver_->setup_b(
                &ws_->which, timepoints[it], model_, ws_->xB_, ws_->dxB_,
                ws_->xQB_
            );
        // for the initial datapoint, diagnosis needs to be stored after setup
        // as it is not called in handleDataPointB
        solver_->store_diagnosis_b(ws_->which);
        --it;
    } else {
        // no data point at t_start, only set up the solver if needed
        // (e.g., during presimulation)
        if (integrate)
            solver_->setup_b(
                &ws_->which, t_start, model_, ws_->xB_, ws_->dxB_, ws_->xQB_
            );
    }

    // DAEs: make the algebraic adjoint states and the derivatives of the
    // differential adjoint states consistent
    auto const calc_ic_b = [&]() {
        auto const tout1
            = (it >= 0 || !ws_->discs_.empty()) ? get_next_t(it) : t_end;
        if (tout1 < t_)
            solver_->calc_ic_b(ws_->which, tout1);
    };
    calc_ic_b();

    while (it >= 0 || !ws_->discs_.empty()) {
        // check if next timepoint is a discontinuity or a data-point
        double const tnext = get_next_t(it);

        if (tnext < t_) {
            solver_->run_b(tnext);
            solver_->write_solution_b(
                t_, ws_->xB_, ws_->dxB_, ws_->xQB_, ws_->which
            );
        }

        // handle data-point
        if (it >= 0 && tnext == timepoints[it]) {
            handle_datapoint_b(it, dJydx);
            it--;
        }

        // handle discontinuity
        if (!ws_->discs_.empty() && tnext == ws_->discs_.back().time) {
            handle_event_b(ws_->discs_.back(), dJzdx);
            ws_->discs_.pop_back();
        }

        // reinitialize state
        if (integrate) {
            solver_->reinit_b(ws_->which, t_, ws_->xB_, ws_->dxB_);
            solver_->reinit_quad_b(ws_->which, ws_->xQB_);
            calc_ic_b();
        }
    }

    // we still need to integrate from first datapoint to t_start
    if (t_ > t_end) {
        solver_->run_b(t_end);
        solver_->write_solution_b(
            t_, ws_->xB_, ws_->dxB_, ws_->xQB_, ws_->which
        );
    }
}

/**
 * @brief Compute the backward quadratures, which contribute to the
 * gradient (xQB) from the quadrature over the backward state itself (xQ)
 * @param model Model instance.
 * @param yQ vector to be multiplied with dxdotdp
 * @param yQB resulting vector after multiplication
 * @param state Simulation state for which to compute dxdot/dp
 */
void computeQBfromQ(
    Model& model, AmiVector const& yQ, AmiVector& yQB, DEStateView const& state
) {
    // Compute the quadrature as the inner product: `yQB = dxdotdp * yQ`.

    // set to zero first, as multiplication adds to existing value
    yQB.zero();
    // yQB += dxdotdp * yQ
    // fill dxdotdp with current values
    auto const& plist = model.get_parameter_list();
    model.fdxdotdp(state.t, state.x, state.dx);
    model.get_dxdotdp_full().multiply(
        yQB.get_nvector(), yQ.get_nvector(), plist, true
    );
}

SteadyStateBackwardProblem::SteadyStateBackwardProblem(
    Solver const& solver, Model& model, SolutionState& final_state,
    gsl::not_null<BwdSimWorkspace*> ws
)
    : xQ_(model.nJ * model.nx_solver, solver.get_sun_context())
    , final_state_(final_state)
    , newton_solver_(NewtonSolver(
          model, solver.get_linear_solver(), solver.get_sun_context()
      ))
    , newton_step_conv_(solver.get_newton_step_steady_state_check())
    , model_(&model)
    , solver_(&solver)
    , ws_(ws) {}

void SteadyStateBackwardProblem::run(realtype const t0) {
    newton_solver_.reinitialize();

    // initialize quadratures
    xQ_.zero();
    ws_->xQB_.zero();

    // Compute quadratures, track computation time
    CpuTimer const cpu_timer;

    compute_steady_state_quadrature(t0);
    cpu_time_b_ = cpu_timer.elapsed_milliseconds();
}

AmiVector const& SteadyStateBackwardProblem::get_adjoint_state() const {
    return ws_->xB_;
}

AmiVector const& SteadyStateBackwardProblem::get_adjoint_quadrature() const {
    return ws_->xQB_;
}

void SteadyStateBackwardProblem::compute_steady_state_quadrature(
    realtype const t0
) {
    // This routine computes the quadratures:
    //     xQB = Integral[ xB(x(t), t, p) * dxdot/dp(x(t), t, p) | dt ]
    // As we're in steady state, we have x(t) = x_ss (x_steadystate), hence
    //     xQB = Integral[ xB(x_ss, t, p) | dt ] * dxdot/dp(x_ss, t, p)
    // We therefore compute the integral over xB first and then do a
    // matrix-vector multiplication.

    auto const sensitivityMode = model_->get_steady_state_sensitivity_mode();

    // Try to compute the analytical solution for quadrature algebraically
    if (sensitivityMode == SteadyStateSensitivityMode::newtonOnly
        || sensitivityMode
               == SteadyStateSensitivityMode::integrateIfNewtonFails)
        compute_quadrature_by_lin_solve();

    // Perform simulation if necessary
    if (sensitivityMode == SteadyStateSensitivityMode::integrationOnly
        || (sensitivityMode
                == SteadyStateSensitivityMode::integrateIfNewtonFails
            && !has_quadrature()))
        compute_quadrature_by_simulation(t0);

    // If the analytic solution and integration did not work, throw
    if (!has_quadrature())
        throw AmiException(
            "Steady state backward computation failed: Linear "
            "system could not be solved (possibly due to singular Jacobian), "
            "and numerical integration did not equilibrate within maxsteps"
        );
}

void SteadyStateBackwardProblem::compute_quadrature_by_lin_solve() {
    // Computes the integral over the adjoint state xB:
    // If the Jacobian has full rank, this has an analytical solution, since
    //   d/dt[ xB(t) ] = JB^T(x(t), p) xB(t) = JB^T(x_ss, p) xB(t)
    // This linear ODE system with time-constant matrix has the solution
    //   xB(t) = exp( t * JB^T(x_ss, p) ) * xB(0)
    // This integral xQ over xB is given as the solution of
    //   JB^T(x_ss, p) * xQ = xB(0)
    // So we first try to solve the linear system, if possible.

    // copy content of xB into vector with integral
    xQ_.copy(ws_->xB_);

    // try to solve the linear system
    try {
        // compute integral over xB and write to xQ
        newton_solver_.prepare_linear_system_b(
            *model_, {final_state_.t, final_state_.x, final_state_.dx}
        );
        newton_solver_.solve_linear_system(xQ_);
        // Compute the quadrature as the inner product xQ * dxdotdp
        computeQBfromQ(
            *model_, xQ_, ws_->xQB_,
            {final_state_.t, final_state_.x, final_state_.dx}
        );
        has_quadrature_ = true;

        // Finalize by setting adjoint state to zero (its steady state)
        ws_->xB_.zero();
    } catch (NewtonFailure const&) {
        has_quadrature_ = false;
    }
}

void SteadyStateBackwardProblem::compute_quadrature_by_simulation(
    realtype const t0
) {
    // If the Jacobian is singular, the integral over xB must be computed
    // by usual integration over time, but simplifications can be applied:
    // x is not time-dependent, no forward trajectory is needed.

    // Set starting timepoint for the simulation solver
    final_state_.t = t0;
    // xQ was written in getQuadratureByLinSolve() -> set to zero
    xQ_.zero();

    auto const sim_solver = std::unique_ptr<Solver>(solver_->clone());
    sim_solver->set_logger(solver_->get_logger());
    sim_solver->set_sensitivity_method(SensitivityMethod::none);
    sim_solver->set_sensitivity_order(SensitivityOrder::none);
    sim_solver->setup(
        t0, model_, ws_->xB_, ws_->dxB_, final_state_.sx, final_state_.sx
    );
    sim_solver->setup_steady_state(
        t0, model_, final_state_.x, final_state_.dx, ws_->xB_, ws_->dxB_, xQ_
    );

    // perform integration and quadrature
    try {
        run_simulation(*sim_solver);
        has_quadrature_ = true;
    } catch (NewtonFailure const&) {
        has_quadrature_ = false;
    }
}

void SteadyStateBackwardProblem::run_simulation(Solver const& solver) {
    if (model_->nx_solver == 0)
        return;

    if (newton_step_conv_) {
        throw NewtonFailure(
            AMICI_NOT_IMPLEMENTED,
            "Newton type convergence check is not implemented for adjoint "
            "steady state computations. Stopping."
        );
    }

    int& sim_steps = num_steps_b_;

    // WRMS computer for xQB
    WRMSComputer wrms_computer_xQB_(
        model_->nplist(), solver.get_sun_context(),
        solver.get_absolute_tolerance_quadratures(),
        solver.get_relative_tolerance_quadratures(), AmiVector()
    );

    // time-derivative of quadrature state vector
    AmiVector xQBdot(model_->nplist(), solver.get_sun_context());

    int const convergence_check_frequency = newton_step_conv_ ? 25 : 1;
    auto const max_steps = (solver.get_max_steps_backward_problem() > 0)
                               ? solver.get_max_steps_backward_problem()
                               : solver.get_max_steps() * 100;

    while (true) {
        if (sim_steps % convergence_check_frequency == 0) {
            // Check for convergence (already before simulation, since we might
            // start in steady state)

            // In the adjoint case, only xQB contributes to the gradient, the
            // exact steadystate is less important, as xB = xQdot may even not
            // converge to zero at all. So we need xQBdot, hence compute xQB
            // first.
            computeQBfromQ(
                *model_, xQ_, ws_->xQB_,
                {final_state_.t, final_state_.x, final_state_.dx}
            );
            computeQBfromQ(
                *model_, ws_->xB_, xQBdot,
                {final_state_.t, final_state_.x, final_state_.dx}
            );
            auto wrms = wrms_computer_xQB_.wrms(xQBdot, ws_->xQB_);
            if (wrms < conv_thresh) {
                break; // converged
            }
        }

        // check for maxsteps
        if (sim_steps >= max_steps) {
            throw IntegrationFailureB(AMICI_TOO_MUCH_WORK, final_state_.t);
        }

        // increase counter
        sim_steps++;

        // One step of ODE integration
        // Reason for tout specification:
        // * max with 1 ensures the correct direction
        //  (any positive value would do)
        // * multiplication with 10 ensures nonzero difference and should
        //   ensure stable computation.
        // The value is not important for AMICI_ONE_STEP mode, only the
        // direction w.r.t. current t.
        solver.step(std::max(final_state_.t, 1.0) * 10);

        solver.write_solution(
            final_state_.t, ws_->xB_, final_state_.dx, final_state_.sx, xQ_
        );
    }
}

} // namespace amici
