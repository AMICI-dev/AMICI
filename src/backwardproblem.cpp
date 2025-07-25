#include "amici/backwardproblem.h"

#include "amici/edata.h"
#include "amici/exception.h"
#include "amici/forwardproblem.h"
#include "amici/model.h"
#include "amici/solver.h"
#include "amici/steadystateproblem.h"

namespace amici {
constexpr realtype conv_thresh = 1.0;

BackwardProblem::BackwardProblem(ForwardProblem& fwd)
    : model_(fwd.model)
    , solver_(fwd.solver)
    , edata_(fwd.edata)
    , t_(fwd.getFinalTime())
    , discs_main_(fwd.getDiscontinuities())
    , dJydx_(fwd.getAdjointUpdates(*model_, *edata_))
    , dJzdx_(fwd.getDJzdx())
    , preeq_problem_(fwd.getPreequilibrationProblem())
    , posteq_problem_(fwd.getPostequilibrationProblem())
    , presim_result(fwd.get_presimulation_result())
    , ws_(model_, solver_)
    , simulator_(model_, solver_, &ws_) {}

void BackwardProblem::workBackwardProblem() {

    if (model_->nx_solver <= 0
        || solver_->getSensitivityOrder() < SensitivityOrder::first
        || solver_->getSensitivityMethod() != SensitivityMethod::adjoint
        || model_->nplist() == 0) {
        return;
    }

    handlePostequilibration();

    // handle main simulation

    // If we have posteq, infinity timepoints were already treated
    int it = model_->nt() - 1;
    while (it >= 0 && std::isinf(model_->getTimepoint(it))) {
        --it;
    }

    ws_.discs_ = discs_main_;
    ws_.nroots_ = compute_nroots(discs_main_, model_->ne, model_->nMaxEvent());
    simulator_.run(
        t_, model_->t0(), it, model_->getTimepoints(), &dJydx_, &dJzdx_
    );

    // handle presimulation
    if (edata_ && edata_->t_presim > 0) {
        ConditionContext cc(
            model_, edata_, FixedParameterContext::presimulation
        );
        ws_.discs_ = presim_result.discs;
        ws_.nroots_
            = compute_nroots(ws_.discs_, model_->ne, model_->nMaxEvent());
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
    if (preeq_problem_
        && preeq_problem_->get_solver()->getSensitivityMethodPreequilibration()
               == SensitivityMethod::adjoint) {
        auto preeq_solver = preeq_problem_->get_solver();

        ConditionContext cc2(
            model_, edata_, FixedParameterContext::preequilibration
        );

        // If we need to reinitialize solver states, this won't work yet.
        if (model_->nx_reinit() > 0)
            throw NewtonFailure(
                AMICI_NOT_IMPLEMENTED,
                "Adjoint preequilibration with reinitialization of "
                "non-constant states is not yet implemented. Stopping."
            );

        auto const t0
            = std::isnan(model_->t0Preeq()) ? model_->t0() : model_->t0Preeq();

        auto const& preeq_result = preeq_problem_->get_result();

        // If there were no discontinuities or no simulation was performed,
        // we can use the steady-state shortcuts.
        // If not we need to do the regular backward integration.
        if (preeq_problem_->getSteadyStateStatus()[1]
                == SteadyStateStatus::not_run
            || preeq_result.discs.empty()) {
            preeq_solver->updateAndReinitStatesAndSensitivities(model_);

            auto preeq_final_fwd_state = preeq_result.final_state_;
            preeq_problem_bwd_.emplace(
                *solver_, *model_, preeq_final_fwd_state.sol, &ws_
            );
            preeq_problem_bwd_->run(t0);
        } else if (preeq_problem_->getSteadyStateStatus()[1]
                   != SteadyStateStatus::not_run) {
            // backward integration of the pre-equilibration problem
            // (only if preequilibration was done via simulation)
            ws_.discs_ = preeq_result.discs;
            ws_.nroots_
                = compute_nroots(ws_.discs_, model_->ne, model_->nMaxEvent());
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
}

void BackwardProblem::handlePostequilibration() {
    if (!posteq_problem_) {
        return;
    }

    // initialize xB - only process the post-equilibration timepoints
    for (int it = 0; it < model_->nt(); it++) {
        if (std::isinf(model_->getTimepoint(it))) {
            for (int ix = 0; ix < model_->nxtrue_solver; ix++)
                ws_.xB_[ix] += dJydx_[ix + it * model_->nx_solver];
        }
    }

    auto const& posteq_result = posteq_problem_->get_result();

    // If there were no discontinuities or no simulation was performed,
    // we can use the steady-state shortcuts.
    // If not we need to do the regular backward integration.
    if (posteq_problem_->getSteadyStateStatus()[1] == SteadyStateStatus::not_run
        || posteq_result.discs.empty()) {

        auto final_state = posteq_problem_->getFinalSimulationState();
        posteq_problem_bwd_.emplace(*solver_, *model_, final_state.sol, &ws_);
        posteq_problem_bwd_->run(model_->t0());

        // re-initialize state vectors for main simulation
        model_->initializeB(ws_.xB_, ws_.dxB_, ws_.xQB_, true);
    } else if (posteq_problem_->getSteadyStateStatus()[1]
               != SteadyStateStatus::not_run) {
        // backward integration of the post-equilibration problem
        // (only if post-equilibration was done via simulation)
        ws_.discs_ = posteq_result.discs;
        ws_.nroots_
            = compute_nroots(ws_.discs_, model_->ne, model_->nMaxEvent());
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

void EventHandlingBwdSimulator::handleEventB(
    Discontinuity const& disc, std::vector<realtype> const* dJzdx
) {
    for (int ie = 0; ie < model_->ne; ie++) {

        if (disc.root_info[ie] != 1) {
            continue;
        }

        model_->addAdjointQuadratureEventUpdate(
            ws_->xQB_, ie, t_, disc.x_post, ws_->xB_, disc.xdot_post,
            disc.xdot_pre, disc.x_pre, disc.dx_post
        );
        model_->addAdjointStateEventUpdate(
            ws_->xB_, ie, t_, disc.x_post, disc.xdot_post, disc.xdot_pre,
            disc.x_pre, disc.dx_post
        );

        if (model_->nz > 0) {
            Expects(ws_->nroots_[ie] >= 0);
            for (int ix = 0; ix < model_->nxtrue_solver; ++ix) {
                for (int iJ = 0; iJ < model_->nJ; ++iJ) {
                    ws_->xB_[ix + iJ * model_->nxtrue_solver] += (*dJzdx
                    )[iJ
                      + (ix + ws_->nroots_[ie] * model_->nx_solver)
                            * model_->nJ];
                }
            }
        }

        ws_->nroots_[ie]--;
    }

    // apply pre-event state
    auto state = model_->getModelState();
    state.h = disc.h_pre;
    state.total_cl = disc.total_cl_pre;
    model_->setModelState(state);
}

void EventHandlingBwdSimulator::handleDataPointB(
    int const it, std::vector<realtype> const* dJydx
) {
    // The solver wasn't reset yet, as xB_ is necessary for solver setup.
    // For initial time point (we are integrating backwards!), diagnosis needs
    // to be stored outside this function.
    if (it < model_->nt() - 1)
        solver_->storeDiagnosisB(ws_->which);

    Expects(dJydx != nullptr);
    for (int ix = 0; ix < model_->nxtrue_solver; ix++) {
        for (int iJ = 0; iJ < model_->nJ; iJ++)
            // we only need the 1:nxtrue_solver (not the nx_true) slice here!
            ws_->xB_[ix + iJ * model_->nxtrue_solver]
                += (*dJydx)[iJ + (ix + it * model_->nx_solver) * model_->nJ];
    }
}

realtype EventHandlingBwdSimulator::getTnext(int const it) {
    if (it < 0 && ws_->discs_.empty()) {
        throw AmiException(
            "No more timepoints (it=%d, ie=%d) available at %f. This should "
            "not happen, please report a bug including this stacktrace at "
            "https://github.com/AMICI-dev/AMICI/issues/new/choose",
            it, ws_->discs_.size(), this->t_
        );
    }

    if (!ws_->discs_.empty()
        && (it < 0 || ws_->discs_.back().time > model_->getTimepoint(it))) {
        double const tdisc = ws_->discs_.back().time;
        return tdisc;
    }

    return model_->getTimepoint(it);
}

BwdSimWorkspace::BwdSimWorkspace(
    gsl::not_null<Model*> model, gsl::not_null<Solver const*> solver
)
    : model_(model)
    , xB_(model_->nx_solver, solver->getSunContext())
    , dxB_(model_->nx_solver, solver->getSunContext())
    , xQB_(model_->nJ * model_->nplist(), solver->getSunContext()) {}

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

    // datapoint at t_start?
    if (it >= 0 && timepoints[it] == t_start) {
        handleDataPointB(it, dJydx);
        solver_->setupB(
            &ws_->which, timepoints[it], model_, ws_->xB_, ws_->dxB_, ws_->xQB_
        );
        // for the initial datapoint, diagnosis needs to be stored after setup
        // as it is not called in handleDataPointB
        solver_->storeDiagnosisB(ws_->which);
        --it;
    } else {
        // no data points, only discontinuities, just set up the solver
        // (e.g., during presimulation)
        solver_->setupB(
            &ws_->which, t_start, model_, ws_->xB_, ws_->dxB_, ws_->xQB_
        );
    }

    while (it >= 0 || !ws_->discs_.empty()) {
        // check if next timepoint is a discontinuity or a data-point
        double const tnext = getTnext(it);

        if (tnext < t_) {
            solver_->runB(tnext);
            solver_->writeSolutionB(
                t_, ws_->xB_, ws_->dxB_, ws_->xQB_, ws_->which
            );
        }

        // handle data-point
        if (it >= 0 && tnext == timepoints[it]) {
            handleDataPointB(it, dJydx);
            it--;
        }

        // handle discontinuity
        if (!ws_->discs_.empty() && tnext == ws_->discs_.back().time) {
            handleEventB(ws_->discs_.back(), dJzdx);
            ws_->discs_.pop_back();
        }

        // reinitialize state
        solver_->reInitB(ws_->which, t_, ws_->xB_, ws_->dxB_);
        solver_->quadReInitB(ws_->which, ws_->xQB_);
    }

    // we still need to integrate from first datapoint to t_start
    if (t_ > t_end) {
        solver_->runB(t_end);
        solver_->writeSolutionB(t_, ws_->xB_, ws_->dxB_, ws_->xQB_, ws_->which);
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
    if (model.pythonGenerated) {
        // fill dxdotdp with current values
        auto const& plist = model.getParameterList();
        model.fdxdotdp(state.t, state.x, state.dx);
        model.get_dxdotdp_full().multiply(
            yQB.getNVector(), yQ.getNVector(), plist, true
        );
    } else {
        for (int ip = 0; ip < model.nplist(); ++ip)
            yQB[ip] = dotProd(yQ, model.get_dxdotdp()[ip]);
    }
}

SteadyStateBackwardProblem::SteadyStateBackwardProblem(
    Solver const& solver, Model& model, SolutionState& final_state,
    gsl::not_null<BwdSimWorkspace*> ws
)
    : xQ_(model.nJ * model.nx_solver, solver.getSunContext())
    , final_state_(final_state)
    , newton_solver_(
          NewtonSolver(model, solver.getLinearSolver(), solver.getSunContext())
      )
    , newton_step_conv_(solver.getNewtonStepSteadyStateCheck())
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
    cpu_timeB_ = cpu_timer.elapsed_milliseconds();
}

AmiVector const& SteadyStateBackwardProblem::getAdjointState() const {
    return ws_->xB_;
}

AmiVector const& SteadyStateBackwardProblem::getAdjointQuadrature() const {
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

    auto const sensitivityMode = model_->getSteadyStateSensitivityMode();

    // Try to compute the analytical solution for quadrature algebraically
    if (sensitivityMode == SteadyStateSensitivityMode::newtonOnly
        || sensitivityMode
               == SteadyStateSensitivityMode::integrateIfNewtonFails)
        compute_quadrature_by_lin_solve();

    // Perform simulation if necessary
    if (sensitivityMode == SteadyStateSensitivityMode::integrationOnly
        || (sensitivityMode
                == SteadyStateSensitivityMode::integrateIfNewtonFails
            && !hasQuadrature()))
        compute_quadrature_by_simulation(t0);

    // If the analytic solution and integration did not work, throw
    if (!hasQuadrature())
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
        newton_solver_.prepareLinearSystemB(
            *model_, {final_state_.t, final_state_.x, final_state_.dx}
        );
        newton_solver_.solveLinearSystem(xQ_);
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
    sim_solver->logger = solver_->logger;
    sim_solver->setSensitivityMethod(SensitivityMethod::none);
    sim_solver->setSensitivityOrder(SensitivityOrder::none);
    sim_solver->setup(
        t0, model_, ws_->xB_, ws_->dxB_, final_state_.sx, final_state_.sx
    );
    sim_solver->setupSteadystate(
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

    int& sim_steps = numstepsB_;

    // WRMS computer for xQB
    WRMSComputer wrms_computer_xQB_(
        model_->nplist(), solver.getSunContext(),
        solver.getAbsoluteToleranceQuadratures(),
        solver.getRelativeToleranceQuadratures(), AmiVector()
    );

    // time-derivative of quadrature state vector
    AmiVector xQBdot(model_->nplist(), solver.getSunContext());

    int const convergence_check_frequency = newton_step_conv_ ? 25 : 1;
    auto const max_steps = (solver.getMaxStepsBackwardProblem() > 0)
                               ? solver.getMaxStepsBackwardProblem()
                               : solver.getMaxSteps() * 100;

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

        solver.writeSolution(
            final_state_.t, ws_->xB_, final_state_.dx, final_state_.sx, xQ_
        );
    }
}

} // namespace amici
