#include "amici/forwardproblem.h"

#include "amici/edata.h"
#include "amici/exception.h"
#include "amici/misc.h"
#include "amici/model.h"
#include "amici/solver.h"
#include "amici/steadystateproblem.h"

#include <algorithm>
#include <cmath>
#include <ranges>

namespace amici {

constexpr realtype conv_thresh = 1.0;

/**
 * @brief Check if the next timepoint is too close to the current timepoint.
 *
 * Based on CVODES' `cvHin`.
 * @param cur_t Current time.
 * @param t_next Next stop time.
 * @return True if too close, false otherwise.
 */
bool is_next_t_too_close(realtype cur_t, realtype t_next) {
    auto tdiff = t_next - cur_t;
    if (tdiff == 0.0)
        return true;

    auto tdist = std::fabs(tdiff);
    auto tround = std::numeric_limits<realtype>::epsilon()
                  * std::max(std::fabs(cur_t), std::fabs(t_next));
    if (tdist < 2.0 * tround)
        return true;

    return false;
}

ForwardProblem::ForwardProblem(
    ExpData const* edata, gsl::not_null<Model*> model,
    gsl::not_null<Solver*> solver
)
    : model(model)
    , solver(solver)
    , edata(edata)
    , dJzdx_(model->nJ * model->nx_solver * model->n_max_event(), 0.0)
    , uses_presimulation_(edata && edata->t_presim > 0)
    , ws_(model, solver)
    , main_simulator_(model, solver, &ws_, &dJzdx_)
    , pre_simulator_(model, solver, &ws_, &dJzdx_) {}

void EventHandlingSimulator::run(
    realtype const t0, ExpData const* edata,
    std::vector<realtype> const& timepoints, bool store_diagnosis
) {
    std::ranges::fill(ws_->nroots, 0);

    t0_ = t0;
    ws_->sol.t = t0;
    ws_->tlastroot = t0;

    handle_initial_events(edata);

    // store initial state and sensitivity
    result.initial_state_ = get_simulation_state();
    // store root information at t0
    model_->froot(ws_->sol.t, ws_->sol.x, ws_->sol.dx, ws_->rootvals);

    // get list of trigger timepoints for fixed-time triggered events
    // and filter for timepoints that are within the simulation time range
    auto trigger_timepoints_tmp = model_->get_trigger_timepoints();
    auto trigger_timepoints = std::ranges::views::filter(
        trigger_timepoints_tmp, [this, timepoints](auto t) {
            return t > ws_->sol.t && !timepoints.empty()
                   && t <= timepoints.at(timepoints.size() - 1);
        }
    );
    auto it_trigger_timepoints = trigger_timepoints.begin();

    // loop over timepoints
    for (it_ = 0; it_ < gsl::narrow<int>(timepoints.size()); it_++) {
        // next output time-point
        auto next_t_out = timepoints[it_];

        if (std::isinf(next_t_out)) {
            // post-equilibration is handled elsewhere
            break;
        }

        if (next_t_out > t0) {
            // Solve for next output timepoint
            while (ws_->sol.t < next_t_out) {
                if (is_next_t_too_close(ws_->sol.t, next_t_out)) {
                    // the next timepoint is too close to the current timepoint.
                    // we use the state of the current timepoint.
                    break;
                }

                // next stop time is next output timepoint or next
                // time-triggered event
                auto next_t_event
                    = it_trigger_timepoints != trigger_timepoints.end()
                          ? *it_trigger_timepoints
                          : std::numeric_limits<realtype>::infinity();
                auto const next_t_stop = std::min(next_t_out, next_t_event);

                int const status = solver_->run(next_t_stop);
                // sx will be copied from solver on demand if sensitivities
                // are computed
                solver_->write_solution(ws_->sol);

                if (status == AMICI_ILL_INPUT) {
                    // clustering of roots => turn off root-finding
                    solver_->turn_off_root_finding();
                } else if (status == AMICI_ROOT_RETURN
                           || ws_->sol.t == next_t_event) {
                    // solver-tracked or time-triggered event
                    solver_->get_root_info(ws_->roots_found.data());

                    // check if we are at a trigger timepoint.
                    // if so, set the root-found flag
                    if (ws_->sol.t == next_t_event) {
                        for (auto const ie :
                             model_->get_explicit_roots().at(ws_->sol.t)) {
                            // determine the direction of root crossing from
                            // root function value at the previous event
                            ws_->roots_found[ie]
                                = std::copysign(1, -ws_->rootvals[ie]);
                        }
                        ++it_trigger_timepoints;
                    }

                    handle_events(false, edata);
                }
            }
        }
        handle_datapoint(next_t_out);
        if (store_diagnosis) {
            // store diagnosis information for debugging
            solver_->store_diagnosis();
        }
    }

    // fill events
    if (model_->nz > 0 && model_->nt() > 0) {
        fill_events(model_->n_max_event(), edata);
    }

    result.final_state_ = {.sol = ws_->sol, .mod = model_->get_model_state()};
}

void EventHandlingSimulator::run_steady_state(
    std::function<bool(bool)> check_convergence,
    int convergence_check_frequency, int& sim_steps
) {
    // NOTE: initial events are assumed to be handled elsewhere
    // so we do not reinitialize this->result here
    ws_->tlastroot = ws_->sol.t;

    // get list of trigger timepoints for fixed-time triggered events
    // and filter for timepoints that are within the simulation time range
    auto trigger_timepoints_tmp = model_->get_trigger_timepoints();
    auto trigger_timepoints
        = std::ranges::views::filter(trigger_timepoints_tmp, [this](auto t) {
              return t > ws_->sol.t;
          });
    auto it_trigger_timepoints = trigger_timepoints.begin();
    auto next_t_event = it_trigger_timepoints != trigger_timepoints.end()
                            ? *it_trigger_timepoints
                            : std::numeric_limits<realtype>::infinity();

    while (true) {
        if (sim_steps % convergence_check_frequency == 0) {
            // Check for convergence (already before simulation, since we might
            // start in steady state)
            if (check_convergence(sim_steps > 0))
                break;
        }

        // check for maxsteps
        if (sim_steps >= solver_->get_max_steps()) {
            throw IntegrationFailure(AMICI_TOO_MUCH_WORK, ws_->sol.t);
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
        auto tout = std::isfinite(next_t_event)
                        ? next_t_event
                        : std::max(ws_->sol.t, 1.0) * 10;

        if(!std::isfinite(tout)) {
            // tout overflowed
            throw IntegrationFailure(AMICI_T_OVERFLOW, tout);
        }

        auto status = solver_->step(tout);
        solver_->write_solution(ws_->sol);

        if (status < 0) {
            throw IntegrationFailure(status, ws_->sol.t);
        } else if (status == AMICI_ROOT_RETURN || ws_->sol.t >= next_t_event) {
            if (ws_->sol.t >= next_t_event) {
                // Solver::step will over-step next_t_event
                solver_->write_solution(next_t_event, ws_->sol);
            }

            // solver-tracked or time-triggered event
            solver_->get_root_info(ws_->roots_found.data());

            // check if we are at a trigger timepoint.
            // if so, set the root-found flag
            if (ws_->sol.t == next_t_event) {
                for (auto const ie :
                     model_->get_explicit_roots().at(ws_->sol.t)) {
                    // determine the direction of root crossing from
                    // root function value at the previous event
                    ws_->roots_found[ie] = std::copysign(1, -ws_->rootvals[ie]);
                }
                ++it_trigger_timepoints;
                next_t_event = it_trigger_timepoints != trigger_timepoints.end()
                                   ? *it_trigger_timepoints
                                   : std::numeric_limits<realtype>::infinity();
            }

            handle_events(false, nullptr);
        }
    }

    result.final_state_ = {.sol = ws_->sol, .mod = model_->get_model_state()};
}

void ForwardProblem::run() {
    handle_preequilibration();

    {
        FinalStateStorer fss(this);

        handle_presimulation();
        handle_main_simulation();
    }

    handle_postequilibration();
}

void ForwardProblem::handle_preequilibration() {
    if (!edata || edata->fixed_parameters_pre_equilibration.empty()) {
        return;
    }

    ConditionContext cc2(model, edata, FixedParameterContext::preequilibration);

    preeq_problem_.emplace(&ws_, *solver, *model, true);
    auto t0 = std::isnan(model->t0_preeq()) ? model->t0() : model->t0_preeq();

    // The solver was not run before, set up everything.
    model->initialize(
        t0, ws_.sol.x, ws_.sol.dx, ws_.sol.sx, ws_.sdx,
        solver->get_sensitivity_order() >= SensitivityOrder::first,
        ws_.roots_found
    );
    if (solver->get_sensitivity_order() >= SensitivityOrder::first
        and solver->get_sensitivity_method() != SensitivityMethod::none) {
        model->initialize_state_sensitivities(t0, ws_.sol.sx, ws_.sol.x);
    }
    solver->setup(t0, model, ws_.sol.x, ws_.sol.dx, ws_.sol.sx, ws_.sdx);

    preeq_problem_->run(*solver, -1, t0);

    ws_.sol.x = preeq_problem_->get_state();
    ws_.sol.sx = preeq_problem_->get_state_sensitivity();
    preequilibrated_ = true;
}

void ForwardProblem::handle_presimulation() {
    if (!uses_presimulation_)
        return;

    // Are there dedicated condition preequilibration parameters provided?
    ConditionContext cond(model, edata, FixedParameterContext::presimulation);

    // If we need to reinitialize solver states, this won't work yet.
    if (model->nx_reinit() > 0)
        throw AmiException(
            "Adjoint presimulation with reinitialization of "
            "non-constant states is not yet implemented. Stopping."
        );

    // compute initial time and setup solver for (pre-)simulation
    ws_.sol.t = model->t0() - edata->t_presim;

    // if preequilibration was done, model was already initialized
    if (!preequilibrated_) {
        model->initialize(
            ws_.sol.t, ws_.sol.x, ws_.sol.dx, ws_.sol.sx, ws_.sdx,
            solver->get_sensitivity_order() >= SensitivityOrder::first,
            ws_.roots_found
        );
    } else if (model->ne) {
        // copy, since model state will be updated in reinit_events
        auto h_old = model->get_model_state().h;
        model->reinit_events(
            ws_.sol.t, ws_.sol.x, ws_.sol.dx, h_old, ws_.roots_found
        );
    }
    solver->setup(ws_.sol.t, model, ws_.sol.x, ws_.sol.dx, ws_.sol.sx, ws_.sdx);
    solver->update_and_reinit_states_and_sensitivities(model);

    std::vector<realtype> const timepoints{model->t0()};
    pre_simulator_.run(ws_.sol.t, edata, timepoints, false);
    solver->write_solution(ws_.sol);
}

void ForwardProblem::handle_main_simulation() {
    // When computing adjoint sensitivity analysis with presimulation,
    // we need to store sx after the reinitialization after preequilibration
    // but before reinitialization after presimulation. As presimulation with
    // ASA will not update sx, we can simply extract the values here.
    if (solver->computing_asa() && uses_presimulation_)
        ws_.sol.sx = solver->get_state_sensitivity(model->t0());

    if (!preequilibrated_ && !uses_presimulation_) {
        // if preequilibration or presimulation was done, the model was already
        // initialized
        model->initialize(
            model->t0(), ws_.sol.x, ws_.sol.dx, ws_.sol.sx, ws_.sdx,
            solver->get_sensitivity_order() >= SensitivityOrder::first,
            ws_.roots_found
        );
    }

    ws_.sol.t = model->t0();

    // in case of presimulation, the solver was set up already
    if (!uses_presimulation_) {
        solver->setup(
            ws_.sol.t, model, ws_.sol.x, ws_.sol.dx, ws_.sol.sx, ws_.sdx
        );
    }

    if (preequilibrated_ || uses_presimulation_) {
        // Reset the time and re-initialize events for the main simulation
        solver->update_and_reinit_states_and_sensitivities(model);
        if (model->ne) {
            // copy, since model state will be updated in reinit_events
            auto h_old = model->get_model_state().h;
            model->reinit_events(
                ws_.sol.t, ws_.sol.x, ws_.sol.dx, h_old, ws_.roots_found
            );
        }
    }

    // update x0 after computing consistence IC/reinitialization
    ws_.sol.x = solver->get_state(model->t0());
    // When computing forward sensitivities, we generally want to update sx
    // after presimulation/preequilibration, and if we didn't do either this
    // also won't harm. when computing ASA, we only want to update here if we
    // didn't update before presimulation (if applicable).
    if (solver->computing_fsa()
        || (solver->computing_asa() && !uses_presimulation_))
        ws_.sol.sx = solver->get_state_sensitivity(model->t0());

    main_simulator_.run(ws_.sol.t, edata, model->get_timepoints(), true);
    it_ = main_simulator_.it_;
}

void ForwardProblem::handle_postequilibration() {
    if (get_current_time_iteration() < model->nt()) {
        posteq_problem_.emplace(&ws_, *solver, *model, false);
        auto it = get_current_time_iteration();
        auto t0 = it < 1 ? model->t0() : model->get_timepoint(it - 1);

        // The solver was run before, extract current state from solver.
        solver->write_solution(ws_.sol);
        Expects(t0 == ws_.sol.t);
        posteq_problem_->run(*solver, it, t0);
    }
}

void EventHandlingSimulator::handle_events(
    bool const initial_event, ExpData const* edata
) {
    // Some event triggered. This may be due to some discontinuity, a bolus to
    // be applied, or an event observable to process.

    if (!initial_event && ws_->sol.t == ws_->tlastroot) {
        throw AmiException(
            "AMICI is stuck in an event at time %g, as the initial "
            "step-size after the event is too small. "
            "To fix this, increase absolute and relative "
            "tolerances!",
            ws_->sol.t
        );
    }
    ws_->tlastroot = ws_->sol.t;

    // store the event info and pre-event simulation state
    // whenever a new event is triggered
    auto store_pre_event_info
        = [this, initial_event, edata](bool const seflag) {
              // store Heaviside information at event occurrence
              model_->froot(ws_->sol.t, ws_->sol.x, ws_->sol.dx, ws_->rootvals);

              // store timepoint at which the event occurred, the root function
              // values, and the direction of any zero crossings of the root
              // function
              result.discs.emplace_back(ws_->sol.t, ws_->roots_found);
              ws_->rval_tmp = ws_->rootvals;

              if (model_->nz > 0)
                  store_event(edata);

              store_pre_event_state(seflag, initial_event);
          };

    // store post-event information that is to be saved
    //  not after processing every single event, but after processing all events
    //  that did not trigger a secondary event
    auto store_post_event_info = [this]() {
        if (solver_->computing_asa()) {
            // store updated x to compute jump in discontinuity
            result.discs.back().x_post = ws_->sol.x;
            result.discs.back().dx_post = ws_->sol.dx;
            // Update xdot after the state update
            model_->fxdot(ws_->sol.t, ws_->sol.x, ws_->sol.dx, ws_->xdot);
            result.discs.back().xdot_post = ws_->xdot;
        }
    };

    store_pre_event_info(false);

    if (!initial_event) {
        model_->update_heaviside(ws_->roots_found);
    }

    // Collect all triggered events waiting for execution
    for (int ie = 0; ie < model_->ne; ie++) {
        // only consider transitions false -> true
        if (ws_->roots_found.at(ie) == 1) {
            auto const& event = model_->get_event(ie);
            ws_->pending_events.push(
                {.event = event,
                 .idx = ie,
                 .state_old
                 = (event.uses_values_from_trigger_time()
                        ? std::optional<SimulationState>(get_simulation_state())
                        : std::nullopt)}
            );
        }
    }

    while (!ws_->pending_events.empty()) {
        // get the next event to be handled
        auto const& pending_event = ws_->pending_events.pop();
        auto const ie = pending_event.idx;
        auto const& state_old = pending_event.state_old;

        gsl_Assert(
            // storing the old state is not always necessary,
            // (e.g., if there is only 1 single event and there are no delays)
            // but for now, that's the assumption
            state_old.has_value()
            || !pending_event.event.uses_values_from_trigger_time()
        );

        // TODO: if this is not the first event, check the "persistent"
        // attribute of the event trigger, re-evaluate the trigger if necessary
        // and process or just remove the event from the queue

        // Execute the event
        // Apply bolus to the state and the sensitivities
        model_->add_state_event_update(
            ws_->sol.x, ie, ws_->sol.t, ws_->xdot, ws_->xdot_old,
            state_old.has_value() ? state_old->sol.x : ws_->sol.x,
            state_old.has_value() ? state_old->mod : model_->get_model_state()
        );
        if (solver_->computing_fsa()) {
            // compute the new xdot
            model_->fxdot(ws_->sol.t, ws_->sol.x, ws_->sol.dx, ws_->xdot);
            model_->add_state_sensitivity_event_update(
                ws_->sol.sx, ie, ws_->sol.t, ws_->sol.x, ws_->x_old, ws_->xdot,
                ws_->xdot_old,
                state_old.has_value() ? state_old->sol.sx : ws_->sol.sx,
                ws_->stau
            );
        }

        // check if the event assignment triggered another event
        // and add it to the list of pending events if necessary
        if (detect_secondary_events()) {
            store_post_event_info();

            store_pre_event_info(true);

            model_->update_heaviside(ws_->roots_found);
        }
    }
    store_post_event_info();

    // reinitialize the solver after all events have been processed
    solver_->reinit(ws_->sol.t, ws_->sol.x, ws_->sol.dx);
    if (solver_->computing_fsa()) {
        solver_->sens_reinit(ws_->sol.sx, ws_->sdx);
    }
}

void EventHandlingSimulator::handle_initial_events(ExpData const* edata) {
    if (model_->ne && std::ranges::any_of(ws_->roots_found, [](int const rf) {
            return rf == 1;
        })) {
        handle_events(true, edata);
    }
}

void EventHandlingSimulator::store_event(ExpData const* edata) {
    bool const is_last_timepoint
        = (ws_->sol.t == model_->get_timepoint(model_->nt() - 1));

    if (is_last_timepoint) {
        // call from fillEvent at last timepoint
        model_->froot(ws_->sol.t, ws_->sol.x, ws_->sol.dx, ws_->rootvals);
        for (int ie = 0; ie < model_->ne; ie++) {
            ws_->roots_found.at(ie)
                = (ws_->nroots.at(ie) < model_->n_max_event()) ? 1 : 0;
        }
        result.discs.back().root_info = ws_->roots_found;
    }

    if (get_root_counter() < get_event_counter()) {
        // update stored state (sensi)
        result.event_states_.at(get_root_counter()) = get_simulation_state();
    } else {
        // add stored state (sensi)
        result.event_states_.push_back(get_simulation_state());
    }

    // EVENT OUTPUT
    for (int ie = 0; ie < model_->ne; ie++) {
        // only look for roots of the root function, not discontinuities
        if (ws_->nroots.at(ie) >= model_->n_max_event())
            continue;

        // only consider transitions false -> true or event filling
        if (ws_->roots_found.at(ie) != 1 && !is_last_timepoint) {
            continue;
        }

        if (edata && solver_->computing_asa()) {
            Expects(dJzdx_ != nullptr);
            model_->get_adjoint_state_event_update(
                slice(
                    *dJzdx_, ws_->nroots.at(ie), model_->nx_solver * model_->nJ
                ),
                ie, ws_->nroots.at(ie), ws_->sol.t, ws_->sol.x, *edata
            );
        }
        ws_->nroots.at(ie)++;
    }

    if (is_last_timepoint) {
        // call from fillEvent at last timepoint
        // loop until all events are filled
        fill_events(model_->n_max_event(), edata);
    }
}

void EventHandlingSimulator::store_pre_event_state(
    bool seflag, bool const initial_event
) {
    // If we need to do forward sensitivities later on, we need to store the old
    // x and the old xdot.
    if (solver_->get_sensitivity_order() >= SensitivityOrder::first) {
        // store x and xdot to compute jump in sensitivities
        ws_->x_old.copy(ws_->sol.x);
        model_->fxdot(ws_->sol.t, ws_->sol.x, ws_->sol.dx, ws_->xdot);
        ws_->xdot_old.copy(ws_->xdot);
    }
    if (solver_->computing_fsa()) {
        // compute event-time derivative only for primary events, we get
        // into trouble with multiple simultaneously firing events here (but
        // is this really well-defined then?), in that case, just use the
        // last ie and hope for the best.
        if (!seflag && !initial_event) {
            for (int ie = 0; ie < model_->ne; ie++) {
                // only consider transitions false -> true
                if (ws_->roots_found.at(ie) == 1) {
                    model_->get_event_time_sensitivity(
                        ws_->stau, ws_->sol.t, ie, ws_->sol.x, ws_->sol.sx,
                        ws_->sol.dx
                    );
                }
            }
        }
        if (initial_event) {
            // t0 has no parameter dependency
            std::ranges::fill(ws_->stau, 0.0);
        }
    } else if (solver_->computing_asa()) {
        result.discs.back().xdot_pre = ws_->xdot_old;
        result.discs.back().x_pre = ws_->x_old;
        result.discs.back().h_pre = model_->get_model_state().h;
        result.discs.back().total_cl_pre = model_->get_model_state().total_cl;
    }
}

int EventHandlingSimulator::detect_secondary_events() {
    int secondevent = 0;

    // check whether we need to fire a secondary event
    model_->froot(ws_->sol.t, ws_->sol.x, ws_->sol.dx, ws_->rootvals);
    for (int ie = 0; ie < model_->ne; ie++) {
        // the same event should not trigger itself
        if (ws_->roots_found.at(ie) == 0) {
            // check whether the value of the Heaviside function changed
            if (heaviside_differs(ws_->rval_tmp.at(ie), ws_->rootvals.at(ie))) {
                if (ws_->rval_tmp.at(ie) < ws_->rootvals.at(ie)) {
                    ws_->roots_found.at(ie) = 1;
                    auto const& event = model_->get_event(ie);
                    ws_->pending_events.push(
                        {.event = event,
                         .idx = ie,
                         .state_old
                         = (event.uses_values_from_trigger_time()
                                ? std::optional<SimulationState>(
                                      get_simulation_state()
                                  )
                                : std::nullopt)}
                    );
                } else {
                    ws_->roots_found.at(ie) = -1;
                }
                secondevent++;
            } else {
                ws_->roots_found.at(ie) = 0;
            }
        } else {
            // don't fire the same event again
            ws_->roots_found.at(ie) = 0;
        }
    }

    // fire the secondary event?
    if (secondevent > 0) {
        // Secondary events may result in wrong forward sensitivities
        // if the secondary event has a bolus...
        if (solver_->computing_fsa() && solver_->get_logger())
            solver_->get_logger()->log(
                LogSeverity::warning, "SECONDARY_EVENT",
                "Secondary event was triggered. Depending on "
                "the bolus of the secondary event, forward "
                "sensitivities can be incorrect."
            );
    }

    return secondevent;
}

void EventHandlingSimulator::handle_datapoint(realtype t) {
    // We only store the simulation state if it's not the initial state, as the
    // initial state is stored anyway, and we want to avoid storing it twice
    if (t != t0_ && !result.timepoint_states_.contains(t))
        result.timepoint_states_[t] = get_simulation_state();
}

std::vector<realtype>
ForwardProblem::get_adjoint_updates(Model& model, ExpData const& edata) {
    std::vector<realtype> dJydx(model.nJ * model.nx_solver * model.nt(), 0.0);

    for (int it = 0; it < model.nt(); it++) {
        if (std::isinf(model.get_timepoint(it)))
            break;
        model.get_adjoint_state_observable_update(
            slice(dJydx, it, model.nx_solver * model.nJ), it,
            get_simulation_state_timepoint(it).sol.x, edata
        );
    }

    if (posteq_problem_.has_value()) {
        // Complement dJydx from postequilibration. This shouldn't overwrite
        // anything but only fill in previously 0 values, as only non-inf
        // timepoints were filled above.
        auto const& x = posteq_problem_->get_state();
        for (int it = 0; it < model.nt(); it++) {
            if (std::isinf(model.get_timepoint(it))) {
                model.get_adjoint_state_observable_update(
                    slice(dJydx, it, model.nx_solver * model.nJ), it, x, edata
                );
            }
        }
    }
    return dJydx;
}

SimulationState EventHandlingSimulator::get_simulation_state() {
    return {.sol = ws_->sol, .mod = model_->get_model_state()};
}

std::vector<int> compute_nroots(
    std::vector<Discontinuity> const& discs, int ne, int nmaxevents
) {
    auto size = gsl::narrow<std::vector<int>::size_type>(ne);
    std::vector<int> nroots(size, 0);
    for (auto const& disc : discs) {
        for (std::vector<int>::size_type i = 0; i < size; ++i) {
            if (disc.root_info[i] == 1) {
                nroots[i]++;
            }
        }
    }

    for (auto& n : nroots) {
        if (n > nmaxevents) {
            n = nmaxevents;
        }
    }
    return nroots;
}

/**
 * @brief Assemble the error message to be thrown according to steady state
 * computation status.
 * @param errorString The error string to append to.
 * @param status Status of the steady state computation.
 */
static void
writeErrorString(std::string& errorString, SteadyStateStatus status) {
    switch (status) {
    case SteadyStateStatus::failed_too_long_simulation:
        errorString.append(": System could not be equilibrated.");
        break;
    case SteadyStateStatus::failed_damping:
        errorString.append(": Damping factor reached lower bound.");
        break;
    case SteadyStateStatus::failed_factorization:
        errorString.append(": RHS could not be factorized.");
        break;
    case SteadyStateStatus::failed_convergence:
        errorString.append(": No convergence was achieved.");
        break;
    default:
        errorString.append(".");
        break;
    }
}

SteadyStateProblem::SteadyStateProblem(
    FwdSimWorkspace* ws, Solver const& solver, Model& model, bool const is_preeq
)
    : is_preeq_(is_preeq)
    , ws_(ws)
    , wrms_computer_x_(
          model.nx_solver, solver.get_sun_context(),
          solver.get_absolute_tolerance_steady_state(),
          solver.get_relative_tolerance_steady_state(),
          AmiVector(model.get_steadystate_mask(), solver.get_sun_context())
      )
    , newton_solver_(NewtonSolver(
          model, solver.get_linear_solver(), solver.get_sun_context()
      ))
    , newtons_method_(
          &model, solver.get_sun_context(), &newton_solver_,
          solver.get_newton_damping_factor_mode(),
          solver.get_newton_damping_factor_lower_bound(),
          solver.get_newton_max_steps(),
          solver.get_newton_step_steady_state_check()
      )
    , newton_step_conv_(solver.get_newton_step_steady_state_check())
    , solver_(&solver)
    , model_(&model) {
    // Check for compatibility of options
    if (solver.get_sensitivity_method() == SensitivityMethod::forward
        && solver.get_sensitivity_method_pre_equilibration()
               == SensitivityMethod::adjoint
        && solver.get_sensitivity_order() > SensitivityOrder::none)
        throw AmiException(
            "Preequilibration using adjoint sensitivities "
            "is not compatible with using forward "
            "sensitivities during simulation"
        );
    if (solver.get_sensitivity_method() == SensitivityMethod::forward
        && model.get_steady_state_computation_mode()
               == SteadyStateComputationMode::newtonOnly
        && model.get_steady_state_sensitivity_mode()
               == SteadyStateSensitivityMode::integrationOnly)
        throw AmiException(
            "For forward sensitivity analysis steady-state "
            "computation mode 'newtonOnly' and steady-state "
            "sensitivity mode 'integrationOnly' are not "
            "compatible as numerical integration of the model "
            "ODEs and corresponding forward sensitivities ODEs "
            "is coupled"
        );
}

void SteadyStateProblem::run(Solver& solver, int it, realtype t0) {
    // Compute steady state, track computation time
    CpuTimer cpu_timer;
    ws_->sol.t = t0;

    // store final results at exit, independent of success or failure
    auto const _ = gsl::finally([&]() {
        cpu_time_ = cpu_timer.elapsed_milliseconds();
        period_result_.final_state_.mod = model_->get_model_state();
        period_result_.final_state_.sol = ws_->sol;
    });

    flag_updated_state();
    newton_solver_.reinitialize();
    find_steady_state(it, t0);

    // Check whether state sensitivities still need to be computed.
    // Did we already compute forward sensitivities?
    bool const forwardSensisAlreadyComputed
        = solver.get_sensitivity_order() >= SensitivityOrder::first
          && steady_state_status_[1] == SteadyStateStatus::success
          && (model_->get_steady_state_sensitivity_mode()
                  == SteadyStateSensitivityMode::integrationOnly
              || model_->get_steady_state_sensitivity_mode()
                     == SteadyStateSensitivityMode::integrateIfNewtonFails);
    bool const simulationStartedInSteadystate
        = steady_state_status_[0] == SteadyStateStatus::success
          && numsteps_[0] == 0;
    // Do we need forward sensis for postequilibration?
    // Do we need forward sensitivities for pre- or post-equilibration?
    bool const needForwardSensisPosteq
        = !is_preeq_ && !forwardSensisAlreadyComputed
          && solver.get_sensitivity_order() >= SensitivityOrder::first
          && solver.get_sensitivity_method() == SensitivityMethod::forward;
    bool const needForwardSensisPreeq
        = is_preeq_ && !forwardSensisAlreadyComputed
          && solver.get_sensitivity_method_pre_equilibration()
                 == SensitivityMethod::forward
          && solver.get_sensitivity_order() >= SensitivityOrder::first;

    // Do we need to do the linear system solve to get forward sensitivities?
    if ((needForwardSensisPreeq || needForwardSensisPosteq)
        && !simulationStartedInSteadystate) {
        try {
            // This might still fail, if the Jacobian is singular and
            // simulation did not find a steady state.
            newton_solver_.compute_newton_sensis(
                ws_->sol.sx, *model_, {ws_->sol}
            );
        } catch (NewtonFailure const&) {
            throw AmiException(
                "Steady state sensitivity computation failed due "
                "to unsuccessful factorization of RHS Jacobian"
            );
        }
    }
}

SteadyStateStatus SteadyStateProblem::find_steady_state_by_simulation(
    int const /*it*/, realtype const t0
) {
    if (simulator_.has_value()) {
        // There was already a (failed) Newton attempt
        // reset the workspace to its initial state and restart with a fresh
        // simulator. (Depending on the combination of methods, the old one
        // might not have stored all required information for sensitivity
        // computation.)
        auto const& init_sim_state = simulator_->result.initial_state_;
        model_->set_model_state(init_sim_state.mod);
        ws_->sol = init_sim_state.sol;
    }

    try {
        // Preequilibration -> Create a new solver instance for simulation
        // Postequilibration -> Solver was already created, use that one
        if (is_preeq_) {
            // TODO(performance): We should be able to avoid this clone
            //   if we aren't using ASA in combination with events.
            auto main_solver = solver_;
            auto new_solver = std::unique_ptr<Solver>(main_solver->clone());
            new_solver->set_logger(main_solver->get_logger());

            // do we need sensitivities?
            if (new_solver->get_sensitivity_order() >= SensitivityOrder::first
                && (model_->get_steady_state_sensitivity_mode()
                        == SteadyStateSensitivityMode::integrationOnly
                    || model_->get_steady_state_sensitivity_mode()
                           == SteadyStateSensitivityMode::
                               integrateIfNewtonFails)) {
                // need FSA to compute sx0 for the pre/main simulation,
                // or enable ASA for backward integration of pre-equibration
                new_solver->set_sensitivity_method(
                    new_solver->get_sensitivity_method_pre_equilibration()
                );
            } else {
                new_solver->set_sensitivity_method(SensitivityMethod::none);
                new_solver->set_sensitivity_order(SensitivityOrder::none);
            }
            new_solver->setup(
                t0, model_, ws_->sol.x, ws_->sol.dx, ws_->sol.sx, ws_->sdx
            );
            preeq_solver_unique_ptr_ = std::move(new_solver);
            solver_ = preeq_solver_unique_ptr_.get();
        }

        run_steadystate_simulation_fwd();

        return SteadyStateStatus::success;
    } catch (IntegrationFailure const& ex) {
        switch (ex.error_code) {
        case AMICI_TOO_MUCH_WORK:
            if (model_->get_logger())
                model_->get_logger()->log(
                    LogSeverity::debug, "EQUILIBRATION_FAILURE",
                    "AMICI equilibration exceeded maximum number of"
                    " integration steps at t=%g.",
                    ex.time
                );
            return SteadyStateStatus::failed_convergence;
        case amici::AMICI_T_OVERFLOW:
            if (model_->get_logger())
                model_->get_logger()->log(
                    LogSeverity::debug, "EQUILIBRATION_FAILURE",
                    "AMICI equilibration was stopped after exceedingly"
                    " long simulation time at t=%g.",
                    ex.time
                );
            return SteadyStateStatus::failed_too_long_simulation;
        default:
            if (model_->get_logger())
                model_->get_logger()->log(
                    LogSeverity::debug, "OTHER",
                    "AMICI equilibration failed at t=%g.", ex.time
                );
            return SteadyStateStatus::failed;
        }
    } catch (AmiException const& ex) {
        if (model_->get_logger())
            model_->get_logger()->log(
                LogSeverity::debug, "OTHER", "AMICI equilibration failed: %s",
                ex.what()
            );
        return SteadyStateStatus::failed;
    }
}

[[noreturn]] void SteadyStateProblem::handle_steady_state_failure(
    bool tried_newton_1, bool tried_simulation, bool tried_newton_2
) const {
    // Throw error message according to error codes
    std::string errorString = "Steady state computation failed.";
    if (tried_newton_1) {
        errorString.append(" First run of Newton solver failed");
        writeErrorString(errorString, steady_state_status_[0]);
    }
    if (tried_simulation) {
        errorString.append(" Simulation to steady state failed");
        writeErrorString(errorString, steady_state_status_[1]);
    }
    if (tried_newton_2) {
        errorString.append(" Second run of Newton solver failed");
        writeErrorString(errorString, steady_state_status_[2]);
    }
    throw AmiException(errorString.c_str());
}

realtype SteadyStateProblem::get_wrms_fsa(WRMSComputer& wrms_computer_sx) {
    // Forward sensitivities: Compute weighted error norm for their RHS
    realtype wrms = 0.0;

    // we don't need to call prepareLinearSystem in this function, since it was
    // already computed in the preceding getWrms call and both equations have
    // the same Jacobian.

    xdot_updated_ = false;
    for (int ip = 0; ip < model_->nplist(); ++ip) {
        model_->fsxdot(
            ws_->sol.t, ws_->sol.x, ws_->sol.dx, ip, ws_->sol.sx[ip],
            ws_->sol.dx, ws_->xdot
        );
        if (newton_step_conv_) {
            newton_solver_.solve_linear_system(ws_->xdot);
        }
        wrms = wrms_computer_sx.wrms(ws_->xdot, ws_->sol.sx[ip]);
        // ideally this function would report the maximum of all wrms over
        // all ip, but for practical purposes we can just report the wrms for
        // the first ip where we know that the convergence threshold is not
        // satisfied.
        if (wrms > conv_thresh)
            break;
    }
    // Just report the parameter for the last `ip`. The value doesn't matter -
    // it's only important that all of them satisfy the convergence threshold.
    return wrms;
}

bool SteadyStateProblem::check_steady_state_success() const {
    // Did one of the attempts yield a steady state?
    return std::ranges::any_of(
        steady_state_status_, [](SteadyStateStatus status) {
            return status == SteadyStateStatus::success;
        }
    );
}

void SteadyStateProblem::run_steadystate_simulation_fwd() {
    if (model_->nx_solver == 0)
        return;

    // Do we also have to check for convergence of sensitivities?
    auto sensitivity_method = SensitivityMethod::none;
    if (solver_->get_sensitivity_order() > SensitivityOrder::none
        && solver_->get_sensitivity_method() == SensitivityMethod::forward) {
        sensitivity_method = SensitivityMethod::forward;
    }
    // If forward sensitivity computation by simulation is disabled,
    // disable forward sensitivity integration in the solver.
    // Sensitivities will be computed by newtonsolver.computeNewtonSensis then.
    if (model_->get_steady_state_sensitivity_mode()
        == SteadyStateSensitivityMode::newtonOnly) {
        solver_->switch_forward_sensis_off();
        sensitivity_method = SensitivityMethod::none;
    }

    // function for sensitivity convergence check or dummy
    std::function<bool()> sensi_converged;
    if (solver_->get_sensi_steady_state_check()
        && sensitivity_method == SensitivityMethod::forward) {
        sensi_converged =
            [&,
             wrms_computer_sx = WRMSComputer(
                 model_->nx_solver, solver_->get_sun_context(),
                 solver_->get_absolute_tolerance_steady_state_sensi(),
                 solver_->get_relative_tolerance_steady_state_sensi(),
                 AmiVector(
                     model_->get_steadystate_mask(), solver_->get_sun_context()
                 )
             )]() mutable -> bool {
            update_sensi_simulation();
            // getWrms needs to be called before getWrmsFSA
            // such that the linear system is prepared for newton-type
            // convergence check
            return get_wrms_fsa(wrms_computer_sx) < conv_thresh;
        };
    } else {
        sensi_converged = []() { return true; };
    }

    // Returns the WRMS for the current state
    auto get_wrms_state = [&]() {
        if (newton_step_conv_) {
            newtons_method_.compute_step(ws_->xdot, {ws_->sol});
            return wrms_computer_x_.wrms(
                newtons_method_.get_delta(), ws_->sol.x
            );
        }

        return wrms_computer_x_.wrms(ws_->xdot, ws_->sol.x);
    };

    // Checks for convergence of the steady state solution.
    // Checks state and delta/xdot depending on the options.
    auto check_convergence = [&](bool state_changed) {
        if (state_changed) {
            // update current solution
            flag_updated_state();
            update_rhs();
        }
        wrms_ = get_wrms_state();
        if (wrms_ < conv_thresh && sensi_converged()) {
            return true;
        }
        return false;
    };

    // Create a new simulator that uses the (potentially changed) solver
    // with the right sensitivity settings.
    simulator_.emplace(model_, solver_, ws_, nullptr);
    simulator_->result.initial_state_.mod = model_->get_model_state();
    simulator_->result.initial_state_.sol = ws_->sol;

    // store results on (any) exit
    auto _ = gsl::finally([&]() { period_result_ = simulator_->result; });

    simulator_->handle_initial_events(nullptr);

    update_rhs();

    int const convergence_check_frequency = newton_step_conv_ ? 25 : 1;

    simulator_->run_steady_state(
        check_convergence, convergence_check_frequency, numsteps_.at(1)
    );

    if (numsteps_.at(1)) {
        flag_updated_state();
    }

    // if check_sensi_conv_ is deactivated, we still have to update sensis
    if (sensitivity_method == SensitivityMethod::forward)
        update_sensi_simulation();
}

void SteadyStateProblem::find_steady_state(int it, realtype t0) {
    steady_state_status_.resize(3, SteadyStateStatus::not_run);
    // Turn off Newton's method if 'integrationOnly' approach is chosen for
    // steady-state computation or newton_maxsteps is set to 0 or
    // if 'integrationOnly' approach is chosen for sensitivity computation
    // in combination with forward sensitivities approach. The latter is
    // necessary as numerical integration of the model ODEs and corresponding
    // forward sensitivities ODEs is coupled. If 'integrationOnly' approach is
    // chosen for sensitivity computation it is enforced that steady state is
    // computed only by numerical integration as well.
    bool const turnOffNewton
        = model_->get_steady_state_computation_mode()
              == SteadyStateComputationMode::integrationOnly
          || solver_->get_newton_max_steps() == 0
          || (solver_->get_sensitivity_order() >= SensitivityOrder::first
              && model_->get_steady_state_sensitivity_mode()
                     == SteadyStateSensitivityMode::integrationOnly
              && ((is_preeq_
                   && solver_->get_sensitivity_method_pre_equilibration()
                          == SensitivityMethod::forward)
                  || solver_->get_sensitivity_method()
                         == SensitivityMethod::forward));

    bool const turnOffSimulation = model_->get_steady_state_computation_mode()
                                   == SteadyStateComputationMode::newtonOnly;

    // First, try to run the Newton solver.
    if (!turnOffNewton)
        find_steady_state_by_newtons_method(false);

    // Newton solver didn't work, so try to simulate to steady state.
    if (!turnOffSimulation && !check_steady_state_success())
        steady_state_status_[1] = find_steady_state_by_simulation(it, t0);

    /* Simulation didn't work, retry the Newton solver from last sim state. */
    if (!turnOffNewton && !turnOffSimulation && !check_steady_state_success())
        find_steady_state_by_newtons_method(true);

    // Nothing worked, throw an as informative error as possible.
    if (!check_steady_state_success())
        handle_steady_state_failure(
            !turnOffNewton, !turnOffSimulation,
            !turnOffNewton && !turnOffSimulation
        );
}

void SteadyStateProblem::find_steady_state_by_newtons_method(
    bool newton_retry
) {
    // store initial state and handle initial events, unless that was already
    // done during a previous (failed) simulation
    if (!newton_retry) {
        simulator_.emplace(model_, solver_, ws_, nullptr);
        // store initial state
        simulator_->result.initial_state_.mod = model_->get_model_state();
        simulator_->result.initial_state_.sol = ws_->sol;

        simulator_->handle_initial_events(nullptr);
        period_result_ = simulator_->result;
        // set time to infinity after processing initial events,
        //  so we'll get warnings in case of explicit time dependencies
        //  which are incompatible with Newton's method
        ws_->sol.t = INFINITY;
    }

    int const stage = newton_retry ? 2 : 0;
    try {
        update_rhs();
        newtons_method_.run(ws_->xdot, {ws_->sol}, wrms_computer_x_);
        steady_state_status_[stage] = SteadyStateStatus::success;
        // store final state
        period_result_.final_state_.mod = model_->get_model_state();
        period_result_.final_state_.sol = ws_->sol;
    } catch (NewtonFailure const& ex) {
        switch (ex.error_code) {
        case AMICI_TOO_MUCH_WORK:
            steady_state_status_[stage] = SteadyStateStatus::failed_convergence;
            break;
        case AMICI_NO_STEADY_STATE:
            steady_state_status_[stage]
                = SteadyStateStatus::failed_too_long_simulation;
            break;
        case AMICI_SINGULAR_JACOBIAN:
            steady_state_status_[stage]
                = SteadyStateStatus::failed_factorization;
            break;
        case AMICI_DAMPING_FACTOR_ERROR:
            steady_state_status_[stage] = SteadyStateStatus::failed_damping;
            break;
        default:
            steady_state_status_[stage] = SteadyStateStatus::failed;
            break;
        }
    }
    numsteps_.at(stage) = newtons_method_.get_num_steps();
    wrms_ = newtons_method_.get_wrms();
    flag_updated_state();
}

void SteadyStateProblem::flag_updated_state() {
    xdot_updated_ = false;
    sensis_updated_ = false;
}

void SteadyStateProblem::update_sensi_simulation() {
    if (sensis_updated_)
        return;
    ws_->sol.sx = solver_->get_state_sensitivity(ws_->sol.t);
    sensis_updated_ = true;
}

void SteadyStateProblem::update_rhs() {
    if (xdot_updated_)
        return;
    model_->fxdot(ws_->sol.t, ws_->sol.x, ws_->sol.dx, ws_->xdot);
    xdot_updated_ = true;
}

} // namespace amici
