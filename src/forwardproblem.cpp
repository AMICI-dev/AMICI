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
    , dJzdx_(model->nJ * model->nx_solver * model->nMaxEvent(), 0.0)
    , t_(model->t0())
    , uses_presimulation_(edata && edata->t_presim > 0)
    , ws_(model, solver)
    , main_simulator_(model, solver, &ws_, &dJzdx_)
    , pre_simulator_(model, solver, &ws_, &dJzdx_) {}

void EventHandlingSimulator::run(
    realtype const t0, ExpData const* edata,
    std::vector<realtype> const& timepoints
) {
    std::ranges::fill(ws_->nroots, 0);

    t0_ = t0;
    t_ = t0;
    ws_->tlastroot = t0;

    // handle initial events
    if (model_->ne && std::ranges::any_of(ws_->roots_found, [](int rf) {
            return rf == 1;
        })) {
        handle_event(true, edata);
    }

    // store initial state and sensitivity
    result.initial_state_ = get_simulation_state();
    // store root information at t0
    model_->froot(t_, ws_->x, ws_->dx, ws_->rootvals);

    // get list of trigger timepoints for fixed-time triggered events
    // and filter for timepoints that are within the simulation time range
    auto trigger_timepoints_tmp = model_->get_trigger_timepoints();
    auto trigger_timepoints = std::ranges::views::filter(
        trigger_timepoints_tmp,
        [this, timepoints](auto t) {
            return t > t_ && t <= timepoints.at(timepoints.size() - 1);
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
            while (t_ < next_t_out) {
                if (is_next_t_too_close(t_, next_t_out)) {
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
                solver_->writeSolution(&t_, ws_->x, ws_->dx, ws_->sx, ws_->dx);

                if (status == AMICI_ILL_INPUT) {
                    // clustering of roots => turn off root-finding
                    solver_->turnOffRootFinding();
                } else if (status == AMICI_ROOT_RETURN || t_ == next_t_event) {
                    // solver-tracked or time-triggered event
                    solver_->getRootInfo(ws_->roots_found.data());

                    // check if we are at a trigger timepoint.
                    // if so, set the root-found flag
                    if (t_ == next_t_event) {
                        for (auto const ie :
                             model_->state_independent_events_[t_]) {
                            // determine the direction of root crossing from
                            // root function value at the previous event
                            ws_->roots_found[ie]
                                = std::copysign(1, -ws_->rootvals[ie]);
                        }
                        ++it_trigger_timepoints;
                    }

                    handle_event(false, edata);
                }
            }
        }
        handle_datapoint(next_t_out);
    }

    // fill events
    if (model_->nz > 0 && model_->nt() > 0) {
        fill_events(model_->nMaxEvent(), edata);
    }

    result.nroots = ws_->nroots;
}

void ForwardProblem::workForwardProblem() {
    handlePreequilibration();

    {
        FinalStateStorer fss(this);

        handlePresimulation();
        handleMainSimulation();
    }

    handlePostequilibration();
}

void ForwardProblem::handlePreequilibration() {
    if (!edata || edata->fixedParametersPreequilibration.empty()) {
        return;
    }

    ConditionContext cc2(model, edata, FixedParameterContext::preequilibration);

    preeq_problem_.emplace(*solver, *model);
    preeq_problem_->workSteadyStateProblem(*solver, *model, -1);

    ws_.x = preeq_problem_->getState();
    ws_.sx = preeq_problem_->getStateSensitivity();
    preequilibrated_ = true;
}

void ForwardProblem::handlePresimulation() {
    if (!uses_presimulation_)
        return;

    if (solver->computingASA()) {
        throw AmiException(
            "Presimulation with adjoint sensitivities"
            " is currently not implemented."
        );
    }

    {
        // Are there dedicated condition preequilibration parameters provided?
        ConditionContext cond(
            model, edata, FixedParameterContext::presimulation
        );

        // compute initial time and setup solver for (pre-)simulation
        t_ = model->t0() - edata->t_presim;

        // if preequilibration was done, model was already initialized
        if (!preequilibrated_) {
            model->initialize(
                t_, ws_.x, ws_.dx, ws_.sx, ws_.sdx,
                solver->getSensitivityOrder() >= SensitivityOrder::first,
                ws_.roots_found
            );
        } else if (model->ne) {
            model->initEvents(t_, ws_.x, ws_.dx, ws_.roots_found);
        }
        solver->setup(t_, model, ws_.x, ws_.dx, ws_.sx, ws_.sdx);
        solver->updateAndReinitStatesAndSensitivities(model);

        std::vector<realtype> const timepoints{model->t0()};
        pre_simulator_.run(t_, edata, timepoints);
        solver->writeSolution(&t_, ws_.x, ws_.dx, ws_.sx, ws_.dx);
    }
}

void ForwardProblem::handleMainSimulation() {
    // When computing adjoint sensitivity analysis with presimulation,
    // we need to store sx after the reinitialization after preequilibration
    // but before reinitialization after presimulation. As presimulation with
    // ASA will not update sx, we can simply extract the values here.
    if (solver->computingASA() && uses_presimulation_)
        ws_.sx = solver->getStateSensitivity(model->t0());

    if (!preequilibrated_ && !uses_presimulation_) {
        // if preequilibration or presimulation was done, the model was already
        // initialized
        model->initialize(
            model->t0(), ws_.x, ws_.dx, ws_.sx, ws_.sdx,
            solver->getSensitivityOrder() >= SensitivityOrder::first,
            ws_.roots_found
        );
    }

    t_ = model->t0();

    solver->setup(t_, model, ws_.x, ws_.dx, ws_.sx, ws_.sdx);

    if (preequilibrated_ || uses_presimulation_) {
        // Reset the time and re-initialize events for the main simulation
        solver->updateAndReinitStatesAndSensitivities(model);
        if (model->ne) {
            model->initEvents(model->t0(), ws_.x, ws_.dx, ws_.roots_found);
        }
    }

    // update x0 after computing consistence IC/reinitialization
    ws_.x = solver->getState(model->t0());
    // When computing forward sensitivities, we generally want to update sx
    // after presimulation/preequilibration, and if we didn't do either this
    // also won't harm. when computing ASA, we only want to update here if we
    // didn't update before presimulation (if applicable).
    if (solver->computingFSA()
        || (solver->computingASA() && !uses_presimulation_))
        ws_.sx = solver->getStateSensitivity(model->t0());

    main_simulator_.run(t_, edata, model->getTimepoints());
    t_ = main_simulator_.t_;
    it_ = main_simulator_.it_;
}

void ForwardProblem::handlePostequilibration() {
    if (getCurrentTimeIteration() < model->nt()) {
        posteq_problem_.emplace(*solver, *model);
        posteq_problem_->workSteadyStateProblem(
            *solver, *model, getCurrentTimeIteration()
        );
    }
}

void EventHandlingSimulator::handle_event(
    bool const initial_event, amici::ExpData const* edata
) {
    // Some event triggered. This may be due to some discontinuity, a bolus to
    // be applied, or an event observable to process.

    if (!initial_event && t_ == ws_->tlastroot) {
        throw AmiException(
            "AMICI is stuck in an event, as the initial "
            "step-size after the event is too small. "
            "To fix this, increase absolute and relative "
            "tolerances!"
        );
    }
    ws_->tlastroot = t_;

    // store the event info and pre-event simulation state
    // whenever a new event is triggered
    auto store_pre_event_info
        = [this, initial_event, edata](bool const seflag) {
              // store Heaviside information at event occurrence
              model_->froot(t_, ws_->x, ws_->dx, ws_->rootvals);

              // store timepoint at which the event occurred, the root function
              // values, and the direction of any zero crossings of the root
              // function
              result.discs.emplace_back(t_, ws_->roots_found);
              ws_->rval_tmp = ws_->rootvals;

              if (model_->nz > 0)
                  store_event(edata);

              store_pre_event_state(seflag, initial_event);

              if (!initial_event) {
                  model_->updateHeaviside(ws_->roots_found);
              }
          };

    // store post-event information that is to be saved
    //  not after processing every single event, but after processing all events
    //  that did not trigger a secondary event
    auto store_post_event_info = [this]() {
        if (solver_->computingASA()) {
            // store updated x to compute jump in discontinuity
            result.discs.back().x_post = ws_->x;
            result.discs.back().xdot_post = ws_->xdot;
        }
    };

    store_pre_event_info(false);

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
        model_->addStateEventUpdate(
            ws_->x, ie, t_, ws_->xdot, ws_->xdot_old,
            state_old.has_value() ? state_old->x : ws_->x,
            state_old.has_value() ? state_old->state : model_->getModelState()
        );
        if (solver_->computingFSA()) {
            // compute the new xdot
            model_->fxdot(t_, ws_->x, ws_->dx, ws_->xdot);
            model_->addStateSensitivityEventUpdate(
                ws_->sx, ie, t_, ws_->x, ws_->x_old, ws_->xdot, ws_->xdot_old,
                state_old.has_value() ? state_old->sx : ws_->sx, ws_->stau
            );
        }

        // check if the event assignment triggered another event
        // and add it to the list of pending events if necessary
        if (detect_secondary_events()) {
            store_post_event_info();
            store_pre_event_info(true);
        }
    }
    store_post_event_info();

    // reinitialize the solver after all events have been processed
    solver_->reInit(t_, ws_->x, ws_->dx);
    if (solver_->computingFSA()) {
        solver_->sensReInit(ws_->sx, ws_->sdx);
    }
}

void EventHandlingSimulator::store_event(ExpData const* edata) {
    bool const is_last_timepoint
        = (t_ == model_->getTimepoint(model_->nt() - 1));

    if (is_last_timepoint) {
        // call from fillEvent at last timepoint
        model_->froot(t_, ws_->x, ws_->dx, ws_->rootvals);
        for (int ie = 0; ie < model_->ne; ie++) {
            ws_->roots_found.at(ie)
                = (ws_->nroots.at(ie) < model_->nMaxEvent()) ? 1 : 0;
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
        if (ws_->nroots.at(ie) >= model_->nMaxEvent())
            continue;

        // only consider transitions false -> true or event filling
        if (ws_->roots_found.at(ie) != 1 && !is_last_timepoint) {
            continue;
        }

        if (edata && solver_->computingASA())
            model_->getAdjointStateEventUpdate(
                slice(
                    *dJzdx_, ws_->nroots.at(ie), model_->nx_solver * model_->nJ
                ),
                ie, ws_->nroots.at(ie), t_, ws_->x, *edata
            );

        ws_->nroots.at(ie)++;
    }

    if (is_last_timepoint) {
        // call from fillEvent at last timepoint
        // loop until all events are filled
        fill_events(model_->nMaxEvent(), edata);
    }
}

void EventHandlingSimulator::store_pre_event_state(
    bool seflag, bool const initial_event
) {
    // If we need to do forward sensitivities later on, we need to store the old
    // x and the old xdot.
    if (solver_->getSensitivityOrder() >= SensitivityOrder::first) {
        // store x and xdot to compute jump in sensitivities
        ws_->x_old.copy(ws_->x);
        model_->fxdot(t_, ws_->x, ws_->dx, ws_->xdot);
        ws_->xdot_old.copy(ws_->xdot);
    }
    if (solver_->computingFSA()) {
        // compute event-time derivative only for primary events, we get
        // into trouble with multiple simultaneously firing events here (but
        // is this really well-defined then?), in that case, just use the
        // last ie and hope for the best.
        if (!seflag && !initial_event) {
            for (int ie = 0; ie < model_->ne; ie++) {
                // only consider transitions false -> true
                if (ws_->roots_found.at(ie) == 1) {
                    model_->getEventTimeSensitivity(
                        ws_->stau, t_, ie, ws_->x, ws_->sx
                    );
                }
            }
        }
        if (initial_event) {
            // t0 has no parameter dependency
            std::ranges::fill(ws_->stau, 0.0);
        }
    } else if (solver_->computingASA()) {
        result.discs.back().xdot_pre = ws_->xdot_old;
    }
}

int EventHandlingSimulator::detect_secondary_events() {
    int secondevent = 0;

    // check whether we need to fire a secondary event
    model_->froot(t_, ws_->x, ws_->dx, ws_->rootvals);
    for (int ie = 0; ie < model_->ne; ie++) {
        // the same event should not trigger itself
        if (ws_->roots_found.at(ie) == 0) {
            // check whether there was a zero-crossing
            if (0 > ws_->rval_tmp.at(ie) * ws_->rootvals.at(ie)) {
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
        if (solver_->computingFSA() && solver_->logger)
            solver_->logger->log(
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
    // store diagnosis information for debugging
    solver_->storeDiagnosis();
}

std::vector<realtype>
ForwardProblem::getAdjointUpdates(Model& model, ExpData const& edata) {
    std::vector<realtype> dJydx(model.nJ * model.nx_solver * model.nt(), 0.0);

    for (int it = 0; it < model.nt(); it++) {
        if (std::isinf(model.getTimepoint(it)))
            break;
        model.getAdjointStateObservableUpdate(
            slice(dJydx, it, model.nx_solver * model.nJ), it,
            getSimulationStateTimepoint(it).x, edata
        );
    }

    if (posteq_problem_.has_value()) {
        // Complement dJydx from postequilibration. This shouldn't overwrite
        // anything but only fill in previously 0 values, as only non-inf
        // timepoints were filled above.
        posteq_problem_->getAdjointUpdates(model, edata, dJydx);
    }
    return dJydx;
}

SimulationState EventHandlingSimulator::get_simulation_state() {
    if (std::isfinite(solver_->gett())) {
        solver_->writeSolution(&t_, ws_->x, ws_->dx, ws_->sx, ws_->dx);
    }
    auto state = SimulationState();
    state.t = t_;
    state.x = ws_->x;
    state.dx = ws_->dx;
    if (solver_->computingFSA() || t_ == t0_)
        state.sx = ws_->sx;
    state.state = model_->getModelState();
    return state;
}

} // namespace amici
