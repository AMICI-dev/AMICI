#include "amici/forwardproblem.h"

#include "amici/edata.h"
#include "amici/exception.h"
#include "amici/misc.h"
#include "amici/model.h"
#include "amici/solver.h"
#include "amici/steadystateproblem.h"

#include <algorithm>
#include <cmath>

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
    , nroots_(gsl::narrow<decltype(nroots_)::size_type>(model->ne), 0)
    , rootvals_(gsl::narrow<decltype(rootvals_)::size_type>(model->ne), 0.0)
    , rval_tmp_(gsl::narrow<decltype(rval_tmp_)::size_type>(model->ne), 0.0)
    , dJzdx_(model->nJ * model->nx_solver * model->nMaxEvent(), 0.0)
    , t_(model->t0())
    , roots_found_(model->ne, 0)
    , x_(model->nx_solver, solver->getSunContext())
    , x_old_(model->nx_solver, solver->getSunContext())
    , dx_(model->nx_solver, solver->getSunContext())
    , xdot_(model->nx_solver, solver->getSunContext())
    , xdot_old_(model->nx_solver, solver->getSunContext())
    , sx_(model->nx_solver, model->nplist(), solver->getSunContext())
    , sdx_(model->nx_solver, model->nplist(), solver->getSunContext())
    , stau_(model->nplist())
    , uses_presimulation_(edata && edata->t_presim > 0) {}

void ForwardProblem::workForwardProblem() {
    handlePreequilibration();

    {
        FinalStateStorer fss(this);

        initialize();
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

    x_ = preeq_problem_->getState();
    sx_ = preeq_problem_->getStateSensitivity();
    preequilibrated_ = true;
}

void ForwardProblem::initialize() {
    // if preequilibration was done, model was already initialized
    if (!preequilibrated_)
        model->initialize(
            x_, dx_, sx_, sdx_,
            solver->getSensitivityOrder() >= SensitivityOrder::first,
            roots_found_
        );
    else if (model->ne) {
        model->initEvents(x_, dx_, roots_found_);
    }

    // compute initial time and setup solver for (pre-)simulation
    auto t0 = model->t0();
    if (uses_presimulation_)
        t0 -= edata->t_presim;
    solver->setup(t0, model, x_, dx_, sx_, sdx_);

    if (model->ne
        && std::ranges::any_of(roots_found_, [](int rf) { return rf == 1; })) {
        handleEvent(t0, true);
    }
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

    if (model->ne > 0) {
        solver->logger->log(
            LogSeverity::warning, "PRESIMULATION",
            "Presimulation with events is not supported. "
            "Events will be ignored during pre- and post-equilibration. "
            "This is subject to change."
        );
    }

    {
        // Are there dedicated condition preequilibration parameters provided?
        ConditionContext cond(
            model, edata, FixedParameterContext::presimulation
        );
        solver->updateAndReinitStatesAndSensitivities(model);

        solver->run(model->t0());
        solver->writeSolution(&t_, x_, dx_, sx_, dx_);
    }

    // Reset the time and re-initialize events for the main simulation
    t_ = model->t0();
    if (model->ne) {
        model->initEvents(x_, dx_, roots_found_);
        if (std::ranges::any_of(roots_found_, [](int rf) { return rf == 1; })) {
            auto t0 = model->t0();
            handleEvent(t0, true);
        }
    }
}

void ForwardProblem::handleMainSimulation() {
    // When computing adjoint sensitivity analysis with presimulation,
    // we need to store sx after the reinitialization after preequilibration
    // but before reinitialization after presimulation. As presimulation with
    // ASA will not update sx, we can simply extract the values here.
    if (solver->computingASA() && uses_presimulation_)
        sx_ = solver->getStateSensitivity(model->t0());

    if (uses_presimulation_ || preequilibrated_)
        solver->updateAndReinitStatesAndSensitivities(model);

    // update x0 after computing consistence IC/reinitialization
    x_ = solver->getState(model->t0());
    // When computing forward sensitivities, we generally want to update sx
    // after presimulation/preequilibration, and if we didn't do either this
    // also wont harm. when computing ASA, we only want to update here, if we
    // didn't update before presimulation (if applicable).
    if (solver->computingFSA()
        || (solver->computingASA() && !uses_presimulation_))
        sx_ = solver->getStateSensitivity(model->t0());

    // store initial state and sensitivity
    initial_state_ = getSimulationState();
    // store root information at t0
    model->froot(t_, x_, dx_, rootvals_);

    // get list of trigger timepoints for fixed-time triggered events
    auto trigger_timepoints = model->get_trigger_timepoints();
    auto it_trigger_timepoints
        = std::ranges::find_if(trigger_timepoints, [this](auto t) {
              return t > this->t_;
          });

    // loop over timepoints
    for (it_ = 0; it_ < model->nt(); it_++) {
        // next output time-point
        auto next_t_out = model->getTimepoint(it_);

        if (std::isinf(next_t_out))
            break;

        if (next_t_out > model->t0()) {
            // Solve for next output timepoint
            while (t_ < next_t_out) {
                if (is_next_t_too_close(t_, next_t_out)) {
                    // next timepoint is too close to current timepoint.
                    // we use the state of the current timepoint.
                    break;
                }

                // next stop time is next output timepoint or next
                // time-triggered event
                auto next_t_event
                    = it_trigger_timepoints != trigger_timepoints.end()
                          ? *it_trigger_timepoints
                          : std::numeric_limits<realtype>::infinity();
                auto next_t_stop = std::min(next_t_out, next_t_event);

                int status = solver->run(next_t_stop);
                // sx will be copied from solver on demand if sensitivities
                // are computed
                solver->writeSolution(&t_, x_, dx_, sx_, dx_);

                if (status == AMICI_ILL_INPUT) {
                    // clustering of roots => turn off root-finding
                    solver->turnOffRootFinding();
                } else if (status == AMICI_ROOT_RETURN || t_ == next_t_event) {
                    // solver-tracked or time-triggered event
                    solver->getRootInfo(roots_found_.data());

                    // check if we are at a trigger timepoint.
                    // if so, set the root-found flag
                    if (t_ == next_t_event) {
                        for (auto ie : model->state_independent_events_[t_]) {
                            // determine direction of root crossing from
                            // root function value at the previous event
                            roots_found_[ie] = std::copysign(1, -rootvals_[ie]);
                        }
                        ++it_trigger_timepoints;
                    }

                    handleEvent(tlastroot_, false);
                }
            }
        }
        handleDataPoint(next_t_out);
    }

    // fill events
    if (model->nz > 0 && model->nt() > 0) {
        fillEvents(model->nMaxEvent());
    }
}

void ForwardProblem::handlePostequilibration() {
    if (getCurrentTimeIteration() < model->nt()) {
        posteq_problem_.emplace(*solver, *model);
        posteq_problem_->workSteadyStateProblem(
            *solver, *model, getCurrentTimeIteration()
        );
    }
}

void ForwardProblem::handleEvent(
    realtype& tlastroot, bool const initial_event
) {
    // Some event triggered. This may be due to some discontinuity, a bolus to
    // be applied, or an event observable to process.

    if (!initial_event) {
        if (t_ == tlastroot) {
            throw AmiException(
                "AMICI is stuck in an event, as the initial "
                "step-size after the event is too small. "
                "To fix this, increase absolute and relative "
                "tolerances!"
            );
        }
        tlastroot = t_;
    }

    // store the event info and pre-event simulation state
    // whenever a new event is triggered
    auto store_pre_event_info = [this, initial_event](bool seflag) {
        // store Heaviside information at event occurrence
        model->froot(t_, x_, dx_, rootvals_);

        // store timepoint at which the event occurred, the root function
        // values, and the direction of any zero crossings of the root function
        discs_.emplace_back(t_, roots_found_);
        rval_tmp_ = rootvals_;

        if (model->nz > 0)
            storeEvent();

        store_pre_event_state(seflag, initial_event);

        if (!initial_event) {
            model->updateHeaviside(roots_found_);
        }
    };

    // store post-event information that is to be saved
    //  not after processing every single event, but after processing all events
    //  that did not trigger a secondary event
    auto store_post_event_info = [this]() {
        if (solver->computingASA()) {
            // store updated x to compute jump in discontinuity
            discs_.back().x_post = x_;
            discs_.back().xdot_post = xdot_;
        }
    };

    store_pre_event_info(false);

    // Collect all triggered events waiting for execution
    for (int ie = 0; ie < model->ne; ie++) {
        // only consider transitions false -> true
        if (roots_found_.at(ie) == 1) {
            auto const& event = model->get_event(ie);
            pending_events_.push(
                {.event = event,
                 .idx = ie,
                 .state_old
                 = (event.uses_values_from_trigger_time()
                        ? std::optional<SimulationState>(getSimulationState())
                        : std::nullopt)}
            );
        }
    }

    while (!pending_events_.empty()) {
        // get the next event to be handled
        auto const& pending_event = pending_events_.pop();
        auto ie = pending_event.idx;
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
        model->addStateEventUpdate(
            x_, ie, t_, xdot_, xdot_old_,
            state_old.has_value() ? state_old->x : x_,
            state_old.has_value() ? state_old->state : model->getModelState()
        );
        if (solver->computingFSA()) {
            // compute the new xdot
            model->fxdot(t_, x_, dx_, xdot_);
            model->addStateSensitivityEventUpdate(
                sx_, ie, t_, x_, x_old_, xdot_, xdot_old_,
                state_old.has_value() ? state_old->sx : sx_, stau_
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
    solver->reInit(t_, x_, dx_);
    if (solver->computingFSA()) {
        solver->sensReInit(sx_, sdx_);
    }
}

void ForwardProblem::storeEvent() {
    bool is_last_timepoint = (t_ == model->getTimepoint(model->nt() - 1));

    if (is_last_timepoint) {
        // call from fillEvent at last timepoint
        model->froot(t_, x_, dx_, rootvals_);
        for (int ie = 0; ie < model->ne; ie++) {
            roots_found_.at(ie) = (nroots_.at(ie) < model->nMaxEvent()) ? 1 : 0;
        }
        discs_.back().root_info = roots_found_;
    }

    if (getRootCounter() < getEventCounter()) {
        // update stored state (sensi)
        event_states_.at(getRootCounter()) = getSimulationState();
    } else {
        // add stored state (sensi)
        event_states_.push_back(getSimulationState());
    }

    // EVENT OUTPUT
    for (int ie = 0; ie < model->ne; ie++) {
        // only look for roots of the rootfunction not discontinuities
        if (nroots_.at(ie) >= model->nMaxEvent())
            continue;

        // only consider transitions false -> true or event filling
        if (roots_found_.at(ie) != 1 && !is_last_timepoint) {
            continue;
        }

        if (edata && solver->computingASA())
            model->getAdjointStateEventUpdate(
                slice(dJzdx_, nroots_.at(ie), model->nx_solver * model->nJ), ie,
                nroots_.at(ie), t_, x_, *edata
            );

        nroots_.at(ie)++;
    }

    if (is_last_timepoint) {
        // call from fillEvent at last timepoint
        // loop until all events are filled
        fillEvents(model->nMaxEvent());
    }
}

void ForwardProblem::store_pre_event_state(bool seflag, bool initial_event) {
    // If we need to do forward sensitivities later on we need to store the old
    // x and the old xdot.
    if (solver->getSensitivityOrder() >= SensitivityOrder::first) {
        // store x and xdot to compute jump in sensitivities
        x_old_.copy(x_);
        model->fxdot(t_, x_, dx_, xdot_);
        xdot_old_.copy(xdot_);
    }
    if (solver->computingFSA()) {
        // compute event-time derivative only for primary events, we get
        // into trouble with multiple simultaneously firing events here (but
        // is this really well defined then?), in that case just use the
        // last ie and hope for the best.
        if (!seflag && !initial_event) {
            for (int ie = 0; ie < model->ne; ie++) {
                // only consider transitions false -> true
                if (roots_found_.at(ie) == 1) {
                    model->getEventTimeSensitivity(stau_, t_, ie, x_, sx_);
                }
            }
        }
        if (initial_event) {
            // t0 has no parameter dependency
            std::ranges::fill(stau_, 0.0);
        }
    } else if (solver->computingASA()) {
        discs_.back().xdot_pre = xdot_old_;
    }
}

int ForwardProblem::detect_secondary_events() {
    int secondevent = 0;

    // check whether we need to fire a secondary event
    model->froot(t_, x_, dx_, rootvals_);
    for (int ie = 0; ie < model->ne; ie++) {
        // the same event should not trigger itself
        if (roots_found_.at(ie) == 0) {
            // check whether there was a zero-crossing
            if (0 > rval_tmp_.at(ie) * rootvals_.at(ie)) {
                if (rval_tmp_.at(ie) < rootvals_.at(ie)) {
                    roots_found_.at(ie) = 1;
                    auto const& event = model->get_event(ie);
                    pending_events_.push(
                        {.event = event,
                         .idx = ie,
                         .state_old
                         = (event.uses_values_from_trigger_time()
                                ? std::optional<SimulationState>(
                                      getSimulationState()
                                  )
                                : std::nullopt)}
                    );
                } else {
                    roots_found_.at(ie) = -1;
                }
                secondevent++;
            } else {
                roots_found_.at(ie) = 0;
            }
        } else {
            // don't fire the same event again
            roots_found_.at(ie) = 0;
        }
    }

    // fire the secondary event?
    if (secondevent > 0) {
        // Secondary events may result in wrong forward sensitivities,
        // if the secondary event has a bolus...
        if (solver->computingFSA() && solver->logger)
            solver->logger->log(
                LogSeverity::warning, "SECONDARY_EVENT",
                "Secondary event was triggered. Depending on "
                "the bolus of the secondary event, forward "
                "sensitivities can be incorrect."
            );
    }

    return secondevent;
}

void ForwardProblem::handleDataPoint(realtype t) {
    // We only store the simulation state if it's not the initial state, as the
    // initial state is stored anyway and we want to avoid storing it twice
    if (t != model->t0() && timepoint_states_.count(t) == 0)
        timepoint_states_[t] = getSimulationState();
    // store diagnosis information for debugging
    solver->storeDiagnosis();
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

SimulationState ForwardProblem::getSimulationState() {
    if (std::isfinite(solver->gett())) {
        solver->writeSolution(&t_, x_, dx_, sx_, dx_);
    }
    auto state = SimulationState();
    state.t = t_;
    state.x = x_;
    state.dx = dx_;
    if (solver->computingFSA() || t_ == model->t0())
        state.sx = sx_;
    state.state = model->getModelState();
    return state;
}

} // namespace amici
