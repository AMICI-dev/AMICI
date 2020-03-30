#include "amici/forwardproblem.h"

#include "amici/cblas.h"
#include "amici/misc.h"
#include "amici/model.h"
#include "amici/solver.h"
#include "amici/exception.h"
#include "amici/edata.h"
#include "amici/rdata.h"
#include "amici/steadystateproblem.h"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace amici {

ForwardProblem::ForwardProblem(ReturnData *rdata, const ExpData *edata,
                               Model *model, Solver *solver,
                               const SteadystateProblem *preeq)
    : model(model),
      rdata(rdata),
      solver(solver),
      edata(edata),
      rootidx(static_cast<decltype (rootidx)::size_type>(model->ne * model->ne * model->nMaxEvent()), 0),
      nroots(static_cast<decltype (nroots)::size_type>(model->ne), 0),
      rootvals(static_cast<decltype (rootvals)::size_type>(model->ne), 0.0),
      rvaltmp(static_cast<decltype (rvaltmp)::size_type>(model->ne), 0.0),
      discs(static_cast<decltype (discs)::size_type>(model->nMaxEvent() * model->ne), 0.0),
      irdiscs(model->nMaxEvent() * model->ne, 0.0),
      x_disc(model->nx_solver, model->nMaxEvent()*model->ne),
      xdot_disc(model->nx_solver, model->nMaxEvent()*model->ne),
      xdot_old_disc(model->nx_solver, model->nMaxEvent()*model->ne),
      dJydx(model->nJ * model->nx_solver * model->nt(), 0.0),
      dJzdx(model->nJ * model->nx_solver * model->nMaxEvent(), 0.0),
      t(model->t0()),
      rootsfound(model->ne, 0),
      Jtmp(SUNMatrixWrapper(model->nx_solver,model->nx_solver)),
      
      x0(model->nx_solver),
      sx0(model->nx_solver,model->nplist()),
      x(model->nx_solver),
      x_old(model->nx_solver),
      dx(model->nx_solver),
      dx_old(model->nx_solver),
      xdot(model->nx_solver),
      xdot_old(model->nx_solver),
      sx(model->nx_solver,model->nplist()),
      sdx(model->nx_solver,model->nplist()),
      stau(model->nplist())
{
    if (preeq) {
        preeq->writeSolution(&t, x, sx);
        preequilibrated = true;
    }
}


void ForwardProblem::workForwardProblem() {
    if(model->nx_solver == 0){
        return;
    }

    /* if preequilibration was done, solver was already set up */
    if (!preequilibrated) {
        model->initialize(x, dx, sx, sdx,
                          solver->getSensitivityOrder() >=
                          SensitivityOrder::first);
        // handlePreequilibration will setup the solver on its own.
        solver->setup(model->t0(), model, x, dx, sx, sdx);
        // update x0 after computing consistence IC, only important for DAEs
        x.copy(solver->getState(model->t0()));
    }

    auto presimulated = edata && edata->t_presim > 0;
    /* perform presimulation if necessary */
    if (presimulated) {
        handlePresimulation();
        t = model->t0();
        solver->updateAndReinitStatesAndSensitivities(model);
    }
    
    /* store initial sensitivity and */
    x0.copy(x);
    if(solver->getSensitivityOrder() >= SensitivityOrder::first &&
       (!presimulated || solver->getSensitivityMethod() == SensitivityMethod::forward))
        sx0.copy(sx);

    /* loop over timepoints */
    for (it = 0; it < model->nt(); it++) {
        auto nextTimepoint = model->getTimepoint(it);
        
        if (std::isinf(nextTimepoint))
            break;

        if (nextTimepoint > model->t0()) {
            // Solve for nextTimepoint
            while (t < nextTimepoint) {
                int status = solver->run(nextTimepoint);
                solver->writeSolution(&t, x, dx, sx);
                /* sx will be copied from solver on demand if sensitivities
                 are computed */
                if (status == AMICI_ILL_INPUT) {
                    /* clustering of roots => turn off rootfinding */
                    solver->turnOffRootFinding();
                } else if (status == AMICI_ROOT_RETURN) {
                    handleEvent(&tlastroot, false);
                }
            }
        }
        handleDataPoint(it);
    }

    /* fill events */
    if (model->nz > 0 & model->nt() > 0) {
        getEventOutput();
    }

    storeJacobianAndDerivativeInReturnData();
    rdata->cpu_time = solver->getCpuTime();
}

void ForwardProblem::handlePresimulation()
{
    // Are there dedicated condition preequilibration parameters provided?
    auto cond = ConditionContext(model, edata,
                                 FixedParameterContext::presimulation);
    solver->updateAndReinitStatesAndSensitivities(model);

    solver->run(model->t0());
    solver->writeSolution(&t, x, dx, sx);
}


void ForwardProblem::handleEvent(realtype *tlastroot, const bool seflag) {
    /* store heaviside information at event occurence */
    model->froot(t, x, dx, rootvals);

    if (!seflag) {
        solver->getRootInfo(rootsfound.data());
    }

    if (iroot < model->nMaxEvent() * model->ne) {
        std::copy(rootsfound.begin(), rootsfound.end(),
                  &rootidx[iroot * model->ne]);
    }

    rvaltmp = rootvals;

    if (!seflag) {
        /* only extract in the first event fired */
        if (solver->getSensitivityOrder() >= SensitivityOrder::first &&
            solver->getSensitivityMethod() == SensitivityMethod::forward) {
            sx.copy(solver->getStateSensitivity(t));
        }

        /* only check this in the first event fired, otherwise this will always
         * be true */
        if (t == *tlastroot) {
            throw AmiException("AMICI is stuck in an event, as the initial"
                               "step-size after the event is too small. To fix "
                               "this, increase absolute and relative tolerances!");
        }
        *tlastroot = t;
    }

    if(model->nz>0)
        getEventOutput();

    /* if we need to do forward sensitivities later on we need to store the old
     * x and the old xdot */
    if (solver->getSensitivityOrder() >= SensitivityOrder::first) {
        /* store x and xdot to compute jump in sensitivities */
        x_old = solver->getState(t);
        if (solver->getSensitivityMethod() == SensitivityMethod::forward) {
            model->fxdot(t, x, dx, xdot);
            xdot_old = xdot;
            dx_old = dx;

            /* compute event-time derivative only for primary events, we get
             * into trouble with multiple simultaneously firing events here (but
             * is this really well defined then?), in that case just use the
             * last ie and hope for the best. */
            if (!seflag) {
                for (int ie = 0; ie < model->ne; ie++) {
                    if (rootsfound.at(ie) == 1) {
                        /* only consider transitions false -> true */
                        model->getEventTimeSensitivity(stau, t, ie, x, sx);
                    }
                }
            }
        } else if (solver->getSensitivityMethod() == SensitivityMethod::adjoint) {
            /* store x to compute jump in discontinuity */
            if (iroot < model->nMaxEvent() * model->ne) {
                x_disc[iroot] = x;
                xdot_disc[iroot] = xdot;
                xdot_old_disc[iroot] = xdot_old;
            }
        }
    }

    model->updateHeaviside(rootsfound);

    applyEventBolus();

    if (iroot < model->nMaxEvent() * model->ne) {
        discs[iroot] = t;
        ++iroot;
    } else {
        solver->app->warning(
                    "AMICI:TOO_MUCH_EVENT",
                        "Event was recorded but not reported as the number of "
                        "occured events exceeded (nmaxevents)*(number of "
                        "events in model definition)!");
        /* reinitialise so that we can continue in peace */
        solver->reInit(t, x, dx);
        return;
    }

    if (solver->getSensitivityOrder() >= SensitivityOrder::first
            && solver->getSensitivityMethod() == SensitivityMethod::forward) {
        /* compute the new xdot  */
        model->fxdot(t, x, dx, xdot);
        applyEventSensiBolusFSA();
    }

    int secondevent = 0;

    /* check whether we need to fire a secondary event */
    model->froot(t, x, dx, rootvals);
    for (int ie = 0; ie < model->ne; ie++) {
        /* the same event should not trigger itself */
        if (rootsfound.at(ie) == 0) {
            /* check whether there was a zero-crossing */
            if (0 > rvaltmp.at(ie) * rootvals.at(ie)) {
                if (rvaltmp.at(ie) < rootvals.at(ie)) {
                    rootsfound.at(ie) = 1;
                } else {
                    rootsfound.at(ie) = -1;
                }
                secondevent++;
            } else {
                rootsfound.at(ie) = 0;
            }
        } else {
            /* don't fire the same event again */
            rootsfound.at(ie) = 0;
        }
    }
    /* fire the secondary event */
    if (secondevent > 0) {
        handleEvent(tlastroot, TRUE);
    }

    /* only reinitialise in the first event fired */
    if (!seflag) {
        solver->reInit(t, x, dx);

        if (solver->getSensitivityOrder() >= SensitivityOrder::first) {
            if (solver->getSensitivityMethod() == SensitivityMethod::forward) {
                solver->sensReInit(sx, sdx);
            }
        }
    }
}

void ForwardProblem::storeJacobianAndDerivativeInReturnData() {
    model->fxdot(t, x, dx, xdot);
    rdata->xdot = xdot.getVector();

    model->fJ(t, 0.0, x, dx, xdot, Jtmp.get());
    // CVODES uses colmajor, so we need to transform to rowmajor
    for (int ix = 0; ix < model->nx_solver; ix++) {
        for (int jx = 0; jx < model->nx_solver; jx++) {
            rdata->J[ix * model->nx_solver + jx] =
                Jtmp.data()[ix + model->nx_solver * jx];
        }
    }
}

void ForwardProblem::getEventOutput() {
    if (t == model->getTimepoint(model->nt() - 1)) {
        // call from fillEvent at last timepoint
        model->froot(t, x, dx, rootvals);
        discs[iroot] = t;
        ++iroot;
    }
    
    if (iroot > x_events.size()) {
        /* update stored state (sensi) */
        x_events.at(iroot-1) = x;
        if (solver->getSensitivityOrder() >= SensitivityOrder::first &&
        solver->getSensitivityMethod() == SensitivityMethod::forward)
            sx_events.at(iroot-1) = sx;
    } else {
        /* add stored state (sensi) */
        x_events.push_back(x);
        if (solver->getSensitivityOrder() >= SensitivityOrder::first &&
            solver->getSensitivityMethod() == SensitivityMethod::forward)
            sx_events.push_back(sx);
    }

    /* EVENT OUTPUT */
    for (int ie = 0; ie < model->ne; ie++) {
        /* only look for roots of the rootfunction not discontinuities */
        if (nroots.at(ie) >= model->nMaxEvent())
            continue;

        /* only consider transitions false -> true or event filling */
        if (rootsfound.at(ie) != 1 &&
            t != model->getTimepoint(model->nt() - 1)) {
            continue;
        }
        
        if (edata && solver->getSensitivityOrder() >= SensitivityOrder::first &&
            solver->getSensitivityMethod() == SensitivityMethod::adjoint)
            model->getAdjointStateEventUpdate(slice(dJzdx, nroots.at(ie),
                                                    model->nx_solver * model->nJ),
                                              ie, nroots.at(ie), t, x, *edata);

        nroots.at(ie)++;
    }

    if (t == model->getTimepoint(model->nt() - 1)) {
        // call from fillEvent at last timepoint
        // loop until all events are filled
        if (std::any_of(nroots.cbegin(), nroots.cend(), [&](int curNRoots) {
                return curNRoots < model->nMaxEvent();
            }))
            getEventOutput();
    }
}

void ForwardProblem::handleDataPoint(int it) {
    
    x_timepoints.push_back(x);
    
    if (solver->getSensitivityMethod() == SensitivityMethod::forward &&
        solver->getSensitivityOrder() >= SensitivityOrder::first) {
        sx_timepoints.push_back(sx);
    }
    
    if (edata && solver->getSensitivityOrder() >= SensitivityOrder::first &&
        solver->getSensitivityMethod() == SensitivityMethod::adjoint) {
        model->getAdjointStateObservableUpdate(
            slice(dJydx, it, model->nx_solver * model->nJ), it, x, *edata
        );
    }
    
    if (model->getTimepoint(it) > model->t0()) {
        solver->getDiagnosis(it, rdata);
    }
}

void ForwardProblem::applyEventBolus() {
    for (int ie = 0; ie < model->ne; ie++)
        if (rootsfound.at(ie) == 1) // only consider transitions false -> true
            model->addStateEventUpdate(x, ie, t, xdot, xdot_old);
}

void ForwardProblem::applyEventSensiBolusFSA() {
    for (int ie = 0; ie < model->ne; ie++)
        if (rootsfound.at(ie) == 1) // only consider transitions false -> true
            /*  */
            model->addStateSensitivityEventUpdate(sx, ie, t, x_old, xdot,
                                                  xdot_old, stau);
}

} // namespace amici
