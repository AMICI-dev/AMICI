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
                               Model *model, Solver *solver)
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
      x(model->nx_solver),
      x_rdata(model->nx_rdata),
      x_old(model->nx_solver),
      dx(model->nx_solver),
      dx_old(model->nx_solver),
      xdot(model->nx_solver),
      xdot_old(model->nx_solver),
      sx(model->nx_solver,model->nplist()),
      sx_rdata(model->nx_rdata,model->nplist()),
      sdx(model->nx_solver,model->nplist()),
      stau(model->nplist())
{
}


void ForwardProblem::workForwardProblem() {
    bool computeSensitivities =
        solver->getSensitivityOrder() >= SensitivityOrder::first &&
        model->nx_solver > 0;

    model->initialize(x, dx, sx, sdx, computeSensitivities);
    solver->setup(model->t0(), model, x, dx, sx, sdx);
    // update x0 after computing consistence IC, only important for DAEs
    x.copy(solver->getState(model->t0()));

    model->fx_rdata(x_rdata, x);
    if(solver->getSensitivityOrder() >= SensitivityOrder::first) {
        model->fsx_rdata(sx_rdata, sx);
    }

    if(edata){
        rdata->initializeObjectiveFunction();
    }

    /* if preequilibration is necessary, start Newton solver */
    if (solver->getPreequilibration() || (edata && !edata->fixedParametersPreequilibration.empty())) {
        handlePreequilibration();
    } else {
        model->fx_rdata(x_rdata, x);
        rdata->x0 = x_rdata.getVector();
        if (solver->getSensitivityMethod() == SensitivityMethod::forward &&
                solver->getSensitivityOrder() >= SensitivityOrder::first) {
            model->fsx_rdata(sx_rdata, sx);
            for (int ix = 0; ix < rdata->nx; ix++) {
                for (int ip = 0; ip < model->nplist(); ip++)
                    rdata->sx0[ip*rdata->nx + ix] = sx_rdata.at(ix,ip);
            }
        }
    }

    /* perform presimulation if necessary */
    if (edata && edata->t_presim > 0)
        handlePresimulation();

    /* loop over timepoints */
    for (int it = 0; it < model->nt(); it++) {
        auto nextTimepoint = model->getTimepoint(it);

        if (nextTimepoint > model->t0()) {
            if (model->nx_solver == 0) {
                t = nextTimepoint;
                break;
            }

            // Solve for nextTimepoint
            while (t < nextTimepoint) {
                if (std::isinf(nextTimepoint)) {
                    SteadystateProblem sstate = SteadystateProblem(solver, x);
                    sstate.workSteadyStateProblem(rdata, solver, model, it);
                    sstate.writeSolution(&t, x, sx);
                } else {
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
        }
        handleDataPoint(it);
    }

    /* fill events */
    if (model->nz > 0) {
        getEventOutput();
    }

    // set likelihood
    if (!edata) {
        rdata->invalidateLLH();
        rdata->invalidateSLLH();
    }

    storeJacobianAndDerivativeInReturnData();
    rdata->cpu_time = solver->getCpuTime();
}

void ForwardProblem::handlePreequilibration() {
    // Are there dedicated condition preequilibration parameters provided?
    bool overrideFixedParameters =
        edata && !edata->fixedParametersPreequilibration.empty();

    std::vector<realtype>
        originalFixedParameters; // to restore after pre-equilibration

    if (overrideFixedParameters) {
        if (edata->fixedParametersPreequilibration.size() !=
            (unsigned)model->nk())
            throw AmiException(
                "Number of fixed parameters (%d) in model does not match "
                "preequilibration parameters in ExpData (%zd).",
                model->nk(), edata->fixedParametersPreequilibration.size());
        originalFixedParameters = model->getFixedParameters();
        model->setFixedParameters(edata->fixedParametersPreequilibration);
        model->initialize(x, dx, sx, sdx,
                          solver->getSensitivityOrder() >=
                              SensitivityOrder::first);
    }

    // pre-equilibrate
    SteadystateProblem sstate = SteadystateProblem(solver, x);

    sstate.workSteadyStateProblem(rdata, solver, model, -1);
    sstate.writeSolution(&t, x, sx);

    if (overrideFixedParameters) {
        // Restore
        model->setFixedParameters(originalFixedParameters);
    }

    updateAndReinitStatesAndSensitivities(true);
}

void ForwardProblem::updateAndReinitStatesAndSensitivities(bool isSteadystate) {

    if (isSteadystate) {
        model->fx_rdata(x_rdata, x);
        rdata->x_ss = x_rdata.getVector();
    }

    model->fx0_fixedParameters(x);
    solver->reInit(t, x, dx);
    model->fx_rdata(x_rdata, x);

    rdata->x0 = x_rdata.getVector();
    if (solver->getSensitivityOrder() >= SensitivityOrder::first) {
        if (isSteadystate) {
            model->fsx_rdata(sx_rdata, sx);
            for (int ip = 0; ip < model->nplist(); ip++)
                std::copy_n(sx_rdata.data(ip), rdata->nx,
                            &rdata->sx_ss.at(ip * rdata->nx));
        }

        model->fsx0_fixedParameters(sx, x);
        model->fsx_rdata(sx_rdata, sx);

        for (int ip = 0; ip < model->nplist(); ip++)
            std::copy_n(sx_rdata.data(ip), rdata->nx,
                        &rdata->sx0.at(ip * rdata->nx));

        if (solver->getSensitivityMethod() == SensitivityMethod::forward)
            solver->sensReInit(sx, sdx);
    }
}

void ForwardProblem::handlePresimulation()
{
    // Are there dedicated condition preequilibration parameters provided?
    bool overrideFixedParameters = edata && !edata->fixedParametersPresimulation.empty();

    std::vector<realtype> originalFixedParameters; // to restore after pre-equilibration

    if(overrideFixedParameters) {
        if(edata->fixedParametersPresimulation.size() != (unsigned) model->nk())
            throw AmiException("Number of fixed parameters (%d) in model does not match presimulation parameters in ExpData (%zd).",
                           model->nk(), edata->fixedParametersPresimulation.size());
        originalFixedParameters = model->getFixedParameters();
        model->setFixedParameters(edata->fixedParametersPresimulation);
    }
    t = model->t0() - edata->t_presim;
    updateAndReinitStatesAndSensitivities(false);

    solver->run(model->t0());
    solver->writeSolution(&t, x, dx, sx);
    if(overrideFixedParameters) {
        model->setFixedParameters(originalFixedParameters);
    }
    t = model->t0();
    updateAndReinitStatesAndSensitivities(false);
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

        /* get event output */
        model->getEvent(slice(rdata->z, nroots.at(ie), rdata->nz), ie, t, x);
        /* if called from fillEvent at last timepoint,
         then also get the root function value */
        if (t == model->getTimepoint(model->nt() - 1))
            model->getEventRegularization(slice(rdata->rz, nroots.at(ie),
                                                rdata->nz), ie, t, x);

        if (edata) {
            model->getEventSigma(slice(rdata->sigmaz, nroots.at(ie), rdata->nz),
                                 ie, nroots.at(ie), t, edata);
            model->addEventObjective(rdata->llh, ie, nroots.at(ie), t, x,
                                     *edata);

            /* if called from fillEvent at last timepoint,
               add regularization based on rz */
            if (t == model->getTimepoint(model->nt() - 1))
                model->addEventObjectiveRegularization(
                    rdata->llh, ie, nroots.at(ie), t, x, *edata);
        }

        if (solver->getSensitivityOrder() >= SensitivityOrder::first) {
            if (solver->getSensitivityMethod() == SensitivityMethod::forward) {
                getEventSensisFSA(ie);
            } else {
                if (edata) {
                model->getAdjointStateEventUpdate(slice(dJzdx, nroots.at(ie),
                                                        model->nx_solver * model->nJ),
                                                  ie, nroots.at(ie), t, x, *edata);
                model->addPartialEventObjectiveSensitivity(rdata->sllh,
                                                           rdata->s2llh,
                                                           ie, nroots.at(ie),
                                                           t, x, *edata);
                }
            }
        }

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

void ForwardProblem::getEventSensisFSA(int ie) {
    if (t == model->getTimepoint(model->nt() - 1)) {
        // call from fillEvent at last timepoint
        model->getUnobservedEventSensitivity(slice(rdata->sz, nroots.at(ie),
                                                   rdata->nz * rdata->nplist),
                                             ie);
        model->getEventRegularizationSensitivity(slice(rdata->srz,
                                                       nroots.at(ie),
                                                       rdata->nz * rdata->nplist),
                                                 ie, t, x, sx);
    } else {
        model->getEventSensitivity(slice(rdata->sz, nroots.at(ie),
                                         rdata->nz * rdata->nplist),
                                   ie, t, x, sx);
    }

    if (edata) {
        model->addEventObjectiveSensitivity(rdata->sllh, rdata->s2llh, ie,
                                            nroots.at(ie), t, x, sx, *edata);
    }
}

void ForwardProblem::handleDataPoint(int it) {
    model->fx_rdata(x_rdata, x);
    std::copy_n(x_rdata.data(), rdata->nx, &rdata->x.at(it * rdata->nx));
    model->getExpression(slice(rdata->w, it, model->nw),
                model->getTimepoint(it), x);
    if (model->getTimepoint(it) > model->t0()) {
        solver->getDiagnosis(it, rdata);
    }

    getDataOutput(it);
}

void ForwardProblem::getDataOutput(int it) {
    model->getObservable(slice(rdata->y, it, rdata->ny), rdata->ts[it], x);
    model->getObservableSigma(slice(rdata->sigmay, it, rdata->ny), it,
                              edata);
    if (edata) {
        model->addObservableObjective(rdata->llh, it, x, *edata);
        rdata->fres(it, *edata);
        rdata->fchi2(it);
    }

    if (solver->getSensitivityOrder() >= SensitivityOrder::first &&
        model->nplist() > 0) {

        model->getObservableSigmaSensitivity(slice(rdata->ssigmay, it,
                                                   model->nplist() * model->ny),
                                             it, edata);

        if (solver->getSensitivityMethod() == SensitivityMethod::forward) {
            getDataSensisFSA(it);
        } else {
            if (edata) {
                model->getAdjointStateObservableUpdate(slice(dJydx, it,
                                                             model->nx_solver * model->nJ),
                                                       it, x, *edata);
                model->addPartialObservableObjectiveSensitivity(rdata->sllh,
                                                                rdata->s2llh,
                                                                it, x, *edata);
            }
        }
    }
}

void ForwardProblem::getDataSensisFSA(int it) {
    model->fsx_rdata(sx_rdata, sx);
    for (int ix = 0; ix < rdata->nx; ix++) {
        for (int ip = 0; ip < model->nplist(); ip++) {
            rdata->sx[(it * model->nplist() + ip) * rdata->nx + ix] =
                sx_rdata.at(ix, ip);
        }
    }

    model->getObservableSensitivity(slice(rdata->sy, it,
                                          model->nplist() * model->ny),
                                    t, x, sx);

    if (edata) {
        model->addObservableObjectiveSensitivity(rdata->sllh, rdata->s2llh,
                                                 it, x, sx, *edata);
        rdata->fsres(it, *edata);
        rdata->fFIM(it);
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
