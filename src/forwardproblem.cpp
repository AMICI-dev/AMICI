#include "amici/cblas.h"
#include "amici/model.h"
#include "amici/forwardproblem.h"
#include "amici/solver.h"
#include "amici/exception.h"
#include "amici/edata.h"
#include "amici/rdata.h"
#include "amici/steadystateproblem.h"

#include <cvodes/cvodes.h> // return/option codes

#include <cmath>
#include <cstring>

namespace amici {

// Ensure AMICI options are in sync with Sundials options
static_assert((int)InternalSensitivityMethod::simultaneous == CV_SIMULTANEOUS, "");
static_assert((int)InternalSensitivityMethod::staggered == CV_STAGGERED, "");
static_assert((int)InternalSensitivityMethod::staggered1 == CV_STAGGERED1, "");

static_assert((int)InterpolationType::hermite == CV_HERMITE, "");
static_assert((int)InterpolationType::polynomial == CV_POLYNOMIAL, "");

static_assert((int)LinearMultistepMethod::adams == CV_ADAMS, "");
static_assert((int)LinearMultistepMethod::BDF == CV_BDF, "");

static_assert(AMICI_ROOT_RETURN == CV_ROOT_RETURN, "");

static_assert((int)NonlinearSolverIteration::functional == CV_FUNCTIONAL, "");
static_assert((int)NonlinearSolverIteration::newton == CV_NEWTON, "");

extern msgIdAndTxtFp warnMsgIdAndTxt;


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
      Jtmp(NewDenseMat(model->nx_solver,model->nx_solver)),
      x(model->nx_solver),
      x_rdata(model->nx_rdata),
      x_old(model->nx_solver),
      dx(model->nx_solver),
      dx_old(model->nx_solver),
      xdot(model->nx_solver),
      xdot_old(model->nx_solver),
      sx(model->nx_solver,model->nplist()),
      sx_rdata(model->nx_rdata,model->nplist()),
      sdx(model->nx_solver,model->nplist())
{
}


ForwardProblem::~ForwardProblem() {
    DestroyMat(Jtmp);
}


void ForwardProblem::workForwardProblem() {

    try {
        solver->setup(&x, &dx, &sx, &sdx, model);
    } catch (std::exception const& ex) {
        throw AmiException("AMICI setup failed:\n(%s)",ex.what());
    } catch (...) {
        throw AmiException("AMICI setup failed due to an unknown error");
    }

    model->fx_rdata(&x_rdata, &x);
    if(solver->getSensitivityOrder() >= SensitivityOrder::first) {
        model->fsx_rdata(&sx_rdata, &sx);
    }

    if(edata){
        rdata->initializeObjectiveFunction();
    }

    /* if preequilibration is necessary, start Newton solver */
    if (solver->getNewtonPreequilibration() || (edata && !edata->fixedParametersPreequilibration.empty())) {
        handlePreequilibration();
    } else {
        model->fx_rdata(&x_rdata, &x);
        rdata->x0 = std::move(x_rdata.getVector());
        if (solver->getSensitivityMethod() == SensitivityMethod::forward &&
            solver->getSensitivityOrder() >= SensitivityOrder::first) {
            model->fsx_rdata(&sx_rdata, &sx);
            for (int ix = 0; ix < rdata->nx; ix++) {
                for (int ip = 0; ip < model->nplist(); ip++)
                    rdata->sx0[ip*rdata->nx + ix] = sx_rdata.at(ix,ip);
            }
        }
    }

    int ncheck = 0; /* the number of (internal) checkpoints stored so far */

    /* perform presimulation if necessary */
    if (edata && edata->t_presim > 0)
        handlePresimulation(&ncheck);

    /* loop over timepoints */
    for (int it = 0; it < model->nt(); it++) {
        auto nextTimepoint = model->t(it);

        if (solver->getSensitivityMethod() == SensitivityMethod::forward &&
            solver->getSensitivityOrder() >= SensitivityOrder::first) {
            solver->setStopTime(nextTimepoint);
        }

        if (nextTimepoint > model->t0()) {
            while (t < nextTimepoint) {
                if (model->nx_solver == 0) {
                    t = nextTimepoint;
                    continue;
                }

                if (std::isinf(nextTimepoint)) {
                    SteadystateProblem sstate = SteadystateProblem(&t, &x, &sx);
                    sstate.workSteadyStateProblem(rdata, solver, model, it);
                } else {
                    int status;
                    if (solver->getSensitivityMethod() == SensitivityMethod::adjoint &&
                            solver->getSensitivityOrder() >= SensitivityOrder::first) {
                        status = solver->solveF(RCONST(nextTimepoint), &x, &dx,
                                                   &(t), AMICI_NORMAL, &ncheck);
                    } else {
                        status = solver->solve(RCONST(nextTimepoint), &x, &dx,
                                                  &(t), AMICI_NORMAL);
                    }

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
}

void ForwardProblem::handlePreequilibration()
{
    // Are there dedicated condition preequilibration parameters provided?
    bool overrideFixedParameters = edata && !edata->fixedParametersPreequilibration.empty();

    std::vector<realtype> originalFixedParameters; // to restore after pre-equilibration

    if(overrideFixedParameters) {
        if(edata->fixedParametersPreequilibration.size() != (unsigned) model->nk())
            throw AmiException("Number of fixed parameters (%d) in model does not match preequilibration parameters in ExpData (%zd).",
                               model->nk(), edata->fixedParametersPreequilibration.size());
        originalFixedParameters = model->getFixedParameters();
        model->setFixedParameters(edata->fixedParametersPreequilibration);
    }

    // pre-equilibrate
    SteadystateProblem sstate = SteadystateProblem(&t, &x, &sx);

    sstate.workSteadyStateProblem(rdata, solver, model, -1);

    if(overrideFixedParameters) { // Restore
        model->setFixedParameters(originalFixedParameters);
    }

    updateAndReinitStatesAndSensitivities();
}

void ForwardProblem::updateAndReinitStatesAndSensitivities() {

    rdata->x_ss = std::move(x_rdata.getVector());

    model->fx0_fixedParameters(&x);
    solver->reInit(t, &x, &dx);
    model->fx_rdata(&x_rdata, &x);

    rdata->x0 = std::move(x_rdata.getVector());
    if (solver->getSensitivityOrder() >= SensitivityOrder::first) {
        for (int ip = 0; ip < model->nplist(); ip++)
            std::copy_n(sx_rdata.data(ip), rdata->nx,
                        &rdata->sx_ss.at(ip * rdata->nx));

        model->fsx0_fixedParameters(&sx, &x);
        model->fsx_rdata(&sx_rdata, &sx);

        for (int ip = 0; ip < model->nplist(); ip++)
            std::copy_n(sx_rdata.data(ip), rdata->nx,
                        &rdata->sx0.at(ip * rdata->nx));

        if (solver->getSensitivityMethod() == SensitivityMethod::forward)
            solver->sensReInit(&sx, &sdx);
    }
}

void ForwardProblem::handlePresimulation(int *ncheck)
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
    updateAndReinitStatesAndSensitivities();

    if (solver->getSensitivityMethod() == SensitivityMethod::adjoint &&
        solver->getSensitivityOrder() >= SensitivityOrder::first) {
        solver->solveF(RCONST(model->t0()), &x, &dx,
                       &(t), AMICI_NORMAL, ncheck);
    } else {
        solver->solve(RCONST(model->t0()), &x, &dx,
                      &(t), AMICI_NORMAL);
    }

    if(overrideFixedParameters) {
        model->setFixedParameters(originalFixedParameters);
    }
    t = model->t0();
    updateAndReinitStatesAndSensitivities();
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

void ForwardProblem::handleEvent(realtype *tlastroot, const bool seflag) {
    /**
     * handleEvent executes everything necessary for the handling of events
     *
     * @param tlastroot pointer to the timepoint of the last event
     */

    /* store heaviside information at event occurence */
    model->froot(t, &x, &dx, rootvals.data());

    if (!seflag) {
        solver->getRootInfo(rootsfound.data());
    }

    if (iroot < model->nMaxEvent() * model->ne) {
        std::copy(rootsfound.begin(), rootsfound.end(), &rootidx[iroot * model->ne]);
    }

    rvaltmp = rootvals;

    if (!seflag) {
        /* only extract in the first event fired */
        if (solver->getSensitivityOrder() >= SensitivityOrder::first &&
            solver->getSensitivityMethod() == SensitivityMethod::forward) {
            solver->getSens(&(t), &sx);
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
        x_old = x;
        if (solver->getSensitivityMethod() == SensitivityMethod::forward) {
            model->fxdot(t, &x, &dx, &xdot);
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
                        model->fstau(t, ie, &x, &sx);
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
        warnMsgIdAndTxt("AMICI:mex:TOO_MUCH_EVENT",
                        "Event was recorded but not reported as the number of "
                        "occured events exceeded (nmaxevents)*(number of "
                        "events in model definition)!");
        solver->reInit(t, &x, &dx); /* reinitialise so that we can continue in peace */
        return;
    }

    if (solver->getSensitivityOrder() >= SensitivityOrder::first
            && solver->getSensitivityMethod() == SensitivityMethod::forward) {
        /* compute the new xdot  */
        model->fxdot(t, &x, &dx, &xdot);
        applyEventSensiBolusFSA();
    }

    int secondevent = 0;

    /* check whether we need to fire a secondary event */
    model->froot(t, &x, &dx, rootvals.data());
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
        solver->reInit(t, &x, &dx);

        /* make time derivative consistent */
        solver->calcIC(t, &x, &dx);

        if (solver->getSensitivityOrder() >= SensitivityOrder::first) {
            if (solver->getSensitivityMethod() == SensitivityMethod::forward) {
                solver->sensReInit(&sx, &sdx);
            }
        }
    }
}

void ForwardProblem::storeJacobianAndDerivativeInReturnData() {
    model->fxdot(t, &x, &dx, &xdot);
    rdata->xdot = xdot.getVector();

    model->fJ(t, 0.0, &x, &dx, &xdot, Jtmp);
    // CVODES uses colmajor, so we need to transform to rowmajor
    for (int ix = 0; ix < model->nx_solver; ix++) {
        for (int jx = 0; jx < model->nx_solver; jx++) {
            rdata->J[ix*model->nx_solver + jx] = Jtmp->data[ix + model->nx_solver*jx];
        }
    }
}

void ForwardProblem::getEventOutput() {
    /**
     * getEventOutput extracts output information for events
     */

    if (t == model->gett(model->nt() - 1,rdata)) {
        // call from fillEvent at last timepoint
        model->froot(t, &x, &dx, rootvals.data());
    }

    /* EVENT OUTPUT */
    for (int ie = 0; ie < model->ne; ie++) {
        /* only look for roots of the rootfunction not discontinuities */
        if (nroots.at(ie) >= model->nMaxEvent())
            continue;

        /* only consider transitions false -> true or event filling */
        if (rootsfound.at(ie) != 1 && t != model->gett(model->nt() - 1, rdata)) {
            continue;
        }

        model->fz(nroots.at(ie), ie, t, &x, rdata);

        if (edata) {
            model->fsigmaz(t, ie, nroots.data(), rdata, edata);
            model->fJz(nroots.at(ie), rdata, edata);

            if (t == model->gett(model->nt() - 1,rdata)) {
                // call from fillEvent at last
                // timepoint, add regularization
                // based on rz
                model->frz(nroots.at(ie), ie, t, &x, rdata);
                model->fJrz(nroots.at(ie), rdata, edata);
            }
        }

        if (solver->getSensitivityOrder() >= SensitivityOrder::first) {
            prepEventSensis(ie);
            if (solver->getSensitivityMethod() == SensitivityMethod::forward) {
                getEventSensisFSA(ie);
            }
        }

        nroots.at(ie)++;
    }

    if (t == model->gett(model->nt() - 1, rdata)) {
        // call from fillEvent at last timepoint
        // loop until all events are filled
        if(std::any_of(nroots.cbegin(), nroots.cend(), [&](int curNRoots){ return curNRoots < model->nMaxEvent(); }))
            getEventOutput();
    }
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

void ForwardProblem::prepEventSensis(int ie) {
    /**
     * prepEventSensis preprocesses the provided experimental data to compute
     * event sensitivities via adjoint or forward methods later on
     *
     * @param ie index of current event
     */

    if(!edata)
        return;

    for (int iz = 0; iz < model->nztrue; iz++) {
        if(model->z2event[iz] - 1 != ie)
            continue;

        model->fdzdp(t, ie, &x);

        model->fdzdx(t, ie, &x);

        if (t == model->gett(model->nt() - 1,rdata)) {
            model->fdrzdp(t, ie, &x);
            model->fdrzdx(t, ie, &x);
        }
        /* extract the value for the standard deviation, if the data
           value is NaN, use the parameter value. Store this value in the return
           struct */
    }
    model->fsigmaz(t, ie, nroots.data(), rdata, edata);
    model->fdsigmazdp(t, ie, nroots.data(), rdata, edata);
    model->fdJzdz(nroots.at(ie), rdata, edata);
    model->fdJzdsigma(nroots.at(ie), rdata, edata);

    if (t == model->t(model->nt() - 1)) {
        model->fdJrzdz(nroots.at(ie), rdata, edata);
        model->fdJrzdsigma(nroots.at(ie), rdata, edata);
    }

    model->fdJzdx(&dJzdx, nroots.at(ie), t, edata, rdata);
    model->fdJzdp(nroots.at(ie), t, edata, rdata);

    if (solver->getSensitivityMethod() == SensitivityMethod::adjoint && model->nz > 0) {
        amici_daxpy(model->nplist(), -1.0, model->dJzdp.data(), 1, rdata->sllh.data(), 1);
        amici_daxpy(model->nplist(), -1.0, &model->dJzdp[1], model->nJ, rdata->s2llh.data(), model->nJ - 1);
    }

}

void ForwardProblem::getEventSensisFSA(int ie) {
    /**
     * getEventSensisFSA extracts event information for forward sensitivity
     * analysis
     *
     * @param ie index of event type
     */
    if (t == model->t(model->nt() - 1)) {
        // call from fillEvent at last timepoint
        model->fsz_tf(nroots.data(),ie, rdata);
        model->fsrz(nroots.at(ie),ie,t,&x,&sx,rdata);
    } else {
        model->fsz(nroots.at(ie),ie,t,&x,&sx,rdata);
    }

    if (edata) {
        model->fsJz(nroots.at(ie),dJzdx,&sx,rdata);
    }
}

void ForwardProblem::handleDataPoint(int it) {
    /**
     * handleDataPoint executes everything necessary for the handling of data
     * points
     *
     * @param it index of data point
     */
    model->fx_rdata(&x_rdata, &x);
    std::copy_n(x_rdata.data(), rdata->nx, &rdata->x.at(it*rdata->nx));

    if (model->t(it) > model->t0()) {
        solver->getDiagnosis(it, rdata);
    }

    getDataOutput(it);
}

void ForwardProblem::getDataOutput(int it) {
    /**
     * getDataOutput extracts output information for data-points
     *
     * @param it index of current timepoint
     */

    model->fy(rdata->ts[it], it, &x, rdata);
    model->fsigmay(it, rdata, edata);
    model->fJy(it, rdata, edata);
    model->fres(it, rdata, edata);
    model->fchi2(it, rdata);

    if (solver->getSensitivityOrder() >= SensitivityOrder::first && model->nplist() > 0) {
        prepDataSensis(it);
        if (solver->getSensitivityMethod() == SensitivityMethod::forward)
            getDataSensisFSA(it);
    }
}

void ForwardProblem::prepDataSensis(int it) {
    /**
     * prepDataSensis preprocesses the provided experimental data to compute
     * sensitivities via adjoint or forward methods later on
     *
     * @param it index of current timepoint
     */

    model->fdydx(rdata->ts[it], &x);
    model->fdydp(rdata->ts[it], &x);

    if (!edata)
        return;

    model->fdsigmaydp(it, rdata, edata);
    model->fdJydy(it, rdata, edata);
    model->fdJydsigma(it, rdata, edata);
    model->fdJydx(&dJydx, it, edata);
    model->fdJydp(it, rdata, edata);
}

void ForwardProblem::getDataSensisFSA(int it) {
    /**
     * getDataSensisFSA extracts data information for forward sensitivity
     * analysis
     *
     * @param it index of current timepoint
     */

    if (!std::isinf(model->t(it)) && model->t(it) > model->t0()) {
        solver->getSens(&(t), &sx);
    }

    model->fsx_rdata(&sx_rdata, &sx);
    for (int ix = 0; ix < rdata->nx; ix++) {
        for (int ip = 0; ip < model->nplist(); ip++) {
            rdata->sx[(it * model->nplist() + ip) * rdata->nx + ix] =
                    sx_rdata.at(ix,ip);
        }
    }

    model->fdsigmaydp(it, rdata, edata);

    model->fsy(it, &sx, rdata);
    if (edata) {
        model->fsJy(it, dJydx, &sx, rdata);
        model->fsres(it, rdata, edata);
        model->fFIM(it, rdata);
    }
}

void ForwardProblem::applyEventBolus() {
    /**
     * applyEventBolus applies the event bolus to the current state
     *
     * @param model pointer to model specification object
     */

    for (int ie = 0; ie < model->ne; ie++) {
        if (rootsfound.at(ie) == 1) {
            /* only consider transitions false -> true */
            model->fdeltax(ie, t, &x, &xdot, &xdot_old);

            amici_daxpy(model->nx_solver, 1.0, model->deltax.data(), 1, x.data(), 1);
        }
    }
}

void ForwardProblem::applyEventSensiBolusFSA() {
    /**
     * applyEventSensiBolusFSA applies the event bolus to the current
     * sensitivities
     */
    for (int ie = 0; ie < model->ne; ie++) {
        if (rootsfound.at(ie) == 1) {
            /* only consider transitions false -> true */
            model->fdeltasx(ie, t, &x_old, &sx, &xdot, &xdot_old);

            for (int ip = 0; ip < model->nplist(); ip++) {
                amici_daxpy(model->nx_solver, 1.0, &model->deltasx[model->nx_solver * ip], 1, sx.data(ip), 1);
            }
        }
    }
}

} // namespace amici
