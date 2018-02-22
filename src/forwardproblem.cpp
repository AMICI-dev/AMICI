#include "include/amici_model.h"
#include "include/forwardproblem.h"
#include "include/amici_solver.h"
#include "include/amici_exception.h"
#include "include/edata.h"
#include "include/rdata.h"
#include "include/steadystateproblem.h"
#include <cvodes/cvodes.h> // return/option codes

#include <cmath>
#include <cstring>

namespace amici {

// Ensure AMICI options are in sync with Sundials options
static_assert(InternalSensitivityMethod::SIMULTANEOUS == CV_SIMULTANEOUS, "");
static_assert(InternalSensitivityMethod::STAGGERED == CV_STAGGERED, "");
static_assert(InternalSensitivityMethod::STAGGERED1 == CV_STAGGERED1, "");

static_assert(InterpolationType::HERMITE == CV_HERMITE, "");
static_assert(InterpolationType::POLYNOMIAL == CV_POLYNOMIAL, "");

static_assert(LinearMultistepMethod::ADAMS == CV_ADAMS, "");
static_assert(LinearMultistepMethod::BDF == CV_BDF, "");

static_assert(AMICI_ROOT_RETURN == CV_ROOT_RETURN, "");
    
static_assert(NonlinearSolverIteration::FUNCTIONAL == CV_FUNCTIONAL, "");
static_assert(NonlinearSolverIteration::NEWTON == CV_NEWTON, "");
    
    /**
     * default constructor
     * @param rdata pointer to ReturnData instance
     * @param edata pointer to ExpData instance
     * @param model pointer to Model instance
     * @param solver pointer to Solver instance
     *
     */
    ForwardProblem::ForwardProblem(ReturnData *rdata, const ExpData *edata,
                   Model *model, Solver *solver) :
    rootidx(model->ne*model->ne*model->ne*rdata->nmaxevent, 0),
    nroots(model->ne, 0),
    rootvals(model->ne, 0.0),
    rvaltmp(model->ne, 0.0),
    discs(rdata->nmaxevent * model->ne, 0.0),
    irdiscs(rdata->nmaxevent * model->ne, 0.0),
    x_disc(model->nx,model->nMaxEvent()*model->ne),
    xdot_disc(model->nx,model->nMaxEvent()*model->ne),
    xdot_old_disc(model->nx,model->nMaxEvent()*model->ne),
    dJydx(model->nJ * model->nx * model->nt(), 0.0),
    dJzdx(model->nJ * model->nx * model->nMaxEvent(), 0.0),
    rootsfound(model->ne, 0),
    x(model->nx),
    x_old(model->nx),
    dx(model->nx),
    dx_old(model->nx),
    xdot(model->nx),
    xdot_old(model->nx),
    sx(model->nx,model->nplist()),
    sdx(model->nx,model->nplist())
    {
        t = model->t0();
        this->model = model;
        this->solver = solver;
        this->edata = edata;
        this->rdata = rdata;
        Jtmp = NewDenseMat(model->nx,model->nx);
    }

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

void ForwardProblem::workForwardProblem() {
    /**
     * workForwardProblem solves the forward problem. if forward sensitivities
     * are enabled this will also compute sensitivies
     *
     */

    try {
        solver->setupAMI(this, model);
    } catch (std::exception& ex) {
        throw AmiException("AMICI setup failed:\n(%s)",ex.what());
    } catch (...) {
        throw AmiException("AMICI setup failed due to an unknown error");
    }
    int ncheck = 0; /* the number of (internal) checkpoints stored so far */
    realtype tlastroot = 0; /* storage for last found root */

    /* if preequilibration is necessary, start Newton solver */
    if (solver->getNewtonPreequilibration()) {
        SteadystateProblem sstate = SteadystateProblem(&t,&x,&sx);
        sstate.workSteadyStateProblem(rdata,
                                       solver, model, -1);
    }

    /* loop over timepoints */
    for (int it = 0; it < rdata->nt; it++) {
        if (rdata->sensi_meth == AMICI_SENSI_FSA &&
            rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
            solver->AMISetStopTime(rdata->ts[it]);
        }
        if (rdata->ts[it] > model->t0()) {
            while (t < rdata->ts[it]) {
                if (model->nx > 0) {
                    if (std::isinf(rdata->ts[it])) {
                        SteadystateProblem sstate = SteadystateProblem(&t,&x,&sx);
                        sstate.workSteadyStateProblem(rdata, solver, model, it);
                    } else {
                        int status;
                        if (rdata->sensi_meth == AMICI_SENSI_ASA &&
                            rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
                            status = solver->AMISolveF(RCONST(rdata->ts[it]), &x, &dx,
                                                       &(t), AMICI_NORMAL, &ncheck);
                            
                        } else {
                            status = solver->AMISolve(RCONST(rdata->ts[it]), &x, &dx,
                                                      &(t), AMICI_NORMAL);
                        }
                        if (status == AMICI_ILL_INPUT) {
                            /* clustering of roots => turn off rootfinding */
                            solver->turnOffRootFinding();
                        }
                        if (status == AMICI_ROOT_RETURN) {
                            handleEvent(&tlastroot,FALSE);
                        }
                    }
                } else {
                    t = rdata->ts[it];
                }
            }
        }
        handleDataPoint(it);
    }

    /* fill events */
    if (model->ne > 0) {
        getEventOutput();
    }

    // set likelihood
    if (!edata) {
        rdata->invalidateLLH();
    }

    storeJacobianAndDerivativeInReturnData();
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

void ForwardProblem::handleEvent(realtype *tlastroot, const bool seflag) {
    /**
     * handleEvent executes everything necessary for the handling of events
     *
     * @param[out] tlastroot pointer to the timepoint of the last event 
     */

    int ie;
    int secondevent = 0;

    /* store heaviside information at event occurence */
    model->froot(t, &x, &dx, rootvals.data());
    
    if (!seflag) {
        solver->AMIGetRootInfo(rootsfound.data());
    }

    if (iroot < rdata->nmaxevent * model->ne) {
        for (ie = 0; ie < model->ne; ie++) {
            rootidx[iroot * model->ne + ie] =
                rootsfound.at(ie);
        }
    }
    for (ie = 0; ie < model->ne; ie++) {
        rvaltmp.at(ie) = rootvals.at(ie);
    }

    if (!seflag) {
        /* only extract in the first event fired */
        if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST &&
            rdata->sensi_meth == AMICI_SENSI_FSA) {
            solver->AMIGetSens(&(t), &sx);
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

    getEventOutput();

    /* if we need to do forward sensitivities later on we need to store the old
     * x and the old xdot */
    if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        /* store x and xdot to compute jump in sensitivities */
        x_old = x;
        if (rdata->sensi_meth == AMICI_SENSI_FSA) {
            model->fxdot(t, &x, &dx, &xdot);
            xdot_old = xdot;
            dx_old = dx;

            /* compute event-time derivative only for primary events, we get
             * into trouble with multiple simultaneously firing events here (but
             * is this really well defined then?), in that case just use the
             * last ie and hope for the best. */
            if (!seflag) {
                for (ie = 0; ie < model->ne; ie++) {
                    if (rootsfound.at(ie) ==
                        1) { /* only consider transitions false -> true */
                        model->fstau(t, ie, &x, &sx);
                    }
                }
            }
        } else if (rdata->sensi_meth == AMICI_SENSI_ASA) {
            /* store x to compute jump in discontinuity */
            if (iroot < rdata->nmaxevent * model->ne) {
                x_disc[iroot] = x;
                xdot_disc[iroot] = xdot;
                xdot_old_disc[iroot] = xdot_old;
            }
        }
    }

    model->updateHeaviside(rootsfound);

    applyEventBolus();

    if (iroot < rdata->nmaxevent * model->ne) {
        discs[iroot] = t;
        ++iroot;
    } else {
        warnMsgIdAndTxt("AMICI:mex:TOO_MUCH_EVENT",
                        "Event was recorded but not reported as the number of "
                        "occured events exceeded (nmaxevents)*(number of "
                        "events in model definition)!");
        solver->AMIReInit(t, &x, &dx); /* reinitialise so that we can continue in peace */
        return;
    }

    if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        if (rdata->sensi_meth == AMICI_SENSI_FSA) {

            /* compute the new xdot  */
            model->fxdot(t, &x, &dx, &xdot);
            applyEventSensiBolusFSA();
        }
    }

    /* check whether we need to fire a secondary event */
    model->froot(t, &x, &dx, rootvals.data());
    for (ie = 0; ie < model->ne; ie++) {
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
        solver->AMIReInit(t, &x, &dx);

        /* make time derivative consistent */
        solver->AMICalcIC(t, &x, &dx);

        if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
            if (rdata->sensi_meth == AMICI_SENSI_FSA) {
                solver->AMISensReInit(solver->getInternalSensitivityMethod(), &sx, &sdx);
            }
        }
    }
}

void ForwardProblem::storeJacobianAndDerivativeInReturnData() {
    /**
     * evalues the Jacobian and differential equation right hand side, stores it
     * in rdata 
     */

    /* entries in rdata are actually (double) while entries in tdata are
       (realtype)
       we should perform proper casting here. */
    model->fxdot(t, &x, &dx, &xdot);
    rdata->xdot = xdot.getVector();

    model->fJ(t, 0.0, &x, &dx, &xdot, Jtmp);
    rdata->J.assign(Jtmp->data,Jtmp->data + model->nx * model->nx);
}

void ForwardProblem::getEventOutput() {
    /**
     * getEventOutput extracts output information for events
     *
     */

    if (t == model->gett(rdata->nt - 1,rdata)) {
        // call from fillEvent at last timepoint
        model->froot(t, &x, &dx, rootvals.data());
    }

    /* EVENT OUTPUT */
    for (int ie = 0; ie < model->ne;
         ie++) { /* only look for roots of the rootfunction not discontinuities
                    */
        if (nroots.at(ie) >= rdata->nmaxevent)
            continue;

        if (rootsfound.at(ie) == 1 || t == model->gett(rdata->nt - 1,rdata)) {
            /* only consider transitions false
             -> true  or event filling*/

            model->fz(nroots.at(ie), ie, t, &x, rdata);

            if (edata) {
                model->fsigma_z(t, ie, nroots.data(), edata, rdata);

                model->fJz(nroots.at(ie), rdata, edata);

                if (t == model->gett(rdata->nt - 1,rdata)) {
                    // call from fillEvent at last
                    // timepoint, add regularization
                    // based on rz
                    model->frz(nroots.at(ie), ie, t, &x, rdata);
                    model->fJrz(nroots.at(ie), rdata, edata);
                }
            }

            if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
                prepEventSensis(ie);
                if (rdata->sensi_meth == AMICI_SENSI_FSA) {
                    getEventSensisFSA(ie);
                }
            }
            nroots.at(ie)++;
        }
    }
    if (t == model->gett(rdata->nt - 1,rdata)) {
        // call from fillEvent at last timepoint
        // loop until all events are filled
        bool continue_loop = false;
        for (int ie = 0; ie < model->ne;
             ie++) {
            if (nroots.at(ie) < rdata->nmaxevent) {
                continue_loop = true;
                break;
            }
        }
        if(continue_loop)
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
     * @param[in] ie index of current event @type int
     */

    if (edata) {
        for (int iz = 0; iz < model->nztrue; iz++) {
            if (model->z2event[iz] - 1 == ie) {
                if (!isNaN(edata->mz[iz * rdata->nmaxevent + nroots.at(ie)])) {
                    model->fdzdp(t, ie, &x);
                    
                    model->fdzdx(t, ie, &x);
                    
                    if (t == model->gett(rdata->nt - 1,rdata)) {
                        model->fdrzdp(t, ie, &x);
                        model->fdrzdx(t, ie, &x);
                    }
                    /* extract the value for the standard deviation, if the data
                     value is NaN, use
                     the parameter value. Store this value in the return
                     struct */
                    if (isNaN(edata->sigmaz[nroots.at(ie) +
                                               rdata->nmaxevent * iz])) {
                        model->fdsigma_zdp(t);
                    } else {
                        for (int ip = 0; ip < model->nplist(); ip++) {
                            model->dsigmazdp[iz + model->nz * ip] = 0;
                        }
                        model->sigmaz[iz] =
                        edata->sigmaz[nroots.at(ie) +
                                      rdata->nmaxevent * iz];
                    }
                    rdata->sigmaz[nroots.at(ie)*model->nz + iz] =
                    model->sigmaz[iz];
                    for (int ip = 0; ip < model->nplist(); ip++) {
                        rdata->ssigmaz[(nroots.at(ie)*model->np()+ip)*
                                       model->nz + iz] =
                        model->dsigmazdp[iz + model->nz * ip];
                    }
                }
            }
        }
        model->fdJzdz(nroots.at(ie), rdata, edata);
        model->fdJzdsigma(nroots.at(ie), rdata, edata);

        if (t == rdata->ts[rdata->nt - 1]) {
            model->fdJrzdz(nroots.at(ie), rdata, edata);
            model->fdJrzdsigma(nroots.at(ie), rdata, edata);
        }
        model->fdJzdx(&dJzdx, nroots.at(ie), t, edata, rdata);
        model->fdJzdp(nroots.at(ie), t, edata, rdata);
        if (rdata->sensi_meth == AMICI_SENSI_ASA) {
            for (int iJ = 0; iJ < model->nJ; iJ++) {
                for (int ip = 0; ip < model->nplist(); ip++) {
                    if (iJ == 0) {
                        if (model->nz > 0) {
                            rdata->sllh[ip] -= model->dJzdp[ip];
                        }
                    } else {
                        if (model->nz > 0) {
                            rdata->s2llh[(iJ - 1) + ip * (model->nJ - 1)] -=
                                model->dJzdp[iJ + ip * model->nJ];
                        }
                    }
                }
            }
        }
    }
}

void ForwardProblem::getEventSensisFSA(int ie) {
    /**
     * getEventSensisFSA extracts event information for forward sensitivity
     * analysis
     *
     * @param[in] ie index of event type @type int
     */
    if (t == rdata->ts[rdata->nt - 1]) {
        // call from fillEvent at last timepoint
        model->fsz_tf(nroots.at(ie), rdata);
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
     * @param[in] it index of data point @type int
     */

    for (int ix = 0; ix < model->nx; ix++) {
        rdata->x.at(it*model->nx + ix) = x[ix];
    }
    
    if (rdata->ts[it] > model->t0()) {
        solver->getDiagnosis(it, rdata);
    }
    
    getDataOutput(it);
}

void ForwardProblem::getDataOutput(int it) {
    /**
     * getDataOutput extracts output information for data-points
     *
     * @param[in] it index of current timepoint @type int
     */

    model->fy(it, rdata);
    model->fsigma_y(it, edata, rdata);
    model->fJy(it, rdata, edata);
    
    if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        prepDataSensis(it);
        if (rdata->sensi_meth == AMICI_SENSI_FSA) {
            getDataSensisFSA(it);
        }
    }
}

void ForwardProblem::prepDataSensis(int it) {
    /**
     * prepDataSensis preprocesses the provided experimental data to compute
     * sensitivities via adjoint or forward methods later on
     *
     * @param[in] it index of current timepoint @type int
     */

    model->fdydx(it, rdata);
    model->fdydp(it, rdata);

    if (!edata)
        return;

    model->fdsigma_ydp(it, rdata);

    for (int iy = 0; iy < model->nytrue; iy++) {
        if (!isNaN(edata->sigmay[iy * rdata->nt + it])) {
            for (int ip = 0; ip < model->nplist(); ip++) {
                model->dsigmaydp[ip * model->ny + iy] = 0;
            }
        }
        for (int ip = 0; ip < model->nplist(); ip++) {
            rdata->ssigmay[it + rdata->nt * (ip * model->ny + iy)] =
                model->dsigmaydp[ip * model->ny + iy];
        }
    }
    model->fdJydy(it, rdata, edata);
    model->fdJydsigma(it, rdata, edata);
    model->fdJydx(&dJydx, it, edata, rdata);
    model->fdJydp(it, edata, rdata);

    if (rdata->sensi_meth != AMICI_SENSI_ASA)
        return;

    for (int iJ = 0; iJ < model->nJ; iJ++) {
        for (int ip = 0; ip < model->nplist(); ip++) {
            if (iJ == 0) {
                if (model->ny > 0) {
                    rdata->sllh.at(ip) -= model->dJydp[ip * model->nJ];
                }
            } else {
                if (model->ny > 0) {
                    rdata->s2llh.at((iJ - 1) + ip * (model->nJ - 1)) -=
                        model->dJydp[iJ + ip * model->nJ];
                }
            }
        }
    }
}

void ForwardProblem::getDataSensisFSA(int it) {
    /**
     * getDataSensisFSA extracts data information for forward sensitivity
     * analysis
     *
     * @param[in] it index of current timepoint @type int
     */

    if (!(std::isinf(rdata->ts[it]))) {
        if (rdata->ts[it] > model->t0()) {
            solver->AMIGetSens(&(t), &sx);
        }
    }
    
    if (model->nx > 0) {
        for (int ix = 0; ix < model->nx; ix++) {
            for (int ip = 0; ip < model->nplist(); ip++) {
                rdata->sx[(it * model->nplist() + ip) * rdata->nx + ix] =
                sx.at(ix,ip);
            }
        }
    }

    for (int iy = 0; iy < model->nytrue; iy++) {
        if (edata) {
            if (isNaN(edata->sigmay[iy * rdata->nt + it])) {
                model->fdsigma_ydp(it, rdata);
            } else {
                for (int ip = 0; ip < model->nplist(); ip++) {
                    model->dsigmaydp[ip * model->ny + iy] = 0;
                }
            }
            for (int ip = 0; ip < model->nplist(); ip++) {
                rdata->ssigmay[(it * model->nplist() + ip) * model->ny + iy] =
                    model->dsigmaydp[ip * model->ny + iy];
            }
        } else {
            for (int ip = 0; ip < model->nplist(); ip++) {
                rdata->ssigmay[(it * model->nplist() + ip) * model->ny + iy] = 0;
            }
        }
    }
    
    if (model->ny > 0) {
        model->fsy(it, rdata);
        if (edata) {
            model->fsJy(it, dJydx, rdata);
        }
    }
}

void ForwardProblem::applyEventBolus() {
    /**
     * applyEventBolus applies the event bolus to the current state
     *
     * @param[in] model pointer to model specification object @type Model
     */

    for (int ie = 0; ie < model->ne; ie++) {
        if (rootsfound.at(ie) ==
            1) { /* only consider transitions false -> true */
            model->fdeltax(ie, t, &x, &xdot, &xdot_old);

            for (int ix = 0; ix < model->nx; ix++) {
                x[ix] += model->deltax[ix];
            }
        }
    }
}

void ForwardProblem::applyEventSensiBolusFSA() {
    /**
     * applyEventSensiBolusFSA applies the event bolus to the current
     * sensitivities
     *
     */
    for (int ie = 0; ie < model->ne; ie++) {
        if (rootsfound.at(ie) ==
            1) { /* only consider transitions false -> true */
            model->fdeltasx(ie, t, &x_old, &sx, &xdot, &xdot_old);

            for (int ip = 0; ip < model->nplist(); ip++) {
                for (int ix = 0; ix < model->nx; ix++) {
                    sx.at(ix,ip) += model->deltasx[ix + model->nx * ip];
                }
            }
        }
    }
}

} // namespace amici
