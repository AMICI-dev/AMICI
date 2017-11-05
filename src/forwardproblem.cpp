#include "include/forwardproblem.h"
#include "include/amici_model.h"
#include "include/amici_solver.h"
#include "include/amici_exception.h"
#include "include/edata.h"
#include "include/rdata.h"
#include "include/steadystateproblem.h"
#include "include/udata.h"
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
    
    
    ForwardProblem::ForwardProblem(const UserData *udata,
                   ReturnData *rdata, const ExpData *edata,
                   Model *model, Solver *solver) :
    dJzdx(model->nJ * model->nx * udata->nme(), 0.0),
    dJydx(model->nJ * model->nx * udata->nt(), 0.0),
    x(model->nx), x_old(model->nx), dx(model->nx), dx_old(model->nx),
    xdot(model->nx), xdot_old(model->nx),
    sx(model->nx,udata->nplist()), sdx(model->nx,udata->nplist())
    {
        t = udata->tstart();
        this->model = model;
        this->solver = solver;
        this->udata = udata;
        this->edata = edata;
        this->rdata = rdata;
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
        solver->setupAMI(udata, model);
    } catch (std::exception& ex) {
        throw AmiException("AMICI setup failed:\n(%s)",ex.what());
    } catch (...) {
        throw AmiException("AMICI setup failed due to an unknown error");
    }
    int ncheck = 0; /* the number of (internal) checkpoints stored so far */
    realtype tlastroot = 0; /* storage for last found root */

    /* if preequilibration is necessary, start Newton solver */
    if (udata->newton_preeq == 1) {
        SteadystateProblem sstate = SteadystateProblem(model->nx);
        sstate.workSteadyStateProblem(udata, rdata,
                                       solver, model, -1);
    }

    /* loop over timepoints */
    for (int it = 0; it < rdata->nt; it++) {
        if (rdata->sensi_meth == AMICI_SENSI_FSA &&
            rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
            solver->AMISetStopTime(rdata->ts[it]);
        }
        if (rdata->ts[it] > udata->tstart) {
            while (t < rdata->ts[it]) {
                if (model->nx > 0) {
                    if (std::isinf(rdata->ts[it])) {
                        SteadystateProblem sstate = SteadystateProblem(model->nx);
                        sstate.workSteadyStateProblem(udata,
                                                      rdata, solver, model, it);
                    } else {
                        int status;
                        if (rdata->sensi_meth == AMICI_SENSI_ASA &&
                            rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
                            status = solver->AMISolveF(RCONST(rdata->ts[it]), x, dx,
                                                       &(t), AMICI_NORMAL, &ncheck);
                            
                        } else {
                            status = solver->AMISolve(RCONST(rdata->ts[it]), x, dx,
                                                      &(t), AMICI_NORMAL);
                        }
                        if (status == AMICI_ILL_INPUT) {
                            /* clustering of roots => turn off rootfinding */
                            solver->turnOffRootFinding();
                        }
                        if (status == AMICI_ROOT_RETURN) {
                            handleEvent(&tlastroot, udata, rdata, edata,
                                        0, solver, model);
                        }
                    }
                } else {
                    t = rdata->ts[it];
                }
            }
        }
        handleDataPoint(it, udata, rdata, edata, solver, model);
    }

    /* fill events */
    if (model->ne > 0) {
        getEventOutput(udata, rdata, edata, model);
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

void ForwardProblem::handleEvent(realtype *tlastroot) {
    /**
     * handleEvent executes everything necessary for the handling of events
     *
     * @param[out] tlastroot pointer to the timepoint of the last event 
     */

    int ie;
    int secondevent = 0;

    /* store heaviside information at event occurence */
    model->froot(t, x, dx, rootvals, udata);
    
    if (seflag == 0) {
        solver->AMIGetRootInfo(rootsfound);
    }

    if (iroot < rdata->nmaxevent * model->ne) {
        for (ie = 0; ie < model->ne; ie++) {
            rootidx[iroot * model->ne + ie] =
                rootsfound[ie];
        }
    }
    for (ie = 0; ie < model->ne; ie++) {
        rvaltmp[ie] = rootvals[ie];
    }

    if (seflag == 0) {
        /* only extract in the first event fired */
        if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST &&
            rdata->sensi_meth == AMICI_SENSI_FSA) {
            solver->AMIGetSens(&(t), sx);
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

    getEventOutput(udata, rdata, edata, model);

    /* if we need to do forward sensitivities later on we need to store the old
     * x and the old xdot */
    if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        /* store x and xdot to compute jump in sensitivities */
        N_VScale(1.0, x, x_old);
        if (rdata->sensi_meth == AMICI_SENSI_FSA) {
            model->fxdot(t, x, dx, xdot, udata);
            N_VScale(1.0, xdot, xdot_old);
            N_VScale(1.0, dx, dx_old);

            /* compute event-time derivative only for primary events, we get
             * into trouble with multiple simultaneously firing events here (but
             * is this really well defined then?), in that case just use the
             * last ie and hope for the best. */
            if (seflag == 0) {
                for (ie = 0; ie < model->ne; ie++) {
                    if (rootsfound[ie] ==
                        1) { /* only consider transitions false -> true */
                        model->fstau(t, ie, x, sx, udata);
                    }
                }
            }
        } else if (rdata->sensi_meth == AMICI_SENSI_ASA) {
            /* store x to compute jump in discontinuity */
            if (iroot < rdata->nmaxevent * model->ne) {
                N_VScale(1.0, x, x_disc[iroot]);
                N_VScale(1.0, xdot, xdot_disc[iroot]);
                N_VScale(1.0, xdot_old,
                         xdot_old_disc[iroot]);
            }
        }
    }

    updateHeaviside();

    applyEventBolus();

    if (iroot < rdata->nmaxevent * model->ne) {
        discs[iroot] = t;
        ++iroot;
    } else {
        warnMsgIdAndTxt("AMICI:mex:TOO_MUCH_EVENT",
                        "Event was recorded but not reported as the number of "
                        "occured events exceeded (nmaxevents)*(number of "
                        "events in model definition)!");
        solver->AMIReInit(
            t, x,
            dx); /* reinitialise so that we can continue in peace */
        return;
    }

    if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        if (rdata->sensi_meth == AMICI_SENSI_FSA) {

            /* compute the new xdot  */
            model->fxdot(t, x, dx, xdot, udata);
            applyEventSensiBolusFSA();
        }
    }

    /* check whether we need to fire a secondary event */
    model->froot(t, x, dx, rootvals, udata);
    for (ie = 0; ie < model->ne; ie++) {
        /* the same event should not trigger itself */
        if (rootsfound[ie] == 0) {
            /* check whether there was a zero-crossing */
            if (0 > rvaltmp[ie] * rootvals[ie]) {
                if (rvaltmp[ie] < rootvals[ie]) {
                    rootsfound[ie] = 1;
                } else {
                    rootsfound[ie] = -1;
                }
                secondevent++;
            } else {
                rootsfound[ie] = 0;
            }
        } else {
            /* don't fire the same event again */
            rootsfound[ie] = 0;
        }
    }
    /* fire the secondary event */
    if (secondevent > 0) {
        handleEvent(tlastroot, udata, rdata, edata, secondevent,
                             solver, model);
    }

    /* only reinitialise in the first event fired */
    if (seflag == 0) {
        solver->AMIReInit(t, x, dx);

        /* make time derivative consistent */
        solver->AMICalcIC(t,udata);

        if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
            if (rdata->sensi_meth == AMICI_SENSI_FSA) {
                solver->AMISensReInit(udata->ism, sx, sdx);
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
    model->fxdot(t, x, dx, xdot, udata);
    memcpy(rdata->xdot, xdot.data(), model->nx * sizeof(realtype));

    model->fJ(model->nx, t, 0, x, dx, xdot,
                       Jtmp, NULL, NULL, NULL);
    memcpy(rdata->J, Jtmp->data,
           model->nx * model->nx * sizeof(realtype));

}

void ForwardProblem::getEventOutput() {
    /**
     * getEventOutput extracts output information for events
     *
     */

    if (t ==
        rdata->ts[rdata->nt - 1]) { // call from fillEvent at last timepoint
        model->froot(t, x, dx, rootvals, udata);
    }

    /* EVENT OUTPUT */
    for (int ie = 0; ie < model->ne;
         ie++) { /* only look for roots of the rootfunction not discontinuities
                    */
        if (nroots[ie] >= rdata->nmaxevent)
            continue;

        if (rootsfound[ie] == 1 || t == rdata->ts[rdata->nt - 1]) {
            /* only consider transitions false
             -> true  or event filling*/

            model->fz(t, ie, x, rdata);

            if (edata) {
                model->fsigma_z(t,ie,nroots,rdata,udata)

                model->fJz(t, ie, x, edata, rdata);

                if (t == rdata->ts[rdata->nt - 1]) {
                    // call from fillEvent at last
                    // timepoint, add regularization
                    // based on rz
                    model->frz(t, ie, x, rdata);
                    model->fJrz(t, ie, x, edata,
                                         rdata);
                }
            }

            if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
                prepEventSensis(ie);
                if (rdata->sensi_meth == AMICI_SENSI_FSA) {
                    getEventSensisFSA(ie);
                }
            }
            nroots[ie]++;
        }
    }
    if (t == rdata->ts[rdata->nt - 1]) {
        // call from fillEvent at last timepoint
        // loop until all events are filled
        bool continue_loop = false;
        for (int ie = 0; ie < model->ne;
             ie++) {
            if (nroots[ie] < rdata->nmaxevent) {
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
                if (!amiIsNaN(
                              edata->mz[iz * rdata->nmaxevent + nroots[ie]])) {
                    model->fdzdp(t, ie, x, udata);
                    
                    model->fdzdx(t, ie, x, udata);
                    
                    if (t == rdata->ts[rdata->nt - 1]) {
                        model->fdrzdp(t, ie, x, udata);
                        model->fdrzdx(t, ie, x, udata);
                    }
                    /* extract the value for the standard deviation, if the data
                     value is NaN, use
                     the parameter value. Store this value in the return
                     struct */
                    if (amiIsNaN(edata->sigmaz[nroots[ie] +
                                               rdata->nmaxevent * iz])) {
                        model->fdsigma_zdp(t, ie, udata);
                    } else {
                        for (int ip = 0; ip < rdata->nplist; ip++) {
                            model->dsigmazdp[iz + model->nz * ip] = 0;
                        }
                        model->sigmaz[iz] =
                        edata->sigmaz[nroots[ie] +
                                      rdata->nmaxevent * iz];
                    }
                    rdata->sigmaz[nroots[ie] + rdata->nmaxevent * iz] =
                    model->sigmaz[iz];
                    for (int ip = 0; ip < rdata->nplist; ip++) {
                        rdata->ssigmaz[nroots[ie] +
                                       rdata->nmaxevent *
                                       (iz + model->nz * ip)] =
                        model->dsigmazdp[iz + model->nz * ip];
                    }
                }
            }
        }
        model->fdJzdz(t, ie, x, edata, rdata);
        model->fdJzdsigma(t, ie, x, edata, rdata);

        if (t == rdata->ts[rdata->nt - 1]) {
            model->fdJrzdz(t, ie, x, edata, rdata);
            model->fdJrzdsigma(t, ie, x, edata, rdata);
        }
        model->fdJzdx(ie, edata);
        model->fdJzdp(ie, edata, rdata);
        if (rdata->sensi_meth == AMICI_SENSI_ASA) {
            for (int iJ = 0; iJ < model->nJ; iJ++) {
                for (int ip = 0; ip < rdata->nplist; ip++) {
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
    if (t ==
        rdata->ts[rdata->nt - 1]) { // call from fillEvent at last timepoint
        model->fsz_tf(ie, rdata);
        model->fsrz(t, ie, x, sx, rdata);
    } else {
        model->fsz(t, ie, x, sx, rdata);
    }

    if (edata) {
        model->fsJz(ie, rdata);
    }
}

void ForwardProblem::handleDataPoint(int it) {
    /**
     * handleDataPoint executes everything necessary for the handling of data
     * points
     *
     * @param[in] it index of data point @type int
     */

    if (model->nx > 0) {
        for (int ix = 0; ix < model->nx; ix++) {
            rdata->x[it + rdata->nt * ix] = x[ix];
        }

        if (rdata->ts[it] > udata->tstart) {
            solver->getDiagnosis(it, rdata);
        }
    }
    getDataOutput(it, udata, rdata, edata, solver, model);
}

void ForwardProblem::getDataOutput(int it) {
    /**
     * getDataOutput extracts output information for data-points
     *
     * @param[in] it index of current timepoint @type int
     */

    model->fy(it, x, rdata);
    model->fsigma_y(it, edata, rdata, udata);
    model->fJy(it, x, edata, rdata);
    
    if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        prepDataSensis(it, rdata, edata, model);
        if (rdata->sensi_meth == AMICI_SENSI_FSA) {
            getDataSensisFSA(it, udata, rdata, edata, solver, model);
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

    model->fdydx(rdata->ts[it], it, x, udata);

    model->fdydp(rdata->ts[it], it, x, udata);

    if (!edata)
        return;

    model->fdsigma_ydp(t, udata);

    for (int iy = 0; iy < model->nytrue; iy++) {
        if (!amiIsNaN(edata->sigmay[iy * rdata->nt + it])) {
            for (int ip = 0; ip < rdata->nplist; ip++) {
                model->dsigmaydp[ip * model->ny + iy] = 0;
            }
        }
        for (int ip = 0; ip < rdata->nplist; ip++) {
            rdata->ssigmay[it + rdata->nt * (ip * model->ny + iy)] =
                model->dsigmaydp[ip * model->ny + iy];
        }
    }
    model->fdJydy(t, it, x, edata, rdata);
    model->fdJydsigma(t, it, x, edata, rdata);
    model->fdJydx(it, edata);
    model->fdJydp(it, edata, rdata);

    if (rdata->sensi_meth != AMICI_SENSI_ASA)
        return;

    for (int iJ = 0; iJ < model->nJ; iJ++) {
        for (int ip = 0; ip < rdata->nplist; ip++) {
            if (iJ == 0) {
                if (model->ny > 0) {
                    rdata->sllh[ip] -= model->dJydp[ip * model->nJ];
                }
            } else {
                if (model->ny > 0) {
                    rdata->s2llh[(iJ - 1) + ip * (model->nJ - 1)] -=
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
        for (int ip = 0; ip < rdata->nplist; ip++) {
            if (model->nx > 0) {
                if (rdata->ts[it] > udata->tstart) {
                    solver->AMIGetSens(&(t), sx);
                }
                for (int ix = 0; ix < model->nx; ix++) {
                    rdata->sx[(ip * model->nx + ix) * rdata->nt + it] =
                        sx.at(ix,ip);
                }
            }
        }
    }

    for (int iy = 0; iy < model->nytrue; iy++) {
        if (edata) {
            if (amiIsNaN(edata->sigmay[iy * rdata->nt + it])) {
                model->fdsigma_ydp(t, udata);
            } else {
                for (int ip = 0; ip < rdata->nplist; ip++) {
                    model->dsigmaydp[ip * model->ny + iy] = 0;
                }
            }
            for (int ip = 0; ip < rdata->nplist; ip++) {
                rdata->ssigmay[it + rdata->nt * (ip * model->ny + iy)] =
                    model->dsigmaydp[ip * model->ny + iy];
            }
        } else {
            for (int ip = 0; ip < rdata->nplist; ip++) {
                rdata->ssigmay[it + rdata->nt * (ip * model->ny + iy)] = 0;
            }
        }
    }
    model->fsy(it, rdata);
    if (edata) {
        model->fsJy(it, rdata);
    }
}

void ForwardProblem::applyEventBolus() {
    /**
     * applyEventBolus applies the event bolus to the current state
     *
     * @param[in] model pointer to model specification object @type Model
     */

    for (int ie = 0; ie < model->ne; ie++) {
        if (rootsfound[ie] ==
            1) { /* only consider transitions false -> true */
            model->fdeltax(t, ie, x, xdot, xdot_old, udata);

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
        if (rootsfound[ie] ==
            1) { /* only consider transitions false -> true */
            model->fdeltasx(t, ie, x_old, xdot,
                                     xdot_old, sx, udata);

            for (int ip = 0; ip < udata->nplist(); ip++) {
                for (int ix = 0; ix < model->nx; ix++) {
                    sx.at(ix,ip) += model->deltasx[ix + model->nx * ip];
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

void ForwardProblem::updateHeaviside() {
    /**
     * updateHeaviside updates the heaviside variables h on event occurences
     *
     * @param[in] ne number of events
     */

    /* rootsfound provides the direction of the zero-crossing, so adding
       it will give
         the right update to the heaviside variables */

    for (int ie = 0; ie < model->ne(); ie++) {
        h[ie] += rootsfound[ie];
    }
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

ForwardProblem::ForwardProblem() {}

} // namespace amici
