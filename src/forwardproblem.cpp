#include "include/forwardproblem.h"
#include "include/amici_model.h"
#include "include/amici_solver.h"
#include "include/edata.h"
#include "include/rdata.h"
#include "include/steadystateproblem.h"
#include "include/tdata.h"
#include "include/udata.h"
#include <cstring>

// Ensure AMICI options are in sync with Sundials options
static_assert(InternalSensitivityMethod::SIMULTANEOUS == CV_SIMULTANEOUS, "");
static_assert(InternalSensitivityMethod::STAGGERED == CV_STAGGERED, "");
static_assert(InternalSensitivityMethod::STAGGERED1 == CV_STAGGERED1, "");

static_assert(InterpolationType::HERMITE == CV_HERMITE, "");
static_assert(InterpolationType::POLYNOMIAL == CV_POLYNOMIAL, "");

static_assert(LinearMultistepMethod::ADAMS == CV_ADAMS, "");
static_assert(LinearMultistepMethod::BDF == CV_BDF, "");

static_assert(NonlinearSolverIteration::FUNCTIONAL == CV_FUNCTIONAL, "");
static_assert(NonlinearSolverIteration::NEWTON == CV_NEWTON, "");

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

int ForwardProblem::workForwardProblem(const UserData *udata, TempData *tdata,
                                       ReturnData *rdata, const ExpData *edata,
                                       Model *model) {
    /**
     * workForwardProblem solves the forward problem. if forward sensitivities
     * are enabled this will also compute sensitivies
     *
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[in] tdata pointer to the temporary data struct @type TempData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[out] edata pointer to the experimental data struct @type ExpData
     * @param[in] model pointer to model specification object @type Model
     * @return int status flag indicating success of execution @type int
     */

    Solver *solver = tdata->solver;

    int status = solver->setupAMI(udata, tdata, model);

    if (status != AMICI_SUCCESS)
        return status;

    int ncheck = 0; /* the number of (internal) checkpoints stored so far */
    realtype *x_tmp;
    realtype tlastroot = 0; /* storage for last found root */

    /* if preequilibration is necessary, start Newton solver */
    if (udata->newton_preeq == 1) {
        status = SteadystateProblem::workSteadyStateProblem(udata, tdata, rdata,
                                                            solver, model, -1);
        if (status != AMICI_SUCCESS)
            goto freturn;
    }

    /* loop over timepoints */
    for (int it = 0; it < rdata->nt; it++) {
        if (rdata->sensi_meth == AMICI_SENSI_FSA &&
            rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
            status = solver->AMISetStopTime(rdata->ts[it]);
        }
        if (status == AMICI_SUCCESS) {
            /* only integrate if no errors occured */
            if (rdata->ts[it] > udata->tstart) {
                while (tdata->t < rdata->ts[it]) {
                    if (rdata->sensi_meth == AMICI_SENSI_ASA &&
                        rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
                        if (model->nx > 0) {
                            status = solver->AMISolveF(
                                RCONST(rdata->ts[it]), tdata->x, tdata->dx,
                                &(tdata->t), AMICI_NORMAL, &ncheck);
                        } else {
                            tdata->t = rdata->ts[it];
                        }
                    } else {
                        if (model->nx > 0) {
                            if (std::isinf(rdata->ts[it])) {
                                status =
                                    SteadystateProblem::workSteadyStateProblem(
                                        udata, tdata, rdata, solver, model, it);
                            } else {
                                status = solver->AMISolve(
                                    RCONST(rdata->ts[it]), tdata->x, tdata->dx,
                                    &(tdata->t), AMICI_NORMAL);
                            }
                        } else {
                            tdata->t = rdata->ts[it];
                        }
                    }
                    if (model->nx > 0) {
                        x_tmp = NV_DATA_S(tdata->x);
                        if (!x_tmp)
                            return AMICI_ERROR_SIMULATION;
                        if (status == -22) {
                            /* clustering of roots => turn off rootfinding */
                            solver->turnOffRootFinding();
                            status = AMICI_SUCCESS;
                        }
                        if (status == AMICI_ROOT_RETURN) {
                            status =
                                handleEvent(&tlastroot, udata, rdata, edata,
                                            tdata, 0, solver, model);
                            if (status != AMICI_SUCCESS)
                                goto freturn;
                        }
                        /* integration error occured */
                        if (status != AMICI_SUCCESS)
                            goto freturn;
                    }
                }
            }
            status =
                handleDataPoint(it, udata, rdata, edata, tdata, solver, model);
            if (status != AMICI_SUCCESS)
                goto freturn;
        } else {
            for (int ix = 0; ix < model->nx; ix++)
                rdata->x[ix * rdata->nt + it] = amiGetNaN();
        }
    }

    /* fill events */
    if (model->ne > 0) {
        getEventOutput(udata, rdata, edata, tdata, model);
    }

    // set likelihood
    if (edata) {
        *rdata->llh = -tdata->Jy[0] - tdata->Jz[0];
    } else {
        *rdata->llh = amiGetNaN();
    }

freturn:
    storeJacobianAndDerivativeInReturnData(tdata, rdata, model);
    return status;
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

int ForwardProblem::handleEvent(realtype *tlastroot, const UserData *udata,
                                ReturnData *rdata, const ExpData *edata,
                                TempData *tdata, int seflag, Solver *solver,
                                Model *model) {
    /**
     * handleEvent executes everything necessary for the handling of events
     *
     * @param[out] tlastroot pointer to the timepoint of the last event @type
     * *realtype
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @param[in] seflag flag indicating whether this is a secondary event @type
     * int
     * @param[in] solver pointer to solver object @type Solver
     * @param[in] model pointer to model specification object @type Model
     * @return status flag indicating success of execution @type int
     */

    int ie;
    int secondevent = 0;
    int status = AMICI_SUCCESS;

    /* store heaviside information at event occurence */
    if (model->froot(tdata->t, tdata->x, tdata->dx, tdata->rootvals, tdata) !=
        AMICI_SUCCESS)
        return AMICI_ERROR_EVENT;

    if (seflag == 0) {
        status = solver->AMIGetRootInfo(tdata->rootsfound);
        if (status != AMICI_SUCCESS)
            return status;
    }

    if (tdata->iroot < rdata->nmaxevent * model->ne) {
        for (ie = 0; ie < model->ne; ie++) {
            tdata->rootidx[tdata->iroot * model->ne + ie] =
                tdata->rootsfound[ie];
        }
    }
    for (ie = 0; ie < model->ne; ie++) {
        tdata->rvaltmp[ie] = tdata->rootvals[ie];
    }

    if (seflag == 0) {
        /* only extract in the first event fired */
        if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST &&
            rdata->sensi_meth == AMICI_SENSI_FSA) {
            status = solver->AMIGetSens(&(tdata->t), tdata->sx);
            if (status != AMICI_SUCCESS)
                return AMICI_ERROR_SA;
        }

        /* only check this in the first event fired, otherwise this will always
         * be true */
        if (tdata->t == *tlastroot) {
            warnMsgIdAndTxt("AMICI:mex:STUCK_EVENT",
                            "AMICI is stuck in an event, as the initial "
                            "step-size after the event is too small. To fix "
                            "this, increase absolute and relative tolerances!");
            return AMICI_ERROR_EVENT;
        }
        *tlastroot = tdata->t;
    }

    status = getEventOutput(udata, rdata, edata, tdata, model);
    if (status != AMICI_SUCCESS)
        return status;

    /* if we need to do forward sensitivities later on we need to store the old
     * x and the old xdot */
    if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        /* store x and xdot to compute jump in sensitivities */
        N_VScale(1.0, tdata->x, tdata->x_old);
        if (rdata->sensi_meth == AMICI_SENSI_FSA) {
            status =
                model->fxdot(tdata->t, tdata->x, tdata->dx, tdata->xdot, tdata);
            N_VScale(1.0, tdata->xdot, tdata->xdot_old);
            N_VScale(1.0, tdata->dx, tdata->dx_old);

            /* compute event-time derivative only for primary events, we get
             * into trouble with multiple simultaneously firing events here (but
             * is this really well defined then?), in that case just use the
             * last ie and hope for the best. */
            if (seflag == 0) {
                for (ie = 0; ie < model->ne; ie++) {
                    if (tdata->rootsfound[ie] ==
                        1) { /* only consider transitions false -> true */
                        model->fstau(tdata->t, ie, tdata->x, tdata->sx, tdata);
                    }
                }
            }
        } else if (rdata->sensi_meth == AMICI_SENSI_ASA) {
            /* store x to compute jump in discontinuity */
            if (tdata->iroot < rdata->nmaxevent * model->ne) {
                N_VScale(1.0, tdata->x, tdata->x_disc[tdata->iroot]);
                N_VScale(1.0, tdata->xdot, tdata->xdot_disc[tdata->iroot]);
                N_VScale(1.0, tdata->xdot_old,
                         tdata->xdot_old_disc[tdata->iroot]);
            }
        }
    }

    status = updateHeaviside(tdata, model->ne);
    if (status != AMICI_SUCCESS)
        return status;

    status = applyEventBolus(tdata, model);
    if (status != AMICI_SUCCESS)
        return status;

    if (tdata->iroot < rdata->nmaxevent * model->ne) {
        tdata->discs[tdata->iroot] = tdata->t;
        ++tdata->iroot;
    } else {
        warnMsgIdAndTxt("AMICI:mex:TOO_MUCH_EVENT",
                        "Event was recorded but not reported as the number of "
                        "occured events exceeded (nmaxevents)*(number of "
                        "events in model definition)!");
        status = solver->AMIReInit(
            tdata->t, tdata->x,
            tdata->dx); /* reinitialise so that we can continue in peace */
        return status;
    }

    if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        if (rdata->sensi_meth == AMICI_SENSI_FSA) {

            /* compute the new xdot  */
            status =
                model->fxdot(tdata->t, tdata->x, tdata->dx, tdata->xdot, tdata);
            if (status != AMICI_SUCCESS)
                return status;

            status = applyEventSensiBolusFSA(tdata, model);
            if (status != AMICI_SUCCESS)
                return status;
        }
    }

    /* check whether we need to fire a secondary event */
    status =
        model->froot(tdata->t, tdata->x, tdata->dx, tdata->rootvals, tdata);
    if (status != AMICI_SUCCESS)
        return status;
    for (ie = 0; ie < model->ne; ie++) {
        /* the same event should not trigger itself */
        if (tdata->rootsfound[ie] == 0) {
            /* check whether there was a zero-crossing */
            if (0 > tdata->rvaltmp[ie] * tdata->rootvals[ie]) {
                if (tdata->rvaltmp[ie] < tdata->rootvals[ie]) {
                    tdata->rootsfound[ie] = 1;
                } else {
                    tdata->rootsfound[ie] = -1;
                }
                secondevent++;
            } else {
                tdata->rootsfound[ie] = 0;
            }
        } else {
            /* don't fire the same event again */
            tdata->rootsfound[ie] = 0;
        }
    }
    /* fire the secondary event */
    if (secondevent > 0) {
        status = handleEvent(tlastroot, udata, rdata, edata, tdata, secondevent,
                             solver, model);
        if (status != AMICI_SUCCESS)
            return status;
    }

    /* only reinitialise in the first event fired */
    if (seflag == 0) {
        status = solver->AMIReInit(tdata->t, tdata->x, tdata->dx);
        if (status != AMICI_SUCCESS)
            return status;

        /* make time derivative consistent */
        status = solver->AMICalcIC(tdata->t);
        if (status != AMICI_SUCCESS)
            return status;

        if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
            if (rdata->sensi_meth == AMICI_SENSI_FSA) {
                status =
                    solver->AMISensReInit(udata->ism, tdata->sx, tdata->sdx);
                if (status != AMICI_SUCCESS)
                    return status;
            }
        }
    }
    return status;
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

int ForwardProblem::storeJacobianAndDerivativeInReturnData(TempData *tdata,
                                                           ReturnData *rdata,
                                                           Model *model) {
    /**
     * evalues the Jacobian and differential equation right hand side, stores it
     * in tdata and
     * and copies it to rdata
     *
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] model pointer to model specification object @type Model
     * @return void
     */

    if (!tdata || model->nx <= 0)
        return AMICI_SUCCESS;

    /* entries in rdata are actually (double) while entries in tdata are
       (realtype)
       we should perform proper casting here. */
    int status =
        model->fxdot(tdata->t, tdata->x, tdata->dx, tdata->xdot, tdata);
    if (status != AMICI_SUCCESS)
        return status;

    realtype *xdot_tmp = NV_DATA_S(tdata->xdot);
    if (!xdot_tmp)
        return AMICI_ERROR_SIMULATION;

    if (rdata->xdot)
        memcpy(rdata->xdot, xdot_tmp, model->nx * sizeof(realtype));

    status = model->fJ(model->nx, tdata->t, 0, tdata->x, tdata->dx, tdata->xdot,
                       tdata->Jtmp, tdata, NULL, NULL, NULL);

    if (status != AMICI_SUCCESS)
        return status;

    if (rdata->J)
        memcpy(rdata->J, tdata->Jtmp->data,
               model->nx * model->nx * sizeof(realtype));

    return status;
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

int ForwardProblem::getEventOutput(const UserData *udata, ReturnData *rdata,
                                   const ExpData *edata, TempData *tdata,
                                   Model *model) {
    /**
     * getEventOutput extracts output information for events
     *
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @param[in] model pointer to model specification object @type Model
     * @return status flag indicating success of execution @type int
     */

    int status = AMICI_SUCCESS;

    if (tdata->t ==
        rdata->ts[rdata->nt - 1]) { // call from fillEvent at last timepoint
        status =
            model->froot(tdata->t, tdata->x, tdata->dx, tdata->rootvals, tdata);
        if (status != AMICI_SUCCESS)
            return status;
    }

    /* EVENT OUTPUT */
    for (int ie = 0; ie < model->ne;
         ie++) { /* only look for roots of the rootfunction not discontinuities
                    */
        if (tdata->nroots[ie] >= rdata->nmaxevent)
            continue;

        if (tdata->rootsfound[ie] == 1 ||
            tdata->t ==
                rdata->ts[rdata->nt - 1]) { /* only consider transitions false
                                               -> true  or event filling*/
            status = model->fz(tdata->t, ie, tdata->x, tdata, rdata);
            if (status != AMICI_SUCCESS)
                return status;

            if (edata) {
                status = model->fsigma_z(tdata->t, ie, tdata);
                if (status != AMICI_SUCCESS)
                    return status;
                for (int iz = 0; iz < model->nztrue; iz++) {
                    if (model->z2event[iz] - 1 == ie) {

                        if (!amiIsNaN(edata->sigmaz[tdata->nroots[ie] +
                                                    rdata->nmaxevent * iz])) {
                            tdata->sigmaz[iz] =
                                edata->sigmaz[tdata->nroots[ie] +
                                              rdata->nmaxevent * iz];
                        }
                        rdata->sigmaz[tdata->nroots[ie] +
                                      rdata->nmaxevent * iz] =
                            tdata->sigmaz[iz];
                    }
                }

                status =
                    model->fJz(tdata->t, ie, tdata->x, tdata, edata, rdata);
                if (status != AMICI_SUCCESS)
                    return status;

                if (tdata->t ==
                    rdata->ts[rdata->nt - 1]) { // call from fillEvent at last
                                                // timepoint, add regularization
                                                // based on rz
                    status = model->frz(tdata->t, ie, tdata->x, tdata, rdata);
                    if (status != AMICI_SUCCESS)
                        return status;

                    status = model->fJrz(tdata->t, ie, tdata->x, tdata, edata,
                                         rdata);
                    if (status != AMICI_SUCCESS)
                        return status;
                }
            }

            if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
                status = prepEventSensis(ie, rdata, edata, tdata, model);
                if (status != AMICI_SUCCESS)
                    return status;
                if (rdata->sensi_meth == AMICI_SENSI_FSA) {
                    status = getEventSensisFSA(ie, rdata, edata, tdata, model);
                    if (status != AMICI_SUCCESS)
                        return status;
                }
            }
            tdata->nroots[ie]++;
        }
    }
    return status;
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

int ForwardProblem::prepEventSensis(int ie, ReturnData *rdata,
                                    const ExpData *edata, TempData *tdata,
                                    Model *model) {
    /**
     * prepEventSensis preprocesses the provided experimental data to compute
     * event sensitivities via adjoint or forward methods later on
     *
     * @param[in] ie index of current event @type int
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @param[in] model pointer to model specification object @type Model
     * @return status flag indicating success of execution @type int
     */

    int status = AMICI_SUCCESS;
    if (edata) {
        for (int iz = 0; iz < model->nztrue; iz++) {
            if (model->z2event[iz] - 1 == ie) {
                if (!amiIsNaN(
                        edata->mz[iz * rdata->nmaxevent + tdata->nroots[ie]])) {
                    status = model->fdzdp(tdata->t, ie, tdata->x, tdata);
                    if (status != AMICI_SUCCESS)
                        return status;

                    status = model->fdzdx(tdata->t, ie, tdata->x, tdata);
                    if (status != AMICI_SUCCESS)
                        return status;

                    if (tdata->t == rdata->ts[rdata->nt - 1]) {
                        status = model->fdrzdp(tdata->t, ie, tdata->x, tdata);
                        if (status != AMICI_SUCCESS)
                            return status;
                        status = model->fdrzdx(tdata->t, ie, tdata->x, tdata);
                        if (status != AMICI_SUCCESS)
                            return status;
                    }
                    /* extract the value for the standard deviation, if the data
                       value is NaN, use
                         the parameter value. Store this value in the return
                       struct */
                    if (amiIsNaN(edata->sigmaz[tdata->nroots[ie] +
                                               rdata->nmaxevent * iz])) {
                        status = model->fdsigma_zdp(tdata->t, ie, tdata);
                        if (status != AMICI_SUCCESS)
                            return status;
                    } else {
                        for (int ip = 0; ip < rdata->nplist; ip++) {
                            tdata->dsigmazdp[iz + model->nz * ip] = 0;
                        }
                        tdata->sigmaz[iz] =
                            edata->sigmaz[tdata->nroots[ie] +
                                          rdata->nmaxevent * iz];
                    }
                    rdata->sigmaz[tdata->nroots[ie] + rdata->nmaxevent * iz] =
                        tdata->sigmaz[iz];
                    for (int ip = 0; ip < rdata->nplist; ip++) {
                        rdata->ssigmaz[tdata->nroots[ie] +
                                       rdata->nmaxevent *
                                           (iz + model->nz * ip)] =
                            tdata->dsigmazdp[iz + model->nz * ip];
                    }
                }
            }
        }
        status = model->fdJzdz(tdata->t, ie, tdata->x, tdata, edata, rdata);
        if (status != AMICI_SUCCESS)
            return status;

        status = model->fdJzdsigma(tdata->t, ie, tdata->x, tdata, edata, rdata);
        if (status != AMICI_SUCCESS)
            return status;

        if (tdata->t == rdata->ts[rdata->nt - 1]) {
            status =
                model->fdJrzdz(tdata->t, ie, tdata->x, tdata, edata, rdata);
            if (status != AMICI_SUCCESS)
                return status;

            status =
                model->fdJrzdsigma(tdata->t, ie, tdata->x, tdata, edata, rdata);
            if (status != AMICI_SUCCESS)
                return status;
        }
        status = model->fdJzdx(ie, tdata, edata);
        if (status != AMICI_SUCCESS)
            return status;
        status = model->fdJzdp(ie, tdata, edata, rdata);
        if (status != AMICI_SUCCESS)
            return status;
        if (rdata->sensi_meth == AMICI_SENSI_ASA) {
            for (int iJ = 0; iJ < model->nJ; iJ++) {
                for (int ip = 0; ip < rdata->nplist; ip++) {
                    if (iJ == 0) {
                        if (model->nz > 0) {
                            rdata->sllh[ip] -= tdata->dJzdp[ip];
                        }
                    } else {
                        if (model->nz > 0) {
                            rdata->s2llh[(iJ - 1) + ip * (model->nJ - 1)] -=
                                tdata->dJzdp[iJ + ip * model->nJ];
                        }
                    }
                }
            }
        }
    }
    return status;
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

int ForwardProblem::getEventSensisFSA(int ie, ReturnData *rdata,
                                      const ExpData *edata, TempData *tdata,
                                      Model *model) {
    /**
     * getEventSensisFSA extracts event information for forward sensitivity
     * analysis
     *
     * @param[in] ie index of event type @type int
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[in] tdata pointer to the temporary data struct @type TempData
     * @param[in] model pointer to model specification object @type Model
     * @return status flag indicating success of execution @type int
     */

    int status = AMICI_SUCCESS;

    if (tdata->t ==
        rdata->ts[rdata->nt - 1]) { // call from fillEvent at last timepoint
        status = model->fsz_tf(ie, tdata, rdata);
        if (status != AMICI_SUCCESS)
            return status;

        status = model->fsrz(tdata->t, ie, tdata->x, tdata->sx, tdata, rdata);
        if (status != AMICI_SUCCESS)
            return status;
    } else {
        status = model->fsz(tdata->t, ie, tdata->x, tdata->sx, tdata, rdata);
        if (status != AMICI_SUCCESS)
            return status;
    }

    if (edata) {
        status = model->fsJz(ie, tdata, rdata);
        if (status != AMICI_SUCCESS)
            return status;
    }
    return AMICI_SUCCESS;
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

int ForwardProblem::handleDataPoint(int it, const UserData *udata,
                                    ReturnData *rdata, const ExpData *edata,
                                    TempData *tdata, Solver *solver,
                                    Model *model) {
    /**
     * handleDataPoint executes everything necessary for the handling of data
     * points
     *
     * @param[in] it index of data point @type int
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @param[in] solver pointer to solver object @type Solver
     * @param[in] model pointer to model specification object @type Model
     * @return status flag indicating success of execution @type int
     */

    if (model->nx > 0) {
        realtype *x_tmp = NV_DATA_S(tdata->x);
        if (!x_tmp)
            return AMICI_ERROR_DATA;
        for (int ix = 0; ix < model->nx; ix++) {
            rdata->x[it + rdata->nt * ix] = x_tmp[ix];
        }

        if (rdata->ts[it] > udata->tstart) {
            int status = solver->getDiagnosis(it, rdata);
            if (status != AMICI_SUCCESS)
                return status;
        }
    }

    return getDataOutput(it, udata, rdata, edata, tdata, solver, model);
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

int ForwardProblem::getDataOutput(int it, const UserData *udata,
                                  ReturnData *rdata, const ExpData *edata,
                                  TempData *tdata, Solver *solver,
                                  Model *model) {
    /**
     * getDataOutput extracts output information for data-points
     *
     * @param[in] it index of current timepoint @type int
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @param[in] solver pointer to solver object @type Solver
     * @param[in] model pointer to model specification object @type Model
     * @return status flag indicating success of execution @type int
     */

    int status = model->fy(rdata->ts[it], it, tdata->x, tdata, rdata);
    if (status != AMICI_SUCCESS)
        return status;

    if (edata) {
        status = model->fsigma_y(tdata->t, tdata);
        if (status != AMICI_SUCCESS)
            return status;
        for (int iy = 0; iy < model->nytrue; iy++) {
            /* extract the value for the standard deviation, if the data value
               is NaN, use
                 the parameter value. Store this value in the return struct */
            if (!amiIsNaN(edata->sigmay[iy * rdata->nt + it])) {
                tdata->sigmay[iy] = edata->sigmay[iy * rdata->nt + it];
            }
            rdata->sigmay[iy * rdata->nt + it] = tdata->sigmay[iy];
        }
        status = model->fJy(rdata->ts[it], it, tdata->x, tdata, edata, rdata);
        if (status != AMICI_SUCCESS)
            return status;
    } else {
        status = model->fsigma_y(tdata->t, tdata);
        if (status != AMICI_SUCCESS)
            return status;
        for (int iy = 0; iy < model->nytrue; iy++) {
            rdata->sigmay[iy * rdata->nt + it] = tdata->sigmay[iy];
        }
    }
    if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        status = prepDataSensis(it, rdata, edata, tdata, model);
        if (status != AMICI_SUCCESS)
            return status;
        if (rdata->sensi_meth == AMICI_SENSI_FSA) {
            status =
                getDataSensisFSA(it, udata, rdata, edata, tdata, solver, model);
            if (status != AMICI_SUCCESS)
                return status;
        }
    }
    return AMICI_SUCCESS;
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

int ForwardProblem::prepDataSensis(int it, ReturnData *rdata,
                                   const ExpData *edata, TempData *tdata,
                                   Model *model) {
    /**
     * prepDataSensis preprocesses the provided experimental data to compute
     * sensitivities via adjoint or forward methods later on
     *
     * @param[in] it index of current timepoint @type int
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @param[in] model pointer to model specification object @type Model
     * @return status flag indicating success of execution @type int
     */

    int status = model->fdydx(rdata->ts[it], it, tdata->x, tdata);
    if (status != AMICI_SUCCESS)
        return status;

    status = model->fdydp(rdata->ts[it], it, tdata->x, tdata);
    if (status != AMICI_SUCCESS)
        return status;

    if (!edata)
        return status;

    status = model->fdsigma_ydp(tdata->t, tdata);
    if (status != AMICI_SUCCESS)
        return status;

    for (int iy = 0; iy < model->nytrue; iy++) {
        if (!amiIsNaN(edata->sigmay[iy * rdata->nt + it])) {
            for (int ip = 0; ip < rdata->nplist; ip++) {
                tdata->dsigmaydp[ip * model->ny + iy] = 0;
            }
        }
        for (int ip = 0; ip < rdata->nplist; ip++) {
            rdata->ssigmay[it + rdata->nt * (ip * model->ny + iy)] =
                tdata->dsigmaydp[ip * model->ny + iy];
        }
    }
    status = model->fdJydy(tdata->t, it, tdata->x, tdata, edata, rdata);
    if (status != AMICI_SUCCESS)
        return status;

    status = model->fdJydsigma(tdata->t, it, tdata->x, tdata, edata, rdata);
    if (status != AMICI_SUCCESS)
        return status;

    status = model->fdJydx(it, tdata, edata);
    if (status != AMICI_SUCCESS)
        return status;

    status = model->fdJydp(it, tdata, edata, rdata);
    if (status != AMICI_SUCCESS)
        return status;

    if (rdata->sensi_meth != AMICI_SENSI_ASA)
        return status;

    for (int iJ = 0; iJ < model->nJ; iJ++) {
        for (int ip = 0; ip < rdata->nplist; ip++) {
            if (iJ == 0) {
                if (model->ny > 0) {
                    rdata->sllh[ip] -= tdata->dJydp[ip * model->nJ];
                }
            } else {
                if (model->ny > 0) {
                    rdata->s2llh[(iJ - 1) + ip * (model->nJ - 1)] -=
                        tdata->dJydp[iJ + ip * model->nJ];
                }
            }
        }
    }

    return status;
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

int ForwardProblem::getDataSensisFSA(int it, const UserData *udata,
                                     ReturnData *rdata, const ExpData *edata,
                                     TempData *tdata, Solver *solver,
                                     Model *model) {
    /**
     * getDataSensisFSA extracts data information for forward sensitivity
     * analysis
     *
     * @param[in] it index of current timepoint @type int
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @param[in] solver pointer to solver object @type Solver
     * @param[in] model pointer to model specification object @type Model
     * @return status flag indicating success of execution @type int
     */

    int status = AMICI_SUCCESS;
    realtype *sx_tmp;

    if (!(std::isinf(rdata->ts[it]))) {
        for (int ip = 0; ip < rdata->nplist; ip++) {
            if (model->nx > 0) {
                if (rdata->ts[it] > udata->tstart) {
                    status = solver->AMIGetSens(&(tdata->t), tdata->sx);
                    if (status != AMICI_SUCCESS)
                        return status;
                }

                sx_tmp = NV_DATA_S(tdata->sx[ip]);
                if (!sx_tmp)
                    return AMICI_ERROR_FSA;
                for (int ix = 0; ix < model->nx; ix++) {
                    rdata->sx[(ip * model->nx + ix) * rdata->nt + it] =
                        sx_tmp[ix];
                }
            }
        }
    }

    for (int iy = 0; iy < model->nytrue; iy++) {
        if (edata) {
            if (amiIsNaN(edata->sigmay[iy * rdata->nt + it])) {
                status = model->fdsigma_ydp(tdata->t, tdata);
                if (status != AMICI_SUCCESS)
                    return status;
            } else {
                for (int ip = 0; ip < rdata->nplist; ip++) {
                    tdata->dsigmaydp[ip * model->ny + iy] = 0;
                }
            }
            for (int ip = 0; ip < rdata->nplist; ip++) {
                rdata->ssigmay[it + rdata->nt * (ip * model->ny + iy)] =
                    tdata->dsigmaydp[ip * model->ny + iy];
            }
        } else {
            for (int ip = 0; ip < rdata->nplist; ip++) {
                rdata->ssigmay[it + rdata->nt * (ip * model->ny + iy)] = 0;
            }
        }
    }
    status = model->fsy(it, tdata, rdata);
    if (status != AMICI_SUCCESS)
        return status;
    if (edata) {
        status = model->fsJy(it, tdata, rdata);
        if (status != AMICI_SUCCESS)
            return status;
    }
    return status;
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

int ForwardProblem::applyEventBolus(TempData *tdata, Model *model) {
    /**
     * applyEventBolus applies the event bolus to the current state
     *
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @param[in] model pointer to model specification object @type Model
     * @return status flag indicating success of execution @type int
     */

    int ix, ie;
    int status = AMICI_SUCCESS;
    realtype *x_tmp;

    for (ie = 0; ie < model->ne; ie++) {
        if (tdata->rootsfound[ie] ==
            1) { /* only consider transitions false -> true */
            status = model->fdeltax(tdata->t, ie, tdata->x, tdata->xdot,
                                    tdata->xdot_old, tdata);
            if (status != AMICI_SUCCESS)
                return status;

            x_tmp = NV_DATA_S(tdata->x);
            if (!x_tmp)
                return AMICI_ERROR_EVENT;
            for (ix = 0; ix < model->nx; ix++) {
                x_tmp[ix] += tdata->deltax[ix];
            }
        }
    }
    return status;
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

int ForwardProblem::applyEventSensiBolusFSA(TempData *tdata, Model *model) {
    /**
     * applyEventSensiBolusFSA applies the event bolus to the current
     * sensitivities
     *
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @param[in] model pointer to model specification object @type Model
     * @return status flag indicating success of execution @type int
     */

    int ix, ip, ie;
    int status = AMICI_SUCCESS;
    realtype *sx_tmp;

    for (ie = 0; ie < model->ne; ie++) {
        if (tdata->rootsfound[ie] ==
            1) { /* only consider transitions false -> true */
            status = model->fdeltasx(tdata->t, ie, tdata->x_old, tdata->xdot,
                                     tdata->xdot_old, tdata->sx, tdata);
            if (status != AMICI_SUCCESS)
                return status;

            for (ip = 0; ip < tdata->rdata->nplist; ip++) {
                sx_tmp = NV_DATA_S(tdata->sx[ip]);
                if (!sx_tmp)
                    return AMICI_ERROR_FSA;
                for (ix = 0; ix < model->nx; ix++) {
                    sx_tmp[ix] += tdata->deltasx[ix + model->nx * ip];
                }
            }
        }
    }
    return status;
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

int ForwardProblem::updateHeaviside(TempData *tdata, const int ne) {
    /**
     * updateHeaviside updates the heaviside variables h on event occurences
     *
     * @param[in] ne number of events
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status = status flag indicating success of execution @type int;
     */

    /* tdata->rootsfound provides the direction of the zero-crossing, so adding
       it will give
         the right update to the heaviside variables */

    for (int ie = 0; ie < ne; ie++) {
        tdata->h[ie] += tdata->rootsfound[ie];
    }
    return AMICI_SUCCESS;
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

ForwardProblem::ForwardProblem() {}
