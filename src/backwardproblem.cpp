#include "../include/backwardproblem.h"
#include "include/amici_model.h"
#include "include/amici_solver.h"
#include "include/edata.h"
#include "include/rdata.h"
#include "include/tdata.h"
#include "include/udata.h"
#include <cstring>

int BackwardProblem::workBackwardProblem(const UserData *udata, TempData *tdata,
                                         ReturnData *rdata, Model *model) {
    /**
     * workBackwardProblem solves the backward problem. if adjoint
     * sensitivities are enabled this will also compute sensitivies
     * workForwardProblem should be called before this is function is called
     *
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[in] tdata pointer to the temporary data struct @type TempData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] model pointer to model specification object @type Model
     * @return int status flag
     */
    int ix, it, ip;
    int status = (int)*rdata->status;
    double tnext;

    if (model->nx <= 0 || rdata->sensi < AMICI_SENSI_ORDER_FIRST ||
        rdata->sensi_meth != AMICI_SENSI_ASA || status != AMICI_SUCCESS) {
        return status;
    }

    Solver *solver = tdata->solver;
    solver->setupAMIB(udata, tdata, model);

    it = rdata->nt - 2;
    --tdata->iroot;
    while (it >= 0 || tdata->iroot >= 0) {

        /* check if next timepoint is a discontinuity or a data-point */
        tnext = getTnext(tdata->discs, tdata->iroot, rdata->ts, it, model);

        if (tnext < tdata->t) {
            status = solver->AMISolveB(tnext, AMICI_NORMAL);
            if (status != AMICI_SUCCESS)
                return status;

            status = solver->AMIGetB(tdata->which, &(tdata->t), tdata->xB,
                                     tdata->dxB);
            if (status != AMICI_SUCCESS)
                return status;
            status = solver->AMIGetQuadB(tdata->which, &(tdata->t), tdata->xQB);
            if (status != AMICI_SUCCESS)
                return status;
        }

        /* handle discontinuity */

        if (model->ne > 0 && rdata->nmaxevent > 0 && tdata->iroot >= 0) {
            if (tnext == tdata->discs[tdata->iroot]) {
                handleEventB(tdata->iroot, tdata, model);
                --tdata->iroot;
            }
        }

        /* handle data-point */
        if (tnext == rdata->ts[it]) {
            handleDataPointB(it, rdata, tdata, solver, model);
            it--;
        }

        /* reinit states */
        status =
            solver->AMIReInitB(tdata->which, tdata->t, tdata->xB, tdata->dxB);
        if (status != AMICI_SUCCESS)
            return status;

        status = solver->AMIQuadReInitB(tdata->which, tdata->xQB);
        if (status != AMICI_SUCCESS)
            return status;

        status =
            solver->AMICalcICB(tdata->which, tdata->t, tdata->xB, tdata->dxB);
        if (status != AMICI_SUCCESS)
            return status;
    }

    /* we still need to integrate from first datapoint to tstart */
    if (tdata->t > udata->tstart) {
        if (status == AMICI_SUCCESS) {
            if (model->nx > 0) {
                /* solve for backward problems */
                status = solver->AMISolveB(udata->tstart, AMICI_NORMAL);
                if (status != AMICI_SUCCESS)
                    return status;

                status =
                    solver->AMIGetQuadB(tdata->which, &(tdata->t), tdata->xQB);
                if (status != AMICI_SUCCESS)
                    return status;
                status = solver->AMIGetB(tdata->which, &(tdata->t), tdata->xB,
                                         tdata->dxB);
                if (status != AMICI_SUCCESS)
                    return status;
            }
        }
    }
    
    realtype *sx_tmp;
    realtype *xB_tmp = NV_DATA_S(tdata->xB);
    if (!xB_tmp)
        return AMICI_ERROR_ASA;

    for (int iJ = 0; iJ < model->nJ; iJ++) {
        if (iJ == 0) {
            if (rdata->sensi <= AMICI_SENSI_ORDER_FIRST) {
                for (ip = 0; ip < rdata->nplist; ++ip) {
                    tdata->llhS0[iJ * rdata->nplist + ip] = 0.0;
                    sx_tmp = NV_DATA_S(tdata->sx[ip]);
                    if (!sx_tmp)
                        return AMICI_ERROR_ASA;
                    for (ix = 0; ix < model->nxtrue; ++ix) {
                        tdata->llhS0[ip] =
                            tdata->llhS0[ip] + xB_tmp[ix] * sx_tmp[ix];
                    }
                }
            } else {
                for (ip = 0; ip < rdata->nplist; ++ip) {
                    for (jp = 0; jp < rdata->nplist; ++jp) {
                        tdata->llhS20[ip * rdata->nplist + jp] = 0.0;
                        sx_tmp = NV_DATA_S(tdata->s2x[ip * rdata->nplist + jp]);
                        if (!sx_tmp)
                            return AMICI_ERROR_ASA;
                        for (ix = 0; ix < model->nxtrue; ++ix) {
                            tdata->llhS20[ip * rdata->nplist + jp] =
                            tdata->llhS20[ip * rdata->nplist + jp] +
                            xB_tmp[ix] * sx_tmp[ix];
                        }
                    }
                }
            }
        } else {
            for (ip = 0; ip < rdata->nplist; ++ip) {
                tdata->llhS0[ip + iJ * rdata->nplist] = 0.0;
                sx_tmp = NV_DATA_S(tdata->sx[ip]);
                if (!sx_tmp)
                    return AMICI_ERROR_ASA;
                for (ix = 0; ix < model->nxtrue; ++ix) {
                    tdata->llhS0[ip + iJ * rdata->nplist] =
                        tdata->llhS0[ip + iJ * rdata->nplist] +
                        xB_tmp[ix + iJ * model->nxtrue] * sx_tmp[ix] +
                        xB_tmp[ix] * sx_tmp[ix + iJ * model->nxtrue];
                }
            }
        }
    }

    realtype *xQB_tmp = NV_DATA_S(tdata->xQB);
    if (!xQB_tmp)
        return AMICI_ERROR_ASA;

    for (int iJ = 0; iJ < model->nJ; iJ++) {
        for (ip = 0; ip < rdata->nplist; ip++) {
            if (iJ == 0) {
                rdata->sllh[ip] -= tdata->llhS0[ip] + xQB_tmp[ip];
            } else {
                rdata->s2llh[iJ - 1 + ip * (model->nJ - 1)] -=
                    tdata->llhS0[ip + iJ * rdata->nplist] +
                    xQB_tmp[ip + iJ * rdata->nplist];
            }
        }
    }

    return status;
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

int BackwardProblem::handleEventB(int iroot, TempData *tdata, Model *model) {
    /**
     * handleEventB executes everything necessary for the handling of events
     * for the backward problem
     *
     * @param[out] iroot index of event @type int
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @param[in] model pointer to model specification object @type Model
     * @return status flag indicating success of execution @type int
     */

    int status = AMICI_SUCCESS;

    /* store current values */
    N_VScale(1.0, tdata->xB, tdata->xB_old);
    N_VScale(1.0, tdata->xQB, tdata->xQB_old);

    realtype *xB_tmp = NV_DATA_S(tdata->xB);
    if (!xB_tmp)
        return AMICI_ERROR_EVENT;
    realtype *xQB_tmp = NV_DATA_S(tdata->xQB);
    if (!xQB_tmp)
        return AMICI_ERROR_DATA;

    for (int ie = 0; ie < model->ne; ie++) {

        if (tdata->rootidx[iroot * model->ne + ie] != 0) {

            status = model->fdeltaqB(tdata->t, ie, tdata->x_disc[iroot],
                                     tdata->xB_old, tdata->xQB_old,
                                     tdata->xdot_disc[iroot],
                                     tdata->xdot_old_disc[iroot], tdata);
            if (status != AMICI_SUCCESS)
                return status;

            status = model->fdeltaxB(tdata->t, ie, tdata->x_disc[iroot],
                                     tdata->xB_old, tdata->xdot_disc[iroot],
                                     tdata->xdot_old_disc[iroot], tdata);
            if (status != AMICI_SUCCESS)
                return status;

            for (int ix = 0; ix < model->nxtrue; ++ix) {
                for (int iJ = 0; iJ < model->nJ; ++iJ) {
                    xB_tmp[ix + iJ * model->nxtrue] +=
                        tdata->deltaxB[ix + iJ * model->nxtrue];
                    if (model->nz > 0) {
                        xB_tmp[ix + iJ * model->nxtrue] +=
                            tdata->dJzdx[tdata->nroots[ie] +
                                         (iJ + ix * model->nJ) *
                                             tdata->rdata->nmaxevent];
                    }
                }
            }

            for (int iJ = 0; iJ < model->nJ; ++iJ) {
                for (int ip = 0; ip < tdata->rdata->nplist; ++ip) {
                    xQB_tmp[ip + iJ * tdata->rdata->nplist] +=
                        tdata->deltaqB[ip + iJ * tdata->rdata->nplist];
                }
            }

            tdata->nroots[ie]--;
        }
    }

    return updateHeavisideB(iroot, tdata, model->ne);
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

int BackwardProblem::handleDataPointB(int it, ReturnData *rdata,
                                      TempData *tdata, Solver *solver,
                                      Model *model) {
    /**
     * handleDataPoint executes everything necessary for the handling of data
     * points for the backward problems
     *
     * @param[in] it index of data point @type int
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @param[in] solver pointer to solver object @type Solver
     * @param[in] model pointer to model specification object @type Model
     * @return status flag indicating success of execution @type int
     */

    realtype *xB_tmp = NV_DATA_S(tdata->xB);
    if (!xB_tmp)
        return AMICI_ERROR_DATA;
    for (int ix = 0; ix < model->nxtrue; ix++) {
        for (int iJ = 0; iJ < model->nJ; iJ++)
            // we only need the 1:nxtrue slice here!
            xB_tmp[ix + iJ * model->nxtrue] +=
                tdata->dJydx[it + (iJ + ix * model->nJ) * rdata->nt];
    }
    return solver->getDiagnosisB(it, rdata, tdata);
}

/* --------------------------------------------------------------------------------
 */
/* --------------------------------------------------------------------------------
 */
/* --------------------------------------------------------------------------------
 */

int BackwardProblem::updateHeavisideB(int iroot, TempData *tdata, int ne) {
    /**
     * updateHeavisideB updates the heaviside variables h on event occurences
     * for the backward problem
     *
     * @param[in] iroot discontinuity occurance index @type int
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @param[in] ne number of events @type int
     * @return status flag indicating success of execution @type int
     */

    /* tdata->rootsfound provides the direction of the zero-crossing, so adding
       it will give
         the right update to the heaviside variables */

    for (int ie = 0; ie < ne; ie++) {
        tdata->h[ie] -= tdata->rootidx[iroot * ne + ie];
    }
    return AMICI_SUCCESS;
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

realtype BackwardProblem::getTnext(realtype *troot, int iroot, realtype *tdata,
                                   int it, Model *model) {
    /**
     * getTnext computes the next timepoint to integrate to. This is the maximum
     * of
     * tdata and troot but also takes into account if it<0 or iroot<0 where
     * these expressions
     * do not necessarily make sense
     *
     * @param[in] troot timepoint of next event @type realtype
     * @param[in] iroot index of next event @type int
     * @param[in] tdata timepoint of next data point @type realtype
     * @param[in] it index of next data point @type int
     * @param[in] model pointer to model specification object @type Model
     * @return tnext next timepoint @type realtype
     */

    realtype tnext;

    if (it < 0) {
        tnext = troot[iroot];
    } else {
        if (iroot < 0) {
            tnext = tdata[it];
        } else {
            if (model->ne > 0) {
                if (troot[iroot] > tdata[it]) {
                    tnext = troot[iroot];
                } else {
                    tnext = tdata[it];
                }
            } else {
                tnext = tdata[it];
            }
        }
    }

    return (tnext);
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

BackwardProblem::BackwardProblem() {
    /**
     * this is a placeholder, nothing needs to be done at initialization.
     */
}
