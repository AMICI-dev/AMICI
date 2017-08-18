/**
 * @file   amici.cpp
 * @brief  core routines for integration
 */

#include <cstdlib>
#include <cstring>
#include <cassert>
#define _USE_MATH_DEFINES /* MS definition of PI and other constants */
#include <cmath>
#ifndef M_PI /* define PI if we still have no definition */
#define M_PI 3.14159265358979323846
#endif

#include <stdio.h>
#include <include/amici.h> /* amici functions */
#include <include/symbolic_functions.h>
#include <include/amici_misc.h>
#include "include/amici_solver.h"
#include "include/amici_model.h"
#include "include/forwardproblem.h"

msgIdAndTxtFp errMsgIdAndTxt = &printErrMsgIdAndTxt;
msgIdAndTxtFp warnMsgIdAndTxt = &printWarnMsgIdAndTxt;

int runAmiciSimulation(UserData *udata, const ExpData *edata, ReturnData *rdata, Model *model, Solver *solver) {
    if(!udata) return AMICI_ERROR_UDATA;
    if(!rdata) return AMICI_ERROR_RDATA;
    
    int status = AMICI_SUCCESS;
    
    if (model->nx <= 0) {
        return AMICI_ERROR_NOTHINGTODO;
    }
    
    TempData *tdata = new TempData(udata, model);
    
    // unscale parameters but keep original
    double *originalParams = NULL;
    if(model->pscale != AMICI_SCALING_NONE) {
        originalParams = (double *) malloc(sizeof(double) * model->np);
        memcpy(originalParams, udata->p, sizeof(double) * model->np);
    }
    status = udata->unscaleParameters(model);

    if (status == AMICI_SUCCESS)
        status = solver->setupAMI(udata, tdata, model);

    if (status != AMICI_SUCCESS)
        goto freturn;

    if (status == AMICI_SUCCESS) status = ForwardProblem::workForwardProblem(udata, tdata, rdata, edata, solver, model);
    if (status == AMICI_SUCCESS) status = workBackwardProblem(udata, tdata, rdata, solver, model);
    
    if (status == AMICI_SUCCESS) status = rdata->applyChainRuleFactorToSimulationResults(udata);
    if (status < AMICI_SUCCESS) rdata->invalidate();
    
    
freturn:
    // reset to original parameters
    if(originalParams) {
        memcpy(udata->p, originalParams, sizeof(double) * model->np);
        free(originalParams);
    }

    delete tdata;

    return status;
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int prepDataSensis(int it, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata, Model *model) {
    /**
     * prepDataSensis preprocesses the provided experimental data to compute sensitivities via adjoint or forward methods later on
     *
     * @param[in] it index of current timepoint @type int
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */
        
    int status = model->fdydx(rdata->ts[it],it,tdata->x,tdata);
    if(status != AMICI_SUCCESS) return status;

    status = model->fdydp(rdata->ts[it],it,tdata->x,tdata);
    if(status != AMICI_SUCCESS) return status;

    if (!edata)
        return status;

    status = model->fdsigma_ydp(tdata->t,tdata);
    if(status != AMICI_SUCCESS) return status;

    for (int iy=0; iy<model->nytrue; iy++) {
        if (!amiIsNaN(edata->sigmay[iy*udata->nt+it])) {
            for (int ip=0; ip<udata->nplist; ip++) {
                tdata->dsigmaydp[ip*model->ny+iy] = 0;
            }
        }
        for (int ip=0; ip<udata->nplist; ip++) {
            rdata->ssigmay[it + udata->nt*(ip*model->ny+iy)] = tdata->dsigmaydp[ip*model->ny+iy];
        }
    }
    status = model->fdJydy(tdata->t,it,tdata->x,tdata,edata,rdata);
    if(status != AMICI_SUCCESS) return status;

    status = model->fdJydsigma(tdata->t,it,tdata->x,tdata,edata,rdata);
    if(status != AMICI_SUCCESS) return status;

    status = model->fdJydx(it,udata,tdata,edata);
    if(status != AMICI_SUCCESS) return status;

    status = model->fdJydp(it,udata,tdata,edata,rdata);
    if(status != AMICI_SUCCESS) return status;

    if (udata->sensi_meth != AMICI_SENSI_ASA)
        return status;

    for(int iJ=0; iJ<model->nJ; iJ++) {
        for(int ip=0; ip < udata->nplist; ip++) {
            if (iJ==0) {
                if (model->ny>0) {
                    rdata->sllh[ip] -= tdata->dJydp[ip * model->nJ];
                }
            } else {
                if (model->ny>0) {
                    rdata->s2llh[(iJ - 1) + ip * (model->nJ-1) ] -= tdata->dJydp[iJ + ip * model->nJ];
                }
            }
        }
    }

    return status;
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int getDataOutput(int it, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata, Solver *solver, Model *model) {
    /**
     * getDataOutput extracts output information for data-points
     *
     * @param[in] it index of current timepoint @type int
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */
       
    int status = model->fy(rdata->ts[it],it,tdata->x,tdata,rdata);
    if(status != AMICI_SUCCESS) return status;
    
    if (edata) {
        status = model->fsigma_y(tdata->t,tdata);
        if(status != AMICI_SUCCESS) return status;
        for (int iy=0; iy<model->nytrue; iy++) {
            /* extract the value for the standard deviation, if the data value is NaN, use
             the parameter value. Store this value in the return struct */
            if (!amiIsNaN(edata->sigmay[iy*rdata->nt+it])) {
                tdata->sigmay[iy] = edata->sigmay[iy*rdata->nt+it];
            }
            rdata->sigmay[iy*rdata->nt+it] = tdata->sigmay[iy];
        }
        status = model->fJy(rdata->ts[it],it,tdata->x,tdata,edata,rdata);
        if(status != AMICI_SUCCESS) return status;
    } else {
        status = model->fsigma_y(tdata->t,tdata);
        if(status != AMICI_SUCCESS) return status;
        for (int iy=0; iy<model->nytrue; iy++) {
            rdata->sigmay[iy*rdata->nt+it] = tdata->sigmay[iy];
        }
    }
    if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        status = prepDataSensis(it, udata, rdata, edata, tdata, model);
        if(status != AMICI_SUCCESS) return status;
        if (rdata->sensi_meth == AMICI_SENSI_FSA) {
            status = getDataSensisFSA(it, udata, rdata, edata, tdata, solver, model);
            if(status != AMICI_SUCCESS) return status;
        }
    }
    return AMICI_SUCCESS;
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int getDataSensisFSA(int it, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata, Solver *solver, Model *model) {
    /**
     * getDataSensisFSA extracts data information for forward sensitivity analysis
     *
     * @param[out] status flag indicating success of execution @type int
     * @param[in] it index of current timepoint @type int
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return void
     */
    
    int status = AMICI_SUCCESS;
    realtype *sx_tmp;
    
    for(int ip=0; ip < udata->nplist; ip++) {
        if (model->nx>0) {
            if (rdata->ts[it] > udata->tstart) {
                status = solver->AMIGetSens(&(tdata->t), tdata->sx);
                if (status != AMICI_SUCCESS) return status;
            }
            
            sx_tmp = NV_DATA_S(tdata->sx[ip]);
            if(!sx_tmp) return AMICI_ERROR_FSA;
            for(int ix=0; ix < model->nx; ix++) {
                rdata->sx[(ip*model->nx + ix)*udata->nt + it] = sx_tmp[ix];
            }
        }
    }
    
    for (int iy=0; iy<model->nytrue; iy++) {
        if (edata){
            if (amiIsNaN(edata->sigmay[iy*udata->nt+it])) {
                status = model->fdsigma_ydp(tdata->t,tdata);
                if(status != AMICI_SUCCESS) return status;
            } else {
                for (int ip=0; ip<udata->nplist; ip++) {
                    tdata->dsigmaydp[ip*model->ny+iy] = 0;
                }
            }
            for (int ip=0; ip<udata->nplist; ip++) {
                rdata->ssigmay[it + udata->nt*(ip*model->ny+iy)] = tdata->dsigmaydp[ip*model->ny+iy];
            }
        } else {
            for (int ip=0; ip<udata->nplist; ip++) {
                rdata->ssigmay[it + udata->nt*(ip*model->ny+iy)] = 0;
            }
        }
    }
    status = model->fsy(it,udata,tdata,rdata);
    if(status != AMICI_SUCCESS) return status;
    if (edata) {
        status = model->fsJy(it,udata,tdata,edata,rdata);
        if(status != AMICI_SUCCESS) return status;
    }
    return status;
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */


/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */



/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int handleDataPointB(int it, UserData *udata, ReturnData *rdata, TempData *tdata, Solver *solver, Model *model) {
    /**
     * handleDataPoint executes everything necessary for the handling of data points for the backward problems
     *
     * @param[in] it index of data point @type int
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */
    
    realtype *xB_tmp = NV_DATA_S(tdata->xB);
    if(!xB_tmp) return AMICI_ERROR_DATA;
    for (int ix=0; ix<model->nxtrue; ix++) {
        for(int iJ=0; iJ<model->nJ; iJ++)
            // we only need the 1:nxtrue slice here!
            xB_tmp[ix + iJ * model->nxtrue] += tdata->dJydx[it + (iJ + ix * model->nJ) * udata->nt];
    }
    return solver->getDiagnosisB(it,udata,rdata,tdata);
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int handleEventB(int iroot, UserData *udata, TempData *tdata, Model *model) {
    /**
     * handleEventB executes everything necessary for the handling of events for the backward problem
     *
     * @param[out] iroot index of event @type int
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */
    
    int status = AMICI_SUCCESS;
    
    /* store current values */
    N_VScale(1.0,tdata->xB,tdata->xB_old);
    N_VScale(1.0,tdata->xQB,tdata->xQB_old);
    
    realtype *xB_tmp = NV_DATA_S(tdata->xB);
    if(!xB_tmp) return AMICI_ERROR_EVENT;
    realtype *xQB_tmp = NV_DATA_S(tdata->xQB);
    if(!xQB_tmp) return AMICI_ERROR_DATA;
    
    for (int ie=0; ie<model->ne; ie++) {
        
        if (tdata->rootidx[iroot*model->ne + ie] != 0) {
            
            status = model->fdeltaqB(tdata->t,ie,tdata->x_disc[iroot],tdata->xB_old,tdata->xQB_old,tdata->xdot_disc[iroot],tdata->xdot_old_disc[iroot],tdata);
            if (status != AMICI_SUCCESS) return status;

            status = model->fdeltaxB(tdata->t,ie,tdata->x_disc[iroot],tdata->xB_old,tdata->xdot_disc[iroot],tdata->xdot_old_disc[iroot],tdata);
            if (status != AMICI_SUCCESS) return status;
            
            for (int ix=0; ix<model->nxtrue; ++ix) {
                for (int iJ = 0; iJ < model->nJ; ++iJ) {
                    xB_tmp[ix + iJ*model->nxtrue] += tdata->deltaxB[ix + iJ*model->nxtrue];
                    if (model->nz>0) {
                        xB_tmp[ix + iJ*model->nxtrue] += tdata->dJzdx[tdata->nroots[ie] + (iJ + ix * model->nJ) * udata->nmaxevent];
                    }
                }
            }
            
            for (int iJ=0; iJ<model->nJ; ++iJ) {
                for (int ip=0; ip<udata->nplist; ++ip) {
                    xQB_tmp[ip + iJ*udata->nplist] += tdata->deltaqB[ip + iJ*udata->nplist];
                }
            }
            
            
            tdata->nroots[ie]--;
        }
    }
    
    return updateHeavisideB(iroot, tdata, model->ne);
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

realtype getTnext(realtype *troot, int iroot, realtype *tdata, int it, Model *model) {
    /**
     * getTnext computes the next timepoint to integrate to. This is the maximum of
     * tdata and troot but also takes into account if it<0 or iroot<0 where these expressions
     * do not necessarily make sense
     *
     * @param[in] troot timepoint of next event @type realtype
     * @param[in] iroot index of next event @type int
     * @param[in] tdata timepoint of next data point @type realtype
     * @param[in] it index of next data point @type int
     * @param[in] udata pointer to the user data struct @type UserData
     * @return tnext next timepoint @type realtype
     */
    
    realtype tnext;
    
    
    if (it<0) {
        tnext = troot[iroot];
    } else {
        if (iroot<0) {
            tnext = tdata[it];
        } else {
            if (model->ne>0) {
                if (troot[iroot]>tdata[it]) {
                    tnext = troot[iroot];
                } else {
                    tnext = tdata[it];
                }
            } else {
                tnext = tdata[it];
            }
        }
    }
    
    return(tnext);
    
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */


/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */


/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int updateHeavisideB(int iroot, TempData *tdata, int ne) {
    /**
     * updateHeavisideB updates the heaviside variables h on event occurences for the backward problem
     *
     * @param[in] iroot discontinuity occurance index @type int
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[ne] number of events
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */
        
    /* tdata->rootsfound provides the direction of the zero-crossing, so adding it will give
     the right update to the heaviside variables */
    
    for (int ie = 0; ie<ne; ie++) {
        tdata->h_udata[ie] -= tdata->rootidx[iroot*ne + ie];
    }
    return AMICI_SUCCESS;
}


/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */


/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */


int workBackwardProblem(UserData *udata, TempData *tdata, ReturnData *rdata, Solver *solver, Model *model) {
    /**
     * workBackwardProblem solves the backward problem. if adjoint sensitivities are enabled this will also compute sensitivies
     * workForwardProblem should be called before this is function is called
     *
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[in] tdata pointer to the temporary data struct @type TempData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[out] edata pointer to the experimental data struct @type ExpData
     * @param[in] iroot pointer to the current root index, the value pointed to will be decreased during the forward solve
     * @return int status flag
     */
    int ix, it, ip;
    int status = (int) *rdata->status;
    double tnext;
    
    if (model->nx <= 0
            || udata->sensi < AMICI_SENSI_ORDER_FIRST
            || udata->sensi_meth != AMICI_SENSI_ASA
            || status != AMICI_SUCCESS) {
        return status;
    }

    solver->setupAMIB(udata, tdata, model);

    it = udata->nt-2;
    --tdata->iroot;
    while (it>=0 || tdata->iroot>=0) {

        /* check if next timepoint is a discontinuity or a data-point */
        tnext = getTnext(tdata->discs, tdata->iroot, udata->ts, it, model);

        if (tnext<tdata->t) {
            status = solver->AMISolveB(tnext, AMICI_NORMAL);
            if (status != AMICI_SUCCESS) return status;


            status = solver->AMIGetB(tdata->which, &(tdata->t), tdata->xB, tdata->dxB);
            if (status != AMICI_SUCCESS) return status;
            status = solver->AMIGetQuadB(tdata->which, &(tdata->t), tdata->xQB);
            if (status != AMICI_SUCCESS) return status;
        }

        /* handle discontinuity */

        if (model->ne>0 && udata->nmaxevent>0 && tdata->iroot >=0) {
            if (tnext == tdata->discs[tdata->iroot]) {
                handleEventB(tdata->iroot, udata, tdata, model);
                --tdata->iroot;
            }
        }

        /* handle data-point */
        if (tnext == udata->ts[it]) {
            handleDataPointB(it, udata, rdata, tdata, solver, model);
            it--;
        }

        /* reinit states */
        status = solver->AMIReInitB(tdata->which, tdata->t, tdata->xB, tdata->dxB);
        if (status != AMICI_SUCCESS) return status;

        status = solver->AMIQuadReInitB(tdata->which, tdata->xQB);
        if (status != AMICI_SUCCESS) return status;

        status = solver->AMICalcICB(tdata->which, tdata->t, tdata->xB, tdata->dxB);
        if (status != AMICI_SUCCESS) return status;
    }

    /* we still need to integrate from first datapoint to tstart */
    if (tdata->t>udata->tstart) {
        if (status == AMICI_SUCCESS) {
            if (model->nx>0) {
                /* solve for backward problems */
                status = solver->AMISolveB(udata->tstart, AMICI_NORMAL);
                if (status != AMICI_SUCCESS) return status;

                status = solver->AMIGetQuadB(tdata->which, &(tdata->t), tdata->xQB);
                if (status != AMICI_SUCCESS) return status;
                status = solver->AMIGetB(tdata->which, &(tdata->t), tdata->xB, tdata->dxB);
                if (status != AMICI_SUCCESS) return status;
            }
        }
    }

    status = model->fx0(tdata->x,tdata);
    if (status != AMICI_SUCCESS) return status;
    status = model->fdx0(tdata->x,tdata->dx,tdata);
    if (status != AMICI_SUCCESS) return status;
    status = model->fsx0(tdata->sx, tdata->x, tdata->dx, tdata);
    if (status != AMICI_SUCCESS) return status;

    realtype *xB_tmp = NV_DATA_S(tdata->xB);
    if(!xB_tmp) return AMICI_ERROR_ASA;
    realtype *sx_tmp;

    for (int iJ=0; iJ<model->nJ; iJ++) {
        if (iJ==0) {
            for (ip=0; ip<udata->nplist; ++ip) {
                tdata->llhS0[iJ*udata->nplist + ip] = 0.0;
                sx_tmp = NV_DATA_S(tdata->sx[ip]);
                if(!sx_tmp) return AMICI_ERROR_ASA;
                for (ix = 0; ix < model->nxtrue; ++ix) {
                    tdata->llhS0[ip] = tdata->llhS0[ip] + xB_tmp[ix] * sx_tmp[ix];
                }
            }
        } else {
            for (ip=0; ip<udata->nplist; ++ip) {
                tdata->llhS0[ip + iJ * udata->nplist] = 0.0;
                sx_tmp = NV_DATA_S(tdata->sx[ip]);
                if(!sx_tmp) return AMICI_ERROR_ASA;
                for (ix = 0; ix < model->nxtrue; ++ix) {
                    tdata->llhS0[ip + iJ * udata->nplist] = tdata->llhS0[ip + iJ * udata->nplist]
                            + xB_tmp[ix + iJ * model->nxtrue] * sx_tmp[ix]
                            + xB_tmp[ix] * sx_tmp[ix + iJ * model->nxtrue];
                }
            }
        }
    }

    realtype *xQB_tmp = NV_DATA_S(tdata->xQB);
    if(!xQB_tmp) return AMICI_ERROR_ASA;

    for(int iJ=0; iJ<model->nJ; iJ++) {
        for(ip=0; ip < udata->nplist; ip++) {
            if (iJ==0) {
                rdata->sllh[ip] -=  tdata->llhS0[ip] + xQB_tmp[ip];
            } else {
                rdata->s2llh[iJ-1 + ip*(model->nJ-1)] -= tdata->llhS0[ip + iJ*udata->nplist] + xQB_tmp[ip + iJ*udata->nplist];
            }
        }
    }

    return status;
}


void printErrMsgIdAndTxt(const char * identifier, const char *msg, ...) {
    printf("[Error] %s: %s\n", identifier, msg);
}

void printWarnMsgIdAndTxt(const char * identifier, const char *msg, ...) {
    printf("[Warning] %s: %s\n", identifier, msg);
}
