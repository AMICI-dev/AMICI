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
#include <time.h>
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


int applyNewtonsMethod(UserData *udata, ReturnData *rdata, TempData *tdata, int newton_try, Model *model) {
    /**
     * applyNewtonsMethod applies Newtons method to the current state x to find the steady state
     *
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[in] tdata pointer to the temporary data struct @type TempData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] newton_try integer for the try number of the Newton solver
     * @return void
     */

    int status = AMICI_ERROR_NEWTONSOLVER;
    int i_newtonstep=0;
    int ix=0;    
    double res_abs;
    double res_rel;
    double res_tmp;
    double gamma = 1.0;
    realtype *x_tmp;
    
    N_Vector delta = N_VNew_Serial(model->nx);
    N_Vector rel_x_newton = N_VNew_Serial(model->nx);
    N_Vector x_newton = N_VNew_Serial(model->nx);
    
    /* initialize output von linear solver for Newton step */
    N_VConst(0.0, delta);
    
    /* Check, how fxdot is used exactly within AMICI... */
    model->fxdot(tdata->t, tdata->x, tdata->dx, tdata->xdot, tdata);
    res_abs = sqrt(N_VDotProd(tdata->xdot,tdata->xdot));
    
    /* Check for relative error, but make sure not to divide by 0!
    Ensure positivity of the state */
    N_VScale(1.0, tdata->x, x_newton);
    N_VAbs(x_newton, x_newton);
    x_tmp = N_VGetArrayPointer(x_newton);
    for (ix=0; ix<model->nx; ix++) {
        if (x_tmp[ix] < udata->atol) {
            x_tmp[ix] = udata->atol;
        }
    }
    N_VDiv(tdata->xdot, x_newton, rel_x_newton);
    res_rel = sqrt(N_VDotProd(rel_x_newton, rel_x_newton));
    
    if (res_abs >= udata->atol && res_rel >= udata->rtol) {
        /* If Newton steps are necessary, compute the inital search direction */
        status = getNewtonStep(udata, rdata, tdata, newton_try, i_newtonstep, delta, model);
    
        if (status == AMICI_SUCCESS) {
            /* The linear solver was successful, now the Newton solver needs to be */
            status = AMICI_ERROR_NEWTONSOLVER;
            
            /* Copy the current state to the old one, make up a new vector for JDiag */
            N_VScale(1.0, tdata->x, tdata->x_old);
            N_VScale(1.0, tdata->xdot, tdata->xdot_old);
    
            /* Newton iterations */
            for(i_newtonstep=0; i_newtonstep<udata->newton_maxsteps; i_newtonstep++) {
        
                /* Try a full, undamped Newton step */
                N_VLinearSum(1.0, tdata->x_old, gamma, delta, tdata->x);
        
                /* Ensure positivity of the state */
                x_tmp = N_VGetArrayPointer(tdata->x);
                for (ix=0; ix<model->nx; ix++) {
                    if (x_tmp[ix] < 0.0) {
                        x_tmp[ix] = 0.0;
                    }
                }
        
                /* Compute new xdot */
                model->fxdot(tdata->t, tdata->x, tdata->dx, tdata->xdot, tdata);
        
                /* Check if new residuals are smaller than old ones */
                res_tmp = sqrt(N_VDotProd(tdata->xdot, tdata->xdot));

                if (res_tmp<res_abs) {
                    /* update state */
                    res_abs = res_tmp;
                    N_VScale(1.0, tdata->x, tdata->x_old);
                    N_VScale(1.0, tdata->xdot, tdata->xdot_old);
            
                    /* Check residuals vs tolerances */
                    if (res_abs < udata->atol) {
                        /* Return number of Newton steps */
                        rdata->newton_numsteps[newton_try-1] = i_newtonstep + 1;
                        status = AMICI_SUCCESS;
                        break;
                    }
            
                    if (status != AMICI_SUCCESS) {
                        /* increase dampening factor */
                        gamma = fmax(1.0, 2.0*gamma);
            
                        /* Do another Newton step */
                        status = getNewtonStep(udata, rdata, tdata, newton_try, i_newtonstep, delta, model);
                        if (status == AMICI_SUCCESS) {
                            /* Newton step was successful, now Newtons method still needs to be */
                            status = AMICI_ERROR_NEWTONSOLVER;
                        } else {
                            /* Linear solver errored, go to clean up and return part */
                            rdata->newton_numsteps[newton_try-1] = amiGetNaN();
                            break;
                        }
                    }
                } else {
                    /* Reduce dampening factor */
                    gamma = gamma/4.0;
                }
            }
        
            /* Set return values */
            rdata->newton_numsteps[newton_try-1] = i_newtonstep;
        } else {
            rdata->newton_numsteps[newton_try-1] = amiGetNaN();
        }
        
    } else {
        /* No Newton steps were necessary */
        status = AMICI_SUCCESS;
        
        /* Set return values */
        rdata->newton_numsteps[newton_try-1] = 0.0;
    }

    /* Clean up worksapce */
    N_VDestroy_Serial(delta);
    N_VDestroy_Serial(rel_x_newton);
    N_VDestroy_Serial(x_newton);
    
    return(status);
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int getNewtonStep(UserData *udata, ReturnData *rdata, TempData *tdata, int ntry, int nnewt, N_Vector ns_delta, Model *model) {
    /**
     * getNewtonStep acomputes the Newton Step by solving a linear system
     *
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] tdata pointer to the temporary data struct @type TempData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @param[in] ntry intger number of Newton solver try
     * @param[in] nnewt intger number of Newton steps in the current Newton solver try
     * @param[out] delta N_Vector solution of the linear system
     * @return int status flag indicating success of execution @type int
     */
    
     int status = AMICI_ERROR_NEWTON_LINSOLVER;
     double rho;
     double rho1;
     double alpha;
     double beta;
     double omega;
     double res;
    
     N_Vector ns_p = N_VNew_Serial(model->nx);
     N_Vector ns_h = N_VNew_Serial(model->nx);
     N_Vector ns_t = N_VNew_Serial(model->nx);
     N_Vector ns_s = N_VNew_Serial(model->nx);
     N_Vector ns_r = N_VNew_Serial(model->nx);
     N_Vector ns_rt = N_VNew_Serial(model->nx);
     N_Vector ns_v = N_VNew_Serial(model->nx);
     N_Vector ns_Jv = N_VNew_Serial(model->nx);
     N_Vector ns_tmp = N_VNew_Serial(model->nx);
     N_Vector ns_Jdiag = N_VNew_Serial(model->nx);
    
     N_VScale(-1.0, tdata->xdot, tdata->xdot);
    
    // Get the diagonal of the Jacobian for preconditioning
    model->fJDiag(tdata->t, ns_Jdiag, tdata->x, tdata);
    
    // Ensure positivity of entries in ns_Jdiag
    N_VConst(1.0, ns_p);
    N_VAbs(ns_Jdiag, ns_tmp);
    N_VCompare(1e-15, ns_tmp, ns_tmp);
    N_VLinearSum(-1.0, ns_tmp, 1.0, ns_p, ns_tmp);
    N_VLinearSum(1.0, ns_Jdiag, 1.0, ns_tmp, ns_Jdiag);
    
    
    // Initialize for linear solve
    N_VConst(0.0, ns_p);
    N_VConst(0.0, ns_v);
    N_VConst(0.0, ns_delta);
    N_VConst(0.0, ns_tmp);
    rho = 1.0;
    omega = 1.0;
    alpha = 1.0;
    
    // can be set to 0 at the moment
    model->fJv(ns_delta, ns_Jv, tdata->t, tdata->x, tdata->xdot, tdata, ns_tmp);
    
    // ns_r = xdot - ns_Jv;
    N_VLinearSum(-1.0, ns_Jv, 1.0, tdata->xdot, ns_r);
    N_VDiv(ns_r, ns_Jdiag, ns_r);
    res = sqrt(N_VDotProd(ns_r, ns_r));
    N_VScale(1.0, ns_r, ns_rt);
    
    for (int i_linstep = 0; i_linstep < udata->newton_maxlinsteps; i_linstep++) {
        // Compute factors
        rho1 = rho;
        rho = N_VDotProd(ns_rt, ns_r);
        beta = rho*alpha / (rho1*omega);
        
        // ns_p = ns_r + beta * (ns_p - omega * ns_v);
        N_VLinearSum(1.0, ns_p, -omega, ns_v, ns_p);
        N_VLinearSum(1.0, ns_r, beta, ns_p, ns_p);
        
        // ns_v = J * ns_p
        model->fJv(ns_p, ns_v, tdata->t, tdata->x, tdata->xdot, tdata, ns_tmp);
        N_VDiv(ns_v, ns_Jdiag, ns_v);
        
        // Compute factor
        alpha = rho / N_VDotProd(ns_rt, ns_v);
        
        // ns_h = ns_delta + alpha * ns_p;
        N_VLinearSum(1.0, ns_delta, alpha, ns_p, ns_h);
        // ns_s = ns_r - alpha * ns_v;
        N_VLinearSum(1.0, ns_r, -alpha, ns_v, ns_s);
        
        // ns_t = J * ns_s
        model->fJv(ns_s, ns_t, tdata->t, tdata->x, tdata->xdot, tdata, ns_tmp);
        N_VDiv(ns_t, ns_Jdiag, ns_t);
        
        // Compute factor
        omega = N_VDotProd(ns_t, ns_s) / N_VDotProd(ns_t, ns_t);
        
        // ns_delta = ns_h + omega * ns_s;
        N_VLinearSum(1.0, ns_h, omega, ns_s, ns_delta);
        // ns_r = ns_s - omega * ns_t;
        N_VLinearSum(1.0, ns_s, -omega, ns_t, ns_r);
        
        // Compute the (unscaled) residual
        N_VProd(ns_r, ns_Jdiag, ns_r);
        res = sqrt(N_VDotProd(ns_r, ns_r));
        
        // Test convergence
        if (res < udata->atol) {
            // Write number of steps needed
            rdata->newton_numlinsteps[(ntry-1) * udata->newton_maxsteps + nnewt] = i_linstep + 1;

            // Return success
            N_VScale(-1.0, tdata->xdot, tdata->xdot);
            status = AMICI_SUCCESS;
            break;
        }
        
        // Scale back
        N_VDiv(ns_r, ns_Jdiag, ns_r);
    }
    
    // Clean up workspace
    N_VDestroy_Serial(ns_p);
    N_VDestroy_Serial(ns_h);
    N_VDestroy_Serial(ns_t);
    N_VDestroy_Serial(ns_s);
    N_VDestroy_Serial(ns_r);
    N_VDestroy_Serial(ns_rt);
    N_VDestroy_Serial(ns_v);
    N_VDestroy_Serial(ns_Jv);
    N_VDestroy_Serial(ns_tmp);
    N_VDestroy_Serial(ns_Jdiag);
    
    N_VScale(-1.0, tdata->xdot, tdata->xdot);
    
    // Return
    return(status);
 }

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int getNewtonOutput(TempData *tdata, ReturnData *rdata, int newton_status, double run_time, int nx) {
    /**
     * getNewtonOutput stores the output of the Newton solver run.
     *
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[in] tdata pointer to the temporary data struct @type TempData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] newton_status integer flag indicating the run of the Newton solver
     * @param[in] run_time double computation time of the Newton solver
     * @return int status flag indicating success of execution @type int
     */
    
    realtype *x_tmp;
    
    // Get time for Newton solve
    rdata->newton_time[0] = run_time;
    
    // Since the steady state was found, current time is set to infinity
    tdata->t = INFINITY;
    
    // Write output
    x_tmp = N_VGetArrayPointer(tdata->x);
    for (int ix=0; ix<nx; ix++) {
        rdata->xss[ix] = x_tmp[ix];
    }
    
    // Write flag for the Newton solver
    *rdata->newton_status = (double) newton_status;
    
    return(AMICI_SUCCESS);
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int getNewtonSimulation(UserData *udata, TempData *tdata, ReturnData *rdata, Solver *solver, Model *model) {
    /**
     * getNewtonSimulation solves the forward problem, if the first Newton solver run did not succeed.
     *
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[in] tdata pointer to the temporary data struct @type TempData
     * @param[in] rdata pointer to the return data struct @type ReturnData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @return int status flag indicating success of execution @type int
     */
    
    
    double res_abs;
    double res_rel;
    double sim_time;
    int status = (int) *rdata->status;
    realtype *x_tmp;
    N_Vector rel_x = N_VNew_Serial(model->nx);
    N_Vector tmp_x = N_VNew_Serial(model->nx);
    
    /* Newton solver did not work, so try a simulation */
    if (tdata->t >= 1e6) {
        sim_time = 10.0*(tdata->t);
    } else {
        sim_time = 1e6;
    }
    status = solver->AMISolve(RCONST(sim_time), tdata->x, tdata->dx, &(tdata->t), AMICI_NORMAL);
    
    if (status == AMICI_SUCCESS) {
        /* Check residuals */
        res_abs = sqrt(N_VDotProd(tdata->xdot, tdata->xdot));
    
        /* Ensure positivity for relative residual */
        N_VScale(1.0, tdata->x, tmp_x);
        N_VAbs(tmp_x, tmp_x);
        x_tmp = N_VGetArrayPointer(tmp_x);
        for (int ix=0; ix<model->nx; ix++) {
            if (x_tmp[ix] < udata->atol) {
                x_tmp[ix] = udata->atol;
            }
        }
        N_VDiv(tdata->xdot, tmp_x, rel_x);
        res_rel = sqrt(N_VDotProd(rel_x, rel_x));
    
        /* residuals are small? */
        if (res_abs < udata->atol || res_rel < udata->rtol) {
            return(AMICI_SUCCESS);
        } else {
            return(-1);
        }
    } else {
        return(status);
    }
}

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

int workSteadyStateProblem(UserData *udata, TempData *tdata, ReturnData *rdata, int it, Solver *solver, Model *model) {
    /**
     * tries to determine the steady state of the ODE system by a Newton solver
     * uses forward intergration, if the Newton solver fails
     *
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[in] tdata pointer to the temporary data struct @type UserData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     */
    
    int status = (int) *rdata->status;
    double run_time;
    clock_t starttime;
    
    /* First, try to do Newton steps */
    starttime = clock();
    status = applyNewtonsMethod(udata, rdata, tdata, 1, model);
    
    if (status == AMICI_SUCCESS) {
        /* if the Newton solver found a steady state */
        run_time = (double)((clock() - starttime) * 1000) / CLOCKS_PER_SEC;
        status = getNewtonOutput(tdata, rdata, 1, run_time, model->nx);
    } else {
        /* Newton solver did not find a steady state, so try integration */
        status = getNewtonSimulation(udata, tdata, rdata, solver, model);
        
        if (status == AMICI_SUCCESS) {
            /* if simulation found a steady state */
            run_time = (double)((clock() - starttime) * 1000) / CLOCKS_PER_SEC;
            status = getNewtonOutput(tdata, rdata, 2, run_time, model->nx);
        } else {
            status = applyNewtonsMethod(udata, rdata, tdata, 2, model);
            
            if (status == AMICI_SUCCESS) {
                /* If the second Newton solver found a steady state */
                run_time = (double)((clock() - starttime) * 1000) / CLOCKS_PER_SEC;
                status = getNewtonOutput(tdata, rdata, 3, run_time, model->nx);
            } else {
                /* integration error occured */
                return(status);
            }
        }
    }
    
    /* if this point was reached, the Newton solver was successful */
    return(AMICI_SUCCESS);
}


void printErrMsgIdAndTxt(const char * identifier, const char *msg, ...) {
    printf("[Error] %s: %s\n", identifier, msg);
}

void printWarnMsgIdAndTxt(const char * identifier, const char *msg, ...) {
    printf("[Warning] %s: %s\n", identifier, msg);
}
