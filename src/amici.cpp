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
#include "include/backwardproblem.h"

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
    if (status == AMICI_SUCCESS) status = BackwardProblem::workBackwardProblem(udata, tdata, rdata, solver, model);
    
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


void printErrMsgIdAndTxt(const char * identifier, const char *msg, ...) {
    printf("[Error] %s: %s\n", identifier, msg);
}

void printWarnMsgIdAndTxt(const char * identifier, const char *msg, ...) {
    printf("[Warning] %s: %s\n", identifier, msg);
}
