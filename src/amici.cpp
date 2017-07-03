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
#include "wrapfunctions.h" /* user functions */
#include <include/amici.h> /* amici functions */
#include <include/symbolic_functions.h>

static int fsy(realtype t, int it, realtype *sy, realtype *dydx, realtype *dydp, N_Vector *sx, void *user_data);
static int fsJy(realtype t, int it, realtype *sJy, realtype *s2Jy, realtype *dJydy, realtype *dJydp, realtype *y, realtype *sigma_y, realtype *sy, realtype *dydp, realtype *my, void *user_data);

int runAmiciSimulation(UserData *udata, const ExpData *edata, ReturnData *rdata) {
    if(!udata) return AMICI_ERROR_UDATA;
    if(!rdata) return AMICI_ERROR_RDATA;
    
    int status = AMICI_SUCCESS;
    int iroot = 0;
    
    if (udata->nx <= 0) {
        return AMICI_ERROR_NOTHINGTODO;
    }
    
    TempData *tdata = new TempData(udata);
    
    status = unscaleParameters(udata);
    if (status == AMICI_SUCCESS) udata->initTemporaryFields();
    
    /* pointer to cvodes memory block */
    void *ami_mem = setupAMI(udata, tdata);
    if (ami_mem == NULL){
        status = AMICI_ERROR_SETUP;
        goto freturn;
    }

    if (status == AMICI_SUCCESS) status = workForwardProblem(udata, tdata, rdata, edata, ami_mem, &iroot);
    if (status == AMICI_SUCCESS) status = workBackwardProblem(udata, tdata, rdata, edata, ami_mem, &iroot);
    
    if (status == AMICI_SUCCESS) status = applyChainRuleFactorToSimulationResults(udata, rdata, edata);
    if (status < AMICI_SUCCESS) invalidateReturnData(udata, rdata);

    if (ami_mem) AMIFree(&ami_mem);
    
freturn:
    udata->freeTemporaryFields();
    delete tdata;
    return status;
}

void invalidateReturnData(UserData* udata, ReturnData* rdata) {
    /**
     * @brief performs all necessary actions to reset return data upon integration failure
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     */
    if (rdata->llh)
    *rdata->llh = amiGetNaN();
    
    if (rdata->sllh)
    fillArray(rdata->sllh, udata->nplist, amiGetNaN());
    
    if (rdata->s2llh)
    fillArray(rdata->s2llh, udata->nplist*(udata->nJ-1), amiGetNaN());
}

void *setupAMI(UserData *udata, TempData *tdata) {
    /**
     * @brief setupAMIs initialises the ami memory object
     * @param[out] status flag indicating success of execution @type int
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[in] tdata pointer to the temporary data struct @type TempData
     * @return ami_mem pointer to the cvodes/idas memory block
     */
    void *ami_mem = NULL; /* pointer to ami memory block */
    N_Vector id;
    
    
    tdata->t = udata->tstart;
    
    if (udata->nx > 0) {
        /* initialise states */
        if (tdata->x == NULL) goto freturn;
        if (udata->x0data == NULL) {
            if (fx0(tdata->x, udata) != AMICI_SUCCESS) goto freturn;
        } else {
            int ix;
            realtype *x_tmp = NV_DATA_S(tdata->x);
            for (ix=0; ix < udata->nx; ix++) {
                x_tmp[ix] = (realtype) udata->x0data[ix];
            }
        }
        if (fdx0(tdata->x, tdata->dx, udata) != AMICI_SUCCESS) goto freturn;
        
        /* initialise heaviside variables */
        if (initHeaviside(udata,tdata) != AMICI_SUCCESS) goto freturn;
        
    }
    
    /* Create AMIS object */
    if (udata->lmm != CV_ADAMS && udata->lmm != CV_BDF) {
        errMsgIdAndTxt("AMICI:mex:lmm","Illegal value for lmm!");
        goto freturn;
    }
    if (udata->iter != CV_NEWTON && udata->iter != CV_FUNCTIONAL) {
        errMsgIdAndTxt("AMICI:mex:iter","Illegal value for iter!");
        goto freturn;
    }
    ami_mem = AMICreate(udata->lmm, udata->iter);
    if (ami_mem == NULL) goto freturn;
    
    /* Initialize AMIS solver*/
    if (wrap_init(ami_mem, tdata->x, tdata->dx, udata->tstart) != AMICI_SUCCESS) goto freturn;
    
    /* Specify integration tolerances */
    if (AMISStolerances(ami_mem, RCONST(udata->rtol), RCONST(udata->atol)) != AMICI_SUCCESS) goto freturn;
    
    /* Set optional inputs */
    if (AMISetErrHandlerFn(ami_mem) != AMICI_SUCCESS) goto freturn;
    
    /* attaches userdata*/
    if (AMISetUserData(ami_mem, udata) != AMICI_SUCCESS) goto freturn;
    
    /* specify maximal number of steps */
    if (AMISetMaxNumSteps(ami_mem, udata->maxsteps) != AMICI_SUCCESS) goto freturn;
    
    /* activates stability limit detection */
    if (AMISetStabLimDet(ami_mem, udata->stldet) != AMICI_SUCCESS) goto freturn;
    
    if (udata->ne > 0) {
        if (wrap_RootInit(ami_mem, udata) != AMICI_SUCCESS) goto freturn;
    }
    
    /* Attach linear solver module */
    switch (udata->linsol) {
            
            /* DIRECT SOLVERS */
            
        case AMICI_DENSE:
            if (AMIDense(ami_mem, udata->nx) != AMICI_SUCCESS) goto freturn;
            
            if (wrap_SetDenseJacFn(ami_mem) != AMICI_SUCCESS) goto freturn;
            
            break;
            
        case AMICI_BAND:
            if (AMIBand(ami_mem, udata->nx, udata->ubw, udata->lbw) != AMICI_SUCCESS) goto freturn;
            
            if (wrap_SetBandJacFn(ami_mem) != AMICI_SUCCESS) goto freturn;
            
            break;
            
        case AMICI_LAPACKDENSE:
            errMsgIdAndTxt("AMICI:mex:lapack","Solver currently not supported!");
            /* status = CVLapackDense(ami_mem, nx);
             if (status != AMICI_SUCCESS) return;
             
             status = wrap_SetDenseJacFn(ami_mem);
             if (status != AMICI_SUCCESS) return;
             
             break;*/
            
        case AMICI_LAPACKBAND:
            
            errMsgIdAndTxt("AMICI:mex:lapack","Solver currently not supported!");
            /* status = CVLapackBand(ami_mem, nx);
             if (status != AMICI_SUCCESS) return;
             
             status = wrap_SetBandJacFn(ami_mem);
             if (status != AMICI_SUCCESS) return;
             
             break;*/
            
        case AMICI_DIAG:
            if (AMIDiag(ami_mem) != AMICI_SUCCESS) goto freturn;
            
            break;
            
            /* ITERATIVE SOLVERS */
            
        case AMICI_SPGMR:
            if (AMISpgmr(ami_mem, PREC_NONE, CVSPILS_MAXL) != AMICI_SUCCESS) goto freturn;
            
            if (wrap_SetJacTimesVecFn(ami_mem) != AMICI_SUCCESS) goto freturn;
            
            break;
            
        case AMICI_SPBCG:
            if (AMISpbcg(ami_mem, PREC_NONE, CVSPILS_MAXL) != AMICI_SUCCESS) goto freturn;
            
            if (wrap_SetJacTimesVecFn(ami_mem) != AMICI_SUCCESS) goto freturn;
            
            break;
            
        case AMICI_SPTFQMR:
            if (AMISptfqmr(ami_mem, PREC_NONE, CVSPILS_MAXL) != AMICI_SUCCESS) goto freturn;
            
            if (wrap_SetJacTimesVecFn(ami_mem) != AMICI_SUCCESS) goto freturn;
            
            break;
            
            /* SPARSE SOLVERS */
            
        case AMICI_KLU:
            if (AMIKLU(ami_mem, udata->nx, udata->nnz, CSC_MAT) != AMICI_SUCCESS) goto freturn;
            
            if (wrap_SetSparseJacFn(ami_mem) != AMICI_SUCCESS) goto freturn;
            
            if (AMIKLUSetOrdering(ami_mem, udata->ordering) != AMICI_SUCCESS) goto freturn;
            
            break;
            
        default:
            errMsgIdAndTxt("AMICI:mex:solver","Invalid choice of solver!");
            break;
    }
    
    if (udata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        if (udata->sensi_meth == AMICI_SENSI_FSA) {
            if (udata->nx>0) {
                
                /* initialise sensitivities, this can either be user provided or come from the model definition */
                realtype *sx_tmp;
                
                if (!udata->sx0data) {
                    if (fsx0(tdata->sx, tdata->x, tdata->dx, udata) != AMICI_SUCCESS) goto freturn;
                } else {
                    int ip;
                    for (ip=0; ip<udata->nplist; ip++) {
                        sx_tmp = NV_DATA_S(tdata->sx[ip]);
                        int ix;
                        for (ix=0; ix<udata->nx; ix++) {
                            sx_tmp[ix] = (realtype) udata->sx0data[ix + udata->nx*ip];
                        }
                    }
                }
                
                if (fsdx0(tdata->sdx, tdata->x, tdata->dx, udata) != AMICI_SUCCESS) goto freturn;
                
                /* Activate sensitivity calculations */
                if (wrap_SensInit1(ami_mem, tdata->sx, tdata->sdx, udata) != AMICI_SUCCESS) goto freturn;
                
                /* Set sensitivity analysis optional inputs */
                if (AMISetSensParams(ami_mem, udata->p, udata->pbar, udata->plist) != AMICI_SUCCESS) goto freturn;
                
                if (AMISetSensErrCon(ami_mem, TRUE) != AMICI_SUCCESS) goto freturn;
                
                if (AMISensEEtolerances(ami_mem) != AMICI_SUCCESS) goto freturn;
            }
        }
        
        if (udata->sensi_meth == AMICI_SENSI_ASA) {
            if (udata->nx>0) {
                /* Allocate space for the adjoint computation */
                if (AMIAdjInit(ami_mem, udata->maxsteps, udata->interpType) != AMICI_SUCCESS) goto freturn;
            }
        }
        
    }
    
    id = N_VNew_Serial(udata->nx);
    memcpy(NV_CONTENT_S(id)->data, udata->idlist, udata->nx * sizeof(realtype));
    
    if (AMISetId(ami_mem, id) != AMICI_SUCCESS) goto freturn;
    
    if (AMISetSuppressAlg(ami_mem, TRUE) != AMICI_SUCCESS) goto freturn;
    
    return(ami_mem);
    
freturn:
    if(ami_mem) AMIFree(&ami_mem);
    return NULL;
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int setupAMIB(void *ami_mem, UserData *udata, TempData *tdata) {
    /**
     * setupAMIB initialises the AMI memory object for the backwards problem
     * @param[out] status flag indicating success of execution @type int
     * @param[in] ami_mem pointer to the solver memory object of the forward problem
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[in] tdata pointer to the temporary data struct @type TempData
     * @return ami_mem pointer to the cvodes/idas memory block for the backward problem
     */
    int ix;
    int status = AMICI_SUCCESS;
    
    /* write initial conditions */
    if (tdata->xB == NULL) return AMICI_ERROR_SETUPB;
    realtype *xB_tmp = NV_DATA_S(tdata->xB);
    memset(xB_tmp,0,sizeof(realtype)*udata->nx);
    for (ix=0; ix<udata->nx; ix++) {
        xB_tmp[ix] += tdata->dJydx[udata->nt-1+ix*udata->nt];
    }
    
    if (tdata->dxB == NULL) return AMICI_ERROR_SETUPB;
    memset(NV_DATA_S(tdata->dxB),0,sizeof(realtype)*udata->nx);
    
    if (tdata->xQB == NULL) return AMICI_ERROR_SETUPB;
    memset(NV_DATA_S(tdata->xQB),0,sizeof(realtype)*udata->nJ*udata->nplist);
    
    /* create backward problem */
    if (udata->lmm>2||udata->lmm<1) {
        errMsgIdAndTxt("AMICI:mex:lmm","Illegal value for lmm!");
    }
    if (udata->iter>2||udata->iter<1) {
        errMsgIdAndTxt("AMICI:mex:iter","Illegal value for iter!");
    }
    
    /* allocate memory for the backward problem */
    status = AMICreateB(ami_mem, udata->lmm, udata->iter, &(tdata->which));
    if (status != AMICI_SUCCESS) return status;
    
    
    /* initialise states */
    status = wrap_binit(ami_mem, tdata->which, tdata->xB, tdata->dxB, tdata->t);
    if(status != AMICI_SUCCESS) return status;
    
    /* specify integration tolerances for backward problem */
    status = AMISStolerancesB(ami_mem, tdata->which, RCONST(udata->rtol), RCONST(udata->atol));
    if(status != AMICI_SUCCESS) return status;
    
    /* Attach user data */
    status = AMISetUserDataB(ami_mem, tdata->which, udata);
    if(status != AMICI_SUCCESS) return status;
    
    /* Number of maximal internal steps */
    if (AMISetMaxNumStepsB(ami_mem, tdata->which, 100*udata->maxsteps) != AMICI_SUCCESS) return AMICI_ERROR_SETUPB;
    
    switch (udata->linsol) {
            
            /* DIRECT SOLVERS */
            
        case AMICI_DENSE:
            status = AMIDenseB(ami_mem, tdata->which, udata->nx);
            if(status != AMICI_SUCCESS) return status;
            
            status = wrap_SetDenseJacFnB(ami_mem, tdata->which);
            if(status != AMICI_SUCCESS) return status;
            
            break;
            
        case AMICI_BAND:
            status = AMIBandB(ami_mem, tdata->which, udata->nx, udata->ubw, udata->lbw);
            if(status != AMICI_SUCCESS) return status;
            
            status = wrap_SetBandJacFnB(ami_mem, tdata->which);
            if(status != AMICI_SUCCESS) return status;
            
            break;
            
        case AMICI_LAPACKDENSE:
            
            /* #if SUNDIALS_BLAS_LAPACK
             status = CVLapackDenseB(ami_mem, tdata->which, nx);
             if (status != AMICI_SUCCESS) return;
             
             status = wrap_SetDenseJacFnB(ami_mem, tdata->which);
             if (status != AMICI_SUCCESS) return;
             #else*/
            errMsgIdAndTxt("AMICI:mex:lapack","Solver currently not supported!");
            /* #endif*/
            break;
            
        case AMICI_LAPACKBAND:
            
            
            /* #if SUNDIALS_BLAS_LAPACK
             status = CVLapackBandB(ami_mem, tdata->which, nx, ubw, lbw);
             if (status != AMICI_SUCCESS) return;
             
             status = wrap_SetBandJacFnB(ami_mem, tdata->which);
             if (status != AMICI_SUCCESS) return;
             #else*/
            errMsgIdAndTxt("AMICI:mex:lapack","Solver currently not supported!");
            /* #endif*/
            break;
            
        case AMICI_DIAG:
            status = AMIDiagB(ami_mem, tdata->which);
            if(status != AMICI_SUCCESS) return status;
            
            status = wrap_SetDenseJacFnB(ami_mem, tdata->which);
            if(status != AMICI_SUCCESS) return status;
            
            break;
            
            /* ITERATIVE SOLVERS */
            
        case AMICI_SPGMR:
            status = AMISpgmrB(ami_mem, tdata->which, PREC_NONE, CVSPILS_MAXL);
            if(status != AMICI_SUCCESS) return status;
            
            status = wrap_SetJacTimesVecFnB(ami_mem, tdata->which);
            if(status != AMICI_SUCCESS) return status;
            
            break;
            
        case AMICI_SPBCG:
            status = AMISpbcgB(ami_mem, tdata->which, PREC_NONE, CVSPILS_MAXL);
            if(status != AMICI_SUCCESS) return status;
            
            status = wrap_SetJacTimesVecFnB(ami_mem, tdata->which);
            if(status != AMICI_SUCCESS) return status;
            
            break;
            
        case AMICI_SPTFQMR:
            status = AMISptfqmrB(ami_mem, tdata->which, PREC_NONE, CVSPILS_MAXL);
            if(status != AMICI_SUCCESS) return status;
            
            status = wrap_SetJacTimesVecFnB(ami_mem, tdata->which);
            if(status != AMICI_SUCCESS) return status;
            
            break;
            
            /* SPARSE SOLVERS */
            
        case AMICI_KLU:
            status = AMIKLUB(ami_mem, tdata->which, udata->nx, udata->nnz, CSC_MAT);
            if(status != AMICI_SUCCESS) return status;
            
            status = wrap_SetSparseJacFnB(ami_mem, tdata->which);
            if(status != AMICI_SUCCESS) return status;
            
            status = AMIKLUSetOrderingB(ami_mem, tdata->which, udata->ordering);
            if(status != AMICI_SUCCESS) return status;
            
            break;
            
        default:
            break;
    }
    
    /* Initialise quadrature calculation */
    status = wrap_qbinit(ami_mem, tdata->which, tdata->xQB);
    if(status != AMICI_SUCCESS) return status;
    
    /* Enable Quadrature Error Control */
    status = AMISetQuadErrConB(ami_mem, tdata->which, TRUE);
    if(status != AMICI_SUCCESS) return status;
    
    status = AMIQuadSStolerancesB(ami_mem, tdata->which, RCONST(udata->rtol), RCONST(udata->atol));
    if(status != AMICI_SUCCESS) return status;
    
    status = AMISetStabLimDetB(ami_mem, tdata->which, udata->stldet);
    if(status != AMICI_SUCCESS) return status;
    
    return status;
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int getDataSensisFSA(int it, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata) {
    /**
     * getDataSensisFSA extracts data information for forward sensitivity analysis
     *
     * @param[out] status flag indicating success of execution @type int
     * @param[in] it index of current timepoint @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return void
     */
    
    int ip, iy, ix;
    int status = AMICI_SUCCESS;
    realtype *sx_tmp;
    
    for(ip=0; ip < udata->nplist; ip++) {
        if (udata->nx>0) {
            if (udata->ts[it] > udata->tstart) {
                status = AMIGetSens(ami_mem, &(tdata->t), tdata->sx);
                if (status != AMICI_SUCCESS) return status;
            }
            
            sx_tmp = NV_DATA_S(tdata->sx[ip]);
            for(ix=0; ix < udata->nx; ix++) {
                rdata->sx[(ip*udata->nx + ix)*udata->nt + it] = sx_tmp[ix];
            }
        }
    }
    
    for (iy=0; iy<udata->nytrue; iy++) {
        if (edata){
            if (amiIsNaN(edata->sigmay[iy*udata->nt+it])) {
                status = fdsigma_ydp(tdata->t,tdata->dsigmaydp,udata);
                if(status != AMICI_SUCCESS) return status;
            } else {
                for (ip=0; ip<udata->nplist; ip++) {
                    tdata->dsigmaydp[ip*udata->ny+iy] = 0;
                }
            }
            for (ip=0; ip<udata->nplist; ip++) {
                rdata->ssigmay[it + udata->nt*(ip*udata->ny+iy)] = tdata->dsigmaydp[ip*udata->ny+iy];
            }
        } else {
            for (ip=0; ip<udata->nplist; ip++) {
                rdata->ssigmay[it + udata->nt*(ip*udata->ny+iy)] = 0;
            }
        }
    }
    status = fsy(udata->ts[it],it,rdata->sy,tdata->dydx,tdata->dydp,tdata->sx,udata);
    if(status != AMICI_SUCCESS) return status;
    if (edata) {
        status = fsJy(udata->ts[it],it,rdata->sllh,rdata->s2llh,tdata->dJydy,tdata->dJydp,rdata->y,tdata->sigmay,rdata->sy,tdata->dydp,edata->my,udata);
        if(status != AMICI_SUCCESS) return status;
    }
    return status;
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int prepDataSensis(int it, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata) {
    /**
     * prepDataSensis preprocesses the provided experimental data to compute sensitivities via adjoint or forward methods later on
     *
     * @param[in] it index of current timepoint @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */
    
    int iy,ip,iJ;
    int status = AMICI_SUCCESS;
    
    status = fdydx(udata->ts[it],it,tdata->dydx,tdata->x,udata);
    if(status != AMICI_SUCCESS) return status;
    status = fdydp(udata->ts[it],it,tdata->dydp,tdata->x,udata);
    if(status != AMICI_SUCCESS) return status;
    if (edata) {
        for (iy=0; iy<udata->nytrue; iy++) {
            if (amiIsNaN(edata->sigmay[iy*udata->nt+it])) {
                status = fdsigma_ydp(tdata->t,tdata->dsigmaydp,udata);
                if(status != AMICI_SUCCESS) return status;
            } else {
                for (ip=0; ip<udata->nplist; ip++) {
                    tdata->dsigmaydp[ip*udata->ny+iy] = 0;
                }
            }
            for (ip=0; ip<udata->nplist; ip++) {
                rdata->ssigmay[it + udata->nt*(ip*udata->ny+iy)] = tdata->dsigmaydp[ip*udata->ny+iy];
            }
        }
        status = fdJydp(udata->ts[it],it,tdata->dJydp,rdata->y,tdata->x,tdata->dydp,edata->my,tdata->sigmay,tdata->dsigmaydp,udata);
        if(status != AMICI_SUCCESS) return status;
        
        
        if (udata->sensi_meth == AMICI_SENSI_ASA) {
            for(iJ=0; iJ<udata->nJ; iJ++) {
                for(ip=0; ip < udata->nplist; ip++) {
                    for(iy=0; iy < udata->nytrue; iy++) {
                        if (iJ==0) {
                            if (udata->ny>0) {
                                rdata->sllh[ip] -= tdata->dJydp[iy + ip*udata->nytrue];
                            }
                        } else {
                            if (udata->ny>0) {
                                rdata->s2llh[(iJ-1)*udata->nplist + ip] -= tdata->dJydp[(iJ*udata->nplist + ip)*udata->nytrue + iy];
                            }
                        }
                    }
                }
            }
        }
        status = fdJydy(udata->ts[it],it,tdata->dJydy,rdata->y,tdata->x,edata->my,tdata->sigmay,udata);
        if(status != AMICI_SUCCESS) return status;
        status = fdJydx(udata->ts[it],it,tdata->dJydx,rdata->y,tdata->x,tdata->dydx,edata->my,tdata->sigmay,udata);
        if(status != AMICI_SUCCESS) return status;
    }
    return status;
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int getDataOutput(int it, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata) {
    /**
     * getDataOutput extracts output information for data-points
     *
     * @param[in] it index of current timepoint @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */
    
    int iy;
    int status = AMICI_SUCCESS;
    
    status = fy(udata->ts[it],it,rdata->y,tdata->x,udata);
    if(status != AMICI_SUCCESS) return status;
    
    if (edata) {
        for (iy=0; iy<udata->nytrue; iy++) {
            /* extract the value for the standard deviation, if the data value is NaN, use
             the parameter value. Store this value in the return struct */
            if (amiIsNaN(edata->sigmay[iy*udata->nt+it])) {
                status = fsigma_y(tdata->t,tdata->sigmay,udata);
                if(status != AMICI_SUCCESS) return status;
            } else {
                tdata->sigmay[iy] = edata->sigmay[iy*udata->nt+it];
            }
            rdata->sigmay[iy*udata->nt+it] = tdata->sigmay[iy];
        }
        status = fJy(udata->ts[it],it,tdata->Jy,rdata->y,tdata->x,edata->my,tdata->sigmay,udata);
        if(status != AMICI_SUCCESS) return status;
    } else {
        status = fsigma_y(tdata->t,tdata->sigmay,udata);
        if(status != AMICI_SUCCESS) return status;
        for (iy=0; iy<udata->nytrue; iy++) {
            rdata->sigmay[iy*udata->nt+it] = tdata->sigmay[iy];
        }
    }
    if (udata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        status = prepDataSensis(it, ami_mem, udata, rdata, edata, tdata);
        if(status != AMICI_SUCCESS) return status;
        if (udata->sensi_meth == AMICI_SENSI_FSA) {
            status = getDataSensisFSA(it, ami_mem, udata, rdata, edata, tdata);
            if(status != AMICI_SUCCESS) return status;
        }
    }
    return AMICI_SUCCESS;
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int getEventSensisFSA(int ie, void *ami_mem, UserData *udata, ReturnData *rdata, TempData *tdata) {
    /**
     * getEventSensisFSA extracts event information for forward sensitivity analysis
     *
     * @param[in] ie index of event type @type int
     * @param[in] ami_mem pointer to the solver memory block @type void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */
    
    return fsz(tdata->t,ie,rdata->sz,tdata->x,tdata->sx,udata,tdata);
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int getEventSensisFSA_tf(int ie, void *ami_mem, UserData *udata, ReturnData *rdata, TempData *tdata) {
    /**
     * getEventSensisFSA_tf extracts event information for forward sensitivity
     *     analysis for events that happen at the end of the considered interval
     *
     * @param[in] ie index of event type @type int
     * @param[in] ami_mem pointer to the solver memory block @type void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */
    
    int status = AMICI_SUCCESS;
    
    status = fsz_tf(tdata->t,ie,rdata->sz,tdata->x,tdata->sx,udata,tdata);
    if(status != AMICI_SUCCESS) return status;
    
    status = fsroot(tdata->t,ie,rdata->srz,tdata->x,tdata->sx,udata,tdata);
    if(status != AMICI_SUCCESS) return status;
    
    if (udata->sensi >= AMICI_SENSI_ORDER_SECOND) {
        if (fs2root(tdata->t,ie,rdata->s2rz,tdata->x,tdata->sx,udata,tdata)!= AMICI_SUCCESS) return AMICI_ERROR_FSA;
    }
    
    return status;
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int getEventSensisASA(int ie, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata) {
    /**
     * getEventSensisASA extracts event information for adjoint sensitivity analysis
     *
     * @param[in] ie index of event type @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */
    int ip, iz;
    int status = AMICI_SUCCESS;
    
    for (iz=0; iz<udata->nztrue; iz++) {
        if ( udata->z2event[iz]-1 == ie ){
            if (!amiIsNaN(edata->mz[iz*udata->nmaxevent+tdata->nroots[ie]])) {
                status = fdzdp(tdata->t,ie,tdata->dzdp,tdata->x,udata);
                if(status != AMICI_SUCCESS) return status;
                status = fdzdx(tdata->t,ie,tdata->dzdx,tdata->x,udata);
                if(status != AMICI_SUCCESS) return status;
                /* extract the value for the standard deviation, if the data value is NaN, use
                 the parameter value. Store this value in the return struct */
                if (amiIsNaN(edata->sigmaz[tdata->nroots[ie] + udata->nmaxevent*iz])) {
                    status = fsigma_z(tdata->t,ie,tdata->sigmaz,udata);
                    if(status != AMICI_SUCCESS) return status;
                    status = fdsigma_zdp(tdata->t,ie,tdata->dsigmazdp,udata);
                    if(status != AMICI_SUCCESS) return status;
                } else {
                    for (ip=0; ip<udata->nplist; ip++) {
                        tdata->dsigmazdp[iz+udata->nz*ip] = 0;
                    }
                    tdata->sigmaz[iz] = edata->sigmaz[tdata->nroots[ie] + udata->nmaxevent*iz];
                }
                rdata->sigmaz[tdata->nroots[ie] + udata->nmaxevent*iz] = tdata->sigmaz[iz];
                for (ip=0; ip<udata->nplist; ip++) {
                    rdata->ssigmaz[tdata->nroots[ie] + udata->nmaxevent*(iz+udata->nz*ip)] = tdata->dsigmazdp[iz+udata->nz*ip];
                }
                
                status = fdJzdp(tdata->t,ie,tdata->dJzdp,rdata->z,tdata->x,tdata->dzdp,edata->mz,tdata->sigmaz,tdata->dsigmazdp,udata,tdata);
                if(status != AMICI_SUCCESS) return status;
                status = fdJzdx(tdata->t,ie,tdata->dJzdx,rdata->z,tdata->x,tdata->dzdx,edata->mz,tdata->sigmaz,udata,tdata);
                if(status != AMICI_SUCCESS) return status;
            }
        }
    }
    return AMICI_SUCCESS;
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int getEventSigma(int ie, int iz, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata) {
    /**
     * getEventSigma extracts fills sigmaz either from the user defined function or from user input
     *
     * @param[in] ie event type index @type int
     * @param[in] iz event output index @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */
    
    int status = AMICI_SUCCESS;
    
    /* extract the value for the standard deviation, if the data value is NaN, use
     the parameter value. Store this value in the return struct */
    if (amiIsNaN(edata->sigmaz[tdata->nroots[ie] + udata->nmaxevent*iz])) {
        status = fsigma_z(tdata->t,ie,tdata->sigmaz,udata);
        if(status != AMICI_SUCCESS) return status;
    } else {
        tdata->sigmaz[iz] = edata->sigmaz[tdata->nroots[ie] + udata->nmaxevent*iz];
    }
    rdata->sigmaz[tdata->nroots[ie] + udata->nmaxevent*iz] = tdata->sigmaz[iz];
    
    return status;
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int getEventObjective(int ie, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata) {
    /**
     * getEventObjective updates the objective function on the occurence of an event
     *
     * @param[in] ie event type index @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */
    
    int status = AMICI_SUCCESS;
    
    if (edata) {
        int iz;
        for (iz=0; iz<udata->nztrue; iz++) {
            if (udata->z2event[iz]-1 == ie) {
                status = getEventSigma(ie, iz, ami_mem, udata, rdata, edata, tdata);
                if(status != AMICI_SUCCESS) return status;
                if (!amiIsNaN(edata->mz[iz*udata->nmaxevent+tdata->nroots[ie]])) {
                    tdata->Jz[0] += 0.5*log(2*pi*pow(tdata->sigmaz[tdata->nroots[ie] + udata->nmaxevent*iz],2))
                    + 0.5*pow( ( rdata->z[tdata->nroots[ie] + udata->nmaxevent*iz] - edata->mz[tdata->nroots[ie] + udata->nmaxevent*iz] )/tdata->sigmaz[iz] , 2);
                    *rdata->chi2 += pow( ( rdata->z[tdata->nroots[ie] + udata->nmaxevent*iz] - edata->mz[tdata->nroots[ie] + udata->nmaxevent*iz] )/tdata->sigmaz[iz] , 2);
                }
            }
        }
    }
    return status;
    
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int getEventOutput(realtype *tlastroot, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata) {
    /**
     * getEventOutput extracts output information for events
     *
     * @param[in] tlastroot timepoint of last occured event @type *realtype
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */
    int iz, ie;
    int status = AMICI_SUCCESS;
    
    /* EVENT OUTPUT */
    for (ie=0; ie<udata->ne; ie++){ /* only look for roots of the rootfunction not discontinuities */
        if (tdata->nroots[ie]<udata->nmaxevent) {
            if (tdata->rootsfound[ie] == 1) { /* only consider transitions false -> true */
                status = fz(tdata->t,ie,rdata->z,tdata->x,udata,tdata);
                if(status != AMICI_SUCCESS) return status;
                if (udata->sensi >= AMICI_SENSI_ORDER_FIRST) {
                    if (udata->sensi_meth == AMICI_SENSI_ASA) {
                        status = getEventSensisASA(ie, ami_mem, udata, rdata, edata, tdata);
                        if(status != AMICI_SUCCESS) return status;
                    } else {
                        status = getEventSensisFSA(ie, ami_mem, udata, rdata, tdata);
                        if(status != AMICI_SUCCESS) return status;
                    }
                }
                
                if (edata) {
                    for (iz=0; iz<udata->nztrue; iz++) {
                        if (udata->z2event[iz]-1 == ie) {
                            status = getEventSigma(ie, iz, ami_mem,udata,rdata,edata,tdata);
                            if(status != AMICI_SUCCESS) return status;
                        }
                    }
                    status = getEventObjective(ie, ami_mem, udata, rdata, edata, tdata);
                    if(status != AMICI_SUCCESS) return status;
                }
                
                tdata->nroots[ie]++;
            }
        }
    }
    return status;
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int fillEventOutput(void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata) {
    /**
     * fillEventOutput fills missing roots at last timepoint
     *
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */
    
    int ie,iz;
    int status = AMICI_SUCCESS;
    
    status = froot(tdata->t,tdata->x,tdata->dx,tdata->rootvals,udata);
    if(status != AMICI_SUCCESS) return status;
    
    /* EVENT OUTPUT */
    if (udata->nztrue>0) {
        for (ie=0; ie<udata->ne; ie++){ /* only look for roots of the rootfunction not discontinuities */
            while (tdata->nroots[ie]<udata->nmaxevent) {
                status = fz(tdata->t,ie,rdata->z,tdata->x,udata,tdata);
                if(status != AMICI_SUCCESS) return status;
                
                
                for (iz=0; iz<udata->nztrue; iz++) {
                    if (udata->z2event[iz]-1 == ie) {
                        rdata->rz[tdata->nroots[ie] + udata->nmaxevent*iz] = tdata->rootvals[ie];
                    }
                }
                
                status = getEventObjective(ie, ami_mem, udata, rdata, edata, tdata);
                if(status != AMICI_SUCCESS) return status;
                
                if (udata->sensi >= AMICI_SENSI_ORDER_FIRST) {
                    if (udata->sensi_meth == AMICI_SENSI_ASA) {
                        status = getEventSensisASA(ie, ami_mem, udata, rdata, edata, tdata);
                        if(status != AMICI_SUCCESS) return status;
                    } else {
                        status = getEventSensisFSA_tf(ie, ami_mem, udata, rdata, tdata);
                        if(status != AMICI_SUCCESS) return status;
                    }
                }
                tdata->nroots[ie]++;
            }
        }
    }
    return status;
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int handleDataPoint(int it, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata) {
    /**
     * handleDataPoint executes everything necessary for the handling of data points
     *
     * @param[in] it index of data point @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */
    
    int ix;
    int status = AMICI_SUCCESS;
    
    rdata->ts[it] = udata->ts[it];
    if (udata->nx>0) {
        realtype *x_tmp = NV_DATA_S(tdata->x);
        for (ix=0; ix<udata->nx; ix++) {
            rdata->x[it+udata->nt*ix] = x_tmp[ix];
        }
        
        if (it == udata->nt-1) {
            if (udata->sensi_meth == AMICI_SENSI_SS) {
                status = fdxdotdp(udata->ts[it],rdata->dxdotdp,tdata->x,tdata->dx,udata);
                if(status != AMICI_SUCCESS) return status;
                status = fdydp(udata->ts[it],it,rdata->dydp,tdata->x,udata);
                if(status != AMICI_SUCCESS) return status;
                status = fdydx(udata->ts[it],it,rdata->dydx,tdata->x,udata);
                if(status != AMICI_SUCCESS) return status;
            }
        }
        
        if (udata->ts[it] > udata->tstart) {
            status = getDiagnosis(it, ami_mem, udata, rdata);
            if(status != AMICI_SUCCESS) return status;
        }
    }
    
    return getDataOutput(it, ami_mem, udata, rdata, edata, tdata);
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int handleDataPointB(int it, void *ami_mem, UserData *udata, ReturnData *rdata, TempData *tdata) {
    /**
     * handleDataPoint executes everything necessary for the handling of data points for the backward problems
     *
     * @param[in] it index of data point @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */
    
    int ix;
    
    realtype *xB_tmp = NV_DATA_S(tdata->xB);
    for (ix=0; ix<udata->nx; ix++) {
        xB_tmp[ix] += tdata->dJydx[it+ix*udata->nt];
    }
    return getDiagnosisB(it,ami_mem,udata,rdata,tdata);
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int handleEvent(int *iroot, realtype *tlastroot, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata, int seflag) {
    /**
     * handleEvent executes everything necessary for the handling of events
     *
     * @param[out] iroot index of event @type int
     * @param[out] tlastroot pointer to the timepoint of the last event @type *realtype
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @param[in] seflag flag indicating whether this is a secondary event @type int
     * @return status flag indicating success of execution @type int
     */
    int ie;
    int secondevent = 0;
    int status = AMICI_SUCCESS;
    
    
    /* store heaviside information at event occurence */
    if(froot(tdata->t,tdata->x,tdata->dx,tdata->rootvals,udata) != AMICI_SUCCESS) return AMICI_ERROR_EVENT;
    
    if (seflag == 0) {
        status = AMIGetRootInfo(ami_mem, tdata->rootsfound);
        if(status != AMICI_SUCCESS) return status;
    }
    
    if (*iroot<udata->nmaxevent*udata->ne) {
        for (ie=0; ie<udata->ne; ie++) {
            tdata->rootidx[*iroot*udata->ne + ie] = tdata->rootsfound[ie];
        }
    }
    for (ie = 0; ie<udata->ne; ie++) {
        tdata->h[ie] = tdata->rootvals[ie];
    }
    
    /* only extract in the first event fired */
    if (seflag == 0) {
        if (udata->sensi >= AMICI_SENSI_ORDER_FIRST){
            if (udata->sensi_meth == AMICI_SENSI_FSA) {
                if (AMIGetSens(ami_mem, &(tdata->t), tdata->sx) != AMICI_SUCCESS) return AMICI_ERROR_SA;
            }
        }
    }
    
    /* only check this in the first event fired, otherwise this will always be true */
    if (seflag == 0) {
        if (tdata->t == *tlastroot) {
            warnMsgIdAndTxt("AMICI:mex:STUCK_EVENT","AMICI is stuck in an event, as the initial step-size after the event is too small. To fix this, increase absolute and relative tolerances!");
            return AMICI_ERROR_EVENT;
        }
        *tlastroot = tdata->t;
    }
    
    status = getEventOutput(tlastroot, ami_mem, udata, rdata, edata, tdata);
    if (status != AMICI_SUCCESS) return status;
    
    /* if we need to do forward sensitivities later on we need to store the old x and the old xdot */
    if (udata->sensi >= AMICI_SENSI_ORDER_FIRST){
        /* store x and xdot to compute jump in sensitivities */
        N_VScale(1.0,tdata->x,tdata->x_old);
        if (udata->sensi_meth == AMICI_SENSI_FSA) {
            status = fxdot(tdata->t,tdata->x,tdata->dx,tdata->xdot,udata);
            N_VScale(1.0,tdata->xdot,tdata->xdot_old);
            N_VScale(1.0,tdata->dx,tdata->dx_old);
            
            /* compute event-time derivative only for primary events, we get into trouble with multiple simultaneously firing events here (but is this really well defined then?), in that case just use the last ie and hope for the best. */
            if (seflag == 0) {
                for (ie = 0; ie<udata->ne; ie++) {
                    if (tdata->rootsfound[ie] == 1) { /* only consider transitions false -> true */
                        fstau(tdata->t,ie,udata->stau,tdata->x,tdata->sx,udata);
                    }
                }
            }
        }
        
        if (udata->sensi_meth == AMICI_SENSI_ASA) {
            /* store x to compute jump in discontinuity */
            if (*iroot<udata->nmaxevent*udata->ne) {
                N_VScale(1.0,tdata->x,tdata->x_disc[*iroot]);
                N_VScale(1.0,tdata->xdot,tdata->xdot_disc[*iroot]);
                N_VScale(1.0,tdata->xdot_old,tdata->xdot_old_disc[*iroot]);
            }
        }
    }
    
    status = updateHeaviside(udata, tdata);
    if (status != AMICI_SUCCESS) return status;
    
    status = applyEventBolus(ami_mem, udata, tdata);
    if (status != AMICI_SUCCESS) return status;
    
    if (*iroot<udata->nmaxevent*udata->ne) {
        tdata->discs[*iroot] = tdata->t;
        (*iroot)++;
    } else {
        warnMsgIdAndTxt("AMICI:mex:TOO_MUCH_EVENT","Event was recorded but not reported as the number of occured events exceeded (nmaxevents)*(number of events in model definition)!");
        status = AMIReInit(ami_mem, tdata->t, tdata->x, tdata->dx); /* reinitialise so that we can continue in peace */
        return status;
    }
    
    if (udata->sensi >= AMICI_SENSI_ORDER_FIRST){
        if (udata->sensi_meth == AMICI_SENSI_FSA) {
            
            /* compute the new xdot  */
            status = fxdot(tdata->t,tdata->x,tdata->dx,tdata->xdot,udata);
            if (status != AMICI_SUCCESS) return status;
            
            status = applyEventSensiBolusFSA(ami_mem, udata, tdata);
            if (status != AMICI_SUCCESS) return status;
        }
    }
    
    /* check whether we need to fire a secondary event */
    status = froot(tdata->t,tdata->x,tdata->dx,tdata->rootvals,udata);
    if (status != AMICI_SUCCESS) return status;
    for (ie = 0; ie<udata->ne; ie++) {
        /* the same event should not trigger itself */
        if (tdata->rootsfound[ie] == 0 ) {
            /* check whether there was a zero-crossing */
            if ( 0 > tdata->h[ie]*tdata->rootvals[ie]) {
                if (tdata->h[ie]<tdata->rootvals[ie]) {
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
    if (secondevent>0) {
        status = handleEvent( iroot, tlastroot, ami_mem, udata, rdata, edata, tdata, secondevent);
        if (status != AMICI_SUCCESS) return status;
    }
    
    /* only reinitialise in the first event fired */
    if (seflag == 0) {
        status = AMIReInit(ami_mem, tdata->t, tdata->x, tdata->dx);
        if (status != AMICI_SUCCESS) return status;
        
        /* make time derivative consistent */
        status = AMICalcIC(ami_mem, tdata->t);
        if (status != AMICI_SUCCESS) return status;
    }
    
    if (udata->sensi >= AMICI_SENSI_ORDER_FIRST){
        if (udata->sensi_meth == AMICI_SENSI_FSA) {
            if (seflag == 0) {
                status = AMISensReInit(ami_mem, udata->ism, tdata->sx, tdata->sdx);
                if (status != AMICI_SUCCESS) return status;
            }
        }
    }
    return status;
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int handleEventB(int iroot, void *ami_mem, UserData *udata, TempData *tdata) {
    /**
     * handleEventB executes everything necessary for the handling of events for the backward problem
     *
     * @param[out] status flag indicating success of execution @type int
     * @param[out] iroot index of event @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */
    
    int ie, ix, ip, iJ;
    int status = AMICI_SUCCESS;
    
    
    /* store current values */
    N_VScale(1.0,tdata->xB,tdata->xB_old);
    N_VScale(1.0,tdata->xQB,tdata->xQB_old);
    
    realtype *xB_tmp = NV_DATA_S(tdata->xB);
    realtype *xQB_tmp = NV_DATA_S(tdata->xQB);
    
    for (ie=0; ie<udata->ne; ie++) {
        
        if (tdata->rootidx[iroot*udata->ne + ie] != 0) {
            
            status = fdeltaqB(tdata->t,ie,tdata->deltaqB,tdata->x_disc[iroot],tdata->xB_old,tdata->xQB_old,tdata->xdot_disc[iroot],tdata->xdot_old_disc[iroot],udata);
            if (status != AMICI_SUCCESS) return status;
            status = fdeltaxB(tdata->t,ie,tdata->deltaxB,tdata->x_disc[iroot],tdata->xB_old,tdata->xdot_disc[iroot],tdata->xdot_old_disc[iroot],udata);
            if (status != AMICI_SUCCESS) return status;
            
            for (ix=0; ix<udata->nx; ix++) {
                xB_tmp[ix] += tdata->deltaxB[ix];
                if (udata->nz>0) {
                    xB_tmp[ix] += tdata->dJzdx[tdata->nroots[ie] + udata->nmaxevent*ix];
                }
            }
            
            for (iJ=0; iJ<udata->nJ; iJ++) {
                for (ip=0; ip<udata->nplist; ip++) {
                    xQB_tmp[iJ*udata->nplist+ip] += tdata->deltaqB[iJ*udata->nplist+ip];
                }
            }
            
            
            tdata->nroots[ie]--;
        }
    }
    
    return updateHeavisideB(iroot, udata, tdata);
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

realtype getTnext(realtype *troot, int iroot, realtype *tdata, int it, UserData *udata) {
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
            if (udata->ne>0) {
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

int applyEventBolus( void *ami_mem, UserData *udata, TempData *tdata) {
    /**
     * applyEventBolus applies the event bolus to the current state
     *
     * @param[out] status flag indicating success of execution @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */
    
    int ix, ie;
    int status = AMICI_SUCCESS;
    realtype *x_tmp;
    
    for (ie=0; ie<udata->ne; ie++){
        if (tdata->rootsfound[ie] == 1) { /* only consider transitions false -> true */
            status = fdeltax(tdata->t,ie,tdata->deltax,tdata->x,tdata->xdot,tdata->xdot_old,udata);
            if (status != AMICI_SUCCESS) return status;
            
            x_tmp = NV_DATA_S(tdata->x);
            for (ix=0; ix<udata->nx; ix++) {
                x_tmp[ix] += tdata->deltax[ix];
            }
        }
    }
    return status;
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int applyEventSensiBolusFSA(void *ami_mem, UserData *udata, TempData *tdata) {
    /**
     * applyEventSensiBolusFSA applies the event bolus to the current sensitivities
     *
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */
    
    int ix, ip, ie;
    int status = AMICI_SUCCESS;
    realtype *sx_tmp;
    
    for (ie=0; ie<udata->ne; ie++){
        if (tdata->rootsfound[ie] == 1) { /* only consider transitions false -> true */
            status = fdeltasx(tdata->t,ie,tdata->deltasx,tdata->x_old,tdata->xdot,tdata->xdot_old,tdata->sx,udata);
            if (status != AMICI_SUCCESS) return status;
            
            for (ip=0; ip<udata->nplist; ip++) {
                sx_tmp = NV_DATA_S(tdata->sx[ip]);
                for (ix=0; ix<udata->nx; ix++) {
                    sx_tmp[ix] += tdata->deltasx[ix + udata->nx*ip];
                }
            }
        }
    }
    return status;
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int initHeaviside(UserData *udata, TempData *tdata) {
    /**
     * initHeaviside initialises the heaviside variables h at the intial time t0
     * heaviside variables activate/deactivate on event occurences
     *
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */
    
    int ie;
    int status = AMICI_SUCCESS;
    
    status = froot(tdata->t,tdata->x,tdata->dx,tdata->rootvals,udata);
    if (status != AMICI_SUCCESS) return status;
    
    for (ie = 0; ie<udata->ne; ie++) {
        if (tdata->rootvals[ie]<0) {
            udata->h[ie] = 0.0;
        } else if (tdata->rootvals[ie]==0) {
            errMsgIdAndTxt("AMICI:mex:initHeaviside","Simulation started in an event. This could lead to unexpected results, aborting simulation! Please specify an earlier simulation start via @amimodel.t0");
            return AMICI_ERROR_EVENT;
        } else {
            udata->h[ie] = 1.0;
        }
    }
    return status;
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int updateHeaviside(UserData *udata, TempData *tdata) {
    /**
     * updateHeaviside updates the heaviside variables h on event occurences
     *
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status = status flag indicating success of execution @type int;
     */
    
    int ie;

    /* tdata->rootsfound provides the direction of the zero-crossing, so adding it will give
     the right update to the heaviside variables */
    
    for (ie = 0; ie<udata->ne; ie++) {
        udata->h[ie] += tdata->rootsfound[ie];
    }
    return AMICI_SUCCESS;
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int updateHeavisideB(int iroot, UserData *udata, TempData *tdata) {
    /**
     * updateHeavisideB updates the heaviside variables h on event occurences for the backward problem
     *
     * @param[in] iroot discontinuity occurance index @type int
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */
    
    int ie;
    
    /* tdata->rootsfound provides the direction of the zero-crossing, so adding it will give
     the right update to the heaviside variables */
    
    for (ie = 0; ie<udata->ne; ie++) {
        udata->h[ie] -= tdata->rootidx[iroot*udata->ne + ie];
    }
    return AMICI_SUCCESS;
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int getDiagnosis(int it, void *ami_mem, UserData *udata, ReturnData *rdata) {
    /**
     * getDiagnosis extracts diagnosis information from solver memory block and writes them into the return data struct
     *
     * @param[out]
     * @param[in] it time-point index @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @return status flag indicating success of execution @type int
     */
    long int number;
    int status = AMICI_SUCCESS;
    int order;
    
    
    status = AMIGetNumSteps(ami_mem, &number);
    if (status != AMICI_SUCCESS) return status;
    rdata->numsteps[it] = (double) number;
    
    status = AMIGetNumRhsEvals(ami_mem, &number);
    if (status != AMICI_SUCCESS) return status;
    rdata->numrhsevals[it] = (double) number;
    
    status = AMIGetNumErrTestFails(ami_mem, &number);
    if (status != AMICI_SUCCESS) return status;
    rdata->numerrtestfails[it] = (double) number;
    
    status = AMIGetNumNonlinSolvConvFails(ami_mem, &number);
    if (status != AMICI_SUCCESS) return status;
    rdata->numnonlinsolvconvfails[it] = (double) number;
    
    status = AMIGetLastOrder(ami_mem, &order);
    if (status != AMICI_SUCCESS) return status;
    rdata->order[it] = (double) order;
    
    return status;
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

int getDiagnosisB(int it, void *ami_mem, UserData *udata, ReturnData *rdata, TempData *tdata) {
    /**
     * getDiagnosisB extracts diagnosis information from solver memory block and writes them into the return data struct for the backward problem
     *
     * @param[in] it time-point index @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */
    long int number;
    int status = AMICI_SUCCESS;
    
    void *ami_memB;
    
    ami_memB = AMIGetAdjBmem(ami_mem, tdata->which);
    
    status = AMIGetNumSteps(ami_memB, &number);
    if (status != AMICI_SUCCESS) return status;
    rdata->numstepsB[it] = (double) number;
    
    status = AMIGetNumRhsEvals(ami_memB, &number);
    if (status != AMICI_SUCCESS) return status;
    rdata->numrhsevalsB[it] = (double) number;
    
    status = AMIGetNumErrTestFails(ami_memB, &number);
    if (status != AMICI_SUCCESS) return status;
    rdata->numerrtestfailsB[it] = (double) number;
    
    status = AMIGetNumNonlinSolvConvFails(ami_memB, &number);
    if (status != AMICI_SUCCESS) return status;
    rdata->numnonlinsolvconvfailsB[it] = (double) number;
    
    return status;
}


int workForwardProblem(UserData *udata, TempData *tdata, ReturnData *rdata, const ExpData *edata, void *ami_mem, int *iroot) {
    /**
     * workForwardProblem solves the forward problem. if forward sensitivities are enabled this will also compute sensitivies
     *
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[in] tdata pointer to the temporary data struct @type TempData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[out] edata pointer to the experimental data struct @type ExpData
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] iroot pointer to the current root index, the value pointed to will be increased during the forward solve
     * @return int status flag indicating success of execution @type int
     */
    
    
    /*******************/
    /* FORWARD PROBLEM */
    /*******************/
    int ix, it;
    int ncheck = 0; /* the number of (internal) checkpoints stored so far */
    realtype *x_tmp;
    realtype tlastroot = 0; /* storage for last found root */
    int status = AMICI_SUCCESS;
    
    /* loop over timepoints */
    for (it=0; it < udata->nt; it++) {
        if (udata->sensi_meth == AMICI_SENSI_FSA && udata->sensi >= AMICI_SENSI_ORDER_FIRST) {
            status = AMISetStopTime(ami_mem, udata->ts[it]);
        }
        if (status == AMICI_SUCCESS) {
            /* only integrate if no errors occured */
            if (udata->ts[it] > udata->tstart) {
                while (tdata->t<udata->ts[it]) {
                    if (udata->sensi_meth == AMICI_SENSI_ASA && udata->sensi >= AMICI_SENSI_ORDER_FIRST) {
                        if (udata->nx>0) {
                            status = AMISolveF(ami_mem, RCONST(udata->ts[it]), tdata->x, tdata->dx, &(tdata->t), AMICI_NORMAL, &ncheck);
                        } else {
                            tdata->t = udata->ts[it];
                        }
                    } else {
                        if (udata->nx>0) {
                            status = AMISolve(ami_mem, RCONST(udata->ts[it]), tdata->x, tdata->dx, &(tdata->t), AMICI_NORMAL);
                        } else {
                            tdata->t = udata->ts[it];
                        }
                    }
                    if (udata->nx>0) {
                        x_tmp = NV_DATA_S(tdata->x);
                        if (status == -22) {
                            /* clustering of roots => turn off rootfinding */
                            AMIRootInit(ami_mem, 0, NULL);
                            status = AMICI_SUCCESS;
                        }
                        /* integration error occured */
                        if (status<AMICI_SUCCESS) {
                            return status;
                        }
                        if (status==AMICI_ROOT_RETURN) {
                            status = handleEvent(iroot, &tlastroot, ami_mem, udata, rdata, edata, tdata, 0);
                            if (status != AMICI_SUCCESS) return status;
                        }
                    }
                }
            }
            status = handleDataPoint(it, ami_mem, udata, rdata, edata, tdata);
            if (status != AMICI_SUCCESS) return status;
        } else {
            for(ix=0; ix < udata->nx; ix++) rdata->x[ix*udata->nt+it] = amiGetNaN();
        }
    }
    
    /* fill events */
    if (udata->ne>0) {
        fillEventOutput(ami_mem, udata, rdata, edata, tdata);
    }
    
    storeJacobianAndDerivativeInReturnData(udata, tdata, rdata);
    return AMICI_SUCCESS;
}

int workBackwardProblem(UserData *udata, TempData *tdata, ReturnData *rdata, const ExpData *edata, void *ami_mem, int *iroot) {
    /**
     * workBackwardProblem solves the backward problem. if adjoint sensitivities are enabled this will also compute sensitivies
     * workForwardProblem should be called before this is function is called
     *
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[in] tdata pointer to the temporary data struct @type TempData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[out] edata pointer to the experimental data struct @type ExpData
     * @param[out] status flag indicating success of execution @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] iroot pointer to the current root index, the value pointed to will be decreased during the forward solve
     * @return int status flag
     */
    int ix, it, ip;
    int status = AMICI_SUCCESS;
    double tnext;
    
    if (udata->nx>0) {
        if (udata->sensi >= AMICI_SENSI_ORDER_FIRST) {
            if (udata->sensi_meth == AMICI_SENSI_ASA) {
                if (status == AMICI_SUCCESS) {
                    setupAMIB(ami_mem, udata, tdata);
                    
                    it = udata->nt-2;
                    (*iroot)--;
                    while (it>=0 || *iroot>=0) {
                        
                        /* check if next timepoint is a discontinuity or a data-point */
                        tnext = getTnext(tdata->discs, *iroot, udata->ts, it, udata);
                        
                        if (tnext<tdata->t) {
                            status = AMISolveB(ami_mem, tnext, AMICI_NORMAL);
                            if (status != AMICI_SUCCESS) return status;
                            
                            
                            status = AMIGetB(ami_mem, tdata->which, &(tdata->t), tdata->xB, tdata->dxB);
                            if (status != AMICI_SUCCESS) return status;
                            status = AMIGetQuadB(ami_mem, tdata->which, &(tdata->t), tdata->xQB);
                            if (status != AMICI_SUCCESS) return status;
                        }
                        
                        /* handle discontinuity */
                        
                        if (udata->ne>0){
                            if (udata->nmaxevent>0){
                                if ((*iroot)>=0){
                                    if (tnext == tdata->discs[*iroot]) {
                                        handleEventB(*iroot, ami_mem, udata, tdata);
                                        (*iroot)--;
                                    }
                                }
                            }
                        }
                        
                        /* handle data-point */
                        
                        if (tnext == udata->ts[it]) {
                            handleDataPointB(it, ami_mem, udata, rdata, tdata);
                            it--;
                        }
                        
                        /* reinit states */
                        status = AMIReInitB(ami_mem, tdata->which, tdata->t, tdata->xB, tdata->dxB);
                        if (status != AMICI_SUCCESS) return status;
                        
                        status = AMIQuadReInitB(ami_mem, tdata->which, tdata->xQB);
                        if (status != AMICI_SUCCESS) return status;
                        
                        status = AMICalcICB(ami_mem, tdata->which, tdata->t, tdata->xB, tdata->dxB);
                        if (status != AMICI_SUCCESS) return status;
                    }
                    
                    /* we still need to integrate from first datapoint to tstart */
                    if (tdata->t>udata->tstart) {
                        if (status == AMICI_SUCCESS) {
                            if (udata->nx>0) {
                                /* solve for backward problems */
                                status = AMISolveB(ami_mem, udata->tstart, AMICI_NORMAL);
                                if (status != AMICI_SUCCESS) return status;
                                
                                status = AMIGetQuadB(ami_mem, tdata->which, &(tdata->t), tdata->xQB);
                                if (status != AMICI_SUCCESS) return status;
                                status = AMIGetB(ami_mem, tdata->which, &(tdata->t), tdata->xB, tdata->dxB);
                                if (status != AMICI_SUCCESS) return status;
                            }
                        }
                    }
                    
                    status = fx0(tdata->x,udata);
                    if (status != AMICI_SUCCESS) return status;
                    status = fdx0(tdata->x,tdata->dx,udata);
                    if (status != AMICI_SUCCESS) return status;
                    status = fsx0(tdata->sx, tdata->x, tdata->dx, udata);
                    if (status != AMICI_SUCCESS) return status;
                    
                    if (status == AMICI_SUCCESS) {
                        
                        realtype *xB_tmp = NV_DATA_S(tdata->xB);
                        realtype *sx_tmp;
                        
                        int iJ;
                        for (iJ=0; iJ<udata->nJ; iJ++) {
                            if (iJ==0) {
                                for (ip=0; ip<udata->nplist; ip++) {
                                    tdata->llhS0[iJ*udata->nplist + ip] = 0.0;
                                    sx_tmp = NV_DATA_S(tdata->sx[ip]);
                                    for (ix = 0; ix < udata->nxtrue; ix++) {
                                        tdata->llhS0[ip] = tdata->llhS0[ip] + xB_tmp[ix] * sx_tmp[ix];
                                    }
                                }
                            } else {
                                for (ip=0; ip<udata->nplist; ip++) {
                                    tdata->llhS0[iJ*udata->nplist + ip] = 0.0;
                                    sx_tmp = NV_DATA_S(tdata->sx[ip]);
                                    for (ix = 0; ix < udata->nxtrue; ix++) {
                                        tdata->llhS0[iJ*udata->nplist + ip] = tdata->llhS0[iJ*udata->nplist + ip]
                                        + xB_tmp[iJ*udata->nxtrue + ix] * sx_tmp[ix]
                                        + xB_tmp[ix] * sx_tmp[iJ*udata->nxtrue + ix];
                                    }
                                }
                            }
                        }
                        
                        realtype *xQB_tmp = NV_DATA_S(tdata->xQB);
                        
                        for(iJ=0; iJ<udata->nJ; iJ++) {
                            for(ip=0; ip < udata->nplist; ip++) {
                                if (iJ==0) {
                                    rdata->sllh[ip] -=  tdata->llhS0[ip] + xQB_tmp[ip];
                                    if (udata->nz>0) {
                                        rdata->sllh[ip] -= tdata->dJzdp[ip];
                                    }
                                } else {
                                    rdata->s2llh[(iJ-1)*udata->nplist + ip] -= tdata->llhS0[iJ*udata->nplist + ip] + xQB_tmp[iJ*udata->nplist + ip];
                                    if (udata->nz>0) {
                                        rdata->s2llh[(iJ-1)*udata->nplist + ip] -= tdata->dJzdp[iJ*udata->nplist + ip];
                                    }
                                }
                            }
                        }
                        
                    } else {
                        int iJ;
                        for(iJ=0; iJ<udata->nJ; iJ++) {
                            for(ip=0; ip < udata->nplist; ip++) {
                                if (iJ==0) {
                                    rdata->sllh[ip] = amiGetNaN();
                                } else {
                                    rdata->s2llh[(iJ-1)*udata->nplist + ip] = amiGetNaN();
                                }
                            }
                        }
                    }
                } else {
                    int iJ;
                    for(iJ=0; iJ<udata->nJ; iJ++) {
                        for(ip=0; ip < udata->nplist; ip++) {
                            if (iJ==0) {
                                rdata->sllh[ip] = amiGetNaN();
                            } else {
                                rdata->s2llh[(iJ-1)*udata->nplist + ip] = amiGetNaN();
                            }
                        }
                    }
                }
            }
        }
    }
    
    /* evaluate likelihood */
    if (edata) {
        *rdata->llh = - tdata->Jy[0] - tdata->Jz[0];
    } else {
        *rdata->llh = amiGetNaN();
    }
    
    return AMICI_SUCCESS;
}

int storeJacobianAndDerivativeInReturnData(UserData *udata, TempData *tdata,  ReturnData *rdata) {
    /**
     * evalues the Jacobian and differential equation right hand side, stores it in tdata and
     and copys it to rdata
     *
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @return void
     */
    
    int status = AMICI_SUCCESS;
    
    /* store current Jacobian and derivative */
    if (udata) {
        if (tdata) {
            if (udata->nx>0){
                status = fxdot(tdata->t,tdata->x,tdata->dx,tdata->xdot,udata);
                if (status != AMICI_SUCCESS) return status;
                realtype *xdot_tmp = NV_DATA_S(tdata->xdot);
                if (rdata->xdot)
                if (xdot_tmp)
                memcpy(rdata->xdot,xdot_tmp,udata->nx*sizeof(realtype));
            }
        }
    }
    if (udata) {
        if (udata->nx>0) {
            status = fJ(udata->nx,tdata->t,0,tdata->x,tdata->dx,tdata->xdot,tdata->Jtmp,udata,NULL,NULL,NULL);
            if (status != AMICI_SUCCESS) return status;
            if (rdata->J)
            if (tdata->Jtmp->data)
            memcpy(rdata->J,tdata->Jtmp->data,udata->nx*udata->nx*sizeof(realtype));
        }
    }
    return AMICI_SUCCESS;
}

int unscaleParameters(UserData *udata) {
    switch(udata->pscale) {
        case AMICI_SCALING_LOG10:
            for(int ip = 0; ip < udata->np; ++ip) {
                udata->p[ip] = pow(10, udata->p[ip]);
            }
            break;
        case AMICI_SCALING_LN:
            for(int ip = 0; ip < udata->np; ++ip)
            udata->p[ip] = exp(udata->p[ip]);
            break;
        case AMICI_SCALING_NONE:
            //this should never be reached
            break;
    }
    return AMICI_SUCCESS;
}

int applyChainRuleFactorToSimulationResults(const UserData *udata, ReturnData *rdata, const ExpData *edata)
{
    if (udata->pscale == AMICI_SCALING_NONE)
    return AMICI_SUCCESS;
    
    // chain-rule factor: multiplier for am_p
    realtype coefficient;
    realtype *pcoefficient, *augcoefficient;
    
    pcoefficient = new realtype[udata->nplist]();
    augcoefficient = new realtype[udata->np]();
    
    switch(udata->pscale) {
        case AMICI_SCALING_LOG10:
            coefficient = log(10.0);
            for(int ip = 0; ip < udata->nplist; ++ip)
            pcoefficient[ip] = udata->p[udata->plist[ip]]*log(10);
            if (udata->sensi == 2)
            if (udata->o2mode == AMICI_O2MODE_FULL)
            for(int ip = 0; ip < udata->np; ++ip)
            augcoefficient[ip] = udata->p[ip]*log(10);
            break;
        case AMICI_SCALING_LN:
            coefficient = 1.0;
            for(int ip = 0; ip < udata->nplist; ++ip)
            pcoefficient[ip] = udata->p[udata->plist[ip]];
            if (udata->sensi == 2)
            if (udata->o2mode == AMICI_O2MODE_FULL)
            for(int ip = 0; ip < udata->np; ++ip)
            augcoefficient[ip] = udata->p[ip];
            break;
        case AMICI_SCALING_NONE:
            //this should never be reached
            break;
    }
    
    if (udata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        // recover first order sensitivies from states for adjoint sensitivity analysis
        if (udata->sensi == AMICI_SENSI_ORDER_SECOND){
            if (udata->sensi_meth == AMICI_SENSI_ASA){
                if (rdata->x)
                if (rdata->sx)
                for(int ip = 0; ip < udata->nplist; ++ip)
                for(int ix = 0; ix < udata->nxtrue; ++ix)
                for(int it = 0; it < udata->nt; ++it)
                rdata->sx[(ip*udata->nxtrue + ix)*udata->nt + it] = rdata->x[(udata->nxtrue + ip*udata->nxtrue + ix)*udata->nt + it];
                
                if (rdata->y)
                if (rdata->sy)
                for(int ip = 0; ip < udata->nplist; ++ip)
                for(int iy = 0; iy < udata->nytrue; ++iy)
                for(int it = 0; it < udata->nt; ++it)
                rdata->sy[(ip*udata->nytrue + iy)*udata->nt + it] = rdata->y[(udata->nytrue + ip*udata->nytrue + iy)*udata->nt + it];
                
                if (rdata->z)
                if (rdata->sz)
                for(int ip = 0; ip < udata->nplist; ++ip)
                for(int iz = 0; iz < udata->nztrue; ++iz)
                for(int it = 0; it < udata->nt; ++it)
                rdata->sy[(ip * udata->nztrue + iz)*udata->nt + it] = rdata->z[(udata->nztrue + ip*udata->nztrue + iz)*udata->nt + it];
                
            }
        }
        
        if (edata) {
            if (rdata->sllh)
            for(int ip = 0; ip < udata->nplist; ++ip)
            rdata->sllh[ip] *= pcoefficient[ip];
        }
        
#define chainRule(QUANT,IND1,N1T,N1,IND2,N2) \
if (rdata->s ## QUANT ) \
for(int ip = 0; ip < udata->nplist; ++ip) \
for(int IND1 = 0; IND1 < N1T; ++IND1) \
for(int IND2 = 0; IND2 < N2; ++IND2){ \
rdata->s ## QUANT [(ip * N1 + IND1) * N2 + IND2] *= pcoefficient[ip];} \

        chainRule(x,ix,udata->nxtrue,udata->nx,it,udata->nt)
        chainRule(y,iy,udata->nytrue,udata->ny,it,udata->nt)
        chainRule(sigmay,iy,udata->nytrue,udata->ny,it,udata->nt)
        chainRule(z,iz,udata->nztrue,udata->nz,ie,udata->nmaxevent)
        chainRule(sigmaz,iz,udata->nztrue,udata->nz,ie,udata->nmaxevent)
        chainRule(rz,iz,udata->nztrue,udata->nz,ie,udata->nmaxevent)
    }
    if (udata->sensi_meth == AMICI_SENSI_SS) {
        if (rdata->dxdotdp)
        for(int ip = 0; ip < udata->nplist; ++ip)
        for(int ix = 0; ix < udata->nx; ++ix)
        rdata->dxdotdp[ip*udata->nxtrue + ix] *= pcoefficient[ip];
        
        if (rdata->dydp)
        for(int ip = 0; ip < udata->nplist; ++ip)
        for(int iy = 0; iy < udata->ny; ++iy)
        rdata->dydp[ip*udata->nxtrue + iy] *= pcoefficient[ip];
    }
    if (udata->o2mode == AMICI_O2MODE_FULL) { //full
        if (edata){
            if (rdata->s2llh) {
                if (rdata->sllh) {
                    for(int ip = 0; ip < udata->nplist; ++ip) {
                        for(int iJ = 1; iJ < udata->nJ; ++iJ) {
                            rdata->s2llh[ip*udata->nplist+(iJ-1)] *= pcoefficient[ip]*augcoefficient[iJ-1];
                            if (udata->plist[ip] == iJ-1)
                            rdata->s2llh[ip*udata->nplist+(iJ-1)] += rdata->sllh[ip]*coefficient;
                        }
                    }
                }
            }
        }
        
#define s2ChainRule(QUANT,IND1,N1T,N1,IND2,N2) \
if (rdata->s ## QUANT ) \
for(int ip = 0; ip < udata->nplist; ++ip) \
for(int iJ = 1; iJ < udata->nJ; ++iJ) \
for(int IND1 = 0; IND1 < N1T; ++IND1) \
for(int IND2 = 0; IND2 < N2; ++IND2){ \
rdata->s ## QUANT [(ip*N1 + iJ*N1T + IND1)*N2 + IND2] *= pcoefficient[ip]*augcoefficient[iJ-1]; \
if (udata->plist[ip]==iJ-1) \
rdata->s  ## QUANT [(ip*N1 + iJ*N1T + IND1)*N2 + IND2] += rdata->s ## QUANT [(ip*N1 + IND1)*N2 + IND2]*coefficient;}
        
        s2ChainRule(x,ix,udata->nxtrue,udata->nx,it,udata->nt)
        s2ChainRule(y,iy,udata->nytrue,udata->ny,it,udata->nt)
        s2ChainRule(sigmay,iy,udata->nytrue,udata->ny,it,udata->nt)
        s2ChainRule(z,iz,udata->nztrue,udata->nz,ie,udata->nmaxevent)
        s2ChainRule(sigmaz,iz,udata->nztrue,udata->nz,ie,udata->nmaxevent)
        s2ChainRule(rz,iz,udata->nztrue,udata->nz,ie,udata->nmaxevent)
    }
    
    if (udata->o2mode == AMICI_O2MODE_DIR) { //directional
        if (rdata->s2llh) {
            if (rdata->sllh) {
                for(int ip = 0; ip < udata->nplist; ++ip) {
                    rdata->s2llh[ip] *= pcoefficient[ip];
                    rdata->s2llh[ip] += udata->k[udata->nk-udata->nplist+ip]*rdata->sllh[ip]/udata->p[udata->plist[ip]];
                }
            }
        }
        
#define s2vecChainRule(QUANT,IND1,N1T,N1,IND2,N2) \
if (rdata->s ## QUANT ) \
for(int ip = 0; ip < udata->nplist; ++ip) \
for(int IND1 = 0; IND1 < N1T; ++IND1) \
for(int IND2 = 0; IND2 < N2; ++IND2){ \
rdata->s ## QUANT [(ip*N1 + N1T + IND1)*N2 + IND2] *= pcoefficient[ip]; \
rdata->s ## QUANT [(ip*N1 + N1T + IND1)*N2 + IND2] += udata->k[udata->nk-udata->nplist+ip]*rdata->s ## QUANT [(ip*N1 + IND1)*N2 + IND2]/udata->p[udata->plist[ip]];}
        
        s2vecChainRule(x,ix,udata->nxtrue,udata->nx,it,udata->nt)
        s2vecChainRule(y,iy,udata->nytrue,udata->ny,it,udata->nt)
        s2vecChainRule(sigmay,iy,udata->nytrue,udata->ny,it,udata->nt)
        s2vecChainRule(z,iz,udata->nztrue,udata->nz,ie,udata->nmaxevent)
        s2vecChainRule(sigmaz,iz,udata->nztrue,udata->nz,ie,udata->nmaxevent)
        s2vecChainRule(rz,iz,udata->nztrue,udata->nz,ie,udata->nmaxevent)
    }
    
    delete[] pcoefficient;
    delete[] augcoefficient;
    return AMICI_SUCCESS;
}

int fsy(realtype t, int it, realtype *sy, realtype *dydx, realtype *dydp, N_Vector *sx, void *user_data){
    // Compute sy = dydx * sx + dydp
    
    int status = AMICI_SUCCESS;
    UserData *udata = (UserData*) user_data;
    
    for(int ip = 0; ip < udata->nplist; ++ip) {
        for(int iy = 0; iy < udata->ny; ++iy)
        // copy dydp to sy
        sy[ip * udata->nt * udata->ny + iy * udata->nt + it] = dydp[iy + ip * udata->ny];
        
        realtype *sxTmp = N_VGetArrayPointer(sx[ip]);
        
        // compute sy = 1.0*dydx*sx + 1.0*sy
        amici_dgemv(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, udata->ny, udata->nx,
                    1.0, dydx, udata->ny, sxTmp, 1,
                    1.0, &sy[ip * udata->nt * udata->ny + it], udata->nt);
    }
    
    return status;
}

int fsJy(realtype t, int it, realtype *sJy, realtype *s2Jy, realtype *dJydy, realtype *dJydp, realtype *y, realtype *sigmay, realtype *sy, realtype *dydp, realtype *my, void *user_data) {
    int status = AMICI_SUCCESS;
    UserData *udata = (UserData*) user_data;
    
    // Compute sy-dydp for current 'it'
    // dydp         ny x nplist
    // sy           nt x ny x nṕlist
    // dydp part needs to be substracted as it is already contained in dJydp
    // we only need to account for sensitivities here
    realtype *diff = new realtype[udata->ny * udata->nplist];
    for(int iy = 0; iy < udata->ny; ++iy)
    for(int ip = 0; ip < udata->nplist; ++ip)
    diff[iy + ip * udata->ny] = sy[ip * udata->nt * udata->ny + iy * udata->nt + it] - dydp[iy + ip * udata->ny];
    
    // sJy          nplist x ng
    // dJydp=dJydp   nytrue x nplist x ng
    // dJydy=dJydy   nytrue x ng x ny
    
    realtype *dJydyTmp = new realtype[udata->nJ * udata->ny];
    realtype *multResult = new realtype[udata->nplist * udata->nJ];
    
    for(int iyt = 0; iyt < udata->nytrue; ++iyt) {
        if (amiIsNaN(my[udata->nt * iyt + it]))
        continue;
        
        // copy current (iyt) dJydy slice
        // dJydyTmp     ng x ny
        for(int iJ = 0; iJ < udata->nJ; ++iJ)
            for(int iy = 0; iy < udata->ny; ++iy)
                dJydyTmp[iJ + iy * udata->nJ] = dJydy[iyt + iJ * udata->nytrue + iy * udata->nytrue * udata->nJ];
        
        // compute multResult = (dJydyTmp * diff)' + dJydp == diff' * dJydyTmp' + dJydp
        // copy dJydp slice (iyt) to result
        for(int ip = 0; ip < udata->nplist; ++ip)
            for(int iJ = 0; iJ < udata->nJ; ++iJ)
                multResult[ip + udata->np * iJ] = dJydp[iyt + ip * udata->nytrue + iJ * udata->nytrue * udata->nplist];
        
        // C := alpha*op(A)*op(B) + beta*C,
        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_Trans, AMICI_BLAS_Trans,
                    udata->nplist, udata->nJ, udata->ny,
                    1.0, diff, udata->ny,
                    dJydyTmp, udata->nJ,
                    1.0, multResult, udata->nplist);
        
        // sJy += multResult
        for(int iJ = 0; iJ < udata->nJ; ++iJ) {
            if (iJ == 0)
            for(int ip = 0; ip < udata->nplist; ++ip)
                sJy[ip] -= multResult[ip];
            else
                for(int ip = 0; ip < udata->nplist; ++ip)
                    s2Jy[ip + udata->nplist * (iJ - 1)] -= multResult[ip+ udata->nplist * iJ];
        }
        
        
    }
    delete[] dJydyTmp;
    delete[] multResult;
    delete[] diff;
    
    return(status);
}

