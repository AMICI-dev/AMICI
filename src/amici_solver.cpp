#include "include/amici_solver.h"
#include "include/amici.h"
#include <cstdio>
#include <cstring>
#include <include/amici_model.h>

void *Solver::setupAMI(UserData *udata, TempData *tdata, Model *model)
{
    int status;
    void *ami_mem = NULL; /* pointer to ami memory block */
    N_Vector id = NULL;


    tdata->t = udata->tstart;

    if(model->initialize(udata, tdata) != AMICI_SUCCESS) goto freturn;


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

    status = setLinearSolver(udata, ami_mem);
    if(status != AMICI_SUCCESS)
        goto freturn;

    if (udata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        if (udata->sensi_meth == AMICI_SENSI_FSA) {
            if (udata->nx>0) {

                /* initialise sensitivities, this can either be user provided or come from the model definition */
                realtype *sx_tmp;

                if (!udata->sx0data) {
                    if (model->fsx0(tdata->sx, tdata->x, tdata->dx, udata) != AMICI_SUCCESS) goto freturn;
                } else {
                    int ip;
                    for (ip=0; ip<udata->nplist; ip++) {
                        sx_tmp = NV_DATA_S(tdata->sx[ip]);
                        if(!sx_tmp) goto freturn;
                        int ix;
                        for (ix=0; ix<udata->nx; ix++) {
                            sx_tmp[ix] = (realtype) udata->sx0data[ix + udata->nx*ip];
                        }
                    }
                }

                if (model->fsdx0(tdata->sdx, tdata->x, tdata->dx, udata) != AMICI_SUCCESS) goto freturn;

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
    if(!id) goto freturn;
    if(!udata->idlist) goto freturn;
    memcpy(NV_CONTENT_S(id)->data, udata->idlist, udata->nx * sizeof(double));

    if (AMISetId(ami_mem, id) != AMICI_SUCCESS) goto freturn;
    if(id) N_VDestroy_Serial(id);

    if (AMISetSuppressAlg(ami_mem, TRUE) != AMICI_SUCCESS) goto freturn;

    return(ami_mem);

freturn:
    if(id) N_VDestroy_Serial(id);
    if(ami_mem) AMIFree(&ami_mem);
    return NULL;

}

int Solver::setupAMIB(void *ami_mem, UserData *udata, TempData *tdata) {
    int status = AMICI_SUCCESS;

    /* write initial conditions */
    if(!tdata->xB) return AMICI_ERROR_SETUPB;
    realtype *xB_tmp = NV_DATA_S(tdata->xB);
    if(!xB_tmp) return AMICI_ERROR_SETUPB;
    memset(xB_tmp,0,sizeof(realtype)*udata->nxtrue*udata->nJ);
    for (int ix=0; ix<udata->nxtrue; ++ix)
        for (int iJ=0; iJ<udata->nJ; ++iJ)
            xB_tmp[ix + iJ * udata->nxtrue] += tdata->dJydx[udata->nt-1 + (iJ + ix * udata->nJ) * udata->nt];

    if(!tdata->dxB) return AMICI_ERROR_SETUPB;
    if(!NV_DATA_S(tdata->dxB)) return AMICI_ERROR_SETUPB;
    memset(NV_DATA_S(tdata->dxB),0,sizeof(realtype)*udata->nx);

    if(!tdata->xQB) return AMICI_ERROR_SETUPB;
    if(!NV_DATA_S(tdata->xQB)) return AMICI_ERROR_SETUPB;
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

void Solver::wrap_ErrHandlerFn(int error_code, const char *module, const char *function, char *msg, void *eh_data)
{
    char buffer[250];
    char buffid[250];
    sprintf(buffer,"AMICI ERROR: in module %s in function %s : %s ",module,function,msg);
    switch (error_code) {
    case 99:
        sprintf(buffid,"AMICI:mex:%s:%s:CV_WARNING",module,function);
        break;

    case -1:
        sprintf(buffid,"AMICI:mex:%s:%s:CV_TOO_MUCH_WORK",module,function);
        break;

    case -2:
        sprintf(buffid,"AMICI:mex:%s:%s:CV_TOO_MUCH_ACC",module,function);
        break;

    case -3:
        sprintf(buffid,"AMICI:mex:%s:%s:CV_ERR_FAILURE",module,function);
        break;

    case -4:
        sprintf(buffid,"AMICI:mex:%s:%s:CV_CONV_FAILURE",module,function);
        break;

    default:
        sprintf(buffid,"AMICI:mex:%s:%s:CV_OTHER",module,function);
        break;
    }

    warnMsgIdAndTxt(buffid,buffer);
}

int Solver::getDiagnosis(int it, void *ami_mem, ReturnData *rdata)
{
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

int Solver::getDiagnosisB(int it, void *ami_mem, UserData *udata, ReturnData *rdata, TempData *tdata)
{
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

int Solver::setLinearSolver(const UserData *udata, void *ami_mem) {
    int status;
    /* Attach linear solver module */

    switch (udata->linsol) {

    /* DIRECT SOLVERS */

    case AMICI_DENSE:

        status = AMIDense(ami_mem, udata->nx);
        if (status != AMICI_SUCCESS) return status;

        return wrap_SetDenseJacFn(ami_mem);

    case AMICI_BAND:
        status = AMIBand(ami_mem, udata->nx, udata->ubw, udata->lbw);
        if (status != AMICI_SUCCESS) return status;

        return wrap_SetBandJacFn(ami_mem);

    case AMICI_LAPACKDENSE:
        errMsgIdAndTxt("AMICI:mex:lapack","Solver currently not supported!");
        /* status = CVLapackDense(ami_mem, nx);
             if (status != AMICI_SUCCESS) return;

             status = wrap_SetDenseJacFn(ami_mem);
             if (status != AMICI_SUCCESS) return;
        */
        return AMICI_ERROR_NOT_IMPLEMENTED;

    case AMICI_LAPACKBAND:

        errMsgIdAndTxt("AMICI:mex:lapack","Solver currently not supported!");
        /* status = CVLapackBand(ami_mem, nx);
             if (status != AMICI_SUCCESS) return;

             status = wrap_SetBandJacFn(ami_mem);
             if (status != AMICI_SUCCESS) return;
        */
        return AMICI_ERROR_NOT_IMPLEMENTED;

    case AMICI_DIAG:
        return AMIDiag(ami_mem);

        /* ITERATIVE SOLVERS */

    case AMICI_SPGMR:
        status = AMISpgmr(ami_mem, PREC_NONE, CVSPILS_MAXL);
        if (status != AMICI_SUCCESS) return status;

        return wrap_SetJacTimesVecFn(ami_mem);

    case AMICI_SPBCG:
        status = AMISpbcg(ami_mem, PREC_NONE, CVSPILS_MAXL);
        if (status != AMICI_SUCCESS) return status;

        return wrap_SetJacTimesVecFn(ami_mem);

    case AMICI_SPTFQMR:
        status = AMISptfqmr(ami_mem, PREC_NONE, CVSPILS_MAXL);
        if (status != AMICI_SUCCESS) return status;

        return wrap_SetJacTimesVecFn(ami_mem);

        /* SPARSE SOLVERS */

    case AMICI_KLU:
        status = AMIKLU(ami_mem, udata->nx, udata->nnz, CSC_MAT);
        if (status != AMICI_SUCCESS) return status;

        status = wrap_SetSparseJacFn(ami_mem);
        if (status != AMICI_SUCCESS) return status;

        return AMIKLUSetOrdering(ami_mem, udata->ordering);
    }

    errMsgIdAndTxt("AMICI:mex:solver","Invalid choice of solver!");

    return AMICI_ERROR_OTHER;
}
