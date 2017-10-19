#include "include/amici_solver.h"
#include "include/amici.h"
#include <cstdio>
#include <cstring>
#include <include/amici_model.h>
#include <include/rdata.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <sundials/sundials_spgmr.h>
// TODO: don't use cvodes includes here
#include <cvodes/cvodes_spils.h>

namespace amici {


/**
 * @brief setupAMIs initialises the ami memory object
 * @param[in] udata pointer to the user data object @type UserData
 * @param[in] tdata pointer to the temporary data object @type TempData
 * @param[in] model pointer to the model object @type Model
 * @return status flag indicating successful execution
 */
int Solver::setupAMI(const UserData *udata, TempData *tdata, Model *model) {
    int status;
    tdata->t = udata->tstart;

    if (model->initialize(udata, tdata) != AMICI_SUCCESS)
        goto freturn;

    /* Create AMIS object */
    if (udata->lmm != CV_ADAMS && udata->lmm != CV_BDF) {
        errMsgIdAndTxt("AMICI:mex:lmm", "Illegal value for lmm!");
        goto freturn;
    }
    if (udata->iter != CV_NEWTON && udata->iter != CV_FUNCTIONAL) {
        errMsgIdAndTxt("AMICI:mex:iter", "Illegal value for iter!");
        goto freturn;
    }
    ami_mem = AMICreate(udata->lmm, udata->iter);
    if (ami_mem == NULL)
        goto freturn;

    /* Initialize AMIS solver*/
    if (init(tdata->x, tdata->dx, udata->tstart) != AMICI_SUCCESS)
        goto freturn;

    /* Specify integration tolerances */
    if (AMISStolerances(RCONST(udata->rtol), RCONST(udata->atol)) !=
        AMICI_SUCCESS)
        goto freturn;

    /* Set optional inputs */
    if (AMISetErrHandlerFn() != AMICI_SUCCESS)
        goto freturn;

    /* attaches userdata*/
    if (AMISetUserData(tdata) != AMICI_SUCCESS)
        goto freturn;

    /* specify maximal number of steps */
    if (AMISetMaxNumSteps(udata->maxsteps) != AMICI_SUCCESS)
        goto freturn;

    /* activates stability limit detection */
    if (AMISetStabLimDet(udata->stldet) != AMICI_SUCCESS)
        goto freturn;

    if (rootInit(model->ne) != AMICI_SUCCESS)
        goto freturn;

    status = setLinearSolver(udata, model);
    if (status != AMICI_SUCCESS)
        goto freturn;

    if (udata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        
        if (model->nx > 0) {

            /* initialise sensitivities, this can either be user provided or
             * come from the model definition */
            realtype *sx_tmp;

            if (!udata->sx0data) {
                if (model->fsx0(tdata->sx, tdata->x, tdata->dx, tdata) !=
                    AMICI_SUCCESS)
                    goto freturn;
            } else {
                for (int ip = 0; ip < udata->nplist; ip++) {
                    sx_tmp = NV_DATA_S(tdata->sx[ip]);
                    if (!sx_tmp)
                        goto freturn;
                    for (int ix = 0; ix < model->nx; ix++) {
                        sx_tmp[ix] =
                            (realtype)udata->sx0data[ix + model->nx * ip];
                    }
                }
            }

            if (model->fsdx0(tdata->sdx, tdata->x, tdata->dx, tdata) !=
                AMICI_SUCCESS)
                goto freturn;
            
            if (udata->sensi_meth == AMICI_SENSI_FSA) {
                
                /* Activate sensitivity calculations */
                if (sensInit1(tdata->sx, tdata->sdx, udata) != AMICI_SUCCESS)
                    goto freturn;

                /* Set sensitivity analysis optional inputs */
                if (AMISetSensParams(tdata->p, udata->pbar, udata->plist) !=
                    AMICI_SUCCESS)
                    goto freturn;

                if (AMISetSensErrCon(TRUE) != AMICI_SUCCESS)
                    goto freturn;

                if (AMISensEEtolerances() != AMICI_SUCCESS)
                    goto freturn;
            }
        }

        if (udata->sensi_meth == AMICI_SENSI_ASA) {
            if (model->nx > 0) {
                /* Allocate space for the adjoint computation */
                if (AMIAdjInit(udata->maxsteps, udata->interpType) !=
                    AMICI_SUCCESS)
                    goto freturn;
            }
        }
    }

    if (AMISetId(model) != AMICI_SUCCESS)
        goto freturn;

    if (AMISetSuppressAlg(TRUE) != AMICI_SUCCESS)
        goto freturn;
    
    /* calculate consistent DAE initial conditions (no effect for ODE) */
    //if(udata->nt>1){
    //    if (AMICalcIC(udata->ts[1]) != AMICI_SUCCESS)
    //        goto freturn;
    //}

    return AMICI_SUCCESS;

freturn:
    if (ami_mem)
        AMIFree();
    return AMICI_ERROR_SETUP;
}

/**
 * setupAMIB initialises the AMI memory object for the backwards problem
 * @param[in] udata pointer to the user data object @type UserData
 * @param[in] tdata pointer to the temporary data object @type TempData
 * @param[in] model pointer to the model object @type Model
 * @return status flag indicating successful execution
 */
int Solver::setupAMIB(const UserData *udata, TempData *tdata, Model *model) {
    int status = AMICI_SUCCESS;

    /* write initial conditions */
    if (!tdata->xB)
        return AMICI_ERROR_SETUPB;
    realtype *xB_tmp = NV_DATA_S(tdata->xB);
    if (!xB_tmp)
        return AMICI_ERROR_SETUPB;
    memset(xB_tmp, 0, sizeof(realtype) * model->nxtrue * model->nJ);
    for (int ix = 0; ix < model->nxtrue; ++ix)
        for (int iJ = 0; iJ < model->nJ; ++iJ)
            xB_tmp[ix + iJ * model->nxtrue] +=
                tdata->dJydx[tdata->rdata->nt - 1 +
                             (iJ + ix * model->nJ) * tdata->rdata->nt];

    if (!tdata->dxB)
        return AMICI_ERROR_SETUPB;
    if (!NV_DATA_S(tdata->dxB))
        return AMICI_ERROR_SETUPB;
    memset(NV_DATA_S(tdata->dxB), 0, sizeof(realtype) * model->nx);

    if (!tdata->xQB)
        return AMICI_ERROR_SETUPB;
    if (!NV_DATA_S(tdata->xQB))
        return AMICI_ERROR_SETUPB;
    memset(NV_DATA_S(tdata->xQB), 0,
           sizeof(realtype) * model->nJ * tdata->rdata->nplist);

    /* create backward problem */
    if (udata->lmm > 2 || udata->lmm < 1) {
        errMsgIdAndTxt("AMICI:mex:lmm", "Illegal value for lmm!");
    }
    if (udata->iter > 2 || udata->iter < 1) {
        errMsgIdAndTxt("AMICI:mex:iter", "Illegal value for iter!");
    }

    /* allocate memory for the backward problem */
    status = AMICreateB(udata->lmm, udata->iter, &(tdata->which));
    if (status != AMICI_SUCCESS)
        return status;

    /* initialise states */
    status = binit(tdata->which, tdata->xB, tdata->dxB, tdata->t);
    if (status != AMICI_SUCCESS)
        return status;

    /* specify integration tolerances for backward problem */
    status = AMISStolerancesB(tdata->which, RCONST(udata->rtol),
                              RCONST(udata->atol));
    if (status != AMICI_SUCCESS)
        return status;

    /* Attach user data */
    status = AMISetUserDataB(tdata->which, tdata);
    if (status != AMICI_SUCCESS)
        return status;

    /* Number of maximal internal steps */
    if (AMISetMaxNumStepsB(tdata->which, 100 * udata->maxsteps) !=
        AMICI_SUCCESS)
        return AMICI_ERROR_SETUPB;

    switch (udata->linsol) {

    /* DIRECT SOLVERS */

    case AMICI_DENSE:
        status = AMIDenseB(tdata->which, model->nx);
        if (status != AMICI_SUCCESS)
            return status;

        status = setDenseJacFnB(tdata->which);
        if (status != AMICI_SUCCESS)
            return status;

        break;

    case AMICI_BAND:
        status = AMIBandB(tdata->which, model->nx, model->ubw, model->lbw);
        if (status != AMICI_SUCCESS)
            return status;

        status = setBandJacFnB(tdata->which);
        if (status != AMICI_SUCCESS)
            return status;

        break;

    case AMICI_LAPACKDENSE:

        /* #if SUNDIALS_BLAS_LAPACK
                 status = CVLapackDenseB(ami_mem, tdata->which, nx);
                 if (status != AMICI_SUCCESS) return;

                 status = SetDenseJacFnB(ami_mem, tdata->which);
                 if (status != AMICI_SUCCESS) return;
                 #else*/
        errMsgIdAndTxt("AMICI:mex:lapack", "Solver currently not supported!");
        /* #endif*/
        break;

    case AMICI_LAPACKBAND:

        /* #if SUNDIALS_BLAS_LAPACK
                 status = CVLapackBandB(ami_mem, tdata->which, nx, ubw, lbw);
                 if (status != AMICI_SUCCESS) return;

                 status = SetBandJacFnB(ami_mem, tdata->which);
                 if (status != AMICI_SUCCESS) return;
                 #else*/
        errMsgIdAndTxt("AMICI:mex:lapack", "Solver currently not supported!");
        /* #endif*/
        break;

    case AMICI_DIAG:
        status = AMIDiagB(tdata->which);
        if (status != AMICI_SUCCESS)
            return status;

        status = setDenseJacFnB(tdata->which);
        if (status != AMICI_SUCCESS)
            return status;

        break;

    /* ITERATIVE SOLVERS */

    case AMICI_SPGMR:
        status = AMISpgmrB(tdata->which, PREC_NONE, CVSPILS_MAXL);
        if (status != AMICI_SUCCESS)
            return status;

        status = setJacTimesVecFnB(tdata->which);
        if (status != AMICI_SUCCESS)
            return status;

        break;

    case AMICI_SPBCG:
        status = AMISpbcgB(tdata->which, PREC_NONE, CVSPILS_MAXL);
        if (status != AMICI_SUCCESS)
            return status;

        status = setJacTimesVecFnB(tdata->which);
        if (status != AMICI_SUCCESS)
            return status;

        break;

    case AMICI_SPTFQMR:
        status = AMISptfqmrB(tdata->which, PREC_NONE, CVSPILS_MAXL);
        if (status != AMICI_SUCCESS)
            return status;

        status = setJacTimesVecFnB(tdata->which);
        if (status != AMICI_SUCCESS)
            return status;

        break;

    /* SPARSE SOLVERS */

    case AMICI_KLU:
        status = AMIKLUB(tdata->which, model->nx, model->nnz, CSC_MAT);
        if (status != AMICI_SUCCESS)
            return status;

        status = setSparseJacFnB(tdata->which);
        if (status != AMICI_SUCCESS)
            return status;

        status = AMIKLUSetOrderingB(tdata->which, udata->ordering);
        if (status != AMICI_SUCCESS)
            return status;

        break;

    default:
        break;
    }

    /* Initialise quadrature calculation */
    status = qbinit(tdata->which, tdata->xQB);
    if (status != AMICI_SUCCESS)
        return status;

    /* Enable Quadrature Error Control */
    status = AMISetQuadErrConB(tdata->which, TRUE);
    if (status != AMICI_SUCCESS)
        return status;

    status = AMIQuadSStolerancesB(tdata->which, RCONST(udata->rtol),
                                  RCONST(udata->atol));
    if (status != AMICI_SUCCESS)
        return status;

    status = AMISetStabLimDetB(tdata->which, udata->stldet);
    if (status != AMICI_SUCCESS)
        return status;

    return status;
}

/**
 * ErrHandlerFn extracts diagnosis information from solver memory block and
 * writes them into the return data object for the backward problem
 *
 * @param[in] error_code error identifier @type int
 * @param[in] module name of the module in which the error occured @type char
 * @param[in] function name of the function in which the error occured @type
 * char
 * @param[in] msg error message @type char
 * @param[in] eh_data unused input
 */
void Solver::wrapErrHandlerFn(int error_code, const char *module,
                              const char *function, char *msg, void *eh_data) {
    char buffer[250];
    char buffid[250];
    sprintf(buffer, "AMICI ERROR: in module %s in function %s : %s ", module,
            function, msg);
    switch (error_code) {
    case 99:
        sprintf(buffid, "AMICI:mex:%s:%s:CV_WARNING", module, function);
        break;

    case -1:
        sprintf(buffid, "AMICI:mex:%s:%s:CV_TOO_MUCH_WORK", module, function);
        break;

    case -2:
        sprintf(buffid, "AMICI:mex:%s:%s:CV_TOO_MUCH_ACC", module, function);
        break;

    case -3:
        sprintf(buffid, "AMICI:mex:%s:%s:CV_ERR_FAILURE", module, function);
        break;

    case -4:
        sprintf(buffid, "AMICI:mex:%s:%s:CV_CONV_FAILURE", module, function);
        break;

    default:
        sprintf(buffid, "AMICI:mex:%s:%s:CV_OTHER", module, function);
        break;
    }

    warnMsgIdAndTxt(buffid, buffer);
}

/**
 * getDiagnosis extracts diagnosis information from solver memory block and
 * writes them into the return data object
 *
 * @param[in] it time-point index @type int
 * @param[out] rdata pointer to the return data object @type ReturnData
 * @return status flag indicating success of execution @type int
 */
int Solver::getDiagnosis(const int it, ReturnData *rdata) {
    long int number;
    int status = AMICI_SUCCESS;
    int order;

    status = AMIGetNumSteps(ami_mem, &number);
    if (status != AMICI_SUCCESS)
        return status;
    rdata->numsteps[it] = (double)number;

    status = AMIGetNumRhsEvals(ami_mem, &number);
    if (status != AMICI_SUCCESS)
        return status;
    rdata->numrhsevals[it] = (double)number;

    status = AMIGetNumErrTestFails(ami_mem, &number);
    if (status != AMICI_SUCCESS)
        return status;
    rdata->numerrtestfails[it] = (double)number;

    status = AMIGetNumNonlinSolvConvFails(ami_mem, &number);
    if (status != AMICI_SUCCESS)
        return status;
    rdata->numnonlinsolvconvfails[it] = (double)number;

    status = AMIGetLastOrder(ami_mem, &order);
    if (status != AMICI_SUCCESS)
        return status;
    rdata->order[it] = (double)order;

    return status;
}

/**
 * getDiagnosisB extracts diagnosis information from solver memory block and
 * writes them into the return data object for the backward problem
 *
 * @param[in] it time-point index @type int
 * @param[out] rdata pointer to the return data object @type ReturnData
 * @param[out] tdata pointer to the temporary data object @type TempData
 * @return status flag indicating success of execution @type int
 */
int Solver::getDiagnosisB(const int it, ReturnData *rdata,
                          const TempData *tdata) {
    long int number;
    int status = AMICI_SUCCESS;

    void *ami_memB = AMIGetAdjBmem(ami_mem, tdata->which);

    status = AMIGetNumSteps(ami_memB, &number);
    if (status != AMICI_SUCCESS)
        return status;
    rdata->numstepsB[it] = (double)number;

    status = AMIGetNumRhsEvals(ami_memB, &number);
    if (status != AMICI_SUCCESS)
        return status;
    rdata->numrhsevalsB[it] = (double)number;

    status = AMIGetNumErrTestFails(ami_memB, &number);
    if (status != AMICI_SUCCESS)
        return status;
    rdata->numerrtestfailsB[it] = (double)number;

    status = AMIGetNumNonlinSolvConvFails(ami_memB, &number);
    if (status != AMICI_SUCCESS)
        return status;
    rdata->numnonlinsolvconvfailsB[it] = (double)number;

    return status;
}

/**
 * setLinearSolver sets the linear solver
 *
 * @param[out] udata pointer to the user data object @type UserData
 * @param[in] model pointer to the model object @type Model
 * @return status flag indicating success of execution @type int
 */
int Solver::setLinearSolver(const UserData *udata, Model *model) {
    int status;
    /* Attach linear solver module */

    switch (udata->linsol) {

    /* DIRECT SOLVERS */

    case AMICI_DENSE:

        status = AMIDense(model->nx);
        if (status != AMICI_SUCCESS)
            return status;

        return setDenseJacFn();

    case AMICI_BAND:
        status = AMIBand(model->nx, model->ubw, model->lbw);
        if (status != AMICI_SUCCESS)
            return status;

        return setBandJacFn();

    case AMICI_LAPACKDENSE:
        errMsgIdAndTxt("AMICI:mex:lapack", "Solver currently not supported!");
        /* status = CVLapackDense(ami_mem, nx);
             if (status != AMICI_SUCCESS) return;

             status = SetDenseJacFn(ami_mem);
             if (status != AMICI_SUCCESS) return;
        */
        return AMICI_ERROR_NOT_IMPLEMENTED;

    case AMICI_LAPACKBAND:

        errMsgIdAndTxt("AMICI:mex:lapack", "Solver currently not supported!");
        /* status = CVLapackBand(ami_mem, nx);
             if (status != AMICI_SUCCESS) return;

             status = SetBandJacFn(ami_mem);
             if (status != AMICI_SUCCESS) return;
        */
        return AMICI_ERROR_NOT_IMPLEMENTED;

    case AMICI_DIAG:
        return AMIDiag();

    /* ITERATIVE SOLVERS */

    case AMICI_SPGMR:
        status = AMISpgmr(PREC_NONE, CVSPILS_MAXL);
        if (status != AMICI_SUCCESS)
            return status;

        return setJacTimesVecFn();

    case AMICI_SPBCG:
        status = AMISpbcg(PREC_NONE, CVSPILS_MAXL);
        if (status != AMICI_SUCCESS)
            return status;

        return setJacTimesVecFn();

    case AMICI_SPTFQMR:
        status = AMISptfqmr(PREC_NONE, CVSPILS_MAXL);
        if (status != AMICI_SUCCESS)
            return status;

        return setJacTimesVecFn();

    /* SPARSE SOLVERS */

    case AMICI_KLU:
        status = AMIKLU(model->nx, model->nnz, CSC_MAT);
        if (status != AMICI_SUCCESS)
            return status;

        status = setSparseJacFn();
        if (status != AMICI_SUCCESS)
            return status;

        return AMIKLUSetOrdering(udata->ordering);
    }

    errMsgIdAndTxt("AMICI:mex:solver", "Invalid choice of solver!");

    return AMICI_ERROR_OTHER;
}

} // namespace amici
