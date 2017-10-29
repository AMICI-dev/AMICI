#include "include/amici_solver.h"
#include "include/amici.h"
#include "include/amici_exception.h"
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
 */
void Solver::setupAMI(const UserData *udata, TempData *tdata, Model *model) {
    tdata->t = udata->tstart;

    model->initialize(udata, tdata);

    /* Create solver memory object */
    if (udata->lmm != CV_ADAMS && udata->lmm != CV_BDF) {
        throw AmiException("Illegal value for lmm!");
    }
    if (udata->iter != CV_NEWTON && udata->iter != CV_FUNCTIONAL) {
        throw AmiException("Illegal value for iter!");
    }
    ami_mem = AMICreate(udata->lmm, udata->iter);
    if (ami_mem == NULL)
        throw AmiException("Failed to allocated solver memory!");
    try {
    /* Initialize AMIS solver*/
    init(tdata->x, tdata->dx, udata->tstart);
    /* Specify integration tolerances */
    AMISStolerances(RCONST(udata->rtol), RCONST(udata->atol));
    /* Set optional inputs */
    AMISetErrHandlerFn();
    /* attaches userdata*/
    AMISetUserData(tdata);
    /* specify maximal number of steps */
    AMISetMaxNumSteps(udata->maxsteps);
    /* activates stability limit detection */
    AMISetStabLimDet(udata->stldet);
    
    rootInit(model->ne);
    
    setLinearSolver(udata, model);
        
    if (udata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        
        if (model->nx > 0) {

            /* initialise sensitivities, this can either be user provided or
             * come from the model definition */
            realtype *sx_tmp;

            if (!udata->sx0data) {
                model->fsx0(tdata->sx, tdata->x, tdata->dx, tdata);
            } else {
                for (int ip = 0; ip < udata->nplist; ip++) {
                    sx_tmp = NV_DATA_S(tdata->sx[ip]);
                    if (!sx_tmp)
                        throw NullPointerException("sx_tmp");
                    for (int ix = 0; ix < model->nx; ix++) {
                        sx_tmp[ix] =
                            (realtype)udata->sx0data[ix + model->nx * ip];
                    }
                }
            }

            model->fsdx0(tdata->sdx, tdata->x, tdata->dx, tdata);
            
            if (udata->sensi_meth == AMICI_SENSI_FSA) {
                
                /* Activate sensitivity calculations */
                sensInit1(tdata->sx, tdata->sdx, udata);
                /* Set sensitivity analysis optional inputs */
                AMISetSensParams(tdata->p, udata->pbar, udata->plist);
                AMISetSensErrCon(TRUE);
                AMISensEEtolerances();
            }
        }

        if (udata->sensi_meth == AMICI_SENSI_ASA) {
            if (model->nx > 0) {
                /* Allocate space for the adjoint computation */
                AMIAdjInit(udata->maxsteps, udata->interpType);
            }
        }
    }

    AMISetId(model);
    AMISetSuppressAlg(TRUE);
    /* calculate consistent DAE initial conditions (no effect for ODE) */
    //if(udata->nt>1){
    //    if (AMICalcIC(udata->ts[1]) != AMICI_SUCCESS)
    //        goto freturn;
    //}
    } catch (...) {
        AMIFree();
        throw AmiException("setupAMI routine failed!");
    }
}

/**
 * setupAMIB initialises the AMI memory object for the backwards problem
 * @param[in] udata pointer to the user data object @type UserData
 * @param[in] tdata pointer to the temporary data object @type TempData
 * @param[in] model pointer to the model object @type Model
 */
void Solver::setupAMIB(const UserData *udata, TempData *tdata, Model *model) {

    /* write initial conditions */
    if (!tdata->xB)
        throw NullPointerException("tdata->xB");
    realtype *xB_tmp = NV_DATA_S(tdata->xB);
    if (!xB_tmp)
        throw NullPointerException("xB_tmp");
    memset(xB_tmp, 0, sizeof(realtype) * model->nxtrue * model->nJ);
    for (int ix = 0; ix < model->nxtrue; ++ix)
        for (int iJ = 0; iJ < model->nJ; ++iJ)
            xB_tmp[ix + iJ * model->nxtrue] +=
                tdata->dJydx[tdata->rdata->nt - 1 +
                             (iJ + ix * model->nJ) * tdata->rdata->nt];

    if (!tdata->dxB)
        throw NullPointerException("tdata->dxB");
    if (!NV_DATA_S(tdata->dxB))
        throw NullPointerException("xdB_tmp");
    memset(NV_DATA_S(tdata->dxB), 0, sizeof(realtype) * model->nx);

    if (!tdata->xQB)
        throw NullPointerException("tdata->xQB");
    if (!NV_DATA_S(tdata->xQB))
        throw NullPointerException("xQB_tmp");
    memset(NV_DATA_S(tdata->xQB), 0,
           sizeof(realtype) * model->nJ * tdata->rdata->nplist);

    /* create backward problem */
    if (udata->lmm > 2 || udata->lmm < 1) {
        throw AmiException("Illegal value for lmm!");
    }
    if (udata->iter > 2 || udata->iter < 1) {
        throw AmiException("Illegal value for iter!");
    }

    /* allocate memory for the backward problem */
    AMICreateB(udata->lmm, udata->iter, &(tdata->which));

    /* initialise states */
    binit(tdata->which, tdata->xB, tdata->dxB, tdata->t);

    /* specify integration tolerances for backward problem */
    AMISStolerancesB(tdata->which, RCONST(udata->rtol),
                              RCONST(udata->atol));

    /* Attach user data */
    AMISetUserDataB(tdata->which, tdata);

    /* Number of maximal internal steps */
    AMISetMaxNumStepsB(tdata->which, 100 * udata->maxsteps);
    
    setLinearSolverB(udata, model, tdata->which);
    
    /* Initialise quadrature calculation */
    qbinit(tdata->which, tdata->xQB);

    /* Enable Quadrature Error Control */
    AMISetQuadErrConB(tdata->which, TRUE);

    AMIQuadSStolerancesB(tdata->which, RCONST(udata->rtol),
                                  RCONST(udata->atol));

    AMISetStabLimDetB(tdata->which, udata->stldet);
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
 */
void Solver::getDiagnosis(const int it, ReturnData *rdata) {
    long int number;
    int order;

    AMIGetNumSteps(ami_mem, &number);
    rdata->numsteps[it] = (double)number;

    AMIGetNumRhsEvals(ami_mem, &number);
    rdata->numrhsevals[it] = (double)number;

    AMIGetNumErrTestFails(ami_mem, &number);
    rdata->numerrtestfails[it] = (double)number;

    AMIGetNumNonlinSolvConvFails(ami_mem, &number);
    rdata->numnonlinsolvconvfails[it] = (double)number;

    AMIGetLastOrder(ami_mem, &order);
    rdata->order[it] = (double)order;
    
    return;
}

/**
 * getDiagnosisB extracts diagnosis information from solver memory block and
 * writes them into the return data object for the backward problem
 *
 * @param[in] it time-point index @type int
 * @param[out] rdata pointer to the return data object @type ReturnData
 * @param[out] tdata pointer to the temporary data object @type TempData
 */
void Solver::getDiagnosisB(const int it, ReturnData *rdata,
                          const TempData *tdata) {
    long int number;

    void *ami_memB = AMIGetAdjBmem(ami_mem, tdata->which);

    AMIGetNumSteps(ami_memB, &number);
    rdata->numstepsB[it] = (double)number;

    AMIGetNumRhsEvals(ami_memB, &number);
    rdata->numrhsevalsB[it] = (double)number;

    AMIGetNumErrTestFails(ami_memB, &number);
    rdata->numerrtestfailsB[it] = (double)number;

    AMIGetNumNonlinSolvConvFails(ami_memB, &number);
    rdata->numnonlinsolvconvfailsB[it] = (double)number;

    return;
}

/**
 * setLinearSolver sets the linear solver for the forward problem
 *
 * @param[out] udata pointer to the user data object @type UserData
 * @param[in] model pointer to the model object @type Model
 */
void Solver::setLinearSolver(const UserData *udata, Model *model) {
    /* Attach linear solver module */

    switch (udata->linsol) {

            /* DIRECT SOLVERS */
            
        case AMICI_DENSE:
            AMIDense(model->nx);
            setDenseJacFn();
            break;
            
        case AMICI_BAND:
            AMIBand(model->nx, model->ubw, model->lbw);
            setBandJacFn();
            break;
            
        case AMICI_LAPACKDENSE:
            throw AmiException("Solver currently not supported!");
            /* status = CVLapackDense(ami_mem, nx);
             if (status != AMICI_SUCCESS) return;
             
             status = SetDenseJacFn(ami_mem);
             if (status != AMICI_SUCCESS) return;
             */
            
        case AMICI_LAPACKBAND:
            throw AmiException("Solver currently not supported!");
            /* status = CVLapackBand(ami_mem, nx);
             if (status != AMICI_SUCCESS) return;
             
             status = SetBandJacFn(ami_mem);
             if (status != AMICI_SUCCESS) return;
             */
            
        case AMICI_DIAG:
            AMIDiag();
            break;
            
            
            /* ITERATIVE SOLVERS */
            
        case AMICI_SPGMR:
            AMISpgmr(PREC_NONE, CVSPILS_MAXL);
            setJacTimesVecFn();
            break;
            
        case AMICI_SPBCG:
            AMISpbcg(PREC_NONE, CVSPILS_MAXL);
            setJacTimesVecFn();
            break;
            
        case AMICI_SPTFQMR:
            AMISptfqmr(PREC_NONE, CVSPILS_MAXL);
            setJacTimesVecFn();
            break;
            
            /* SPARSE SOLVERS */
            
        case AMICI_KLU:
            AMIKLU(model->nx, model->nnz, CSC_MAT);
            setSparseJacFn();
            AMIKLUSetOrdering(udata->ordering);
            break;
            
        default:
            throw AmiException("Invalid choice of solver!");
            
    }
    return;
}
    
    /**
     * setLinearSolverB sets the linear solver for the backward problem
     *
     * @param[out] udata pointer to the user data object @type UserData
     * @param[in] model pointer to the model object @type Model
     */
void Solver::setLinearSolverB(const UserData *udata, Model *model, int which) {
    switch (udata->linsol) {
            
            /* DIRECT SOLVERS */
            
        case AMICI_DENSE:
            AMIDenseB(which, model->nx);
            setDenseJacFnB(which);
            break;
            
        case AMICI_BAND:
            AMIBandB(which, model->nx, model->ubw, model->lbw);
            setBandJacFnB(which);
            break;
            
        case AMICI_LAPACKDENSE:
            
            /* #if SUNDIALS_BLAS_LAPACK
             status = CVLapackDenseB(ami_mem, tdata->which, nx);
             if (status != AMICI_SUCCESS) return;
             
             status = SetDenseJacFnB(ami_mem, tdata->which);
             if (status != AMICI_SUCCESS) return;
             #else*/
            throw AmiException("Solver currently not supported!");
            /* #endif*/
            break;
            
        case AMICI_LAPACKBAND:
            
            /* #if SUNDIALS_BLAS_LAPACK
             status = CVLapackBandB(ami_mem, tdata->which, nx, ubw, lbw);
             if (status != AMICI_SUCCESS) return;
             
             status = SetBandJacFnB(ami_mem, tdata->which);
             if (status != AMICI_SUCCESS) return;
             #else*/
            throw AmiException("Solver currently not supported!");
            /* #endif*/
            break;
            
        case AMICI_DIAG:
            AMIDiagB(which);
            setDenseJacFnB(which);
            break;
            
            /* ITERATIVE SOLVERS */
            
        case AMICI_SPGMR:
            AMISpgmrB(which, PREC_NONE, CVSPILS_MAXL);
            setJacTimesVecFnB(which);
            break;
            
        case AMICI_SPBCG:
            AMISpbcgB(which, PREC_NONE, CVSPILS_MAXL);
            setJacTimesVecFnB(which);
            break;
            
        case AMICI_SPTFQMR:
            AMISptfqmrB(which, PREC_NONE, CVSPILS_MAXL);
            setJacTimesVecFnB(which);
            break;
            
            /* SPARSE SOLVERS */
            
        case AMICI_KLU:
            AMIKLUB(which, model->nx, model->nnz, CSC_MAT);
            setSparseJacFnB(which);
            AMIKLUSetOrderingB(which, udata->ordering);
            break;
            
        default:
            throw AmiException("Invalid local Solver!");
            break;
    }
}


} // namespace amici
