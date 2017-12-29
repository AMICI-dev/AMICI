#include "include/amici_solver.h"
#include "include/amici.h"
#include "include/amici_exception.h"
#include "include/forwardproblem.h"
#include "include/backwardproblem.h"
#include <cstdio>
#include <cstring>
#include <include/amici_model.h>
#include <include/rdata.h>
#include <include/udata.h>
#include <sundials/sundials_spgmr.h>
// TODO: don't use cvodes includes here
#include <cvodes/cvodes_spils.h>

namespace amici {


/**
 * @brief setupAMIs initialises the ami memory object
 * @param fwd pointer to forward problem
 * @param udata pointer to the user data object
 * @param model pointer to the model object
 */
void Solver::setupAMI(ForwardProblem *fwd, const UserData *udata, Model *model) {
    model->initialize(fwd->getxptr(), fwd->getdxptr(), udata);

    /* Create solver memory object */
    if (udata->getLinearMultistepMethod() != ADAMS && udata->getLinearMultistepMethod() != BDF) {
        throw AmiException("Illegal value for lmm!");
    }
    if (udata->getNonlinearSolverIteration() != NEWTON && udata->getNonlinearSolverIteration() != FUNCTIONAL) {
        throw AmiException("Illegal value for iter!");
    }
    ami_mem = AMICreate(udata->getLinearMultistepMethod(), udata->getNonlinearSolverIteration());
    if (ami_mem == NULL)
        throw AmiException("Failed to allocated solver memory!");
    try {
    /* Initialize AMIS solver*/
    init(fwd->getxptr(), fwd->getdxptr(), udata->t0());
    /* Specify integration tolerances */
    AMISStolerances(RCONST(udata->rtol), RCONST(udata->atol));
    /* Set optional inputs */
    AMISetErrHandlerFn();
    /* attaches userdata*/
    AMISetUserData(model);
    /* specify maximal number of steps */
    AMISetMaxNumSteps(udata->maxsteps);
    /* activates stability limit detection */
    AMISetStabLimDet(udata->getStabilityLimitFlag());
    
    rootInit(model->ne);
    
    setLinearSolver(udata, model);
        
    if (udata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        
        if (model->nx > 0) {
            /* initialise sensitivities, this can either be user provided or
             * come from the model definition */
            std::vector<double> sx0 = udata->getInitialSensitivityStates();
            if (sx0.empty()) {
                model->fsx0(fwd->getsxptr(), fwd->getxptr(), udata);
            } else {
                AmiVectorArray *sx = fwd->getsxptr();
                for (int ip = 0; ip < udata->nplist(); ip++) {
                    for (int ix = 0; ix < model->nx; ix++) {
                        sx->at(ix,ip) =
                            (realtype)sx0.at(ix + model->nx * ip);
                    }
                }
            }

            model->fsdx0();
            
            if (udata->sensmeth() == AMICI_SENSI_FSA) {
                
                /* Activate sensitivity calculations */
                sensInit1(fwd->getsxptr(), fwd->getsdxptr(), udata);
                /* Set sensitivity analysis optional inputs */
                std::vector<int> plist(udata->plist());
                std::vector<double> par;
                par.assign(udata->unp(),udata->unp()+udata->nplist());
                std::vector<double> pbar;
                pbar.assign(udata->getPbar(),udata->getPbar()+udata->nplist());
                AMISetSensParams(par.data(), pbar.data(), plist.data());
                AMISetSensErrCon(TRUE);
                AMISensEEtolerances();
            }
        }

        if (udata->sensmeth() == AMICI_SENSI_ASA) {
            if (model->nx > 0) {
                /* Allocate space for the adjoint computation */
                AMIAdjInit(udata->maxsteps, udata->getInterpolationType());
            }
        }
    }

    AMISetId(model);
    AMISetSuppressAlg(TRUE);
    /* calculate consistent DAE initial conditions (no effect for ODE) */
    if(udata->nt()>1)
        AMICalcIC(udata->t(1), fwd->getxptr(), fwd->getdxptr());
    } catch (...) {
        AMIFree();
        throw AmiException("setupAMI routine failed!");
    }
}

/**
 * setupAMIB initialises the AMI memory object for the backwards problem
 * @param bwd pointer to backward problem
 * @param udata pointer to the user data object
 * @param model pointer to the model object
 */
void Solver::setupAMIB(BackwardProblem *bwd, const UserData *udata, Model *model) {

    /* write initial conditions */
    std::vector<realtype> dJydx = bwd->getdJydx();
    AmiVector *xB = bwd->getxBptr();
    xB->reset();
    for (int ix = 0; ix < model->nxtrue; ++ix)
        for (int iJ = 0; iJ < model->nJ; ++iJ)
            xB->at(ix + iJ * model->nxtrue) +=
                dJydx.at(iJ + ( ix + (udata->nt() - 1)  * model->nx ) * model->nJ);
    bwd->getdxBptr()->reset();
    bwd->getxQBptr()->reset();

    /* create backward problem */
    if (udata->getLinearMultistepMethod() > 2 || udata->getLinearMultistepMethod() < 1) {
        throw AmiException("Illegal value for lmm!");
    }
    if (udata->getNonlinearSolverIteration() > 2 || udata->getNonlinearSolverIteration() < 1) {
        throw AmiException("Illegal value for iter!");
    }

    /* allocate memory for the backward problem */
    AMICreateB(udata->getLinearMultistepMethod(), udata->getNonlinearSolverIteration(), bwd->getwhichptr());

    /* initialise states */
    binit(bwd->getwhich(), bwd->getxBptr(), bwd->getdxBptr(), bwd->gett());

    /* specify integration tolerances for backward problem */
    AMISStolerancesB(bwd->getwhich(), RCONST(udata->rtol),
                              RCONST(udata->atol));

    /* Attach user data */
    AMISetUserDataB(bwd->getwhich(), model);

    /* Number of maximal internal steps */
    AMISetMaxNumStepsB(bwd->getwhich(), 100 * udata->maxsteps);
    
    setLinearSolverB(udata, model, bwd->getwhich());
    
    /* Initialise quadrature calculation */
    qbinit(bwd->getwhich(), bwd->getxQBptr());

    /* Enable Quadrature Error Control */
    AMISetQuadErrConB(bwd->getwhich(), TRUE);

    AMIQuadSStolerancesB(bwd->getwhich(), RCONST(udata->rtol),
                                  RCONST(udata->atol));

    AMISetStabLimDetB(bwd->getwhich(), udata->getStabilityLimitFlag());
}

/**
 * ErrHandlerFn extracts diagnosis information from solver memory block and
 * writes them into the return data object for the backward problem
 *
 * @param error_code error identifier
 * @param module name of the module in which the error occured
 * @param function name of the function in which the error occured @type
 * char
 * @param msg error message
 * @param eh_data unused input
 */
void Solver::wrapErrHandlerFn(int error_code, const char *module,
                              const char *function, char *msg, void *eh_data) {
    char buffer[250];
    char buffid[250];
    sprintf(buffer, "AMICI ERROR: in module %s in function %s : %s ", module,
            function, msg);
    switch (error_code) {
    case 99:
        sprintf(buffid, "AMICI:mex:%s:%s:WARNING", module, function);
        break;

    case -1:
        sprintf(buffid, "AMICI:mex:%s:%s:TOO_MUCH_WORK", module, function);
        break;

    case -2:
        sprintf(buffid, "AMICI:mex:%s:%s:TOO_MUCH_ACC", module, function);
        break;

    case -3:
        sprintf(buffid, "AMICI:mex:%s:%s:ERR_FAILURE", module, function);
        break;

    case -4:
        sprintf(buffid, "AMICI:mex:%s:%s:CONV_FAILURE", module, function);
        break;

    default:
        sprintf(buffid, "AMICI:mex:%s:%s:OTHER", module, function);
        break;
    }

    warnMsgIdAndTxt(buffid, buffer);
}

/**
 * getDiagnosis extracts diagnosis information from solver memory block and
 * writes them into the return data object
 *
 * @param it time-point index
 * @param rdata pointer to the return data object
 */
void Solver::getDiagnosis(const int it, ReturnData *rdata) {
    long int number;
    int order;

    if(solverWasCalled) {
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
    }
}

/**
 * getDiagnosisB extracts diagnosis information from solver memory block and
 * writes them into the return data object for the backward problem
 *
 * @param it time-point index
 * @param rdata pointer to the return data object
 * @param bwd pointer to backward problem
 */
void Solver::getDiagnosisB(const int it, ReturnData *rdata, const BackwardProblem *bwd) {
    long int number;

    void *ami_memB = AMIGetAdjBmem(ami_mem, bwd->getwhich());
    
    if(solverWasCalled && ami_memB) {
        AMIGetNumSteps(ami_memB, &number);
        rdata->numstepsB[it] = (double)number;
        
        AMIGetNumRhsEvals(ami_memB, &number);
        rdata->numrhsevalsB[it] = (double)number;
        
        AMIGetNumErrTestFails(ami_memB, &number);
        rdata->numerrtestfailsB[it] = (double)number;
        
        AMIGetNumNonlinSolvConvFails(ami_memB, &number);
        rdata->numnonlinsolvconvfailsB[it] = (double)number;
    }
}

/**
 * setLinearSolver sets the linear solver for the forward problem
 *
 * @param udata pointer to the user data object
 * @param model pointer to the model object
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
            AMIKLUSetOrdering(udata->getStateOrdering());
            break;
            
        default:
            throw AmiException("Invalid choice of solver!");
            
    }
}
    
    /**
     * setLinearSolverB sets the linear solver for the backward problem
     *
     * @param udata pointer to the user data object
     * @param model pointer to the model object
     * @param which index of the backward problem
     */
void Solver::setLinearSolverB(const UserData *udata, Model *model, const int which) {
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
             status = CVLapackDenseB(ami_mem, bwd->getwhich(), nx);
             if (status != AMICI_SUCCESS) return;
             
             status = SetDenseJacFnB(ami_mem, bwd->getwhich());
             if (status != AMICI_SUCCESS) return;
             #else*/
            throw AmiException("Solver currently not supported!");
            /* #endif*/
            break;
            
        case AMICI_LAPACKBAND:
            
            /* #if SUNDIALS_BLAS_LAPACK
             status = CVLapackBandB(ami_mem, bwd->getwhich(), nx, ubw, lbw);
             if (status != AMICI_SUCCESS) return;
             
             status = SetBandJacFnB(ami_mem, bwd->getwhich());
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
            AMIKLUSetOrderingB(which, udata->getStateOrdering());
            break;
            
        default:
            throw AmiException("Invalid local Solver!");
            break;
    }
}


} // namespace amici
