#include "amici/solver.h"
#include "amici/exception.h"
#include "amici/forwardproblem.h"
#include "amici/backwardproblem.h"
#include "amici/model.h"
#include "amici/rdata.h"

#include <cmath>
#include <cstdio>
#include <cstring>

#include <sundials/sundials_spgmr.h>
#include <cvodes/cvodes_spils.h>

namespace amici {


/**
 * @brief setupAMIs initialises the ami memory object
 * @param fwd pointer to forward problem
 * @param model pointer to the model object
 */
void Solver::setupAMI(ForwardProblem *fwd, Model *model) {
    if (ami_mem) {
        /* This solver was used before. Not sure what we need to do to reuse the allocated
         *  memory, so just free here and then reallocate. */
        AMIFree();
    }

    model->initialize(fwd->getStatePointer(), fwd->getStateDerivativePointer());

    /* Create solver memory object */
    ami_mem = AMICreate(lmm, iter);
    if (ami_mem == NULL)
        throw AmiException("Failed to allocated solver memory!");
    try {
    /* Initialize AMIS solver*/
    init(fwd->getStatePointer(), fwd->getStateDerivativePointer(), model->t0());
    /* Specify integration tolerances */
    AMISStolerances(RCONST(rtol), RCONST(atol));
    /* Set optional inputs */
    AMISetErrHandlerFn();
    /* attaches userdata*/
    AMISetUserData(model);
    /* specify maximal number of steps */
    AMISetMaxNumSteps(maxsteps);
    /* activates stability limit detection */
    AMISetStabLimDet(stldet);
    
    rootInit(model->ne);
    
    initializeLinearSolver(model);
        
    if (sensi >= AMICI_SENSI_ORDER_FIRST) {
        
        if (model->nx > 0) {
            /* initialise sensitivities, this can either be user provided or
             * come from the model definition */
            std::vector<double> sx0 = model->getInitialStateSensitivities();
            if (sx0.empty()) {
                model->fsx0(fwd->getStateSensitivityPointer(), fwd->getStatePointer());
            } else {
                AmiVectorArray *sx = fwd->getStateSensitivityPointer();
                for (int ip = 0; ip < model->nplist(); ip++) {
                    for (int ix = 0; ix < model->nx; ix++) {
                        sx->at(ix,ip) =
                            (realtype)sx0.at(ix + model->nx * ip);
                    }
                }
            }

            model->fsdx0();
            
            auto plist = model->getParameterList();
            
            if (sensi_meth == AMICI_SENSI_FSA && plist.size() > 0) {
                /* Set sensitivity analysis optional inputs */
                auto par = model->getUnscaledParameters();

                /* Activate sensitivity calculations */
                sensInit1(fwd->getStateSensitivityPointer(), fwd->getStateDerivativeSensitivityPointer(), plist.size());
                AMISetSensParams(par.data(), nullptr, plist.data());
                std::vector<realtype> atols(plist.size(),atol);
                AMISensSStolerances( rtol, atols.data());
                AMISetSensErrCon(TRUE);
            }
        }

        if (sensi_meth == AMICI_SENSI_ASA) {
            if (model->nx > 0) {
                if (getNewtonPreequilibration()) {
                } else {
                    /* Allocate space for the adjoint computation */
                    AMIAdjInit(maxsteps, interpType);
                }
            }
        }
    }

    AMISetId(model);
    AMISetSuppressAlg(TRUE);
    /* calculate consistent DAE initial conditions (no effect for ODE) */
    if(model->nt()>1)
        AMICalcIC(model->t(1), fwd->getStatePointer(), fwd->getStateDerivativePointer());
    } catch (...) {
        AMIFree();
        throw AmiException("setupAMI routine failed!");
    }
}

/**
 * setupAMIB initialises the AMI memory object for the backwards problem
 * @param bwd pointer to backward problem
 * @param model pointer to the model object
 */
void Solver::setupAMIB(BackwardProblem *bwd, Model *model) {

    /* write initial conditions */
    std::vector<realtype> dJydx = bwd->getdJydx();
    AmiVector *xB = bwd->getxBptr();
    xB->reset();
    for (int ix = 0; ix < model->nxtrue; ++ix)
        for (int iJ = 0; iJ < model->nJ; ++iJ)
            xB->at(ix + iJ * model->nxtrue) +=
                dJydx.at(iJ + ( ix + (model->nt() - 1)  * model->nx ) * model->nJ);
    bwd->getdxBptr()->reset();
    bwd->getxQBptr()->reset();

    /* allocate memory for the backward problem */
    AMICreateB(lmm, iter, bwd->getwhichptr());

    /* initialise states */
    binit(bwd->getwhich(), bwd->getxBptr(), bwd->getdxBptr(), bwd->gett());

    /* specify integration tolerances for backward problem */
    AMISStolerancesB(bwd->getwhich(), RCONST(rtol), RCONST(atol));

    /* Attach user data */
    AMISetUserDataB(bwd->getwhich(), model);

    /* Number of maximal internal steps */
    AMISetMaxNumStepsB(bwd->getwhich(), (maxstepsB == 0) ? maxsteps * 100 : maxstepsB);
    
    initializeLinearSolverB(model, bwd->getwhich());
    
    /* Initialise quadrature calculation */
    qbinit(bwd->getwhich(), bwd->getxQBptr());
    
    double quad_rtol = isNaN(this->quad_rtol) ? rtol : this->quad_rtol;
    double quad_atol = isNaN(this->quad_atol) ? atol : this->quad_atol;
    
    /* Enable Quadrature Error Control */
    if (std::isinf(quad_atol) || std::isinf(quad_rtol)) {
        AMISetQuadErrConB(bwd->getwhich(), FALSE);
    } else {
        AMISetQuadErrConB(bwd->getwhich(), TRUE);
    }

    AMIQuadSStolerancesB(bwd->getwhich(), RCONST(quad_rtol),
                         RCONST(quad_atol));

    AMISetStabLimDetB(bwd->getwhich(), stldet);
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
        rdata->numsteps[it] = number;
        
        AMIGetNumRhsEvals(ami_mem, &number);
        rdata->numrhsevals[it] = number;
        
        AMIGetNumErrTestFails(ami_mem, &number);
        rdata->numerrtestfails[it] = number;
        
        AMIGetNumNonlinSolvConvFails(ami_mem, &number);
        rdata->numnonlinsolvconvfails[it] = number;
        
        AMIGetLastOrder(ami_mem, &order);
        rdata->order[it] = order;
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
 * @param model pointer to the model object
 */
void Solver::initializeLinearSolver(Model *model) {
    /* Attach linear solver module */

    switch (linsol) {

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
            AMIKLUSetOrdering(getStateOrdering());
            break;
            
        default:
            throw AmiException("Invalid choice of solver!");
            
    }
}
    
    /**
     * Sets the linear solver for the backward problem
     *
     * @param model pointer to the model object
     * @param which index of the backward problem
     */
void Solver::initializeLinearSolverB(Model *model, const int which) {
    switch (linsol) {
            
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
            AMIKLUSetOrderingB(which, getStateOrdering());
            break;
            
        default:
            throw AmiException("Invalid local Solver!");
            break;
    }
}

bool operator ==(const Solver &a, const Solver &b)
{
    if (typeid(a) != typeid(b))
            return false;

    return (a.sensi_meth == b.sensi_meth)
            && (a.sensi_meth == b.sensi_meth)
            && (a.interpType == b.interpType)
            && (a.lmm == b.lmm)
            && (a.iter == b.iter)
            && (a.stldet == b.stldet)
            && (a.ordering == b.ordering)
            && (a.newton_maxsteps == b.newton_maxsteps)
            && (a.newton_maxlinsteps == b.newton_maxlinsteps)
            && (a.newton_preeq == b.newton_preeq)
            && (a.ism == b.ism)
            && (a.linsol == b.linsol)
            && (a.atol == b.atol)
            && (a.rtol == b.rtol)
            && (a.maxsteps == b.maxsteps)
            && (a.quad_atol == b.quad_atol)
            && (a.quad_rtol == b.quad_rtol)
            && (a.maxstepsB == b.maxstepsB)
            && (a.sensi == b.sensi);
}


} // namespace amici
