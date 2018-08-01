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

extern msgIdAndTxtFp warnMsgIdAndTxt;

/**
 * @brief setupAMIs initialises the ami memory object
 * @param fwd pointer to forward problem
 * @param model pointer to the model object
 */
void Solver::setup(ForwardProblem *fwd, Model *model) {
    model->initialize(fwd->getStatePointer(), fwd->getStateDerivativePointer());

    /* Create solver memory object */
    create(lmm, iter);
    if (!solverMemory)
        throw AmiException("Failed to allocated solver memory!");

    /* Initialize AMIS solver*/
    init(fwd->getStatePointer(), fwd->getStateDerivativePointer(), model->t0());
    /* Specify integration tolerances */
    setSStolerances(RCONST(rtol), RCONST(atol));
    /* Set optional inputs */
    setErrHandlerFn();
    /* attaches userdata*/
    setUserData(model);
    /* specify maximal number of steps */
    setMaxNumSteps(maxsteps);
    /* activates stability limit detection */
    setStabLimDet(stldet);
    
    rootInit(model->ne);
    
    initializeLinearSolver(model);

    if (sensi >= AMICI_SENSI_ORDER_FIRST && model->nx > 0) {
        /* initialise sensitivities, this can either be user provided or
         * come from the model definition */
        auto sx0 = model->getInitialStateSensitivities();
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

        if (sensi_meth == AMICI_SENSI_FSA && !plist.empty()) {
            /* Set sensitivity analysis optional inputs */
            auto par = model->getUnscaledParameters();

            /* Activate sensitivity calculations */
            sensInit1(fwd->getStateSensitivityPointer(), fwd->getStateDerivativeSensitivityPointer(), plist.size());
            setSensParams(par.data(), nullptr, plist.data());
            std::vector<realtype> atols(plist.size(),atol);
            setSensSStolerances( rtol, atols.data());
            setSensErrCon(true);
        } else if (sensi_meth == AMICI_SENSI_ASA) {
            /* Allocate space for the adjoint computation */
            adjInit(maxsteps, interpType);
        }
    }

    setId(model);
    setSuppressAlg(true);
    /* calculate consistent DAE initial conditions (no effect for ODE) */
    if(model->nt()>1)
        calcIC(model->t(1), fwd->getStatePointer(), fwd->getStateDerivativePointer());
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
    createB(lmm, iter, bwd->getwhichptr());

    /* initialise states */
    binit(bwd->getwhich(), bwd->getxBptr(), bwd->getdxBptr(), bwd->gett());

    /* specify integration tolerances for backward problem */
    setSStolerancesB(bwd->getwhich(), RCONST(rtol), RCONST(atol));

    /* Attach user data */
    setUserDataB(bwd->getwhich(), model);

    /* Number of maximal internal steps */
    setMaxNumStepsB(bwd->getwhich(), (maxstepsB == 0) ? maxsteps * 100 : maxstepsB);
    
    initializeLinearSolverB(model, bwd->getwhich());
    
    /* Initialise quadrature calculation */
    qbinit(bwd->getwhich(), bwd->getxQBptr());
    
    double quad_rtol = isNaN(this->quad_rtol) ? rtol : this->quad_rtol;
    double quad_atol = isNaN(this->quad_atol) ? atol : this->quad_atol;
    
    /* Enable Quadrature Error Control */
    setQuadErrConB(bwd->getwhich(), !std::isinf(quad_atol) && !std::isinf(quad_rtol));

    quadSStolerancesB(bwd->getwhich(), RCONST(quad_rtol),
                         RCONST(quad_atol));

    setStabLimDetB(bwd->getwhich(), stldet);
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
                              const char *function, char *msg, void * /*eh_data*/) {
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

    if(solverWasCalled) {
        getNumSteps(solverMemory.get(), &number);
        rdata->numsteps[it] = number;
        
        getNumRhsEvals(solverMemory.get(), &number);
        rdata->numrhsevals[it] = number;
        
        getNumErrTestFails(solverMemory.get(), &number);
        rdata->numerrtestfails[it] = number;
        
        getNumNonlinSolvConvFails(solverMemory.get(), &number);
        rdata->numnonlinsolvconvfails[it] = number;
        
        getLastOrder(solverMemory.get(), &rdata->order[it]);
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

    void *ami_memB = getAdjBmem(solverMemory.get(), bwd->getwhich());
    
    if(solverWasCalled && ami_memB) {
        getNumSteps(ami_memB, &number);
        rdata->numstepsB[it] = (double)number;
        
        getNumRhsEvals(ami_memB, &number);
        rdata->numrhsevalsB[it] = (double)number;
        
        getNumErrTestFails(ami_memB, &number);
        rdata->numerrtestfailsB[it] = (double)number;
        
        getNumNonlinSolvConvFails(ami_memB, &number);
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
            dense(model->nx);
            setDenseJacFn();
            break;
            
        case AMICI_BAND:
            band(model->nx, model->ubw, model->lbw);
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
            diag();
            break;
            
            
            /* ITERATIVE SOLVERS */
            
        case AMICI_SPGMR:
            spgmr(PREC_NONE, CVSPILS_MAXL);
            setJacTimesVecFn();
            break;
            
        case AMICI_SPBCG:
            spbcg(PREC_NONE, CVSPILS_MAXL);
            setJacTimesVecFn();
            break;
            
        case AMICI_SPTFQMR:
            sptfqmr(PREC_NONE, CVSPILS_MAXL);
            setJacTimesVecFn();
            break;
            
            /* SPARSE SOLVERS */
            
        case AMICI_KLU:
            klu(model->nx, model->nnz, CSC_MAT);
            setSparseJacFn();
            kluSetOrdering(getStateOrdering());
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
            denseB(which, model->nx);
            setDenseJacFnB(which);
            break;
            
        case AMICI_BAND:
            bandB(which, model->nx, model->ubw, model->lbw);
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
            
        case AMICI_LAPACKBAND:
            
            /* #if SUNDIALS_BLAS_LAPACK
             status = CVLapackBandB(ami_mem, bwd->getwhich(), nx, ubw, lbw);
             if (status != AMICI_SUCCESS) return;
             
             status = SetBandJacFnB(ami_mem, bwd->getwhich());
             if (status != AMICI_SUCCESS) return;
             #else*/
            throw AmiException("Solver currently not supported!");
            /* #endif*/
            
        case AMICI_DIAG:
            diagB(which);
            setDenseJacFnB(which);
            break;
            
            /* ITERATIVE SOLVERS */
            
        case AMICI_SPGMR:
            spgmrB(which, PREC_NONE, CVSPILS_MAXL);
            setJacTimesVecFnB(which);
            break;
            
        case AMICI_SPBCG:
            spbcgB(which, PREC_NONE, CVSPILS_MAXL);
            setJacTimesVecFnB(which);
            break;
            
        case AMICI_SPTFQMR:
            sptfqmrB(which, PREC_NONE, CVSPILS_MAXL);
            setJacTimesVecFnB(which);
            break;
            
            /* SPARSE SOLVERS */
            
        case AMICI_KLU:
            kluB(which, model->nx, model->nnz, CSC_MAT);
            setSparseJacFnB(which);
            kluSetOrderingB(which, getStateOrdering());
            break;
            
        default:
            throw AmiException("Invalid local Solver!");
    }
}

bool operator ==(const Solver &a, const Solver &b)
{
    if (typeid(a) != typeid(b))
            return false;

    return (a.sensi_meth == b.sensi_meth)
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
