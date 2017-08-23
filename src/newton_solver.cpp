#include "include/amici_defines.h"
#include "include/newton_solver.h"
#include "include/steadystateproblem.h"
#include "include/forwardproblem.h"
#include "include/udata.h"
#include "include/rdata.h"
#include "include/tdata.h"
#include "include/edata.h"
#include "include/amici_model.h"
#include "include/amici_solver.h"
#include <cstring>
#include <ctime>
#include "include/newton_solver.h"

NewtonSolver::NewtonSolver(Model *model, ReturnData *rdata, UserData *udata, TempData *tdata):model(model), rdata(rdata), udata(udata), tdata(tdata) {
}

NewtonSolver *NewtonSolver::getSolver(int linsolType, Model *model, ReturnData *rdata, UserData *udata, TempData *tdata, int *status) {
    /**
     * getNewtonStep computes the Newton Step by solving the linear system
     *
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] tdata pointer to the temporary data struct @type TempData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @return int status flag indicating success of execution @type int
     */
    
    NewtonSolver *newtonSolver;
    
    switch (linsolType) {
            /* DIRECT SOLVERS */
            
        case AMICI_DENSE:
            return new NewtonSolverDense(model, rdata, udata, tdata);
            
        case AMICI_BAND:
            errMsgIdAndTxt("AMICI:mex:dense","Solver currently not supported!");
            *status = AMICI_ERROR_NEWTONSOLVER;
            return NULL;
            
        case AMICI_LAPACKDENSE:
            errMsgIdAndTxt("AMICI:mex:lapack","Solver currently not supported!");
            *status = AMICI_ERROR_NEWTONSOLVER;
            return NULL;
            
        case AMICI_LAPACKBAND:
            errMsgIdAndTxt("AMICI:mex:lapack","Solver currently not supported!");
            *status = AMICI_ERROR_NEWTONSOLVER;
            return NULL;
            
        case AMICI_DIAG:
            errMsgIdAndTxt("AMICI:mex:dense","Solver currently not supported!");
            *status = AMICI_ERROR_NEWTONSOLVER;
            return NULL;
            
            /* ITERATIVE SOLVERS */
            
        case AMICI_SPGMR:
            errMsgIdAndTxt("AMICI:mex:spils","Solver currently not supported!");
            *status = AMICI_ERROR_NEWTONSOLVER;
            return NULL;
            
        case AMICI_SPBCG:
            return new NewtonSolverIterative(model, rdata, udata, tdata);
            
        case AMICI_SPTFQMR:
            errMsgIdAndTxt("AMICI:mex:spils","Solver currently not supported!");
            *status = AMICI_ERROR_NEWTONSOLVER;
            return NULL;
            
            /* SPARSE SOLVERS */
            
        case AMICI_KLU:
            newtonSolver = new NewtonSolverSparse(model, rdata, udata, tdata);
            newtonSolver->solverStatus = klu_defaults (&(((NewtonSolverSparse *)newtonSolver)->common));
            if (newtonSolver->solverStatus == 1) {
                return newtonSolver;
            } else {
                *status = AMICI_ERROR_NEWTONSOLVER;
                errMsgIdAndTxt("AMICI:mex:solver","Initialization of KLU solver feiled!");
                return NULL;
            }
        
        
        default:
            errMsgIdAndTxt("AMICI:mex:solver","Invalid choice of solver!");
            *status = AMICI_ERROR_NEWTONSOLVER;
            return NULL;
    }
}

NewtonSolver::~NewtonSolver(){
}



/* derived class for dense linear solver */
NewtonSolverDense::NewtonSolverDense(Model *model, ReturnData *rdata, UserData *udata, TempData *tdata):NewtonSolver(model, rdata, udata, tdata) {
    pivots = NewLintArray(model->nx);
    tmp1 = N_VNew_Serial(model->nx);
    tmp2 = N_VNew_Serial(model->nx);
    tmp3 = N_VNew_Serial(model->nx);
}

int NewtonSolverDense::getStep(int ntry, int nnewt, N_Vector delta) {
    
    int status = model->fJ(model->nx, tdata->t, 0, tdata->x, tdata->dx, tdata->xdot, tdata->Jtmp, tdata, tmp1, tmp2, tmp3);
    realtype *x_tmp;
    
    if (status == AMICI_SUCCESS) {
        /* Compute factorization and pivoting */
        status = DenseGETRF(tdata->Jtmp, pivots);
        if (status == AMICI_SUCCESS) {
            /* Solve linear system by elimiation using pivotiing */
            model->fxdot(tdata->t, tdata->x, tdata->dx, delta, tdata);
            N_VScale(-1.0, delta, delta);
            x_tmp = N_VGetArrayPointer(delta);
            DenseGETRS(tdata->Jtmp, pivots, x_tmp);
        } else {
            return(AMICI_ERROR_NEWTON_LINSOLVER);
        }
    } else {
        return(AMICI_ERROR_NEWTON_LINSOLVER);
    }
    
    return status;
}

NewtonSolverDense::~NewtonSolverDense() {
    N_VDestroy_Serial(tmp1);
    N_VDestroy_Serial(tmp2);
    N_VDestroy_Serial(tmp3);
    DestroyArray(pivots);
}



/* derived class for sparse linear solver */
NewtonSolverSparse::NewtonSolverSparse(Model *model, ReturnData *rdata, UserData *udata, TempData *tdata):NewtonSolver(model, rdata, udata, tdata) {
    
    /* Initialize the KLU solver */
    tmp1 = N_VNew_Serial(model->nx);
    tmp2 = N_VNew_Serial(model->nx);
    tmp3 = N_VNew_Serial(model->nx);
}

int NewtonSolverSparse::getStep(int ntry, int nnewt, N_Vector delta) {
    
    realtype *x_tmp;
    int status = model->fJSparse(tdata->t, tdata->x, tdata->xdot, tdata->J, tdata, tmp1, tmp2, tmp3);
    
    symbolic = klu_analyze (model->nx, (tdata->J)->indexptrs, (tdata->J)->indexvals, &common) ;
    numeric = klu_factor((tdata->J)->indexptrs, (tdata->J)->indexvals, (tdata->J)->data, symbolic, &common) ;
    N_VScale(-1.0, tdata->xdot, delta);
    x_tmp = N_VGetArrayPointer(delta);
    klu_solve(symbolic, numeric, model->nx, 1, x_tmp, &common);
    
    return status;
}

NewtonSolverSparse::~NewtonSolverSparse() {
    N_VDestroy_Serial(tmp1);
    N_VDestroy_Serial(tmp2);
    N_VDestroy_Serial(tmp3);
    klu_free_symbolic(&symbolic, &common);
    klu_free_numeric(&numeric, &common);
}


/* derived class for iterative linear solver */
NewtonSolverIterative::NewtonSolverIterative(Model *model, ReturnData *rdata, UserData *udata, TempData *tdata):NewtonSolver(model, rdata, udata, tdata) {
}

int NewtonSolverIterative::getStep(int ntry, int nnewt, N_Vector delta) {
    return SteadystateProblem::linsolveSPBCG(udata, rdata, tdata, ntry, nnewt, delta, model);
}

NewtonSolverIterative::~NewtonSolverIterative() {};