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
     * @param[out] status flag indicating success of execution @type int
     */
    
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
            return new NewtonSolverSparse(model, rdata, udata, tdata);
        
        default:
            errMsgIdAndTxt("AMICI:mex:solver","Invalid choice of solver!");
            *status = AMICI_ERROR_NEWTONSOLVER;
            return NULL;
    }
}

int NewtonSolver::getStep(int ntry, int nnewt, N_Vector delta) {
    
    int status = this->prepareLinearSystem(ntry, nnewt);
    
    N_VScale(-1.0, tdata->xdot, delta);
    if (status == AMICI_SUCCESS) {
        status = this->solveLinearSystem(delta);
    }
    
    return status;
}

int NewtonSolver::getSensis(int it) {
    
    N_Vector sx_ip = N_VNew_Serial(model->nx);
    realtype *x_tmp;
    int status = this->prepareLinearSystem(0, 0);
    
    if (status == AMICI_SUCCESS) {
        status = model->fdxdotdp(tdata->t, tdata->x, tdata->dx, tdata);
        if (status == AMICI_SUCCESS) {
            for (int ip=0; ip<udata->nplist; ip++) {
                
                /* Copy the content of dxdotdp to sx_ip */
                x_tmp = N_VGetArrayPointer(sx_ip);
                for (int ix=0; ix<model->nx; ix++) {
                    x_tmp[ix] = -tdata->dxdotdp[model->nx * ip + ix];
                }
                status = this->solveLinearSystem(sx_ip);
                
                /* Copy result to return data */
                if (status == AMICI_SUCCESS) {
                    for (int ix=0; ix<model->nx; ix++) {
                        rdata->sx[(ip * model->nx + ix) * rdata->nt + it] = x_tmp[ix];
                        N_VScale(1.0, sx_ip, tdata->sx[ip]);
                    }
                } else {
                    N_VDestroy_Serial(sx_ip);
                    return AMICI_ERROR_SS_SENSIS;
                }
            }
            
            /* Compute sy from sx */
            status = model->fsy(it, tdata, rdata);
            if (status != AMICI_SUCCESS) {
                N_VDestroy_Serial(sx_ip);
                return AMICI_ERROR_SS_SENSIS;
            }
        }
    }
    
    return status;
}

NewtonSolver::~NewtonSolver(){
}



/* ---------------------------------------------------------------------------------- */
/* - Dense linear solver ------------------------------------------------------------ */
/* ---------------------------------------------------------------------------------- */


/* derived class for dense linear solver */
NewtonSolverDense::NewtonSolverDense(Model *model, ReturnData *rdata, UserData *udata, TempData *tdata):NewtonSolver(model, rdata, udata, tdata) {
    pivots = NewLintArray(model->nx);
    tmp1 = N_VNew_Serial(model->nx);
    tmp2 = N_VNew_Serial(model->nx);
    tmp3 = N_VNew_Serial(model->nx);
}

int NewtonSolverDense::prepareLinearSystem(int ntry, int nnewt) {
    
    /* Get Jacobian */
    int status = model->fJ(model->nx, tdata->t, 0, tdata->x, tdata->dx, tdata->xdot, tdata->Jtmp, tdata, tmp1, tmp2, tmp3);
    
    if (status == AMICI_SUCCESS) {
        /* Compute factorization and pivoting */
        status = DenseGETRF(tdata->Jtmp, pivots);
    }
    
    return status;
}

int NewtonSolverDense::solveLinearSystem(N_Vector rhs) {
    
    /* Pass pointer to the linear solver */
    realtype *x_tmp = N_VGetArrayPointer(rhs);
    DenseGETRS(tdata->Jtmp, pivots, x_tmp);
    return AMICI_SUCCESS;
}

NewtonSolverDense::~NewtonSolverDense() {
    N_VDestroy_Serial(tmp1);
    N_VDestroy_Serial(tmp2);
    N_VDestroy_Serial(tmp3);
    DestroyArray(pivots);
}



/* ---------------------------------------------------------------------------------- */
/* - Sparse linear solver ----------------------------------------------------------- */
/* ---------------------------------------------------------------------------------- */

/* Derived class for sparse linear solver */
NewtonSolverSparse::NewtonSolverSparse(Model *model, ReturnData *rdata, UserData *udata, TempData *tdata):NewtonSolver(model, rdata, udata, tdata) {
    
    /* Initialize the KLU solver */
    klu_status = klu_defaults(&common);
    tmp1 = N_VNew_Serial(model->nx);
    tmp2 = N_VNew_Serial(model->nx);
    tmp3 = N_VNew_Serial(model->nx);
}

int NewtonSolverSparse::prepareLinearSystem(int ntry, int nnewt) {
    
    /* Check if KLU was initialized successfully */
    if (klu_status != 1)
        return AMICI_ERROR_NEWTON_LINSOLVER;
    
    /* Get sparse Jacobian */
    int status = model->fJSparse(tdata->t, tdata->x, tdata->xdot, tdata->J, tdata, tmp1, tmp2, tmp3);
    
    /* Get factorization of sparse Jacobian */
    if (status == AMICI_SUCCESS) {
        symbolic = klu_analyze(model->nx, (tdata->J)->indexptrs, (tdata->J)->indexvals, &common);
        if (symbolic != NULL) {
            numeric = klu_factor((tdata->J)->indexptrs, (tdata->J)->indexvals, (tdata->J)->data, symbolic, &common);
            if (numeric != NULL) {
                status = AMICI_SUCCESS;
            } else {
                status = AMICI_ERROR_NEWTON_LINSOLVER;
            }
        } else {
            status = AMICI_ERROR_NEWTON_LINSOLVER;
        }
    }
    
    return status;
}

int NewtonSolverSparse::solveLinearSystem(N_Vector rhs) {
    
    int status = AMICI_ERROR_NEWTON_LINSOLVER;
    realtype *x_tmp = N_VGetArrayPointer(rhs);
    
    /* Pass pointer to the linear solver */
    klu_status = klu_solve(symbolic, numeric, model->nx, 1, x_tmp, &common);
    if (klu_status == 1)
        status = AMICI_SUCCESS;

    return status;
}

NewtonSolverSparse::~NewtonSolverSparse() {
    N_VDestroy_Serial(tmp1);
    N_VDestroy_Serial(tmp2);
    N_VDestroy_Serial(tmp3);
    klu_free_symbolic(&symbolic, &common);
    klu_free_numeric(&numeric, &common);
}



/* ---------------------------------------------------------------------------------- */
/* - Iteartive linear solver -------------------------------------------------------- */
/* ---------------------------------------------------------------------------------- */

/* Derived class for iterative linear solver */
NewtonSolverIterative::NewtonSolverIterative(Model *model, ReturnData *rdata, UserData *udata, TempData *tdata):NewtonSolver(model, rdata, udata, tdata) {
}

int NewtonSolverIterative::getSensis(int it) {
    errMsgIdAndTxt("AMICI:mex:SPILS","Solver does not currently suppport sensitivity calculation!");
    return AMICI_ERROR_NEWTONSOLVER;
}

int NewtonSolverIterative::prepareLinearSystem(int ntry, int nnewt) {
    newton_try = ntry;
    i_newton = nnewt;
    return AMICI_SUCCESS;
}

int NewtonSolverIterative::solveLinearSystem(N_Vector rhs) {
    return SteadystateProblem::linsolveSPBCG(udata, rdata, tdata, newton_try, i_newton, rhs, model);;
}


NewtonSolverIterative::~NewtonSolverIterative() {};