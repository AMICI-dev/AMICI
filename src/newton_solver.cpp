#include "include/newton_solver.h"
#include "include/amici_defines.h"
#include "include/amici_model.h"
#include "include/amici_solver.h"
#include "include/edata.h"
#include "include/forwardproblem.h"
#include "include/newton_solver.h"
#include "include/rdata.h"
#include "include/steadystateproblem.h"
#include "include/tdata.h"
#include "include/udata.h"
#include <cstring>
#include <ctime>

namespace amici {

NewtonSolver::NewtonSolver(Model *model, ReturnData *rdata,
                           const UserData *udata, TempData *tdata)
    : model(model), rdata(rdata), udata(udata), tdata(tdata) {
    /**
     * default constructor, initializes all members with the provided objects
     *
     * @param[in] model pointer to the AMICI model object @type Model
     * @param[in] rdata pointer to the return data object @type ReturnData
     * @param[in] udata pointer to the user data object @type UserData
     * @param[in] tdata pointer to the temporary data object @type TempData
     */
    sx_ip = N_VNew_Serial(model->nx);
}

/* ----------------------------------------------------------------------------------
 */

NewtonSolver *NewtonSolver::getSolver(int linsolType, Model *model,
                                      ReturnData *rdata, const UserData *udata,
                                      TempData *tdata) {
    /**
     * Tries to determine the steady state of the ODE system by a Newton
     * solver, uses forward intergration, if the Newton solver fails,
     * restarts Newton solver, if integration fails.
     * Computes steady state sensitivities
     *
     * @param[in] linsolType integer indicating which linear solver to use
     * @param[in] model pointer to the AMICI model object @type Model
     * @param[in] udata pointer to the user data object @type UserData
     * @param[in,out] tdata pointer to the temporary data object @type TempData
     * @param[out] rdata pointer to the return data object @type ReturnData
     * @return solver NewtonSolver according to the specified linsolType
     */

    switch (linsolType) {

    /* DIRECT SOLVERS */
    case AMICI_DENSE:
        return new NewtonSolverDense(model, rdata, udata, tdata);

    case AMICI_BAND:
        throw NewtonFailure("Solver currently not supported!");

    case AMICI_LAPACKDENSE:
        throw NewtonFailure("Solver currently not supported!");

    case AMICI_LAPACKBAND:
        throw NewtonFailure("Solver currently not supported!");

    case AMICI_DIAG:
        throw NewtonFailure("Solver currently not supported!");

    /* ITERATIVE SOLVERS */
    case AMICI_SPGMR:
        throw NewtonFailure("Solver currently not supported!");

    case AMICI_SPBCG:
        return new NewtonSolverIterative(model, rdata, udata, tdata);

    case AMICI_SPTFQMR:
        throw NewtonFailure("Solver currently not supported!");

    /* SPARSE SOLVERS */
    case AMICI_KLU:
        return new NewtonSolverSparse(model, rdata, udata, tdata);

    default:
        throw NewtonFailure("Invalid Choice of Solver!");
    }
}

/* ----------------------------------------------------------------------------------
 */

void NewtonSolver::getStep(int ntry, int nnewt, N_Vector delta) {
    /**
     * Computes the solution of one Newton iteration
     *
     * @param[in] ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param[in] nnewt integer number of current Newton step
     * @param[in,out] delta containing the RHS of the linear system, will be
     * overwritten by solution to the linear system @type N_Vector
     */

    this->prepareLinearSystem(ntry, nnewt);

    N_VScale(-1.0, tdata->xdot, delta);
    this->solveLinearSystem(delta);
}

/* ----------------------------------------------------------------------------------
 */

void NewtonSolver::getSensis(int it) {
    /**
     * Computes steady state sensitivities
     *
     * @param[in] it integer index of current time step
     * @return stats integer flag indicating success of the method
     */
    realtype *sx_tmp;
    realtype *sx0_tmp;
    this->prepareLinearSystem(0, -1);

    model->fdxdotdp(tdata->t, tdata->x, tdata->dx, tdata);
    for (int ip = 0; ip < udata->nplist; ip++) {
        
        /* Copy the content of dxdotdp to sx_ip */
        sx_tmp = N_VGetArrayPointer(sx_ip);
        for (int ix = 0; ix < model->nx; ix++) {
            sx_tmp[ix] = -tdata->dxdotdp[model->nx * ip + ix];
        }
        this->solveLinearSystem(sx_ip);
        
        /* Copy result to return data */
        if (it == AMICI_PREEQUILIBRATE) {
            /* Case of preequlibration */
            sx0_tmp = N_VGetArrayPointer(tdata->sx[ip]);
            for (int ix = 0; ix < model->nx; ix++) {
                rdata->sx0[ip * model->nx + ix] = sx_tmp[ix];
                sx0_tmp[ix] = sx_tmp[ix];
            }
        } else {
            /* Classical steady state computation */
            for (int ix = 0; ix < model->nx; ix++) {
                rdata->sx[(ip * model->nx + ix) * rdata->nt + it] =
                sx_tmp[ix];
            }
        }
    }
}

/* ----------------------------------------------------------------------------------
 */
/* - Dense linear solver
 * ------------------------------------------------------------ */
/* ----------------------------------------------------------------------------------
 */

/* Derived class for dense linear solver */
NewtonSolverDense::NewtonSolverDense(Model *model, ReturnData *rdata,
                                     const UserData *udata, TempData *tdata)
    : NewtonSolver(model, rdata, udata, tdata) {
    /**
     * default constructor, initializes all members with the provided objects
     * and
     * initializes temporary storage objects
     *
     * @param[in] model pointer to the AMICI model object @type Model
     * @param[in] rdata pointer to the return data object @type ReturnData
     * @param[in] udata pointer to the user data object @type UserData
     * @param[in] tdata pointer to the temporary data object @type TempData
     */
    pivots = NewLintArray(model->nx);
    tmp1 = N_VNew_Serial(model->nx);
    tmp2 = N_VNew_Serial(model->nx);
    tmp3 = N_VNew_Serial(model->nx);
}

/* ----------------------------------------------------------------------------------
 */

void NewtonSolverDense::prepareLinearSystem(int ntry, int nnewt) {
    /**
     * Writes the Jacobian for the Newton iteration and passes it to the linear
     * solver
     *
     * @param[in] ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param[in] nnewt integer number of current Newton step
     */

    /* Get Jacobian */
    model->fJ(model->nx, tdata->t, 0, tdata->x, tdata->dx,
                           tdata->xdot, tdata->Jtmp, tdata, tmp1, tmp2, tmp3);
    int status = DenseGETRF(tdata->Jtmp, pivots);
    if(status != 0)
        throw NewtonFailure("Dense factorization failed!");
}

/* ----------------------------------------------------------------------------------
 */

void NewtonSolverDense::solveLinearSystem(N_Vector rhs) {
    /**
     * Solves the linear system for the Newton step
     *
     * @param[in,out] rhs containing the RHS of the linear system, will be
     * overwritten by solution to the linear system @type N_Vector
     */

    /* Pass pointer to the linear solver */
    realtype *x_tmp = N_VGetArrayPointer(rhs);
    DenseGETRS(tdata->Jtmp, pivots, x_tmp);
}

/* ----------------------------------------------------------------------------------
 */

NewtonSolverDense::~NewtonSolverDense() {
    N_VDestroy_Serial(tmp1);
    N_VDestroy_Serial(tmp2);
    N_VDestroy_Serial(tmp3);
    DestroyArray(pivots);
}

/* ----------------------------------------------------------------------------------
 */
/* - Sparse linear solver
 * ----------------------------------------------------------- */
/* ----------------------------------------------------------------------------------
 */

/* Derived class for sparse linear solver */
NewtonSolverSparse::NewtonSolverSparse(Model *model, ReturnData *rdata,
                                       const UserData *udata, TempData *tdata)
    : NewtonSolver(model, rdata, udata, tdata) {
    /**
     * default constructor, initializes all members with the provided objects,
     * initializes temporary storage objects and the klu solver
     *
     * @param[in] model pointer to the AMICI model object @type Model
     * @param[in] rdata pointer to the return data object @type ReturnData
     * @param[in] udata pointer to the user data object @type UserData
     * @param[in] tdata pointer to the temporary data object @type TempData
     */

    /* Initialize the KLU solver */
    klu_status = klu_defaults(&common);
    tmp1 = N_VNew_Serial(model->nx);
    tmp2 = N_VNew_Serial(model->nx);
    tmp3 = N_VNew_Serial(model->nx);
}

/* ----------------------------------------------------------------------------------
 */

void NewtonSolverSparse::prepareLinearSystem(int ntry, int nnewt) {
    /**
     * Writes the Jacobian for the Newton iteration and passes it to the linear
     * solver
     *
     * @param[in] ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param[in] nnewt integer number of current Newton step
     */

    /* Check if KLU was initialized successfully */
    if (klu_status != 1)
        throw NewtonFailure("KLU was not initialized!");

    /* Get sparse Jacobian */
    model->fJSparse(tdata->t, 0.0, tdata->x, tdata->dx, tdata->xdot, tdata->J,
                                 tdata, tmp1, tmp2, tmp3);

    /* Get factorization of sparse Jacobian */
    if(symbolic) /* if symbolic was already created free first to avoid memory leak */
        klu_free_symbolic(&symbolic, &common);
    symbolic = klu_analyze(model->nx, (tdata->J)->indexptrs,
                           (tdata->J)->indexvals, &common);
    if (symbolic) {
        if(numeric) /* if numeric was already created free first to avoid memory leak */
            klu_free_numeric(&numeric, &common);
        numeric = klu_factor((tdata->J)->indexptrs, (tdata->J)->indexvals,
                             (tdata->J)->data, symbolic, &common);
        if (numeric) {
            return;
        } else {
            throw NewtonFailure("KLU factorization failed!");
        }
    } else {
        throw NewtonFailure("KLU symbolic analisis failed!");
    }
}

/* ----------------------------------------------------------------------------------
 */

void NewtonSolverSparse::solveLinearSystem(N_Vector rhs) {
    /**
     * Solves the linear system for the Newton step
     *
     * @param[in] rhs containing the RHS of the linear system,will be
     * overwritten by solution to the linear system @type N_Vector
     * @return stats integer flag indicating success of the method
     */
    realtype *x_tmp = N_VGetArrayPointer(rhs);

    /* Pass pointer to the linear solver */
    klu_status = klu_solve(symbolic, numeric, model->nx, 1, x_tmp, &common);
    if (klu_status != 1)
        throw NewtonFailure("KLU solver failed");
}

/* ----------------------------------------------------------------------------------
 */

NewtonSolverSparse::~NewtonSolverSparse() {
    N_VDestroy_Serial(tmp1);
    N_VDestroy_Serial(tmp2);
    N_VDestroy_Serial(tmp3);
    if(symbolic)
        klu_free_symbolic(&symbolic, &common);
    if(numeric)
        klu_free_numeric(&numeric, &common);
}

/* ----------------------------------------------------------------------------------
 */
/* - Iterative linear solver
 * -------------------------------------------------------- */
/* ----------------------------------------------------------------------------------
 */

/* Derived class for iterative linear solver */
NewtonSolverIterative::NewtonSolverIterative(Model *model, ReturnData *rdata,
                                             const UserData *udata,
                                             TempData *tdata)
    : NewtonSolver(model, rdata, udata, tdata) {
    /**
     * default constructor, initializes all members with the provided objects
     * @param[in] model pointer to the AMICI model object @type Model
     * @param[in] rdata pointer to the return data object @type ReturnData
     * @param[in] udata pointer to the user data object @type UserData
     * @param[in] tdata pointer to the temporary data object @type TempData
     */
        
        ns_p = N_VNew_Serial(model->nx);
        ns_h = N_VNew_Serial(model->nx);
        ns_t = N_VNew_Serial(model->nx);
        ns_s = N_VNew_Serial(model->nx);
        ns_r = N_VNew_Serial(model->nx);
        ns_rt = N_VNew_Serial(model->nx);
        ns_v = N_VNew_Serial(model->nx);
        ns_Jv = N_VNew_Serial(model->nx);
        ns_tmp = N_VNew_Serial(model->nx);
        ns_Jdiag = N_VNew_Serial(model->nx);
        
        N_VScale(-1.0, tdata->xdot, tdata->xdot);
}

/* ----------------------------------------------------------------------------------
 */

void NewtonSolverIterative::prepareLinearSystem(int ntry, int nnewt) {
    /**
     * Writes the Jacobian for the Newton iteration and passes it to the linear
     * solver.
     * Also wraps around getSensis for iterative linear solver.
     *
     * @param[in] ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param[in] nnewt integer number of current Newton step
     */

    newton_try = ntry;
    i_newton = nnewt;
    if (nnewt == -1) {
        throw NewtonFailure("No setup necessary");
    }
}

/* ----------------------------------------------------------------------------------
 */

void NewtonSolverIterative::solveLinearSystem(N_Vector rhs) {
    /**
     * Solves the linear system for the Newton step by passing it to
     * linsolveSPBCG
     *
     * @param[in,out] rhs containing the RHS of the linear system, will be
     * overwritten by solution to the linear system @type N_Vector
     */

    linsolveSPBCG(newton_try, i_newton, rhs);
}
    
    
void NewtonSolverIterative::linsolveSPBCG(int ntry,int nnewt, N_Vector ns_delta) {
    /**
     * Iterative linear solver created from SPILS BiCG-Stab.
     * Solves the linear system within each Newton step if iterative solver is
     * chosen.
     *
     * @param[in] ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param[in] nnewt integer number of current Newton step
     * @param[in] ns_delta ???
     */
    
    double rho;
    double rho1;
    double alpha;
    double beta;
    double omega;
    double res;
    
    // Get the diagonal of the Jacobian for preconditioning
    model->fJDiag(tdata->t, ns_Jdiag, 0.0, tdata->x, tdata->dx, tdata);
    
    // Ensure positivity of entries in ns_Jdiag
    N_VConst(1.0, ns_p);
    N_VAbs(ns_Jdiag, ns_tmp);
    N_VCompare(1e-15, ns_tmp, ns_tmp);
    N_VLinearSum(-1.0, ns_tmp, 1.0, ns_p, ns_tmp);
    N_VLinearSum(1.0, ns_Jdiag, 1.0, ns_tmp, ns_Jdiag);
    
    // Initialize for linear solve
    N_VConst(0.0, ns_p);
    N_VConst(0.0, ns_v);
    N_VConst(0.0, ns_delta);
    N_VConst(0.0, ns_tmp);
    rho = 1.0;
    omega = 1.0;
    alpha = 1.0;
    
    // can be set to 0 at the moment
    model->fJv(tdata->t, tdata->x, tdata->dx, tdata->xdot, ns_delta, ns_Jv, 0.0,  tdata, ns_tmp, NULL);
    
    // ns_r = xdot - ns_Jv;
    N_VLinearSum(-1.0, ns_Jv, 1.0, tdata->xdot, ns_r);
    N_VDiv(ns_r, ns_Jdiag, ns_r);
    res = sqrt(N_VDotProd(ns_r, ns_r));
    N_VScale(1.0, ns_r, ns_rt);
    
    for (int i_linstep = 0; i_linstep < udata->newton_maxlinsteps;
         i_linstep++) {
        // Compute factors
        rho1 = rho;
        rho = N_VDotProd(ns_rt, ns_r);
        beta = rho * alpha / (rho1 * omega);
        
        // ns_p = ns_r + beta * (ns_p - omega * ns_v);
        N_VLinearSum(1.0, ns_p, -omega, ns_v, ns_p);
        N_VLinearSum(1.0, ns_r, beta, ns_p, ns_p);
        
        // ns_v = J * ns_p
        model->fJv(tdata->t, tdata->x, tdata->dx, tdata->xdot, ns_p, ns_v, 0.0,  tdata, ns_tmp, NULL);
        N_VDiv(ns_v, ns_Jdiag, ns_v);
        
        // Compute factor
        alpha = rho / N_VDotProd(ns_rt, ns_v);
        
        // ns_h = ns_delta + alpha * ns_p;
        N_VLinearSum(1.0, ns_delta, alpha, ns_p, ns_h);
        // ns_s = ns_r - alpha * ns_v;
        N_VLinearSum(1.0, ns_r, -alpha, ns_v, ns_s);
        
        // ns_t = J * ns_s
        model->fJv(tdata->t, tdata->x, tdata->dx, tdata->xdot, ns_s, ns_t, 0.0,  tdata, ns_tmp, NULL);
        N_VDiv(ns_t, ns_Jdiag, ns_t);
        
        // Compute factor
        omega = N_VDotProd(ns_t, ns_s) / N_VDotProd(ns_t, ns_t);
        
        // ns_delta = ns_h + omega * ns_s;
        N_VLinearSum(1.0, ns_h, omega, ns_s, ns_delta);
        // ns_r = ns_s - omega * ns_t;
        N_VLinearSum(1.0, ns_s, -omega, ns_t, ns_r);
        
        // Compute the (unscaled) residual
        N_VProd(ns_r, ns_Jdiag, ns_r);
        res = sqrt(N_VDotProd(ns_r, ns_r));
        
        // Test convergence
        if (res < udata->atol) {
            // Write number of steps needed
            rdata->newton_numlinsteps[(ntry - 1) * udata->newton_maxsteps +
                                      nnewt] = i_linstep + 1;
            
            // Return success
            return;
        }
        
        // Scale back
        N_VDiv(ns_r, ns_Jdiag, ns_r);
    }
    throw NewtonFailure("SPBCG solver failed to converge");
}
    

/* ----------------------------------------------------------------------------------
 */

NewtonSolverIterative::~NewtonSolverIterative(){
    N_VDestroy_Serial(ns_p);
    N_VDestroy_Serial(ns_h);
    N_VDestroy_Serial(ns_t);
    N_VDestroy_Serial(ns_s);
    N_VDestroy_Serial(ns_r);
    N_VDestroy_Serial(ns_rt);
    N_VDestroy_Serial(ns_v);
    N_VDestroy_Serial(ns_Jv);
    N_VDestroy_Serial(ns_tmp);
    N_VDestroy_Serial(ns_Jdiag);
    N_VScale(-1.0, tdata->xdot, tdata->xdot);

};

} // namespace amici
