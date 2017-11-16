#include "include/amici_defines.h"
#include "include/amici_model.h"
#include "include/amici_solver.h"
#include "include/amici_vector.h"
#include "include/steadystateproblem.h"
#include "include/forwardproblem.h"
#include "include/rdata.h"
#include "include/udata.h"
#include "include/edata.h"
#include "include/newton_solver.h"
#include <cstring>
#include <ctime>

namespace amici {

NewtonSolver::NewtonSolver(realtype *t, Model *model, ReturnData *rdata,
                           const UserData *udata)
    : model(model), rdata(rdata), udata(udata),
    sx_ip(model->nx), xdot(model->nx), x(model->nx), dx(model->nx)
    {
    /**
     * default constructor, initializes all members with the provided objects
     *
     * @param t current timepoint
     * @param model pointer to the AMICI model object
     * @param rdata pointer to the return data object
     * @param udata pointer to the user data object
     */
    this->t = t;
}

/* ----------------------------------------------------------------------------------
 */

NewtonSolver *NewtonSolver::getSolver(realtype *t, int linsolType, Model *model,
                                      ReturnData *rdata, const UserData *udata) {
    /**
     * Tries to determine the steady state of the ODE system by a Newton
     * solver, uses forward intergration, if the Newton solver fails,
     * restarts Newton solver, if integration fails.
     * Computes steady state sensitivities
     *
     * @param[in] linsolType integer indicating which linear solver to use
     * @param[in] model pointer to the AMICI model object @type Model
     * @param[in] udata pointer to the user data object @type UserData
     * @param[out] rdata pointer to the return data object @type ReturnData
     * @return solver NewtonSolver according to the specified linsolType
     */

    switch (linsolType) {

    /* DIRECT SOLVERS */
    case AMICI_DENSE:
        return new NewtonSolverDense(t, model, rdata, udata);

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
        return new NewtonSolverIterative(t, model, rdata, udata);

    case AMICI_SPTFQMR:
        throw NewtonFailure("Solver currently not supported!");

    /* SPARSE SOLVERS */
    case AMICI_KLU:
        return new NewtonSolverSparse(t, model, rdata, udata);

    default:
        throw NewtonFailure("Invalid Choice of Solver!");
    }
}

/* ----------------------------------------------------------------------------------
 */

void NewtonSolver::getStep(int ntry, int nnewt, AmiVector *delta) {
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

    *delta = xdot;
    this->solveLinearSystem(delta);
}

/* ----------------------------------------------------------------------------------
 */

void NewtonSolver::getSensis(int it) {
    /**
     * Computes steady state sensitivities
     *
     * @param[in] it integer index of current time step     */
    this->prepareLinearSystem(0, -1);

    model->fdxdotdp(*t, &x, &dx);
    for (int ip = 0; ip < udata->nplist(); ip++) {
        
        for (int ix = 0; ix < model->nx; ix++) {
            sx_ip[ix] = -model->dxdotdp[model->nx * ip + ix];
        }
        this->solveLinearSystem(&sx_ip);
        
        /* Copy result to return data */
        if (it == AMICI_PREEQUILIBRATE) {
            for (int ix = 0; ix < model->nx; ix++) {
                rdata->sx0[ip * model->nx + ix] = sx_ip[ix];
            }
        } else {
            /* Classical steady state computation */
            for (int ix = 0; ix < model->nx; ix++) {
                rdata->sx[(ip * model->nx + ix) * rdata->nt + it] =
                sx_ip[ix];
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
NewtonSolverDense::NewtonSolverDense(realtype *t, Model *model, ReturnData *rdata,
                                     const UserData *udata)
    : NewtonSolver(t, model, rdata, udata) {
    /**
     * default constructor, initializes all members with the provided objects
     * and
     * initializes temporary storage objects
     *
     * @param[in] model pointer to the AMICI model object @type Model
     * @param[in] rdata pointer to the return data object @type ReturnData
     * @param[in] udata pointer to the user data object @type UserData
     */
     pivots = NewLintArray(model->nx);
     Jtmp = NewDenseMat(model->nx,model->nx);
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
    model->fJ(*t, 0.0, &x, &dx, &xdot, Jtmp);
    int status = DenseGETRF(Jtmp, pivots);
    if(status != 0)
        throw NewtonFailure("Dense factorization failed!");
}

/* ----------------------------------------------------------------------------------
 */

void NewtonSolverDense::solveLinearSystem(AmiVector *rhs) {
    /**
     * Solves the linear system for the Newton step
     *
     * @param[in,out] rhs containing the RHS of the linear system, will be
     * overwritten by solution to the linear system @type N_Vector
     */

    /* Pass pointer to the linear solver */
    DenseGETRS(Jtmp, pivots, rhs->data());
}

/* ----------------------------------------------------------------------------------
 */

NewtonSolverDense::~NewtonSolverDense() {
    DestroyMat(Jtmp);
    DestroyArray(pivots);
}

/* ----------------------------------------------------------------------------------
 */
/* - Sparse linear solver
 * ----------------------------------------------------------- */
/* ----------------------------------------------------------------------------------
 */

/* Derived class for sparse linear solver */
NewtonSolverSparse::NewtonSolverSparse(realtype *t, Model *model, ReturnData *rdata,
                                       const UserData *udata)
    : NewtonSolver(t, model, rdata, udata) {
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
    Jtmp = SparseNewMat(model->nx, model->nx, model->nnz, CSC_MAT);
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
    model->fJSparse(*t, 0.0, &x, &dx, &xdot, Jtmp);

    /* Get factorization of sparse Jacobian */
    if(symbolic) /* if symbolic was already created free first to avoid memory leak */
        klu_free_symbolic(&symbolic, &common);
    symbolic = klu_analyze(model->nx, Jtmp->indexptrs,
                           Jtmp->indexvals, &common);
    if (symbolic) {
        if(numeric) /* if numeric was already created free first to avoid memory leak */
            klu_free_numeric(&numeric, &common);
        numeric = klu_factor(Jtmp->indexptrs, Jtmp->indexvals,
                             Jtmp->data, symbolic, &common);
        if (numeric) {
            return;
        } else {
            throw NewtonFailure("KLU factorization failed!");
        }
    } else {
        throw NewtonFailure("KLU symbolic analysis failed!");
    }
}

/* ----------------------------------------------------------------------------------
 */

void NewtonSolverSparse::solveLinearSystem(AmiVector *rhs) {
    /**
     * Solves the linear system for the Newton step
     *
     * @param[in] rhs containing the RHS of the linear system,will be
     * overwritten by solution to the linear system     */

    /* Pass pointer to the linear solver */
    klu_status = klu_solve(symbolic, numeric, model->nx, 1, rhs->data(), &common);
    if (klu_status != 1)
        throw NewtonFailure("KLU solver failed");
}

/* ----------------------------------------------------------------------------------
 */

NewtonSolverSparse::~NewtonSolverSparse() {
    SparseDestroyMat(Jtmp);
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
NewtonSolverIterative::NewtonSolverIterative(realtype *t, Model *model, ReturnData *rdata,
                                             const UserData *udata)
    : NewtonSolver(t, model, rdata, udata), ns_p(model->nx), ns_h(model->nx),
    ns_t(model->nx), ns_s(model->nx), ns_r(model->nx), ns_rt(model->nx), ns_v(model->nx),
    ns_Jv(model->nx), ns_tmp(model->nx), ns_Jdiag(model->nx)
    {
    /**
     * default constructor, initializes all members with the provided objects
     * @param[in] model pointer to the AMICI model object @type Model
     * @param[in] rdata pointer to the return data object @type ReturnData
     * @param[in] udata pointer to the user data object @type UserData
     */
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
        throw NewtonFailure("Linear solver SPBCG does not support sensitivity computation for steady state problems.");
    }
}

/* ----------------------------------------------------------------------------------
 */

void NewtonSolverIterative::solveLinearSystem(AmiVector *rhs) {
    /**
     * Solves the linear system for the Newton step by passing it to
     * linsolveSPBCG
     *
     * @param[in,out] rhs containing the RHS of the linear system, will be
     * overwritten by solution to the linear system @type N_Vector
     */

    linsolveSPBCG(newton_try, i_newton, rhs);
}
    
    
void NewtonSolverIterative::linsolveSPBCG(int ntry,int nnewt, AmiVector *ns_delta) {
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
    model->fJDiag(*t, &ns_Jdiag, 0.0, &x, &dx);
    
    // Ensure positivity of entries in ns_Jdiag
    ns_p.set(1.0);
    N_VAbs(ns_Jdiag.getNVector(), ns_tmp.getNVector());
    N_VCompare(1e-15, ns_tmp.getNVector(), ns_tmp.getNVector());
    N_VLinearSum(-1.0, ns_tmp.getNVector(), 1.0, ns_p.getNVector(), ns_tmp.getNVector());
    N_VLinearSum(1.0, ns_Jdiag.getNVector(), 1.0, ns_tmp.getNVector(), ns_Jdiag.getNVector());
    
    // Initialize for linear solve
    ns_p.reset();
    ns_v.reset();
    ns_delta->reset();
    ns_tmp.reset();
    rho = 1.0;
    omega = 1.0;
    alpha = 1.0;
    
    // can be set to 0 at the moment
    model->fJv(*t, &x, &dx, &xdot, ns_delta, &ns_Jv, 0.0);
    
    // ns_r = xdot - ns_Jv;
    N_VLinearSum(-1.0, ns_Jv.getNVector(), 1.0, xdot.getNVector(), ns_r.getNVector());
    N_VDiv(ns_r.getNVector(), ns_Jdiag.getNVector(), ns_r.getNVector());
    res = sqrt(N_VDotProd(ns_r.getNVector(), ns_r.getNVector()));
    ns_rt = ns_r;
    
    for (int i_linstep = 0; i_linstep < udata->newton_maxlinsteps;
         i_linstep++) {
        // Compute factors
        rho1 = rho;
        rho = N_VDotProd(ns_rt.getNVector(), ns_r.getNVector());
        beta = rho * alpha / (rho1 * omega);
        
        // ns_p = ns_r + beta * (ns_p - omega * ns_v);
        N_VLinearSum(1.0, ns_p.getNVector(), -omega, ns_v.getNVector(), ns_p.getNVector());
        N_VLinearSum(1.0, ns_r.getNVector(), beta, ns_p.getNVector(), ns_p.getNVector());
        
        // ns_v = J * ns_p
        model->fJv(*t, &x, &dx, &xdot, &ns_p, &ns_v, 0.0);
        N_VDiv(ns_v.getNVector(), ns_Jdiag.getNVector(), ns_v.getNVector());
        
        // Compute factor
        alpha = rho / N_VDotProd(ns_rt.getNVector(), ns_v.getNVector());
        
        // ns_h = ns_delta + alpha * ns_p;
        N_VLinearSum(1.0, ns_delta->getNVector(), alpha, ns_p.getNVector(), ns_h.getNVector());
        // ns_s = ns_r - alpha * ns_v;
        N_VLinearSum(1.0, ns_r.getNVector(), -alpha, ns_v.getNVector(), ns_s.getNVector());
        
        // ns_t = J * ns_s
        model->fJv(*t, &x, &dx, &xdot, &ns_s, &ns_t, 0.0);
        N_VDiv(ns_t.getNVector(), ns_Jdiag.getNVector(), ns_t.getNVector());
        
        // Compute factor
        omega = N_VDotProd(ns_t.getNVector(), ns_s.getNVector()) / N_VDotProd(ns_t.getNVector(), ns_t.getNVector());
        
        // ns_delta = ns_h + omega * ns_s;
        N_VLinearSum(1.0, ns_h.getNVector(), omega, ns_s.getNVector(), ns_delta->getNVector());
        // ns_r = ns_s - omega * ns_t;
        N_VLinearSum(1.0, ns_s.getNVector(), -omega, ns_t.getNVector(), ns_r.getNVector());
        
        // Compute the (unscaled) residual
        N_VProd(ns_r.getNVector(), ns_Jdiag.getNVector(), ns_r.getNVector());
        res = sqrt(N_VDotProd(ns_r.getNVector(), ns_r.getNVector()));
        
        // Test convergence
        if (res < udata->atol) {
            // Write number of steps needed
            rdata->newton_numlinsteps[(ntry - 1) * udata->newton_maxsteps +
                                      nnewt] = i_linstep + 1;
            
            // Return success
            return;
        }
        
        // Scale back
        N_VDiv(ns_r.getNVector(), ns_Jdiag.getNVector(), ns_r.getNVector());
    }
    throw NewtonFailure("SPBCG solver failed to converge");
}
    

/* ----------------------------------------------------------------------------------
 */

NewtonSolverIterative::~NewtonSolverIterative(){
};

} // namespace amici
