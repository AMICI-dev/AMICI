#include "amici/newton_solver.h"

#include "amici/defines.h"
#include "amici/model.h"
#include "amici/solver.h"
#include "amici/vector.h"
#include "amici/steadystateproblem.h"
#include "amici/forwardproblem.h"
#include "amici/rdata.h"
#include "amici/edata.h"

#include <cstring>
#include <ctime>
#include <cmath>

namespace amici {

NewtonSolver::NewtonSolver(realtype *t, AmiVector *x, Model *model, ReturnData *rdata)
    : model(model), rdata(rdata), xdot(x->getLength()), dx(x->getLength())
    {
    /**
     * default constructor, initializes all members with the provided objects
     *
     * @param t pointer to time variable
     * @param x pointer to state variables
     * @param model pointer to the AMICI model object
     * @param rdata pointer to the return data object
     */
    this->t = t;
    this->x = x;
}

/* ----------------------------------------------------------------------------------
 */

std::unique_ptr<NewtonSolver> NewtonSolver::getSolver(realtype *t, AmiVector *x, LinearSolver linsolType, Model *model,
                                      ReturnData *rdata, int maxlinsteps, int maxsteps, double atol, double rtol) {
    /**
     * Tries to determine the steady state of the ODE system by a Newton
     * solver, uses forward intergration, if the Newton solver fails,
     * restarts Newton solver, if integration fails.
     * Computes steady state sensitivities
     *
     * @param t pointer to time variable
     * @param x pointer to state variables
     * @param linsolType integer indicating which linear solver to use
     * @param model pointer to the AMICI model object
     * @param rdata pointer to the return data object
     * @param maxlinsteps maximum number of allowed linear steps per Newton step for steady state computation
     * @param maxsteps maximum number of allowed Newton steps for steady state computation
     * @param atol absolute tolerance
     * @param rtol relative tolerance
     * @return solver NewtonSolver according to the specified linsolType
     */

    std::unique_ptr<NewtonSolver> solver;

    switch (linsolType) {

    /* DIRECT SOLVERS */
    case LinearSolver::dense:
        solver.reset(new NewtonSolverDense(t, x, model, rdata));
        break;

    case LinearSolver::band:
        throw NewtonFailure(AMICI_NOT_IMPLEMENTED, "getSolver");

    case LinearSolver::LAPACKDense:
        throw NewtonFailure(AMICI_NOT_IMPLEMENTED, "getSolver");

    case LinearSolver::LAPACKBand:
        throw NewtonFailure(AMICI_NOT_IMPLEMENTED, "getSolver");

    case LinearSolver::diag:
        throw NewtonFailure(AMICI_NOT_IMPLEMENTED, "getSolver");

    /* ITERATIVE SOLVERS */
    case LinearSolver::SPGMR:
        throw NewtonFailure(AMICI_NOT_IMPLEMENTED, "getSolver");

    case LinearSolver::SPBCG:
        solver.reset(new NewtonSolverIterative(t, x, model, rdata));
        break;

    case LinearSolver::SPTFQMR:
        throw NewtonFailure(AMICI_NOT_IMPLEMENTED, "getSolver");

    /* SPARSE SOLVERS */
    case LinearSolver::KLU:
        solver.reset(new NewtonSolverSparse(t, x, model, rdata));
        break;
    default:
        throw NewtonFailure(AMICI_NOT_IMPLEMENTED, "getSolver");
    }

    solver->atol = atol;
    solver->rtol = rtol;
    solver->maxlinsteps = maxlinsteps;
    solver->maxsteps = maxsteps;

    return solver;
}

/* ----------------------------------------------------------------------------------
 */

void NewtonSolver::getStep(int ntry, int nnewt, AmiVector *delta) {
    /**
     * Computes the solution of one Newton iteration
     *
     * @param ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param nnewt integer number of current Newton step
     * @param delta containing the RHS of the linear system, will be
     * overwritten by solution to the linear system
     */

    this->prepareLinearSystem(ntry, nnewt);

    delta->minus();
    this->solveLinearSystem(delta);
}

/* ----------------------------------------------------------------------------------
 */

void NewtonSolver::computeNewtonSensis(AmiVectorArray *sx) {
    /**
     * Computes steady state sensitivities
     *
     * @param sx pointer to state variable sensitivities
     */
    prepareLinearSystem(0, -1);

    model->fdxdotdp(*t, x, &dx);
    for (int ip = 0; ip < model->nplist(); ip++) {
        
        for (int ix = 0; ix < model->nx; ix++) {
            sx->at(ix,ip) = -model->dxdotdp[model->nx * ip + ix];
        }
        solveLinearSystem(&((*sx)[ip]));
    }
}
/* ----------------------------------------------------------------------------------
 */
/* - Dense linear solver
 * ------------------------------------------------------------ */
/* ----------------------------------------------------------------------------------
 */

/* Derived class for dense linear solver */
NewtonSolverDense::NewtonSolverDense(realtype *t, AmiVector *x, Model *model, ReturnData *rdata)
    : NewtonSolver(t, x, model, rdata) {
    /**
     * default constructor, initializes all members with the provided objects
     * and
     * initializes temporary storage objects
     *
     * @param t pointer to time variable
     * @param x pointer to state variables
     * @param model pointer to the AMICI model object
     * @param rdata pointer to the return data object
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
     * @param ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param nnewt integer number of current Newton step
     */

    /* Get Jacobian */
    model->fJ(*t, 0.0, x, &dx, &xdot, Jtmp);
    int status = DenseGETRF(Jtmp, pivots);
    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "DenseGETRF");
}

/* ----------------------------------------------------------------------------------
 */

void NewtonSolverDense::solveLinearSystem(AmiVector *rhs) {
    /**
     * Solves the linear system for the Newton step
     *
     * @param rhs containing the RHS of the linear system, will be
     * overwritten by solution to the linear system
     */

    /* Pass pointer to the linear solver */
    DenseGETRS(Jtmp, pivots, rhs->data());
}

/* ----------------------------------------------------------------------------------
 */

NewtonSolverDense::~NewtonSolverDense() {
    if(Jtmp)
        DestroyMat(Jtmp);
    if(pivots)
        DestroyArray(pivots);
}

/* ----------------------------------------------------------------------------------
 */
/* - Sparse linear solver
 * ----------------------------------------------------------- */
/* ----------------------------------------------------------------------------------
 */

/* Derived class for sparse linear solver */
NewtonSolverSparse::NewtonSolverSparse(realtype *t, AmiVector *x, Model *model, ReturnData *rdata)
    : NewtonSolver(t, x, model, rdata) {
    /**
     * default constructor, initializes all members with the provided objects,
     * initializes temporary storage objects and the klu solver
     *
     * @param t pointer to time variable
     * @param x pointer to state variables
     * @param model pointer to the AMICI model object
     * @param rdata pointer to the return data object
     */

    /* Initialize the KLU solver */
    klu_status = klu_defaults(&common);
    /* Check if KLU was initialized successfully */
    if (klu_status != 1)
        throw NewtonFailure(common.status, "klu_defaults");
    Jtmp = SparseNewMat(model->nx, model->nx, model->nnz, CSC_MAT);
}

/* ----------------------------------------------------------------------------------
 */

void NewtonSolverSparse::prepareLinearSystem(int ntry, int nnewt) {
    /**
     * Writes the Jacobian for the Newton iteration and passes it to the linear
     * solver
     *
     * @param ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param nnewt integer number of current Newton step
     */

    /* Get sparse Jacobian */
    model->fJSparse(*t, 0.0, x, &dx, &xdot, Jtmp);

    /* Get factorization of sparse Jacobian */
    if(symbolic) /* if symbolic was already created free first to avoid memory leak */
        klu_free_symbolic(&symbolic, &common);
    symbolic = klu_analyze(model->nx, Jtmp->indexptrs,
                           Jtmp->indexvals, &common);
    if (!symbolic) {
        throw NewtonFailure(common.status,"klu_analyze");
    }

    if(numeric) /* if numeric was already created free first to avoid memory leak */
        klu_free_numeric(&numeric, &common);
    numeric = klu_factor(Jtmp->indexptrs, Jtmp->indexvals,
                         Jtmp->data, symbolic, &common);
    if (!numeric) {
        throw NewtonFailure(common.status,"klu_factor");
    }
}

/* ----------------------------------------------------------------------------------
 */

void NewtonSolverSparse::solveLinearSystem(AmiVector *rhs) {
    /**
     * Solves the linear system for the Newton step
     *
     * @param rhs containing the RHS of the linear system,will be
     * overwritten by solution to the linear system
     */

    /* Pass pointer to the linear solver */
    klu_status = klu_solve(symbolic, numeric, model->nx, 1, rhs->data(), &common);
    if (klu_status != 1)
        throw NewtonFailure(common.status, "klu_solve");
}

/* ----------------------------------------------------------------------------------
 */

NewtonSolverSparse::~NewtonSolverSparse() {
    if(Jtmp)
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
NewtonSolverIterative::NewtonSolverIterative(realtype *t, AmiVector *x, Model *model, ReturnData *rdata)
    : NewtonSolver(t, x, model, rdata), ns_p(model->nx), ns_h(model->nx),
    ns_t(model->nx), ns_s(model->nx), ns_r(model->nx), ns_rt(model->nx), ns_v(model->nx),
    ns_Jv(model->nx), ns_tmp(model->nx), ns_Jdiag(model->nx)
    {
    /**
     * default constructor, initializes all members with the provided objects
     * @param t pointer to time variable
     * @param x pointer to state variables
     * @param model pointer to the AMICI model object
     * @param rdata pointer to the return data object
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
     * @param ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param nnewt integer number of current Newton step
     */

    newton_try = ntry;
    i_newton = nnewt;
    if (nnewt == -1) {
        throw AmiException("Linear solver SPBCG does not support sensitivity computation for steady state problems.");
    }
}

/* ----------------------------------------------------------------------------------
 */

void NewtonSolverIterative::solveLinearSystem(AmiVector *rhs) {
    /**
     * Solves the linear system for the Newton step by passing it to
     * linsolveSPBCG
     *
     * @param rhs containing the RHS of the linear system, will be
     * overwritten by solution to the linear system
     */

    linsolveSPBCG(newton_try, i_newton, rhs);
}
    
    
void NewtonSolverIterative::linsolveSPBCG(int ntry,int nnewt, AmiVector *ns_delta) {
    /**
     * Iterative linear solver created from SPILS BiCG-Stab.
     * Solves the linear system within each Newton step if iterative solver is
     * chosen.
     *
     * @param ntry integer newton_try integer start number of Newton solver
     * (1 or 2)
     * @param nnewt integer number of current Newton step
     * @param ns_delta Newton step
     */
    
    double rho;
    double alpha;
    double omega;
    double res;
    
    xdot = *ns_delta;
    xdot.minus();
    
    // Get the diagonal of the Jacobian for preconditioning
    model->fJDiag(*t, &ns_Jdiag, 0.0, x, &dx);
    
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
    model->fJv(*t, x, &dx, &xdot, ns_delta, &ns_Jv, 0.0);
    
    // ns_r = xdot - ns_Jv;
    N_VLinearSum(-1.0, ns_Jv.getNVector(), 1.0, xdot.getNVector(), ns_r.getNVector());
    N_VDiv(ns_r.getNVector(), ns_Jdiag.getNVector(), ns_r.getNVector());
    res = sqrt(N_VDotProd(ns_r.getNVector(), ns_r.getNVector()));
    ns_rt = ns_r;
    
    for (int i_linstep = 0; i_linstep < maxlinsteps;
         i_linstep++) {
        // Compute factors
        double rho1 = rho;
        rho = N_VDotProd(ns_rt.getNVector(), ns_r.getNVector());
        double beta = rho * alpha / (rho1 * omega);
        
        // ns_p = ns_r + beta * (ns_p - omega * ns_v);
        N_VLinearSum(1.0, ns_p.getNVector(), -omega, ns_v.getNVector(), ns_p.getNVector());
        N_VLinearSum(1.0, ns_r.getNVector(), beta, ns_p.getNVector(), ns_p.getNVector());
        
        // ns_v = J * ns_p
        model->fJv(*t, x, &dx, &xdot, &ns_p, &ns_v, 0.0);
        N_VDiv(ns_v.getNVector(), ns_Jdiag.getNVector(), ns_v.getNVector());
        
        // Compute factor
        alpha = rho / N_VDotProd(ns_rt.getNVector(), ns_v.getNVector());
        
        // ns_h = ns_delta + alpha * ns_p;
        N_VLinearSum(1.0, ns_delta->getNVector(), alpha, ns_p.getNVector(), ns_h.getNVector());
        // ns_s = ns_r - alpha * ns_v;
        N_VLinearSum(1.0, ns_r.getNVector(), -alpha, ns_v.getNVector(), ns_s.getNVector());
        
        // ns_t = J * ns_s
        model->fJv(*t, x, &dx, &xdot, &ns_s, &ns_t, 0.0);
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
        if (res < atol) {
            // Write number of steps needed
            rdata->newton_numlinsteps[(ntry - 1) * maxsteps +
                                      nnewt] = i_linstep + 1;
            
            // Return success
            return;
        }
        
        // Scale back
        N_VDiv(ns_r.getNVector(), ns_Jdiag.getNVector(), ns_r.getNVector());
    }
    throw NewtonFailure(AMICI_CONV_FAILURE, "linsolveSPBCG");
}
    

} // namespace amici
