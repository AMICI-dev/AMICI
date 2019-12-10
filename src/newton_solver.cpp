#include "amici/newton_solver.h"

#include "amici/model.h"
#include "amici/solver.h"
#include "amici/steadystateproblem.h"
#include "amici/forwardproblem.h"
#include "amici/rdata.h"
#include "amici/edata.h"

#include "sunlinsol/sunlinsol_klu.h" // sparse solver
#include "sunlinsol/sunlinsol_dense.h" // dense solver

#include <cstring>
#include <ctime>
#include <cmath>

namespace amici {

NewtonSolver::NewtonSolver(realtype *t, AmiVector *x, Model *model,
                           ReturnData *rdata)
    : model(model), rdata(rdata), xdot(model->nx_solver), dx(model->nx_solver)
    {
    this->t = t;
    this->x = x;
}

/* ------------------------------------------------------------------------- */

std::unique_ptr<NewtonSolver> NewtonSolver::getSolver(
        realtype *t, AmiVector *x, LinearSolver linsolType, Model *model,
        ReturnData *rdata, int maxlinsteps, int maxsteps, double atol, double rtol,
        NewtonDampingFactorMode dampingFactorMode, double dampingFactorLowerBound) {

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
    case LinearSolver::SuperLUMT:
        throw NewtonFailure(AMICI_NOT_IMPLEMENTED, "getSolver");
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
    solver->dampingFactorMode = dampingFactorMode;
    solver->dampingFactorLowerBound = dampingFactorLowerBound;

    return solver;
}

/* ------------------------------------------------------------------------- */

void NewtonSolver::getStep(int ntry, int nnewt, AmiVector &delta) {
    this->prepareLinearSystem(ntry, nnewt);

    delta.minus();
    this->solveLinearSystem(delta);
}

/* ------------------------------------------------------------------------- */

void NewtonSolver::computeNewtonSensis(AmiVectorArray &sx) {
    prepareLinearSystem(0, -1);
    model->fdxdotdp(*t, *x, dx);

    if (model->pythonGenerated) {
        for (int ip = 0; ip < model->nplist(); ip++) {
            N_VConst(0.0, sx.getNVector(ip));

            // copy explicit version
            if (model->ndxdotdp_explicit > 0) {
                auto col = model->dxdotdp_explicit.indexptrs();
                auto row = model->dxdotdp_explicit.indexvals();
                auto data_ptr = model->dxdotdp_explicit.data();
                for (sunindextype iCol = col[model->plist(ip)];
                     iCol < col[model->plist(ip) + 1]; ++iCol)
                    sx.at(row[iCol], ip) -= data_ptr[iCol];
            }

            // copy implicit version
            if (model->ndxdotdp_implicit > 0) {
                auto col = model->dxdotdp_implicit.indexptrs();
                auto row = model->dxdotdp_implicit.indexvals();
                auto data_ptr = model->dxdotdp_implicit.data();
                for (sunindextype iCol = col[model->plist(ip)];
                     iCol < col[model->plist(ip) + 1]; ++iCol)
                    sx.at(row[iCol], ip) -= data_ptr[iCol];
            }

            solveLinearSystem(sx[ip]);
        }
    } else {
        for (int ip = 0; ip < model->nplist(); ip++) {
            for (int ix = 0; ix < model->nx_solver; ix++)
                sx.at(ix,ip) = -model->dxdotdp.at(ix, ip);

            solveLinearSystem(sx[ip]);
        }
    }
}
/* ------------------------------------------------------------------------- */
/* - Dense linear solver --------------------------------------------------- */
/* ------------------------------------------------------------------------- */

/* Derived class for dense linear solver */
NewtonSolverDense::NewtonSolverDense(realtype *t, AmiVector *x, Model *model,
                                     ReturnData *rdata)
    : NewtonSolver(t, x, model, rdata),
      Jtmp(model->nx_solver, model->nx_solver),
      linsol(SUNLinSol_Dense(x->getNVector(), Jtmp.get()))
{
    int status = SUNLinSolInitialize_Dense(linsol);
    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolInitialize_Dense");
}

/* ------------------------------------------------------------------------- */

void NewtonSolverDense::prepareLinearSystem(int  /*ntry*/, int  /*nnewt*/) {
    model->fJ(*t, 0.0, *x, dx, xdot, Jtmp.get());
    int status = SUNLinSolSetup_Dense(linsol, Jtmp.get());
    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolSetup_Dense");
}

/* ------------------------------------------------------------------------- */

void NewtonSolverDense::solveLinearSystem(AmiVector &rhs) {
    int status = SUNLinSolSolve_Dense(linsol, Jtmp.get(),
                                      rhs.getNVector(), rhs.getNVector(),
                                      0.0);
    // last argument is tolerance and does not have any influence on result

    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolSolve_Dense");
}

/* ------------------------------------------------------------------------- */

NewtonSolverDense::~NewtonSolverDense() {
    if(linsol)
        SUNLinSolFree_Dense(linsol);
}

/* ------------------------------------------------------------------------- */
/* - Sparse linear solver -------------------------------------------------- */
/* ------------------------------------------------------------------------- */

/* Derived class for sparse linear solver */
NewtonSolverSparse::NewtonSolverSparse(realtype *t, AmiVector *x, Model *model,
                                       ReturnData *rdata)
    : NewtonSolver(t, x, model, rdata),
      Jtmp(model->nx_solver, model->nx_solver, model->nnz, CSC_MAT),
      linsol(SUNKLU(x->getNVector(), Jtmp.get()))
{
    int status = SUNLinSolInitialize_KLU(linsol);
    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolInitialize_KLU");
}

/* ------------------------------------------------------------------------- */

void NewtonSolverSparse::prepareLinearSystem(int  /*ntry*/, int  /*nnewt*/) {
    /* Get sparse Jacobian */
    model->fJSparse(*t, 0.0, *x, dx, xdot, Jtmp.get());
    int status = SUNLinSolSetup_KLU(linsol, Jtmp.get());
    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolSetup_KLU");
}

/* ------------------------------------------------------------------------- */

void NewtonSolverSparse::solveLinearSystem(AmiVector &rhs) {
    /* Pass pointer to the linear solver */
    int status = SUNLinSolSolve_KLU(linsol, Jtmp.get(),
                                    rhs.getNVector(), rhs.getNVector(), 0.0);
    // last argument is tolerance and does not have any influence on result

    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolSolve_Dense");
}

/* ------------------------------------------------------------------------- */

NewtonSolverSparse::~NewtonSolverSparse() {
    if(linsol)
        SUNLinSolFree_KLU(linsol);
}

/* ------------------------------------------------------------------------- */
/* - Iterative linear solver------------------------------------------------ */
/* ------------------------------------------------------------------------- */

NewtonSolverIterative::NewtonSolverIterative(realtype *t, AmiVector *x,
                                             Model *model, ReturnData *rdata)
    : NewtonSolver(t, x, model, rdata), ns_p(model->nx_solver),
    ns_h(model->nx_solver), ns_t(model->nx_solver), ns_s(model->nx_solver),
    ns_r(model->nx_solver), ns_rt(model->nx_solver), ns_v(model->nx_solver),
    ns_Jv(model->nx_solver), ns_tmp(model->nx_solver),
    ns_Jdiag(model->nx_solver)
    {
}

/* ------------------------------------------------------------------------- */

void NewtonSolverIterative::prepareLinearSystem(int ntry, int nnewt) {
    newton_try = ntry;
    i_newton = nnewt;
    if (nnewt == -1) {
        throw AmiException("Linear solver SPBCG does not support sensitivity "
                           "computation for steady state problems.");
    }
}

/* ------------------------------------------------------------------------- */

void NewtonSolverIterative::solveLinearSystem(AmiVector &rhs) {
    linsolveSPBCG(newton_try, i_newton, rhs);
    rhs.minus();
}


void NewtonSolverIterative::linsolveSPBCG(int ntry, int nnewt,
                                          AmiVector &ns_delta) {
    xdot = ns_delta;
    xdot.minus();

    // Get the diagonal of the Jacobian for preconditioning
    model->fJDiag(*t, ns_Jdiag, 0.0, *x, dx);

    // Ensure positivity of entries in ns_Jdiag
    ns_p.set(1.0);
    N_VAbs(ns_Jdiag.getNVector(), ns_Jdiag.getNVector());
    N_VCompare(1e-15, ns_Jdiag.getNVector(), ns_tmp.getNVector());
    N_VLinearSum(-1.0, ns_tmp.getNVector(), 1.0, ns_p.getNVector(), ns_tmp.getNVector());
    N_VLinearSum(1.0, ns_Jdiag.getNVector(), 1.0, ns_tmp.getNVector(), ns_Jdiag.getNVector());

    // Initialize for linear solve
    ns_p.reset();
    ns_v.reset();
    ns_delta.reset();
    ns_tmp.reset();
    double rho = 1.0;
    double omega = 1.0;
    double alpha = 1.0;

    // can be set to 0 at the moment
    model->fJv(*t, *x, dx, xdot, ns_delta, ns_Jv, 0.0);

    // ns_r = xdot - ns_Jv;
    N_VLinearSum(-1.0, ns_Jv.getNVector(), 1.0, xdot.getNVector(), ns_r.getNVector());
    N_VDiv(ns_r.getNVector(), ns_Jdiag.getNVector(), ns_r.getNVector());
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
        model->fJv(*t, *x, dx, xdot, ns_p, ns_v, 0.0);
        N_VDiv(ns_v.getNVector(), ns_Jdiag.getNVector(), ns_v.getNVector());

        // Compute factor
        alpha = rho / N_VDotProd(ns_rt.getNVector(), ns_v.getNVector());

        // ns_h = ns_delta + alpha * ns_p;
        N_VLinearSum(1.0, ns_delta.getNVector(), alpha, ns_p.getNVector(),
                     ns_h.getNVector());
        // ns_s = ns_r - alpha * ns_v;
        N_VLinearSum(1.0, ns_r.getNVector(), -alpha, ns_v.getNVector(), ns_s.getNVector());

        // ns_t = J * ns_s
        model->fJv(*t, *x, dx, xdot, ns_s, ns_t, 0.0);
        N_VDiv(ns_t.getNVector(), ns_Jdiag.getNVector(), ns_t.getNVector());

        // Compute factor
        omega = N_VDotProd(ns_t.getNVector(), ns_s.getNVector()) / N_VDotProd(ns_t.getNVector(), ns_t.getNVector());

        // ns_delta = ns_h + omega * ns_s;
        N_VLinearSum(1.0, ns_h.getNVector(), omega, ns_s.getNVector(),
                     ns_delta.getNVector());
        // ns_r = ns_s - omega * ns_t;
        N_VLinearSum(1.0, ns_s.getNVector(), -omega, ns_t.getNVector(), ns_r.getNVector());

        // Compute the (unscaled) residual
        N_VProd(ns_r.getNVector(), ns_Jdiag.getNVector(), ns_r.getNVector());
        double res = sqrt(N_VDotProd(ns_r.getNVector(), ns_r.getNVector()));

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
