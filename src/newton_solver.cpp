#include "amici/newton_solver.h"

#include "amici/model.h"
#include "amici/solver.h"

#include "sunlinsol/sunlinsol_klu.h" // sparse solver
#include "sunlinsol/sunlinsol_dense.h" // dense solver

#include <cstring>
#include <ctime>
#include <cmath>

namespace amici {

NewtonSolver::NewtonSolver(realtype *t, AmiVector *x, Model *model)
    : t_(t), model_(model), xdot_(model->nx_solver), x_(x),
      dx_(model->nx_solver), xB_(model->nx_solver), dxB_(model->nx_solver) {
}

/* ------------------------------------------------------------------------- */

std::unique_ptr<NewtonSolver> NewtonSolver::getSolver(realtype *t, AmiVector *x,
                                                      Solver &simulationSolver,
                                                      Model *model) {

    std::unique_ptr<NewtonSolver> solver;

    switch (simulationSolver.getLinearSolver()) {

    /* DIRECT SOLVERS */
    case LinearSolver::dense:
        solver.reset(new NewtonSolverDense(t, x, model));
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
        solver.reset(new NewtonSolverIterative(t, x, model));
        break;

    case LinearSolver::SPTFQMR:
        throw NewtonFailure(AMICI_NOT_IMPLEMENTED, "getSolver");

    /* SPARSE SOLVERS */
    case LinearSolver::SuperLUMT:
        throw NewtonFailure(AMICI_NOT_IMPLEMENTED, "getSolver");
    case LinearSolver::KLU:
        solver.reset(new NewtonSolverSparse(t, x, model));
        break;
    default:
        throw NewtonFailure(AMICI_NOT_IMPLEMENTED, "getSolver");
    }
    solver->atol_ = simulationSolver.getAbsoluteTolerance();
    solver->rtol_ = simulationSolver.getRelativeTolerance();
    solver->max_lin_steps_ = simulationSolver.getNewtonMaxLinearSteps();
    solver->max_steps = simulationSolver.getNewtonMaxSteps();
    solver->damping_factor_mode_ = simulationSolver.getNewtonDampingFactorMode();
    solver->damping_factor_lower_bound =
        simulationSolver.getNewtonDampingFactorLowerBound();
    if (simulationSolver.getLinearSolver() == LinearSolver::SPBCG)
        solver->num_lin_steps_.resize(simulationSolver.getNewtonMaxSteps(), 0);

    return solver;
}

/* ------------------------------------------------------------------------- */

void NewtonSolver::getStep(int ntry, int nnewt, AmiVector &delta) {
    prepareLinearSystem(ntry, nnewt);

    delta.minus();
    solveLinearSystem(delta);
}

/* ------------------------------------------------------------------------- */

void NewtonSolver::computeNewtonSensis(AmiVectorArray &sx) {
    prepareLinearSystem(0, -1);
    model_->fdxdotdp(*t_, *x_, dx_);

    if (model_->pythonGenerated) {
        for (int ip = 0; ip < model_->nplist(); ip++) {
            N_VConst(0.0, sx.getNVector(ip));
            model_->get_dxdotdp_full().scatter(model_->plist(ip), -1.0, nullptr,
                                               gsl::make_span(sx.getNVector(ip)),
                                               0, nullptr, 0);

            solveLinearSystem(sx[ip]);
        }
    } else {
        for (int ip = 0; ip < model_->nplist(); ip++) {
            for (int ix = 0; ix < model_->nx_solver; ix++)
                sx.at(ix,ip) = -model_->get_dxdotdp().at(ix, ip);

            solveLinearSystem(sx[ip]);
        }
    }
}

/* ------------------------------------------------------------------------- */
/* - Dense linear solver --------------------------------------------------- */
/* ------------------------------------------------------------------------- */

/* Derived class for dense linear solver */
NewtonSolverDense::NewtonSolverDense(realtype *t, AmiVector *x, Model *model)
    : NewtonSolver(t, x, model), Jtmp_(model->nx_solver, model->nx_solver),
      linsol_(SUNLinSol_Dense(x->getNVector(), Jtmp_.get())) {
    int status = SUNLinSolInitialize_Dense(linsol_);
    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolInitialize_Dense");
}

/* ------------------------------------------------------------------------- */

void NewtonSolverDense::prepareLinearSystem(int  /*ntry*/, int  /*nnewt*/) {
    model_->fJ(*t_, 0.0, *x_, dx_, xdot_, Jtmp_.get());
    Jtmp_.refresh();
    int status = SUNLinSolSetup_Dense(linsol_, Jtmp_.get());
    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolSetup_Dense");
}

/* ------------------------------------------------------------------------- */

void NewtonSolverDense::prepareLinearSystemB(int  /*ntry*/, int  /*nnewt*/) {
    model_->fJB(*t_, 0.0, *x_, dx_, xB_, dxB_, xdot_, Jtmp_.get());
    Jtmp_.refresh();
    int status = SUNLinSolSetup_Dense(linsol_, Jtmp_.get());
    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolSetup_Dense");
}


/* ------------------------------------------------------------------------- */

void NewtonSolverDense::solveLinearSystem(AmiVector &rhs) {
    int status = SUNLinSolSolve_Dense(linsol_, Jtmp_.get(),
                                      rhs.getNVector(), rhs.getNVector(),
                                      0.0);
    Jtmp_.refresh();
    // last argument is tolerance and does not have any influence on result

    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolSolve_Dense");
}

/* ------------------------------------------------------------------------- */

NewtonSolverDense::~NewtonSolverDense() {
    if(linsol_)
        SUNLinSolFree_Dense(linsol_);
}

/* ------------------------------------------------------------------------- */
/* - Sparse linear solver -------------------------------------------------- */
/* ------------------------------------------------------------------------- */

/* Derived class for sparse linear solver */
NewtonSolverSparse::NewtonSolverSparse(realtype *t, AmiVector *x, Model *model)
    : NewtonSolver(t, x, model),
      Jtmp_(model->nx_solver, model->nx_solver, model->nnz, CSC_MAT),
      linsol_(SUNKLU(x->getNVector(), Jtmp_.get())) {
    int status = SUNLinSolInitialize_KLU(linsol_);
    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolInitialize_KLU");
}

/* ------------------------------------------------------------------------- */

void NewtonSolverSparse::prepareLinearSystem(int  /*ntry*/, int  /*nnewt*/) {
    /* Get sparse Jacobian */
    model_->fJSparse(*t_, 0.0, *x_, dx_, xdot_, Jtmp_.get());
    Jtmp_.refresh();
    int status = SUNLinSolSetup_KLU(linsol_, Jtmp_.get());
    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolSetup_KLU");
}

/* ------------------------------------------------------------------------- */

void NewtonSolverSparse::prepareLinearSystemB(int  /*ntry*/, int  /*nnewt*/) {
    /* Get sparse Jacobian */
    model_->fJSparseB(*t_, 0.0, *x_, dx_, xB_, dxB_, xdot_, Jtmp_.get());
    Jtmp_.refresh();
    int status = SUNLinSolSetup_KLU(linsol_, Jtmp_.get());
    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolSetup_KLU");
}

/* ------------------------------------------------------------------------- */

void NewtonSolverSparse::solveLinearSystem(AmiVector &rhs) {
    /* Pass pointer to the linear solver */
    int status = SUNLinSolSolve_KLU(linsol_, Jtmp_.get(),
                                    rhs.getNVector(), rhs.getNVector(), 0.0);
    // last argument is tolerance and does not have any influence on result

    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolSolve_KLU");
}

/* ------------------------------------------------------------------------- */

NewtonSolverSparse::~NewtonSolverSparse() {
    if(linsol_)
        SUNLinSolFree_KLU(linsol_);
}

/* ------------------------------------------------------------------------- */
/* - Iterative linear solver------------------------------------------------ */
/* ------------------------------------------------------------------------- */

NewtonSolverIterative::NewtonSolverIterative(realtype *t, AmiVector *x,
                                             Model *model)
    : NewtonSolver(t, x, model), ns_p_(model->nx_solver),
    ns_h_(model->nx_solver), ns_t_(model->nx_solver), ns_s_(model->nx_solver),
    ns_r_(model->nx_solver), ns_rt_(model->nx_solver), ns_v_(model->nx_solver),
    ns_Jv_(model->nx_solver), ns_tmp_(model->nx_solver),
    ns_Jdiag_(model->nx_solver), ns_J_(model->nx_solver, model->nx_solver)
    {
}

/* ------------------------------------------------------------------------- */

void NewtonSolverIterative::prepareLinearSystem(int ntry, int nnewt) {
    newton_try_ = ntry;
    i_newton_ = nnewt;
    if (nnewt == -1) {
        throw NewtonFailure(AMICI_NOT_IMPLEMENTED,
                            "Linear solver SPBCG does not support sensitivity "
                            "computation for steady state problems.");
    }

    // Get the Jacobian and its diagonal for preconditioning
    model_->fJ(*t_, 0.0, *x_, dx_, xdot_, ns_J_.get());
    ns_J_.refresh();
    model_->fJDiag(*t_, ns_Jdiag_, 0.0, *x_, dx_);

    // Ensure positivity of entries in ns_Jdiag
    ns_p_.set(1.0);
    ns_Jdiag_.abs();
    N_VCompare(1e-15, ns_Jdiag_.getNVector(), ns_tmp_.getNVector());
    linearSum(-1.0, ns_tmp_, 1.0, ns_p_, ns_tmp_);
    linearSum(1.0, ns_Jdiag_, 1.0, ns_tmp_, ns_Jdiag_);
}

/* ------------------------------------------------------------------------- */

void NewtonSolverIterative::prepareLinearSystemB(int ntry, int nnewt) {
    newton_try_ = ntry;
    i_newton_ = nnewt;
    if (nnewt == -1) {
        throw AmiException("Linear solver SPBCG does not support sensitivity "
                           "computation for steady state problems.");
    }

    // Get the Jacobian and its diagonal for preconditioning
    model_->fJB(*t_, 0.0, *x_, dx_, xB_, dxB_, xdot_, ns_J_.get());
    ns_J_.refresh();
    // Get the diagonal and ensure negativity of entries is ns_J. Note that diag(JB) = -diag(J).
    model_->fJDiag(*t_, ns_Jdiag_, 0.0, *x_, dx_);

    ns_p_.set(1.0);
    ns_Jdiag_.abs();
    N_VCompare(1e-15, ns_Jdiag_.getNVector(), ns_tmp_.getNVector());
    linearSum(-1.0, ns_tmp_, 1.0, ns_p_, ns_tmp_);
    linearSum(1.0, ns_Jdiag_, 1.0, ns_tmp_, ns_Jdiag_);
    ns_Jdiag_.minus();
}

/* ------------------------------------------------------------------------- */

void NewtonSolverIterative::solveLinearSystem(AmiVector &rhs) {
    linsolveSPBCG(newton_try_, i_newton_, rhs);
    rhs.minus();
}

/* ------------------------------------------------------------------------- */

void NewtonSolverIterative::linsolveSPBCG(int /*ntry*/, int nnewt,
                                          AmiVector &ns_delta) {
    xdot_ = ns_delta;
    xdot_.minus();

    // Initialize for linear solve
    ns_p_.zero();
    ns_v_.zero();
    ns_delta.zero();
    ns_tmp_.zero();
    double rho = 1.0;
    double omega = 1.0;
    double alpha = 1.0;

    ns_J_.multiply(ns_Jv_, ns_delta);

    // ns_r = xdot - ns_Jv;
    linearSum(-1.0, ns_Jv_, 1.0, xdot_, ns_r_);
    ns_r_ /= ns_Jdiag_;
    ns_rt_ = ns_r_;

    for (int i_linstep = 0; i_linstep < max_lin_steps_; i_linstep++) {
        // Compute factors
        double rho1 = rho;
        rho = dotProd(ns_rt_, ns_r_);
        double beta = rho * alpha / (rho1 * omega);

        // ns_p = ns_r + beta * (ns_p - omega * ns_v);
        linearSum(1.0, ns_p_, -omega, ns_v_, ns_p_);
        linearSum(1.0, ns_r_, beta, ns_p_, ns_p_);

        // ns_v = J * ns_p
        ns_v_.zero();
        ns_J_.multiply(ns_v_, ns_p_);
        ns_v_ /= ns_Jdiag_;

        // Compute factor
        alpha = rho / dotProd(ns_rt_, ns_v_);

        // ns_h = ns_delta + alpha * ns_p;
        linearSum(1.0, ns_delta, alpha, ns_p_, ns_h_);
        // ns_s = ns_r - alpha * ns_v;
        linearSum(1.0, ns_r_, -alpha, ns_v_, ns_s_);

        // ns_t = J * ns_s
        ns_t_.zero();
        ns_J_.multiply(ns_t_, ns_s_);
        ns_t_ /= ns_Jdiag_;

        // Compute factor
        omega = dotProd(ns_t_, ns_s_) / dotProd(ns_t_, ns_t_);

        // ns_delta = ns_h + omega * ns_s;
        linearSum(1.0, ns_h_, omega, ns_s_, ns_delta);
        // ns_r = ns_s - omega * ns_t;
        linearSum(1.0, ns_s_, -omega, ns_t_, ns_r_);

        // Compute the (unscaled) residual
        ns_r_ *= ns_Jdiag_;
        double res = sqrt(dotProd(ns_r_, ns_r_));

        // Test convergence
        if (res < atol_) {
            // Write number of steps needed
            num_lin_steps_.at(nnewt) = i_linstep + 1;

            // Return success
            return;
        }

        // Scale back
        ns_r_ /= ns_Jdiag_;
    }
    throw NewtonFailure(AMICI_CONV_FAILURE, "linsolveSPBCG");
}


} // namespace amici
