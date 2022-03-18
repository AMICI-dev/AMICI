#include "amici/newton_solver.h"

#include "amici/model.h"
#include "amici/solver.h"

#include "sunlinsol/sunlinsol_klu.h" // sparse solver
#include "sunlinsol/sunlinsol_dense.h" // dense solver

#include <cstring>
#include <ctime>
#include <cmath>

namespace amici {

NewtonSolver::NewtonSolver(Model *model)
    : xdot_(model->nx_solver), x_(model->nx_solver), dxB_(model->nx_solver) {}

/* ------------------------------------------------------------------------- */

std::unique_ptr<NewtonSolver> NewtonSolver::getSolver(
    const Solver &simulationSolver, Model *model) {

    std::unique_ptr<NewtonSolver> solver;

    switch (simulationSolver.getLinearSolver()) {

    /* DIRECT SOLVERS */
    case LinearSolver::dense:
        solver.reset(new NewtonSolverDense(model));
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
        throw NewtonFailure(AMICI_NOT_IMPLEMENTED, "getSolver");

    case LinearSolver::SPTFQMR:
        throw NewtonFailure(AMICI_NOT_IMPLEMENTED, "getSolver");

    /* SPARSE SOLVERS */
    case LinearSolver::SuperLUMT:
        throw NewtonFailure(AMICI_NOT_IMPLEMENTED, "getSolver");
    case LinearSolver::KLU:
        solver.reset(new NewtonSolverSparse(model));
        break;
    default:
        throw NewtonFailure(AMICI_NOT_IMPLEMENTED, "getSolver");
    }
    solver->max_lin_steps_ = simulationSolver.getNewtonMaxLinearSteps();
    solver->max_steps = simulationSolver.getNewtonMaxSteps();
    solver->damping_factor_mode_ = simulationSolver.getNewtonDampingFactorMode();
    solver->damping_factor_lower_bound =
        simulationSolver.getNewtonDampingFactorLowerBound();

    return solver;
}

/* ------------------------------------------------------------------------- */

void NewtonSolver::getStep(int ntry, int nnewt, AmiVector &delta,
                           Model *model, const SimulationState &state) {
    prepareLinearSystem(ntry, nnewt, model, state);

    delta.minus();
    solveLinearSystem(delta);
}

/* ------------------------------------------------------------------------- */

void NewtonSolver::computeNewtonSensis(AmiVectorArray &sx, Model *model,
                                       const SimulationState &state) {
    prepareLinearSystem(0, -1, model, state);
    model->fdxdotdp(state.t, state.x, state.dx);

    if (model->pythonGenerated) {
        for (int ip = 0; ip < model->nplist(); ip++) {
            N_VConst(0.0, sx.getNVector(ip));
            model->get_dxdotdp_full().scatter(model->plist(ip), -1.0, nullptr,
                                               gsl::make_span(sx.getNVector(ip)),
                                               0, nullptr, 0);

            solveLinearSystem(sx[ip]);
        }
    } else {
        for (int ip = 0; ip < model->nplist(); ip++) {
            for (int ix = 0; ix < model->nx_solver; ix++)
                sx.at(ix,ip) = -model->get_dxdotdp().at(ix, ip);

            solveLinearSystem(sx[ip]);
        }
    }
}

NewtonSolverDense::NewtonSolverDense(Model *model)
    : NewtonSolver(model), Jtmp_(model->nx_solver, model->nx_solver),
      linsol_(SUNLinSol_Dense(x_.getNVector(), Jtmp_.get())) {
    auto status = SUNLinSolInitialize_Dense(linsol_);
    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolInitialize_Dense");
}

void NewtonSolverDense::prepareLinearSystem(int  /*ntry*/, int  /*nnewt*/,
                                            Model *model,
                                            const SimulationState &state) {
    model->fJ(state.t, 0.0, state.x, state.dx, xdot_, Jtmp_.get());
    Jtmp_.refresh();
    auto status = SUNLinSolSetup_Dense(linsol_, Jtmp_.get());
    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolSetup_Dense");
}

void NewtonSolverDense::prepareLinearSystemB(int  /*ntry*/, int  /*nnewt*/,
                                             Model *model,
                                             const SimulationState &state) {
    model->fJB(state.t, 0.0, state.x, state.dx, xB_, dxB_, xdot_, Jtmp_.get());
    Jtmp_.refresh();
    auto status = SUNLinSolSetup_Dense(linsol_, Jtmp_.get());
    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolSetup_Dense");
}

void NewtonSolverDense::solveLinearSystem(AmiVector &rhs) {
    auto status = SUNLinSolSolve_Dense(linsol_, Jtmp_.get(),
                                       rhs.getNVector(), rhs.getNVector(),
                                       0.0);
    Jtmp_.refresh();
    // last argument is tolerance and does not have any influence on result

    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolSolve_Dense");
}

/* does not need reinitialization */
void NewtonSolverDense::reinitialize() {};

NewtonSolverDense::~NewtonSolverDense() {
    if(linsol_)
        SUNLinSolFree_Dense(linsol_);
}

NewtonSolverSparse::NewtonSolverSparse(Model *model)
    : NewtonSolver(model),
      Jtmp_(model->nx_solver, model->nx_solver, model->nnz, CSC_MAT),
      linsol_(SUNKLU(x_.getNVector(), Jtmp_.get())) {
    auto status = SUNLinSolInitialize_KLU(linsol_);
    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolInitialize_KLU");
}

void NewtonSolverSparse::prepareLinearSystem(int  /*ntry*/, int  /*nnewt*/,
                                             Model *model,
                                             const SimulationState &state) {
    /* Get sparse Jacobian */
    model->fJSparse(state.t, 0.0, state.x, state.dx, xdot_, Jtmp_.get());
    Jtmp_.refresh();
    auto status = SUNLinSolSetup_KLU(linsol_, Jtmp_.get());
    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolSetup_KLU");
}

void NewtonSolverSparse::prepareLinearSystemB(int  /*ntry*/, int  /*nnewt*/,
                                              Model *model,
                                              const SimulationState &state) {
    /* Get sparse Jacobian */
    model->fJSparseB(state.t, 0.0, state.x, state.dx, xB_, dxB_, xdot_,
                     Jtmp_.get());
    Jtmp_.refresh();
    auto status = SUNLinSolSetup_KLU(linsol_, Jtmp_.get());
    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolSetup_KLU");
}

void NewtonSolverSparse::solveLinearSystem(AmiVector &rhs) {
    /* Pass pointer to the linear solver */
    auto status = SUNLinSolSolve_KLU(linsol_, Jtmp_.get(),
                                     rhs.getNVector(), rhs.getNVector(), 0.0);
    // last argument is tolerance and does not have any influence on result

    if(status != AMICI_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolSolve_KLU");
}

void NewtonSolverSparse::reinitialize() {
    /* partial reinitialization, don't need to reallocate Jtmp_ */
    auto status = SUNLinSol_KLUReInit(linsol_, Jtmp_.get(), Jtmp_.capacity(),
                                      SUNKLU_REINIT_PARTIAL);
    if(status != SUNLS_SUCCESS)
        throw NewtonFailure(status, "SUNLinSol_KLUReInit");
}

NewtonSolverSparse::~NewtonSolverSparse() {
    if(linsol_)
        SUNLinSolFree_KLU(linsol_);
}


} // namespace amici
