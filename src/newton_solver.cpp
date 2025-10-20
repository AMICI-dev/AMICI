#include "amici/newton_solver.h"

#include <amici/amici.h>
#include <amici/model.h>
#include <amici/solver.h>
#include <amici/vector.h>

#include <sundials/sundials_config.h>  // roundoffs
#include <sunlinsol/sunlinsol_dense.h> // dense solver
#include <sunlinsol/sunlinsol_klu.h>   // sparse solver

namespace amici {

NewtonSolver::NewtonSolver(
    Model const& model, LinearSolver const linsol_type, SUNContext const sunctx
)
    : xdot_(model.nx_solver, sunctx)
    , x_(model.nx_solver, sunctx)
    , xB_(model.nJ * model.nx_solver, sunctx)
    , dxB_(model.nJ * model.nx_solver, sunctx) {
    try {
        // NOTE: when adding new linear solvers, make sure to also add them to
        // NewtonSolver::reinitialize and NewtonSolver::is_singular
        switch (linsol_type) {
        case LinearSolver::dense:
            linsol_.reset(new SUNLinSolDense(x_));
            break;
        case LinearSolver::KLU:
            linsol_.reset(new SUNLinSolKLU(x_, model.nnz, CSC_MAT));
            break;
        default:
            throw NewtonFailure(
                AMICI_NOT_IMPLEMENTED,
                "Invalid solver for steady state simulation"
            );
        }
    } catch (NewtonFailure const&) {
        throw;
    } catch (AmiException const& e) {
        throw NewtonFailure(AMICI_ERROR, e.what());
    }
}

void NewtonSolver::get_step(
    AmiVector& delta, Model& model, DEStateView const& state
) {
    prepare_linear_system(model, state);

    delta.minus();
    solve_linear_system(delta);
}

void NewtonSolver::compute_newton_sensis(
    AmiVectorArray& sx, Model& model, DEStateView const& state
) {
    prepare_linear_system(model, state);
    model.fdxdotdp(state.t, state.x, state.dx);

    if (model.get_logger() && is_singular(model, state)) {
        model.get_logger()->log(
            LogSeverity::warning, "NEWTON_JAC_SINGULAR",
            "Jacobian is singular at steadystate, "
            "sensitivities may be inaccurate."
        );
    }

    for (int ip = 0; ip < model.nplist(); ip++) {
        N_VConst(0.0, sx.get_nvector(ip));
        model.get_dxdotdp_full().scatter(
            model.plist(ip), -1.0, nullptr, gsl::make_span(sx.get_nvector(ip)),
            0, nullptr, 0
        );

        solve_linear_system(sx[ip]);
    }
}

void NewtonSolver::prepare_linear_system(
    Model& model, DEStateView const& state
) {
    auto& J = linsol_->getMatrix();
    if (J.matrix_id() == SUNMATRIX_SPARSE) {
        model.fJSparse(state.t, 0.0, state.x, state.dx, xdot_, J);
    } else {
        model.fJ(state.t, 0.0, state.x, state.dx, xdot_, J);
    }
    J.refresh();
    try {
        linsol_->setup();
    } catch (AmiException const& e) {
        throw NewtonFailure(AMICI_SINGULAR_JACOBIAN, e.what());
    }
}

void NewtonSolver::prepare_linear_system_b(
    Model& model, DEStateView const& state
) {
    auto& J = linsol_->getMatrix();
    if (J.matrix_id() == SUNMATRIX_SPARSE) {
        model.fJSparseB(state.t, 0.0, state.x, state.dx, xB_, dxB_, xdot_, J);
    } else {
        model.fJB(state.t, 0.0, state.x, state.dx, xB_, dxB_, xdot_, J);
    }
    J.refresh();
    try {
        linsol_->setup();
    } catch (AmiException const& e) {
        throw NewtonFailure(AMICI_SINGULAR_JACOBIAN, e.what());
    }
}

void NewtonSolver::solve_linear_system(AmiVector& rhs) {
    // last argument is tolerance and does not have any influence on result
    auto status = linsol_->solve(rhs.get_nvector(), rhs.get_nvector(), 0.0);
    if (status != SUN_SUCCESS)
        throw NewtonFailure(status, "SUNLinSolSolve");
}

void NewtonSolver::reinitialize() {
    // direct solvers do not need reinitialization
    if (dynamic_cast<SUNLinSolDense*>(linsol_.get())) {
        return;
    }
    if (auto const s = dynamic_cast<SUNLinSolKLU*>(linsol_.get())) {
        try {
            s->reinit(s->getMatrix().capacity(), SUNKLU_REINIT_PARTIAL);
            return;
        } catch (AmiException const&) {
            throw NewtonFailure(
                AMICI_ERROR, "Failed to reinitialize KLU solver"
            );
        }
    }
    throw AmiException(
        "Unhandled linear solver type in NewtonSolver::reinitialize."
    );
}

bool NewtonSolver::is_singular(Model& model, DEStateView const& state) const {
    if (auto s = dynamic_cast<SUNLinSolKLU const*>(linsol_.get())) {
        return s->is_singular();
    }

    // dense solver doesn't have any implementation for rcond/condest, so use
    // sparse solver interface, not the most efficient solution, but who is
    // concerned about speed and used the dense solver anyways ¯\_(ツ)_/¯
    NewtonSolver sparse_solver(
        model, LinearSolver::KLU, linsol_->get()->sunctx
    );
    sparse_solver.prepare_linear_system(model, state);
    return sparse_solver.is_singular(model, state);
}

} // namespace amici
