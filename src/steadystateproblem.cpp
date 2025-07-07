#include "amici/steadystateproblem.h"
#include "amici/backwardproblem.h"
#include "amici/defines.h"
#include "amici/edata.h"
#include "amici/model.h"
#include "amici/newton_solver.h"
#include "amici/solver.h"

#include <cmath>
#include <cvodes/cvodes.h>
#include <sundials/sundials_dense.h>

namespace amici {

/**
 * @brief Compute the weighted root-mean-square norm of xdot.
 *
 * The weights are computed according to x:
 * w_i = 1 / ( rtol * x_i + atol )
 * @param x current state (sx[ip] for sensitivities)
 * @param xdot current rhs (sxdot[ip] for sensitivities)
 * @param mask mask for state variables to include in WRMS norm.
 * Positive value: include; non-positive value: exclude; empty: include all.
 * @param atol absolute tolerance
 * @param rtol relative tolerance
 * @param ewt error weight vector
 * @return root-mean-square norm
 */
realtype getWrmsNorm(
    AmiVector const& x, AmiVector const& xdot, AmiVector const& mask,
    realtype atol, realtype rtol, AmiVector& ewt
) {
    // ewt = x
    N_VAbs(const_cast<N_Vector>(x.getNVector()), ewt.getNVector());
    // ewt *= rtol
    N_VScale(rtol, ewt.getNVector(), ewt.getNVector());
    // ewt += atol
    N_VAddConst(ewt.getNVector(), atol, ewt.getNVector());
    // ewt = 1/ewt (ewt = 1/(rtol*x+atol))
    N_VInv(ewt.getNVector(), ewt.getNVector());

    // wrms = sqrt(sum((xdot/ewt)**2)/n) where n = size of state vector
    if (mask.getLength()) {
        return N_VWrmsNormMask(
            const_cast<N_Vector>(xdot.getNVector()), ewt.getNVector(),
            const_cast<N_Vector>(mask.getNVector())
        );
    }
    return N_VWrmsNorm(
        const_cast<N_Vector>(xdot.getNVector()), ewt.getNVector()
    );
}

realtype WRMSComputer::wrms(AmiVector const& x, AmiVector const& x_ref) {
    return getWrmsNorm(x_ref, x, mask_, atol_, rtol_, ewt_);
}

NewtonsMethod::NewtonsMethod(
    gsl::not_null<Model*> model, SUNContext sunctx,
    gsl::not_null<NewtonSolver*> solver,
    NewtonDampingFactorMode damping_factor_mode,
    realtype damping_factor_lower_bound, int max_steps, bool check_delta
)
    : model_(model)
    , max_steps_(max_steps)
    , damping_factor_mode_(damping_factor_mode)
    , damping_factor_lower_bound_(damping_factor_lower_bound)
    , check_delta_(check_delta)
    , solver_(solver)
    , delta_(model->nx_solver, sunctx)
    , delta_old_(model->nx_solver, sunctx)
    , x_old_(model->nx_solver, sunctx) {}

void NewtonsMethod::run(
    AmiVector& xdot, DEStateView const& state, WRMSComputer& wrms_computer
) {
    i_step = 0;

    if (model_->nx_solver == 0) {
        wrms_ = 0.0;
        return;
    }

    wrms_ = INFINITY;
    delta_.zero();

    // The Newton step size.
    double gamma{1.0};
    bool update_direction = true;

    wrms_ = compute_wrms(xdot, state, wrms_computer);
    bool converged = has_converged(xdot, state, wrms_computer);

    // Whether the step was successful
    bool step_successful = true;

    while (!converged && i_step < max_steps_) {
        if (step_successful) {
            // If new residuals are smaller than the old ones, update state
            x_old_.copy(state.x);
        }

        // If Newton steps are necessary, compute the initial search direction
        if (update_direction) {
            // compute the next step if not already done during the previous
            // delta-convergence check
            if (!check_delta_) {
                compute_step(xdot, state);
            }

            // we store delta_ here as later convergence checks may update it
            delta_old_.copy(delta_);
        }

        // Try step with new gamma_/delta_, evaluate rhs
        // x = x_old + delta_[old_] * gamma
        linearSum(
            1.0, x_old_, gamma, update_direction ? delta_ : delta_old_, state.x
        );
        model_->fxdot(state.t, state.x, state.dx, xdot);

        realtype wrms_tmp = compute_wrms(xdot, state, wrms_computer);
        step_successful = wrms_tmp < wrms_;
        if (step_successful) {
            wrms_ = wrms_tmp;
            converged = has_converged(xdot, state, wrms_computer);
        }

        update_direction = update_damping_factor(step_successful, gamma);
        ++i_step;
    }

    if (!converged)
        throw NewtonFailure(AMICI_TOO_MUCH_WORK, "applyNewtonsMethod");
}

void NewtonsMethod::compute_step(
    AmiVector const& xdot, DEStateView const& state
) {
    delta_.copy(xdot);
    solver_->getStep(delta_, *model_, state);
}

bool NewtonsMethod::update_damping_factor(bool step_successful, double& gamma) {
    if (damping_factor_mode_ != NewtonDampingFactorMode::on)
        return true;

    if (step_successful) {
        gamma = fmin(1.0, 2.0 * gamma);
    } else {
        gamma /= 4.0;
    }

    if (gamma < damping_factor_lower_bound_) {
        throw NewtonFailure(
            AMICI_DAMPING_FACTOR_ERROR,
            "Newton solver failed: the damping factor "
            "reached its lower bound"
        );
    }
    return step_successful;
}

realtype NewtonsMethod::compute_wrms(
    AmiVector const& xdot, DEStateView const& state, WRMSComputer& wrms_computer
) {
    if (check_delta_) {
        compute_step(xdot, state);
        return wrms_computer.wrms(delta_, state.x);
    } else {
        return wrms_computer.wrms(xdot, state.x);
    }
}

bool NewtonsMethod::has_converged(
    AmiVector& xdot, DEStateView const& state, WRMSComputer& wrms_computer
) {
    // pre-check convergence
    if (wrms_ >= conv_thresh)
        return false;

    if (!model_->get_any_state_nonnegative()) {
        // no constraints to check for
        return true;
    }

    // Ensure state positivity if requested,
    // and repeat the convergence check if necessary
    auto nonnegative = model_->getStateIsNonNegative();
    Expects(nonnegative.size() == state.x.getVector().size());
    auto state_modified = false;
    for (int ix = 0; ix < state.x.getLength(); ix++) {
        if (state.x[ix] < 0.0 && nonnegative[ix]) {
            state.x[ix] = 0.0;
            state_modified = true;
        }
    }
    if (!state_modified)
        return true;

    model_->fxdot(state.t, state.x, state.dx, xdot);
    wrms_ = compute_wrms(xdot, state, wrms_computer);

    return wrms_ < conv_thresh;
}

} // namespace amici
