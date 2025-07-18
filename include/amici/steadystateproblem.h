#ifndef AMICI_STEADYSTATE_PROBLEM_H
#define AMICI_STEADYSTATE_PROBLEM_H

#include <amici/defines.h>
#include <amici/model_state.h>
#include <amici/newton_solver.h>
#include <amici/vector.h>

namespace amici {

class ExpData;
class Solver;
class Model;
class BackwardProblem;
struct FwdSimWorkspace;
struct BwdSimWorkspace;

/**
 * @brief Computes the weighted root-mean-square norm.
 *
 * This class is used to compute the weighted root-mean-square of the residuals
 * and maintains its work space to avoid reallocation.
 */
class WRMSComputer {
  public:
    /**
     * @brief Constructor.
     * @param n The length of the vectors for which to compute the WRMS.
     * @param sunctx A SUNDIALS context for the NVector.
     * @param atol Absolute tolerance to compute error weights.
     * @param rtol Relative tolerance to compute error weights.
     * @param mask Mask for entries to include in the WRMS norm.
     * Positive value: include; non-positive value: exclude; empty: include all.
     */
    WRMSComputer(
        int const n, SUNContext const sunctx, realtype const atol,
        realtype const rtol, AmiVector mask
    )
        : ewt_(n, sunctx)
        , rtol_(rtol)
        , atol_(atol)
        , mask_(mask) {}

    /**
     * @brief Compute the weighted root-mean-square of the residuals.
     * @param x Vector to compute the WRMS for.
     * @param x_ref The reference vector from which to compute the error
     * weights.
     * @return The WRMS norm.
     */
    realtype wrms(AmiVector const& x, AmiVector const& x_ref);

  private:
    /** Error weights for the residuals. */
    AmiVector ewt_;
    /** Relative tolerance to compute error weights. */
    realtype rtol_;
    /** Absolute tolerance to compute error weights. */
    realtype atol_;
    /**
     * Mask for entries to include in the WRMS norm.
     * Positive value: include; non-positive value: exclude; empty: include all.
     */
    AmiVector mask_;
};

/**
 * @brief Implements Newton's method for finding steady states.
 *
 * See also:
 *  Lines et al. (2019), IFAC-PapersOnLine 52 (26): 32–37.
 *  https://doi.org/10.1016/j.ifacol.2019.12.232
 */
class NewtonsMethod {
  public:
    /**b
     * @brief Constructor.
     * @param model Number of solver states (nx_solver).
     * @param solver NewtonSolver instance to compute the Newton step.
     * Expected to be correctly initialized.
     * @param sunctx A SUNDIALS context for the NVector.
     * @param max_steps
     * @param damping_factor_mode
     * @param damping_factor_lower_bound
     * @param check_delta
     */
    NewtonsMethod(
        gsl::not_null<Model*> model, SUNContext sunctx,
        gsl::not_null<NewtonSolver*> solver,
        NewtonDampingFactorMode damping_factor_mode,
        realtype damping_factor_lower_bound, int max_steps, bool check_delta
    );

    /**
     * @brief Run the Newton solver iterations and checks for convergence
     * to steady state.
     * @param xdot Time derivative of the state vector `state.x`.
     * @param state SimulationState instance containing the current state.
     * @param wrms_computer WRMSComputer instance to compute the WRMS norm.
     */
    void
    run(AmiVector& xdot, DEStateView const& state, WRMSComputer& wrms_computer);

    /**
     * @brief Compute the Newton step for the current state_.x and xdot and
     * store it in delta_.
     * @param xdot Time derivative of the state vector `state.x`.
     * @param state SimulationState instance containing the current state.
     */
    void compute_step(AmiVector const& xdot, DEStateView const& state);

    /**
     * @brief Get the last Newton step.
     * @return Newton step
     */
    [[nodiscard]] AmiVector const& get_delta() const { return delta_; }

    /**
     * @brief Get the number of steps taken in the current iteration.
     * @return Number of steps taken.
     */
    [[nodiscard]] int get_num_steps() const { return i_step; }

    /**
     * @brief Get the current WRMS norm.
     * @return The current WRMS norm.
     */
    [[nodiscard]] realtype get_wrms() const { return wrms_; }

  private:
    /**
     * @brief Update the damping factor gamma that determines step size.
     *
     * @param step_successful flag indicating whether the previous step was
     * successful
     * @param gamma reference to the damping factor that is updated
     * @return boolean flag indicating whether search direction should be
     * updated (true) or the same direction should be retried with the updated
     * dampening (false)
     */

    bool update_damping_factor(bool step_successful, double& gamma);

    /**
     * @brief Compute the weighted root-mean-square of the residuals.
     * @param xdot
     * @param state
     * @param wrms_computer
     * @return WRMS norm.
     */
    realtype compute_wrms(
        AmiVector const& xdot, DEStateView const& state,
        WRMSComputer& wrms_computer
    );

    /**
     * @brief Check for convergence.
     *
     * Check if NewtonsMethod::wrms_ is below the convergence threshold,
     * make the state non-negative if requested, and recompute and check
     * the WRMS norm again.
     *
     * @param xdot
     * @param state
     * @param wrms_computer
     * @return Whether convergence has been reached.
     */
    bool has_converged(
        AmiVector& xdot, DEStateView const& state, WRMSComputer& wrms_computer
    );

    static constexpr realtype conv_thresh = 1.0;

    /** Pointer to the model instance. */
    gsl::not_null<Model*> model_;

    /** Maximum number of iterations. */
    int max_steps_{0};

    /** damping factor flag */
    NewtonDampingFactorMode damping_factor_mode_{NewtonDampingFactorMode::on};

    /** damping factor lower bound */
    realtype damping_factor_lower_bound_{1e-8};

    /**
     * Whether to check the Newton step (delta) or the right-hand side (xdot)
     * during the convergence check.
     */
    bool check_delta_;

    /** Pointer to the Newton solver instance to compute the Newton step. */
    gsl::not_null<NewtonSolver*> solver_;

    /** Newton step (size: nx_solver). */
    AmiVector delta_;

    /** Previous Newton step (size: nx_solver). */
    AmiVector delta_old_;

    /** Newton step (size: nx_solver). */
    AmiVector x_old_;

    /**
     * WRMS norm based on the current state and delta or xdot
     * (depending on `check_delta_`).
     */
    realtype wrms_ = INFINITY;

    /** The current number of Newton iterations. */
    int i_step = 0;
};

} // namespace amici
#endif // AMICI_STEADYSTATE_PROBLEM_H
