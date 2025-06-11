#ifndef AMICI_STEADYSTATE_PROBLEM_H
#define AMICI_STEADYSTATE_PROBLEM_H

#include <amici/defines.h>
#include <amici/model_state.h>
#include <amici/newton_solver.h>
#include <amici/vector.h>

#include <memory>

namespace amici {

class ExpData;
class Solver;
class Model;
class BackwardProblem;

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
        int n, SUNContext sunctx, realtype atol, realtype rtol, AmiVector mask
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
    /**
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
    run(AmiVector& xdot, SimulationState& state, WRMSComputer& wrms_computer);

    /**
     * @brief Compute the Newton step for the current state_.x and xdot and
     * store it in delta_.
     * @param xdot Time derivative of the state vector `state.x`.
     * @param state SimulationState instance containing the current state.
     */
    void compute_step(AmiVector const& xdot, SimulationState const& state);

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
        AmiVector const& xdot, SimulationState const& state,
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
        AmiVector& xdot, SimulationState& state, WRMSComputer& wrms_computer
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

/**
 * @brief The SteadystateProblem class solves a steady-state problem using
 * Newton's method and falls back to integration on failure.
 */
class SteadystateProblem {
  public:
    /**
     * @brief Constructor
     *
     * @param solver Solver instance
     * @param model Model instance
     */
    explicit SteadystateProblem(Solver const& solver, Model& model);

    /**
     * @brief Compute the steady state in the forward case.
     *
     * Tries to determine the steady state of the ODE system and computes
     * steady state sensitivities if requested.
     *
     * @param solver The solver instance
     * @param model The model instance
     * @param it Index of the current output time point.
     */
    void workSteadyStateProblem(Solver const& solver, Model& model, int it);

    /**
     * @brief Compute the gradient via adjoint steady state sensitivities.
     *
     * Integrates over the adjoint state backward in time by solving a linear
     * system of equations, which gives the analytical solution.
     *
     * @param solver The solver instance
     * @param model The model instance
     * @param xB0 Initial adjoint state vector.
     * @param is_preeq Flag indicating whether this is a preequilibration.
     */
    void workSteadyStateBackwardProblem(
        Solver const& solver, Model& model, AmiVector const& xB0, bool is_preeq
    );

    /**
     * @brief Return the stored SimulationState.
     * @return stored SimulationState
     */
    [[nodiscard]] SimulationState const& getFinalSimulationState() const {
        return state_;
    };

    /**
     * @brief Return the quadratures from pre- or postequilibration
     * @return xQB Vector with quadratures
     */
    [[nodiscard]] AmiVector const& getEquilibrationQuadratures() const {
        return xQB_;
    }
    /**
     * @brief Return state at steady state
     * @return x
     */
    [[nodiscard]] AmiVector const& getState() const { return state_.x; };

    /**
     * @brief Return state sensitivity at steady state
     * @return sx
     */
    [[nodiscard]] AmiVectorArray const& getStateSensitivity() const {
        return state_.sx;
    };

    /**
     * @brief Get the CPU time taken to solve the forward problem.
     * @return The CPU time in milliseconds.
     */
    [[nodiscard]] double getCPUTime() const { return cpu_time_; }

    /**
     * @brief Get the CPU time taken to solve the backward problem.
     * @return The CPU time in milliseconds.
     */
    [[nodiscard]] double getCPUTimeB() const { return cpu_timeB_; }

    /**
     * @brief Get the steady state computation status.
     * @return Execution status of the different approaches
     * [newton, simulation, newton].
     */
    [[nodiscard]] std::vector<SteadyStateStatus> const&
    getSteadyStateStatus() const {
        return steady_state_status_;
    }

    /**
     * @brief Get model time at which steady state was found through simulation.
     * @return Time at which steady state was found (model time units).
     */
    [[nodiscard]] realtype getSteadyStateTime() const { return state_.t; }

    /**
     * @brief Get the weighted root mean square of the residuals.
     * @return The weighted root-mean-square of the residuals.
     */
    [[nodiscard]] realtype getResidualNorm() const { return wrms_; }

    /**
     * @brief Get the number of steps taken to find the steady state.
     * @return Number of steps taken to find the steady state as
     * [newton, simulation, newton].
     */
    [[nodiscard]] std::vector<int> const& getNumSteps() const {
        return numsteps_;
    }

    /**
     * @brief Get the number of steps taken to find the steady state in the
     * adjoint case.
     * @return Number of steps.
     */
    [[nodiscard]] int getNumStepsB() const { return numstepsB_; }

    /**
     * @brief Compute adjoint updates dJydx according to the provided model and
     * data.
     * @param model Model instance
     * @param edata Experimental data
     * @param dJydx output argument for dJydx
     */
    void getAdjointUpdates(
        Model& model, ExpData const& edata, std::vector<realtype>& dJydx
    );

    /**
     * @brief Return the adjoint state
     * @return xB adjoint state
     */
    [[nodiscard]] AmiVector const& getAdjointState() const { return xB_; }

    /**
     * @brief Get the adjoint quadratures (xQB).
     * @return xQB
     */
    [[nodiscard]] AmiVector const& getAdjointQuadrature() const { return xQB_; }

    /**
     * @brief Accessor for hasQuadrature_
     * @return hasQuadrature_
     */
    [[nodiscard]] bool hasQuadrature() const { return hasQuadrature_; }

    /**
     * @brief Check, whether any approach to find the steady state was
     * successful.
     * @return Whether any approach to find the steady state was successful.
     */
    [[nodiscard]] bool checkSteadyStateSuccess() const;

  private:
    /**
     * @brief Handle the computation of the steady state.
     *
     * Throws an AmiException if no steady state was found.
     *
     * @param solver Solver instance.
     * @param model Model instance.
     * @param it Index of the current output time point.
     */
    void findSteadyState(Solver const& solver, Model& model, int it);

    /**
     * @brief Try to determine the steady state by using Newton's method.
     * @param model Model instance.
     * @param newton_retry Flag indicating whether Newton's method is being
     * relaunched.
     */
    void findSteadyStateByNewtonsMethod(Model& model, bool newton_retry);

    /**
     * @brief Try to determine the steady state by using forward simulation.
     * @param solver Solver instance.
     * @param model Model instance.
     * @param it Index of the current output time point.
     * @return SteadyStateStatus indicating whether the steady state was found
     * successfully, or if it failed.
     */
    SteadyStateStatus
    findSteadyStateBySimulation(Solver const& solver, Model& model, int it);

    /**
     * @brief Compute quadratures in adjoint mode
     * @param solver Solver instance.
     * @param model Model instance.
     */
    void computeSteadyStateQuadrature(Solver const& solver, Model& model);

    /**
     * @brief Compute the quadrature in steady state backward mode by
     * solving the linear system defined by the backward Jacobian.
     * @param model Model instance.
     */
    void getQuadratureByLinSolve(Model& model);

    /**
     * @brief Computes the quadrature in steady state backward mode by
     * numerical integration of xB forward in time.
     * @param solver Solver instance.
     * @param model Model instance.
     */
    void getQuadratureBySimulation(Solver const& solver, Model& model);

    /**
     * @brief Store state and throw an exception if equilibration failed
     * @param tried_newton_1 Whether any Newton step was attempted before
     * simulation
     * @param tried_simulation Whether simulation was attempted
     * @param tried_newton_2 Whether any Newton step was attempted after
     * simulation
     */
    [[noreturn]] void handleSteadyStateFailure(
        bool tried_newton_1, bool tried_simulation, bool tried_newton_2
    ) const;

    /**
     * @brief Get whether state sensitivities need to be computed.
     *
     * Checks depending on the status of the Newton solver,
     * solver settings, and the model, whether state sensitivities
     * still need to be computed (via a linear system solve or integration).
     * @param model Model instance.
     * @param solver Solver instance.
     * @param it Index of the current output time point.
     * @param context SteadyStateContext giving the situation for the flag
     * @return Whether sensitivities have to be computed.
     */
    bool requires_state_sensitivities(
        Model const& model, Solver const& solver, int it,
        SteadyStateContext context
    ) const;

    /**
     * @brief Checks steady-state convergence for state variables
     * @param model Model instance
     * @return weighted root mean squared residuals of the RHS
     */
    realtype getWrmsState(Model& model);

    /**
     * @brief Checks convergence for state sensitivities
     * @param model Model instance
     * @return weighted root mean squared residuals of the RHS
     */
    realtype getWrmsFSA(Model& model);

    /**
     * @brief Launch simulation if Newton solver or linear system solve
     * fail or are disabled.
     * @param solver Solver instance.
     * @param model Model instance.
     * simulation.
     */
    void runSteadystateSimulationFwd(Solver const& solver, Model& model);

    /**
     * @brief Launch backward simulation if Newton solver or linear system solve
     * fail or are disabled.
     * @param solver Solver instance.
     * @param model Model instance.
     */
    void runSteadystateSimulationBwd(Solver const& solver, Model& model);

    /**
     * @brief Create and initialize a CVodeSolver instance for
     * preequilibration simulation.
     * @param solver Solver instance
     * @param model Model instance.
     * @param forwardSensis flag switching on integration with FSA
     * @param backward flag switching on quadrature computation
     * @return A unique pointer to the created Solver instance.
     */
    std::unique_ptr<Solver> createSteadystateSimSolver(
        Solver const& solver, Model& model, bool forwardSensis, bool backward
    ) const;

    /**
     * @brief Initialize forward computation
     * @param it Index of the current output time point.
     * @param solver pointer to the solver object
     * @param model pointer to the model object
     */
    void initializeForwardProblem(int it, Solver const& solver, Model& model);

    /**
     * @brief Update member variables to indicate that state_.x has been
     * updated and xdot_, delta_, etc. need to be recomputed.
     */
    void flagUpdatedState();

    /**
     * @brief Retrieve simulation sensitivities from the provided solver and
     * set the corresponding flag to indicate they are up to date
     * @param solver simulation solver instance
     */
    void updateSensiSimulation(Solver const& solver);

    /**
     * @brief Compute the right-hand side for the current state_.x and set the
     * corresponding flag to indicate xdot_ is up to date.
     * @param model model instance
     */
    void updateRightHandSide(Model& model);

    /** WRMS computer for x */
    WRMSComputer wrms_computer_x_;
    /** WRMS computer for xQB */
    WRMSComputer wrms_computer_xQB_;
    /** WRMS computer for sx */
    WRMSComputer wrms_computer_sx_;
    /** old state vector */
    AmiVector x_old_;
    /** time derivative state vector */
    AmiVector xdot_;
    /** state differential sensitivities */
    AmiVectorArray sdx_;
    /** adjoint state vector */
    AmiVector xB_;
    /** integral over adjoint state vector */
    AmiVector xQ_;
    /** quadrature state vector */
    AmiVector xQB_;
    /** time-derivative of quadrature state vector */
    AmiVector xQBdot_;

    /** weighted root-mean-square error */
    realtype wrms_{NAN};

    SimulationState state_;

    /** stores diagnostic information about employed number of steps */
    std::vector<int> numsteps_{std::vector<int>(3, 0)};

    /** The employed number of backward steps */
    int numstepsB_{0};

    /** CPU time for solving the forward problem (milliseconds) */
    double cpu_time_{0.0};

    /** CPU time for solving the backward problem (milliseconds) */
    double cpu_timeB_{0.0};

    /** flag indicating whether backward mode was run */
    bool hasQuadrature_{false};

    /**
     * Execution status of the different approaches
     * [newton, simulation, newton] (length = 3)
     */
    std::vector<SteadyStateStatus> steady_state_status_;

    /** Newton solver */
    NewtonSolver newton_solver_;

    /** Newton's method for finding steady states */
    NewtonsMethod newtons_method_;

    /**
     * Whether the Newton step should be used instead of xdot for convergence
     * checks during simulation and Newton's method.
     */
    bool newton_step_conv_{false};
    /**
     * whether sensitivities should be checked for convergence to steady state
     */
    bool check_sensi_conv_{true};

    /** flag indicating whether xdot_ has been computed for the current state */
    bool xdot_updated_{false};
    /**
     * flag indicating whether simulation sensitivities have been retrieved for
     * the current state
     */
    bool sensis_updated_{false};
};

} // namespace amici
#endif // AMICI_STEADYSTATE_PROBLEM_H
