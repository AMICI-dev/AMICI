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
    explicit SteadystateProblem(Solver const& solver, Model const& model);

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
     * @param bwd The backward problem instance
     */
    void workSteadyStateBackwardProblem(
        Solver const& solver, Model& model, BackwardProblem const* bwd
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
     * @brief Accessor for dJydx
     * @return dJydx
     */
    [[nodiscard]] std::vector<realtype> const& getDJydx() const {
        return dJydx_;
    }

    /**
     * @brief Get the CPU time taken to solvethe forward problem.
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
     * @return The weighted root mean square of the residuals.
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
     */
    void getAdjointUpdates(Model& model, ExpData const& edata);

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
     * Throws an AmiException, if no steady state was found.
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
     */
    void
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
     * @brief Checks depending on the status of the Newton solver,
     * solver settings, and the model, whether state sensitivities
     * still need to be computed via a linear system solve or stored
     * @param model Model instance.
     * @param solver Solver instance.
     * @param it Index of the current output time point.
     * @param context SteadyStateContext giving the situation for the flag
     * @return Whether sensitivities have to be computed.
     */
    bool getSensitivityFlag(
        Model const& model, Solver const& solver, int it,
        SteadyStateContext context
    ) const;

    /**
     * @brief Checks convergence for state or adjoint quadratures, depending on
     * sensi method
     * @param model Model instance
     * @param sensi_method sensitivity method
     * @return weighted root mean squared residuals of the RHS
     */
    realtype getWrms(Model& model, SensitivityMethod sensi_method);

    /**
     * @brief Checks convergence for state sensitivities
     * @param model Model instance
     * @return weighted root mean squared residuals of the RHS
     */
    realtype getWrmsFSA(Model& model);

    /**
     * @brief Run the Newton solver iterations and checks for convergence
     * to steady state.
     * @param model Model instance.
     * @param newton_retry flag indicating if Newton solver is rerun
     */
    void applyNewtonsMethod(Model& model, bool newton_retry);

    /**
     * @brief Launch simulation if Newton solver or linear system solve
     * fail or are disabled.
     * @param solver Solver instance.
     * @param model Model instance.
     * @param backward Whether to run the forward (`false`) or backward (`true`)
     * simulation.
     */
    void
    runSteadystateSimulation(Solver const& solver, Model& model, bool backward);

    /**
     * @brief Create and initialize a CVodeSolver instance for
     * preequilibration simulation.
     * @param solver Solver instance
     * @param model Model instance.
     * @param forwardSensis flag switching on integration with FSA
     * @param backward flag switching on quadratures computation
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
     * @brief Initialize backward computation.
     * @param solver Solver instance
     * @param model Model instance.
     * @param bwd pointer to backward problem
     * @return flag indicating whether backward computation to be carried out
     */
    bool initializeBackwardProblem(
        Solver const& solver, Model& model, BackwardProblem const* bwd
    );

    /**
     * @brief Ensure state positivity if requested, and repeat the convergence
     * check if necessary.
     * @param model Model instance.
     */
    bool makePositiveAndCheckConvergence(Model& model);

    /**
     * @brief Update the damping factor gamma that determines step size.
     *
     * @param step_successful flag indicating whether the previous step was
     * successful
     * @return boolean flag indicating whether search direction should be
     * updated (true) or the same direction should be retried with the updated
     * dampening (false)
     */
    bool updateDampingFactor(bool step_successful);

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

    /**
     * @brief Compute the Newton step for the current state_.x and set the
     * corresponding flag to indicate delta_ is up to date.
     * @param model Model instance
     */
    void getNewtonStep(Model& model);

    /** Newton step (size: nx_solver). */
    AmiVector delta_;
    /** previous newton step (size: nx_solver). */
    AmiVector delta_old_;
    /** error weights for solver state, dimension nx_solver */
    AmiVector ewt_;
    /** error weights for backward quadratures, dimension nplist() */
    AmiVector ewtQB_;
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
    /** NVector around Model::steadystate_mask_ */
    AmiVector steadystate_mask_;

    /** maximum number of steps for Newton solver for allocating numlinsteps */
    int max_steps_{0};

    /** weighted root-mean-square error */
    realtype wrms_{NAN};

    /** state derivative of data likelihood
     * (dimension nJ x nx x nt, ordering =?) */
    std::vector<realtype> dJydx_;

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

    /** Stepsize for the the Newton step. */
    double gamma_{1.0};

    /** Execution status of the different approaches
     * [newton, simulation, newton] (length = 3)
     */
    std::vector<SteadyStateStatus> steady_state_status_;

    /** absolute tolerance for convergence check (state)*/
    realtype atol_{NAN};
    /** relative tolerance for convergence check (state)*/
    realtype rtol_{NAN};
    /** absolute tolerance for convergence check (state sensi)*/
    realtype atol_sensi_{NAN};
    /** relative tolerance for convergence check (state sensi)*/
    realtype rtol_sensi_{NAN};
    /** absolute tolerance for convergence check (quadratures)*/
    realtype atol_quad_{NAN};
    /** relative tolerance for convergence check (quadratures)*/
    realtype rtol_quad_{NAN};

    /** Newton solver */
    NewtonSolver newton_solver_;

    /** damping factor flag */
    NewtonDampingFactorMode damping_factor_mode_{NewtonDampingFactorMode::on};
    /** damping factor lower bound */
    realtype damping_factor_lower_bound_{1e-8};
    /** whether newton step should be used for convergence steps */
    bool newton_step_conv_{false};
    /** whether sensitivities should be checked for convergence to steady state
     */
    bool check_sensi_conv_{true};

    /** flag indicating whether xdot_ has been computed for the current state */
    bool xdot_updated_{false};
    /** flag indicating whether delta_ has been computed for the current state
     */
    bool delta_updated_{false};
    /** flag indicating whether simulation sensitivities have been retrieved for
     * the current state */
    bool sensis_updated_{false};
};

} // namespace amici
#endif // AMICI_STEADYSTATE_PROBLEM_H
