#ifndef AMICI_BACKWARDPROBLEM_H
#define AMICI_BACKWARDPROBLEM_H

#include "amici/defines.h"
#include "amici/forwardproblem.h"
#include "amici/vector.h"

#include <optional>
#include <vector>

namespace amici {
class ExpData;
class Solver;
class Model;
class ForwardProblem;
class SteadystateProblem;

/**
 * @brief The BwdSimWorkspace class is used to store temporary simulation
 * state during backward simulations.
 */
struct BwdSimWorkspace {
    /**
     * @brief Constructor
     * @param model The model for which to set up the workspace.
     * @param solver The solver for which to set up this workspace.
     */
    BwdSimWorkspace(
        gsl::not_null<Model*> model, gsl::not_null<Solver const*> solver
    );

    /** The model. */
    Model* model_;

    /** adjoint state vector (size: nx_solver)  */
    AmiVector xB_;
    /** differential adjoint state vector (size: nx_solver) */
    AmiVector dxB_;
    /** quadrature state vector (size: nJ x nplist, col-major) */
    AmiVector xQB_;

    /** array of number of found roots for a certain event type */
    std::vector<int> nroots_;
    /** array containing the time-points of discontinuities*/
    std::vector<Discontinuity> discs_;
    /** index of the backward problem */
    int which = 0;
};

/**
 * @brief The EventHandlingBwdSimulator class runs a backward simulation
 * and processes events and measurements general.
 */
class EventHandlingBwdSimulator {
  public:
    /**
     * @brief EventHandlingBwdSimulator constructor.
     * @param model The model to simulate.
     * @param solver The solver to use for the simulation.
     * @param ws The workspace to use for the simulation.
     */
    EventHandlingBwdSimulator(
        gsl::not_null<Model*> model, gsl::not_null<Solver const*> solver,
        gsl::not_null<BwdSimWorkspace*> ws
    )
        : model_(model)
        , solver_(solver)
        , ws_(ws) {}

    /**
     * @brief Run the simulation.
     *
     * It will run the backward simulation from the initial time of this period
     * to the final timepoint of this period, handling events
     * and data points as they occur.
     *
     * Expects the model and the solver to be set up, and `ws` to be initialized
     * for this period.
     *
     * @param t_start The initial time of this period.
     * @param t_end The final time of this period.
     * @param it The index of the timepoint in `timepoints` to start with.
     * @param timepoints The output timepoints or measurement timepoints of
     * this period. This must contain at least the final timepoint of this
     * period.
     * @param dJydx State-derivative of data likelihood. Must be non-null if
     * there are any data points in this period.
     * @param dJzdx State-derivative of event likelihood. Must be non-null if
     * the model has any event-observables.
     */
    void
    run(realtype t_start, realtype t_end, realtype it,
        std::vector<realtype> const& timepoints,
        std::vector<realtype> const* dJydx, std::vector<realtype> const* dJzdx);

  private:
    /**
     * @brief Execute everything necessary for the handling of events
     * for the backward problem
     * @param disc The discontinuity to handle
     * @param dJzdx State-derivative of event likelihood
     */
    void
    handleEventB(Discontinuity const& disc, std::vector<realtype> const* dJzdx);

    /**
     * @brief Execute everything necessary for the handling of data
     * points for the backward problems
     *
     * @param it index of data point
     * @param dJydx State-derivative of data likelihood
     */
    void handleDataPointB(int it, std::vector<realtype> const* dJydx);

    /**
     * @brief Compute the next timepoint to integrate to.
     *
     * This is the maximum of tdata and troot but also takes into account if
     * it<0 or iroot<0 where these expressions do not necessarily make sense.
     *
     * @param it index of next data point
     * @return tnext next timepoint
     */
    realtype getTnext(int it);

    /** The model to simulate. */
    Model* model_;

    /** The solver to use for the simulation. */
    Solver const* solver_;

    /** The workspace to use for the simulation. */
    gsl::not_null<BwdSimWorkspace*> ws_;

    /** current time */
    realtype t_{0};
};

/**
 * @brief The SteadyStateBackwardProblem class computes the adjoint state and
 * quadratures for pre- or post-equilibration.
 */
class SteadyStateBackwardProblem {
  public:
    /**
     * @brief SteadyStateBackwardProblem ctor
     * @param solver
     * @param model
     * @param final_state Final state from pre/post-equilibration forward
     * problem
     * @param ws Workspace for backward simulation
     */
    SteadyStateBackwardProblem(
        Solver const& solver, Model& model, SolutionState& final_state,
        gsl::not_null<BwdSimWorkspace*> ws
    );

    /**
     * @brief Compute the gradient via adjoint steady-state sensitivities.
     *
     * Integrates over the adjoint state backward in time by solving a linear
     * system of equations, which gives the analytical solution.
     *
     * Expects the workspace `ws_` to be initialized.
     *
     * The results will be written to the workspace `ws_`.
     *
     * @param t0 Initial time for the steady state simulation.
     */
    void run(realtype t0);

    /**
     * @brief Get the CPU time taken to solve the backward problem.
     * @return The CPU time in milliseconds.
     */
    [[nodiscard]] double getCPUTimeB() const { return cpu_timeB_; }

    /**
     * @brief Get the number of steps taken to find the steady state in the
     * adjoint case.
     * @return Number of steps.
     */
    [[nodiscard]] int getNumStepsB() const { return numstepsB_; }

    /**
     * @brief Return the adjoint state
     *
     * Accessible only after run() has been called and before ws_ is used
     * elsewhere.
     *
     * @return xB adjoint state
     */
    [[nodiscard]] AmiVector const& getAdjointState() const;

    /**
     * @brief Get the adjoint quadratures (xQB).
     *
     * Accessible only after run() has been called and before ws_ is used
     * elsewhere.

     * @return xQB
     */
    [[nodiscard]] AmiVector const& getAdjointQuadrature() const;

    /**
     * @brief Accessor for has_quadrature_
     * @return has_quadrature_
     */
    [[nodiscard]] bool hasQuadrature() const { return has_quadrature_; }

  private:
    /**
     * @brief Launch backward simulation if Newton solver or linear system solve
     * fail or are disabled.
     *
     * This does not perform any event-handling.
     * For event-handling, see EventHandlingBwdSimulator.
     *
     * @param solver Solver instance.
     */
    void run_simulation(Solver const& solver);

    /**
     * @brief Compute quadratures in adjoint mode
     * @param t0 Initial time for the steady state simulation.
     */
    void compute_steady_state_quadrature(realtype t0);

    /**
     * @brief Compute the quadrature in steady state backward mode by
     * solving the linear system defined by the backward Jacobian.
     */
    void compute_quadrature_by_lin_solve();

    /**
     * @brief Computes the quadrature in steady state backward mode by
     * numerical integration of xB forward in time.
     * @param t0 Initial time for the steady state simulation.
     */
    void compute_quadrature_by_simulation(realtype t0);

    /** CPU time for solving the backward problem (milliseconds) */
    double cpu_timeB_{0.0};

    /** flag indicating whether backward mode was run */
    bool has_quadrature_{false};

    /** The employed number of backward steps */
    int numstepsB_{0};

    /** integral over adjoint state vector */
    AmiVector xQ_;

    /** Final state from pre/post-equilibration forward problem */
    SolutionState& final_state_;

    /** Newton solver */
    NewtonSolver newton_solver_;

    /**
     * Whether the Newton step should be used instead of xdot for convergence
     * checks during simulation.
     */
    bool newton_step_conv_{false};

    Model* model_{nullptr};
    Solver const* solver_{nullptr};
    BwdSimWorkspace* ws_{nullptr};
};

//!  class to solve backwards problems.
/*!
  solves the backwards problem for adjoint sensitivity analysis and handles
  events and data-points
*/

class BackwardProblem {
  public:
    /**
     * @brief Construct backward problem from forward problem
     * @param fwd pointer to corresponding forward problem
     */
    explicit BackwardProblem(ForwardProblem& fwd);

    /**
     * @brief Solve the backward problem.
     *
     * If adjoint sensitivities are enabled, this will also compute
     * sensitivities. workForwardProblem must be called before this function is
     * called.
     */
    void workBackwardProblem();

    /**
     * @brief The adjoint state vector from before pre-equilibration.
     * @return xB
     */
    [[nodiscard]] AmiVector const& getAdjointStatePrePreeq() const {
        return xB_pre_preeq_;
    }

    /**
     * @brief The quadrature state vector from before pre-equilibration.
     * @return xQB
     */
    [[nodiscard]] AmiVector const& getAdjointQuadraturePrePreeq() const {
        return xQB_pre_preeq_;
    }

    /**
     * @brief The final adjoint state vector
     * @return xB
     */
    [[nodiscard]] AmiVector const& getAdjointState() const { return ws_.xB_; }

    /**
     * @brief The final quadrature state vector.
     * @return xQB
     */
    [[nodiscard]] AmiVector const& getAdjointQuadrature() const {
        return ws_.xQB_;
    }

    /**
     * @brief Return the postequilibration SteadyStateBwdProblem.
     * @return The postequilibration SteadyStateBackwardProblem, if any.
     */
    [[nodiscard]] SteadyStateBackwardProblem const*
    getPostequilibrationBwdProblem() const {
        if (posteq_problem_bwd_.has_value())
            return &*posteq_problem_bwd_;
        return nullptr;
    }

    /**
     * @brief Return the preequilibration SteadyStateBackwardProblem.
     * @return The preequilibration SteadyStateBackwardProblem, if any.
     */
    [[nodiscard]] SteadyStateBackwardProblem const*
    getPreequilibrationBwdProblem() const {
        if (preeq_problem_bwd_.has_value())
            return &*preeq_problem_bwd_;
        return nullptr;
    }

  private:
    void handlePostequilibration();

    Model* model_;
    Solver* solver_;
    ExpData const* edata_;

    /** current time */
    realtype t_;

    /** The discontinuities encountered during the main simulation. */
    std::vector<Discontinuity> discs_main_;

    /**
     * state derivative of data likelihood
     * dimensions: nt * nJ * nx_solver
     */
    std::vector<realtype> dJydx_;
    /**
     * state derivative of event likelihood
     * dimensions: nJ * nx_solver * nmaxevent
     */
    std::vector<realtype> const dJzdx_;

    /** The preequilibration steadystate problem from the forward problem. */
    SteadystateProblem* preeq_problem_;

    /** The postequilibration steadystate problem from the forward problem. */
    SteadystateProblem* posteq_problem_;

    /** Presimulation results */
    PeriodResult presim_result;

    BwdSimWorkspace ws_;

    EventHandlingBwdSimulator simulator_;

    /** The preequilibration steady-state backward problem, if any. */
    std::optional<SteadyStateBackwardProblem> preeq_problem_bwd_;

    /** The postequilibration steady-state backward problem, if any. */
    std::optional<SteadyStateBackwardProblem> posteq_problem_bwd_;

    /**
     * The adjoint state vector before pre-equilibration.
     */
    AmiVector xB_pre_preeq_;

    /**
     * The quadrature state vector before pre-equilibration.
     */
    AmiVector xQB_pre_preeq_;
};

} // namespace amici

#endif // AMICI_BACKWARDPROBLEM_H
