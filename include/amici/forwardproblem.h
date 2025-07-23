#ifndef AMICI_FORWARDPROBLEM_H
#define AMICI_FORWARDPROBLEM_H

#include "amici/defines.h"
#include "amici/edata.h"
#include "amici/misc.h"
#include "amici/model.h"
#include "amici/solver.h"
#include "amici/steadystateproblem.h"
#include "amici/vector.h"

#include <optional>
#include <vector>
namespace amici {
class ExpData;
class Solver;
class SteadystateProblem;
class FinalStateStorer;

/**
 * @brief Data structure to store some state of a simulation at a discontinuity.
 */
struct Discontinuity {
    /**
     * @brief Constructor.
     * @param time
     * @param root_info
     * @param x_pre
     * @param xdot_pre
     * @param x_post
     * @param xdot_post
     * @param h_pre
     * @param total_cl_pre
     */
    explicit Discontinuity(
        realtype const time,
        std::vector<int> const& root_info = std::vector<int>(),
        AmiVector const& x_pre = AmiVector(),
        AmiVector const& xdot_pre = AmiVector(),
        AmiVector const& x_post = AmiVector(),
        AmiVector const& xdot_post = AmiVector(),
        std::vector<realtype> const& h_pre = std::vector<realtype>(),
        std::vector<realtype> const& total_cl_pre = std::vector<realtype>(0)
    )
        : time(time)
        , x_post(x_post)
        , x_pre(x_pre)
        , xdot_post(xdot_post)
        , xdot_pre(xdot_pre)
        , root_info(root_info)
        , h_pre(h_pre)
        , total_cl_pre(total_cl_pre) {}

    /** Time of discontinuity. */
    realtype time;

    /** Post-event state vector (dimension nx). */
    AmiVector x_post;

    /** Pre-event state vector (dimension nx). */
    AmiVector x_pre;

    /** Post-event differential state vectors (dimension nx). */
    AmiVector xdot_post;

    /** Pre-event differential state vectors (dimension nx). */
    AmiVector xdot_pre;

    /**
     * @brief Array of flags indicating which root has been found.
     *
     * Array of length nr (ne) with the indices of the user functions gi found
     * to have a root. For i = 0, ..., nr: 1 or -1 if gi has a root, and = 0
     * if not. See CVodeGetRootInfo for details.
     */
    std::vector<int> root_info;

    /**
     * Flag indicating whether a certain Heaviside function should be active or
     * not, pre-event value (dimension: `ne`)
     */
    std::vector<realtype> h_pre;

    /** time derivative of state (DAE only) post-event */
    AmiVector dx_post;

    /** Total abundances for conservation laws
     (dimension: `nx_rdata - nx_solver`) */
    std::vector<realtype> total_cl_pre;
};

/**
 * @brief Compute the number of roots for each root function from a vector of
 * discontinuities.
 * @param discs Encountered discontinuities.
 * @param ne Number of root functions (ne).
 * @param nmaxevents Maximum number of events to track (nmaxevents).
 * @return The number of roots for each root function.
 */
std::vector<int>
compute_nroots(std::vector<Discontinuity> const& discs, int ne, int nmaxevents);

/**
 * @brief The FwdSimWorkspace class is used to store temporary simulation
 * state during forward simulations.
 */
struct FwdSimWorkspace {
    /**
     * @brief Constructor
     * @param model The model for which to set up the workspace.
     * @param solver The solver for which to set up this workspace.
     */
    FwdSimWorkspace(
        gsl::not_null<Model*> const& model, gsl::not_null<Solver*> solver
    )
        : sol(NAN, model->nx_solver, model->nplist(), solver->getSunContext())
        , x_old(model->nx_solver, solver->getSunContext())
        , xdot(model->nx_solver, solver->getSunContext())
        , xdot_old(model->nx_solver, solver->getSunContext())
        , sdx(model->nx_solver, model->nplist(), solver->getSunContext())
        , stau(model->nplist())
        , roots_found(model->ne, 0)
        , rval_tmp(gsl::narrow<decltype(rval_tmp)::size_type>(model->ne), 0.0)
        , nroots(gsl::narrow<decltype(nroots)::size_type>(model->ne), 0)
        , rootvals(gsl::narrow<decltype(rootvals)::size_type>(model->ne), 0.0)

    {}
    /** Current solution state */
    SolutionState sol;

    /** old state vector (dimension: nx_solver) */
    AmiVector x_old;

    /** time derivative state vector (dimension: nx_solver) */
    AmiVector xdot;

    /** old time derivative state vector (dimension: nx_solver) */
    AmiVector xdot_old;

    /** differential sensitivity state vector array
     * (dimension: nx_cl x nplist, row-major) */
    AmiVectorArray sdx;

    /** sensitivity of the event timepoint (dimension: nplist) */
    std::vector<realtype> stau;

    /**
     * @brief Array of flags indicating which root has been found.
     *
     * Array of length nr (ne) with the indices of the user functions gi found
     * to have a root. For i = 0, . . . ,nr 1 or -1 if gi has a root, and = 0
     * if not. See CVodeGetRootInfo for details.
     */
    std::vector<int> roots_found;

    /** Timepoint, at which the last root was found */
    realtype tlastroot{0.0};

    /** Events that are waiting to be handled at the current timepoint. */
    EventQueue pending_events;

    /** temporary rootval storage to check crossing in secondary event
     * (dimension: ne) */
    std::vector<realtype> rval_tmp;

    /** array of number of found roots for a certain event type
     * (dimension: ne) */
    std::vector<int> nroots;

    /** array of values of the root function (dimension: ne) */
    std::vector<realtype> rootvals;
};

/**
 * @brief The PeriodResult class stores the result of a simulation period.
 */
struct PeriodResult {

    /** Discontinuities encountered so far (dimension: dynamic) */
    std::vector<Discontinuity> discs;

    /** simulation states history at timepoints */
    std::map<realtype, SimulationState> timepoint_states_;

    /** simulation state history at events */
    std::vector<SimulationState> event_states_;

    /** simulation state after initialization*/
    SimulationState initial_state_;

    /** simulation state after simulation */
    SimulationState final_state_;
};

/**
 * @brief The EventHandlingSimulator class runs a forward simulation
 * and processes events, measurements or output timepoints in general.
 */
class EventHandlingSimulator {
  public:
    /**
     * @brief EventHandlingSimulator
     * @param model The model to simulate.
     * @param solver The solver to use for the simulation.
     * @param ws The workspace to use for the simulation.
     * @param dJzdx State-derivative of event likelihood
     * (dimension nJ x nx x nMaxEvent, ordering =?)
     */
    EventHandlingSimulator(
        gsl::not_null<Model*> model, gsl::not_null<Solver const*> solver,
        gsl::not_null<FwdSimWorkspace*> ws, std::vector<realtype>* dJzdx
    )
        : model_(model)
        , solver_(solver)
        , ws_(ws)
        , dJzdx_(dJzdx) {}

    /**
     * @brief run Run the simulation.
     *
     * It will run the simulation from the initial time of this period
     * to the final timepoint of this period, handling events
     * and data points as they occur.
     *
     * Expects the model and the solver to be set up, and `ws` to be initialized
     * for this period.
     *
     * @param t The initial time of this period.
     * @param edata Any experimental data associated with this simulation.
     * `nullptr` if no experimental data is associated with this
     * period.
     * @param timepoints The output timepoints or measurement timepoints of
     * this period. This must contain at least the final timepoint of this
     * period.
     * @param store_diagnosis Whether to store diagnosis in the solver at each
     * timepoint in `timepoints`.
     */
    void
    run(realtype t, ExpData const* edata,
        std::vector<realtype> const& timepoints, bool store_diagnosis);

    /**
     * @brief Simulate until steady state and handle events along the way.
     *
     * This will run the simulation until the steady state is reached.
     * Events will be handled as they occur.
     * No data points will be handled during this simulation
     * (not for event-resolved nor for time-resolved observables).
     *
     * @param check_convergence Function to check convergence based on the
     * current values in `ws_`.
     * Input: whether `xdot` has to be recomputed. Output: whether convergence
     * has been reached.
     * @param convergence_check_frequency Frequency of convergence checks.
     * @param sim_steps Number of simulation steps taken until convergence or
     * error.
     */
    void run_steady_state(
        std::function<bool(bool)> check_convergence,
        int convergence_check_frequency, int& sim_steps
    );

    /**
     * @brief Returns maximal event index for which the timepoint is available
     * @return index
     */
    [[nodiscard]] int get_root_counter() const {
        return gsl::narrow<int>(result.discs.size()) - 1;
    }

    /**
     * @brief Returns maximal event index for which simulations are available
     * @return index
     */
    [[nodiscard]] int get_event_counter() const {
        return gsl::narrow<int>(result.event_states_.size()) - 1;
    }

    /**
     * @brief Creates a carbon copy of the current simulation state variables
     * @return state
     */
    SimulationState get_simulation_state();

    /** Results for the current simulation period. */
    PeriodResult result;

    /** Time index in the current list of timepoints */
    int it_ = -1;

    /**
     * @brief Execute everything necessary for the handling of events.
     *
     * Assumes that at least one event triggered.
     *
     * @param initial_event initial event flag
     * @param edata experimental data
     */
    void handle_events(bool initial_event, ExpData const* edata);

    /**
     * @brief Handle initial events.
     *
     * Check if there are any initial events that need to be handled
     * and execute the necessary steps.
     *
     * @param edata experimental data
     */
    void handle_initial_events(ExpData const* edata);

  private:
    /**
     * @brief Store pre-event model state
     *
     * @param seflag Secondary event flag
     * @param initial_event initial event flag
     */
    void store_pre_event_state(bool seflag, bool initial_event);

    /**
     * @brief Check for, and if applicable, handle any secondary events
     * @return the number of secondary events found
     */
    int detect_secondary_events();

    /**
     * @brief Extract output information for events
     */
    void store_event(ExpData const* edata);

    /**
     * @brief Execute everything necessary for the handling of data points
     *
     * @param t measurement timepoint
     */
    void handle_datapoint(realtype t);

    /**
     * @brief fills events at final timepoint if necessary
     *
     * @param nmaxevent maximal number of events
     * @param edata experimental data
     */
    void fill_events(int nmaxevent, ExpData const* edata) {
        if (!std::ranges::any_of(ws_->nroots, [nmaxevent](int const curNRoots) {
                return curNRoots < nmaxevent;
            }))
            return;

        result.discs.emplace_back(ws_->sol.t);
        store_event(edata);
    }

    /** Initial time of the current period. */
    realtype t0_{NAN};

    /** The model to simulate. */
    Model* model_{nullptr};

    /** The solver to use for the simulation. */
    Solver const* solver_{nullptr};

    /** The workspace to use for the simulation. */
    gsl::not_null<FwdSimWorkspace*> ws_;

    /** state derivative of event likelihood
     * (dimension nJ x nx x nMaxEvent, ordering =?) */
    std::vector<realtype>* dJzdx_ = nullptr;
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
     * @param ws Workspace for forward simulation
     * @param solver Solver instance
     * @param model Model instance
     * @param is_preeq Whether this is a pre-equilibration (`true`) or
     * post-equilibration problem (`false`).
     */
    explicit SteadystateProblem(
        FwdSimWorkspace* ws, Solver const& solver, Model& model, bool is_preeq
    );

    /**
     * @brief Compute the steady state in the forward case.
     *
     * Tries to determine the steady state of the ODE system and computes
     * steady state sensitivities if requested.
     * Expects that solver, model, and ws_ are already initialized.
     *
     * @param solver The solver instance
     * @param it Index of the current output time point.
     * @param t0 Initial time for the steady state simulation.
     */
    void workSteadyStateProblem(Solver& solver, int it, realtype t0);

    /**
     * @brief Return the stored SimulationState.
     * @return stored SimulationState
     */
    [[nodiscard]] SimulationState const& getFinalSimulationState() const {
        return period_result_.final_state_;
    }

    /**
     * @brief Return state at steady state
     * @return x
     */
    [[nodiscard]] AmiVector const& getState() const {
        return period_result_.final_state_.sol.x;
    }

    /**
     * @brief Return state sensitivity at steady state
     * @return sx
     */
    [[nodiscard]] AmiVectorArray const& getStateSensitivity() const {
        return period_result_.final_state_.sol.sx;
    }

    /**
     * @brief Get the CPU time taken to solve the forward problem.
     * @return The CPU time in milliseconds.
     */
    [[nodiscard]] double getCPUTime() const { return cpu_time_; }

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
    [[nodiscard]] realtype getSteadyStateTime() const {
        return period_result_.final_state_.sol.t;
    }

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
     * @brief Check, whether any approach to find the steady state was
     * successful.
     * @return Whether any approach to find the steady state was successful.
     */
    [[nodiscard]] bool checkSteadyStateSuccess() const;

    /**
     * @brief Get the pre-equilibration solver.
     * @return The preequilibration solver.
     */
    [[nodiscard]] Solver const* get_solver() const { return solver_; }

    /**
     * @brief Get the preequilibration result.
     * @return
     */
    [[nodiscard]] PeriodResult const& get_result() const {
        return period_result_;
    }

  private:
    /**
     * @brief Handle the computation of the steady state.
     *
     * Throws an AmiException if no steady state was found.
     *
     * @param it Index of the current output time point.
     * @param t0 Initial time for the steady state simulation.
     */
    void findSteadyState(int it, realtype t0);

    /**
     * @brief Try to determine the steady state by using Newton's method.
     * @param newton_retry Flag indicating whether Newton's method is being
     * relaunched.
     */
    void findSteadyStateByNewtonsMethod(bool newton_retry);

    /**
     * @brief Try to determine the steady state by using forward simulation.
     * @param it Index of the current output time point.
     * @param t0 Initial time for the steady state simulation.
     * @return SteadyStateStatus indicating whether the steady state was found
     * successfully, or if it failed.
     */
    SteadyStateStatus findSteadyStateBySimulation(int it, realtype t0);

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
     * @brief Checks convergence for state sensitivities
     * @param wrms_computer_sx WRMSComputer instance for state sensitivities
     * @return weighted root mean squared residuals of the RHS
     */
    realtype getWrmsFSA(WRMSComputer& wrms_computer_sx);

    /**
     * @brief Launch simulation if Newton solver or linear system solve
     * fail or are disabled.
     * simulation.
     */
    void runSteadystateSimulationFwd();

    /**
     * @brief Update member variables to indicate that state_.x has been
     * updated and xdot_, delta_, etc. need to be recomputed.
     */
    void flagUpdatedState();

    /**
     * @brief Retrieve simulation sensitivities from the provided solver and
     * set the corresponding flag to indicate they are up to date
     */
    void updateSensiSimulation();

    /**
     * @brief Compute the right-hand side for the current state_.x and set the
     * corresponding flag to indicate xdot_ is up to date.
     */
    void updateRightHandSide();

    /** Whether this is a pre- or post-equilibration problem */
    bool is_preeq_;

    /** Workspace for forward simulation */
    FwdSimWorkspace* ws_;

    /** WRMS computer for x */
    WRMSComputer wrms_computer_x_;

    /** weighted root-mean-square error */
    realtype wrms_{NAN};

    /** Results for this period. */
    PeriodResult period_result_;

    /** stores diagnostic information about employed number of steps */
    std::vector<int> numsteps_{std::vector<int>(3, 0)};

    /** CPU time for solving the forward problem (milliseconds) */
    double cpu_time_{0.0};

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

    /** flag indicating whether xdot_ has been computed for the current state */
    bool xdot_updated_{false};
    /**
     * flag indicating whether simulation sensitivities have been retrieved for
     * the current state
     */
    bool sensis_updated_{false};

    /**
     * A dedicated solver for pre-equilibration is required if we need to
     * perform a backward simulation for adjoint sensitivities for models with
     * events.
     *
     * This unique_ptr is used to manage the lifetime of that solver and may
     * be null. Always use the `solver_` pointer to access the solver,
     * which will point to a dedicated pre-equilibration solver if it exists,
     * or to the main solver otherwise.
     */
    std::unique_ptr<Solver> preeq_solver_unique_ptr_;

    /** Pointer to the pre-equilibration solver */
    Solver const* solver_{nullptr};

    /** The model to equilibrate */
    Model* model_{nullptr};

    /** Simulator for event handling */
    std::optional<EventHandlingSimulator> simulator_;
};

/**
 * @brief The ForwardProblem class groups all functions for solving the
 * forward problem.
 */
class ForwardProblem {
  public:
    /**
     * @brief Constructor
     * @param edata pointer to ExpData instance
     * @param model pointer to Model instance
     * @param solver pointer to Solver instance
     * problem, pass nullptr for no initialization
     */
    ForwardProblem(
        ExpData const* edata, gsl::not_null<Model*> model,
        gsl::not_null<Solver*> solver
    );

    ~ForwardProblem() = default;

    /** allow FinalStateStorer to access private members and functions */
    friend ::amici::FinalStateStorer;

    /**
     * @brief Solve the forward problem.
     *
     * If forward sensitivities are enabled this will also compute
     * sensitivities.
     */
    void workForwardProblem();

    /**
     * @brief Computes adjoint updates dJydx according to the provided model
     * and ExpData.
     * @param model Model instance
     * @param edata experimental data
     * @return dJydx
     */
    std::vector<realtype> getAdjointUpdates(Model& model, ExpData const& edata);

    /**
     * @brief Accessor for sx
     * @return sx
     */
    [[nodiscard]] AmiVectorArray const& getStateSensitivity() const {
        return ws_.sol.sx;
    }

    /**
     * @brief Get information on the discontinuities encountered so far.
     * @return The vector of discontinuities.
     */
    [[nodiscard]] std::vector<Discontinuity> const& getDiscontinuities() const {
        return main_simulator_.result.discs;
    }

    /**
     * @brief Accessor for dJzdx
     * @return dJzdx
     */
    [[nodiscard]] std::vector<realtype> const& getDJzdx() const {
        return dJzdx_;
    }

    /**
     * @brief Accessor for it
     * @return it
     */
    [[nodiscard]] int getCurrentTimeIteration() const { return it_; }

    /**
     * @brief Returns final time point for which simulations are available
     * @return time point
     */
    [[nodiscard]] realtype getFinalTime() const {
        return main_simulator_.result.final_state_.sol.t;
    }

    /**
     * @brief Returns maximal event index for which simulations are available
     * @return index
     */
    [[nodiscard]] int getEventCounter() const {
        return main_simulator_.get_event_counter();
    }

    /**
     * @brief Retrieves the carbon copy of the simulation state variables at
     * the specified timepoint index
     * @param it timepoint index
     * @return state
     */
    [[nodiscard]] SimulationState const&
    getSimulationStateTimepoint(int const it) const {
        if (model->getTimepoint(it)
            == main_simulator_.result.initial_state_.sol.t)
            return getInitialSimulationState();
        auto const map_iter = main_simulator_.result.timepoint_states_.find(
            model->getTimepoint(it)
        );
        Ensures(map_iter != main_simulator_.result.timepoint_states_.end());
        return map_iter->second;
    }

    /**
     * @brief Retrieves the carbon copy of the simulation state variables at
     * the specified event index
     * @param iroot event index
     * @return SimulationState
     */
    [[nodiscard]] SimulationState const&
    getSimulationStateEvent(int const iroot) const {
        return main_simulator_.result.event_states_.at(iroot);
    }

    /**
     * @brief Retrieves the carbon copy of the simulation state variables at the
     * initial timepoint
     * @return SimulationState
     */
    [[nodiscard]] SimulationState const& getInitialSimulationState() const {
        return main_simulator_.result.initial_state_;
    }

    /**
     * @brief Retrieves the carbon copy of the simulation state variables at the
     * final timepoint (or when the simulation failed)
     * @return SimulationState
     */
    [[nodiscard]] SimulationState const& getFinalSimulationState() const {
        return main_simulator_.result.final_state_;
    }

    /**
     * @brief Return the preequilibration SteadystateProblem.
     * @return The preequilibration SteadystateProblem, if any.
     */
    [[nodiscard]] SteadystateProblem* getPreequilibrationProblem() {
        if (preeq_problem_.has_value())
            return &*preeq_problem_;
        return nullptr;
    }

    /**
     * @brief Return the preequilibration SteadystateProblem.
     * @return The preequilibration SteadystateProblem, if any.
     */
    [[nodiscard]] SteadystateProblem const* getPreequilibrationProblem() const {
        if (preeq_problem_.has_value())
            return &*preeq_problem_;
        return nullptr;
    }

    /**
     * @brief Return the postequilibration SteadystateProblem.
     * @return The postequilibration SteadystateProblem, if any.
     */
    [[nodiscard]] SteadystateProblem* getPostequilibrationProblem() {
        if (posteq_problem_.has_value())
            return &*posteq_problem_;
        return nullptr;
    }

    /**
     * @brief Return the postequilibration SteadystateProblem.
     * @return The postequilibration SteadystateProblem, if any.
     */
    [[nodiscard]] SteadystateProblem const*
    getPostequilibrationProblem() const {
        if (posteq_problem_.has_value())
            return &*posteq_problem_;
        return nullptr;
    }

    /**
     * @brief Get the presimulation results.
     * @return Presimulation results.
     */
    PeriodResult const& get_presimulation_result() const {
        return pre_simulator_.result;
    }

    /**
     * @brief Whether pre-equilibration was performed successfully.
     * @return
     */
    bool was_preequilibrated() const { return preequilibrated_; }

    /** pointer to model instance */
    Model* model;

    /** pointer to solver instance */
    Solver* solver;

    /** pointer to experimental data instance */
    ExpData const* edata;

  private:
    /**
     * @brief Handle preequilibration if necessary.
     *
     * Preequilibration starts at `Model::t0()`.
     *
     * So far, no event handling takes place during preequilibration.
     */
    void handlePreequilibration();

    /**
     * @brief Handle pre-simulation if required.
     *
     * Pre-simulation starts at `Model::t0() - ExpData::t_presim`.
     *
     * So far, no event handling takes place during presimulation.
     */
    void handlePresimulation();

    /**
     * @brief Handle the main simulation.
     *
     * Simulation starts at `Model::t0()`.
     * During this period, events are processed and data points are
     * handled.
     */
    void handleMainSimulation();

    /**
     * @brief Handle postequilibration if necessary.
     *
     * Postequilibration starts after the last finite output timepoint
     * or the last event timepoint that is known a priori, whichever is later.
     *
     * So far, no event handling takes place during postequilibration.
     * This also includes the processing of event observables.
     */
    void handlePostequilibration();

    /** state derivative of event likelihood
     * (dimension nJ x nx x nMaxEvent, ordering =?) */
    std::vector<realtype> dJzdx_;

    /** flag to indicate whether solver was preeinitialized via preequilibration
     */
    bool preequilibrated_{false};

    /** current iteration number for time index */
    int it_ = 0;

    /** Whether the current model/data requires presimulation. */
    bool uses_presimulation_{false};

    /** The preequilibration steady-state problem, if any. */
    std::optional<SteadystateProblem> preeq_problem_;

    /** The postequilibration steady-state problem, if any. */
    std::optional<SteadystateProblem> posteq_problem_;

    FwdSimWorkspace ws_;
    EventHandlingSimulator main_simulator_;
    EventHandlingSimulator pre_simulator_;
};

/**
 * @brief stores the stimulation state when it goes out of scope
 */
class FinalStateStorer : public ContextManager {
  public:
    /**
     * @brief constructor, attaches problem pointer
     * @param fwd problem from which the simulation state is to be stored
     */
    explicit FinalStateStorer(ForwardProblem* fwd)
        : fwd_(fwd) {}

    FinalStateStorer& operator=(FinalStateStorer const& other) = delete;

    /**
     * @brief destructor, stores simulation state
     */
    ~FinalStateStorer() noexcept(false) {
        if (fwd_) {
            try {
                // This may throw in `CVodeSolver::getSens`
                // due to https://github.com/LLNL/sundials/issues/82.
                // Therefore, this dtor must be `noexcept(false)` to avoid
                // program termination.
                fwd_->main_simulator_.result.final_state_
                    = fwd_->main_simulator_.get_simulation_state();
                // if there is an associated output timepoint, also store it in
                // timepoint_states if it's not present there.
                // this may happen if there is an error just at
                // (or indistinguishably before) an output timepoint
                auto const final_time = fwd_->getFinalTime();
                auto const timepoints = fwd_->model->getTimepoints();
                if (!fwd_->main_simulator_.result.timepoint_states_.contains(
                        final_time
                    )
                    && std::ranges::find(timepoints, final_time)
                           != timepoints.cend()) {
                    fwd_->main_simulator_.result.timepoint_states_[final_time]
                        = fwd_->main_simulator_.result.final_state_;
                }
            } catch (std::exception const&) {
                // We must not throw in case we are already in the stack
                // unwinding phase due to some other active exception, otherwise
                // this will also lead to termination.
                //
                // In case there is another active exception,
                // `fwd_->{final_state_,timepoint_states_}` won't be set,
                // and we assume that they are either not accessed anymore, or
                // that there is appropriate error handling in place.
                if (!std::uncaught_exceptions()) {
                    throw;
                }
            }
        }
    }

  private:
    ForwardProblem* fwd_;
};

} // namespace amici

#endif // AMICI_FORWARDPROBLEM_H
