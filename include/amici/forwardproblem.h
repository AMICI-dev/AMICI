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
     * @param xdot_pre
     * @param x_post
     * @param xdot_post
     */
    explicit Discontinuity(
        realtype time, std::vector<int> const& root_info = std::vector<int>(),
        AmiVector const& xdot_pre = AmiVector(),
        AmiVector const& x_post = AmiVector(),
        AmiVector const& xdot_post = AmiVector()
    )
        : time(time)
        , x_post(x_post)
        , xdot_post(xdot_post)
        , xdot_pre(xdot_pre)
        , root_info(root_info) {}

    /** Time of discontinuity. */
    realtype time;

    /** Post-event state vector (dimension nx). */
    AmiVector x_post;

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
};

/**
 * @brief Compute the number of roots for each root function from a vector of
 * discontinuities.
 * @param discs Encountered discontinuities.
 * @param ne Number of root functions (ne).
 * @param nmaxevents Maximum number of events to track (nmaxevents).
 * @return The number of roots for each root function.
 */
std::vector<int> compute_nroots(std::vector<Discontinuity> const& discs, int ne, int nmaxevents);

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
        : x(model->nx_solver, solver->getSunContext())
        , x_old(model->nx_solver, solver->getSunContext())
        , dx(model->nx_solver, solver->getSunContext())
        , xdot(model->nx_solver, solver->getSunContext())
        , xdot_old(model->nx_solver, solver->getSunContext())
        , sx(model->nx_solver, model->nplist(), solver->getSunContext())
        , sdx(model->nx_solver, model->nplist(), solver->getSunContext())
        , stau(model->nplist())
        , roots_found(model->ne, 0)
        , rval_tmp(gsl::narrow<decltype(rval_tmp)::size_type>(model->ne), 0.0)
        , nroots(gsl::narrow<decltype(nroots)::size_type>(model->ne), 0)
        , rootvals(gsl::narrow<decltype(rootvals)::size_type>(model->ne), 0.0)

    {};

    /** state vector (dimension: nx_solver) */
    AmiVector x;

    /** old state vector (dimension: nx_solver) */
    AmiVector x_old;

    /** differential state vector (dimension: nx_solver) */
    AmiVector dx;

    /** time derivative state vector (dimension: nx_solver) */
    AmiVector xdot;

    /** old time derivative state vector (dimension: nx_solver) */
    AmiVector xdot_old;

    /** sensitivity state vector array (dimension: nx_cl x nplist, row-major) */
    AmiVectorArray sx;

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
    /** array of number of found roots for a certain event type
     * (dimension: ne) */
    std::vector<int> nroots;

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
        gsl::not_null<Model*> model, gsl::not_null<Solver*> solver,
        gsl::not_null<FwdSimWorkspace*> ws,
        gsl::not_null<std::vector<realtype>*> dJzdx
    )
        : model_(model)
        , solver_(solver)
        , ws_(ws)
        , dJzdx_(dJzdx) {};

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
     */
    void
    run(realtype t, ExpData const* edata,
        std::vector<amici::realtype> const& timepoints);

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

    /** The current time. */
    realtype t_;

    /** Time index in the current list of timepoints */
    int it_;

  private:
    /**
     * @brief Execute everything necessary for the handling of events
     *
     * @param initial_event initial event flag
     * @param edata experimental data
     */
    void handle_event(bool initial_event, ExpData const* edata);

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
    void store_event(amici::ExpData const* edata);

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
        if (!std::ranges::any_of(ws_->nroots, [nmaxevent](int curNRoots) {
                return curNRoots < nmaxevent;
            }))
            return;

        result.discs.emplace_back(t_);
        store_event(edata);
    }

    /** Initial time of the current period. */
    realtype t0_;

    /** The model to simulate. */
    Model* model_;

    /** The solver to use for the simulation. */
    Solver* solver_;

    /** The workspace to use for the simulation. */
    gsl::not_null<FwdSimWorkspace*> ws_;

    /** state derivative of event likelihood
     * (dimension nJ x nx x nMaxEvent, ordering =?) */
    std::vector<realtype>* dJzdx_ = nullptr;
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
     * @brief Accessor for t
     * @return t
     */
    [[nodiscard]] realtype getTime() const { return t_; }

    /**
     * @brief Accessor for sx
     * @return sx
     */
    [[nodiscard]] AmiVectorArray const& getStateSensitivity() const {
        return ws_.sx;
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
        return main_simulator_.result.final_state_.t;
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
    getSimulationStateTimepoint(int it) const {
        if (model->getTimepoint(it) == main_simulator_.result.initial_state_.t)
            return getInitialSimulationState();
        auto map_iter = main_simulator_.result.timepoint_states_.find(
            model->getTimepoint(it)
        );
        Ensures(map_iter != main_simulator_.result.timepoint_states_.end());
        return map_iter->second;
    };

    /**
     * @brief Retrieves the carbon copy of the simulation state variables at
     * the specified event index
     * @param iroot event index
     * @return SimulationState
     */
    [[nodiscard]] SimulationState const&
    getSimulationStateEvent(int iroot) const {
        return main_simulator_.result.event_states_.at(iroot);
    };

    /**
     * @brief Retrieves the carbon copy of the simulation state variables at the
     * initial timepoint
     * @return SimulationState
     */
    [[nodiscard]] SimulationState const& getInitialSimulationState() const {
        return main_simulator_.result.initial_state_;
    };

    /**
     * @brief Retrieves the carbon copy of the simulation state variables at the
     * final timepoint (or when the simulation failed)
     * @return SimulationState
     */
    [[nodiscard]] SimulationState const& getFinalSimulationState() const {
        return main_simulator_.result.final_state_;
    };

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

    /** current time */
    realtype t_;

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
                auto final_time = fwd_->getFinalTime();
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
