#ifndef AMICI_FORWARDPROBLEM_H
#define AMICI_FORWARDPROBLEM_H

#include "amici/defines.h"
#include "amici/edata.h"
#include "amici/misc.h"
#include "amici/model.h"
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
     * to have a root. For i = 0, . . . ,nr 1 or -1 if gi has a root, and = 0
     * if not. See CVodeGetRootInfo for details.
     */
    std::vector<int> root_info;
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
    realtype getTime() const { return t_; }

    /**
     * @brief Accessor for x
     * @return x
     */
    AmiVector const& getState() const { return x_; }

    /**
     * @brief Accessor for dx
     * @return dx
     */
    AmiVector const& getStateDerivative() const { return dx_; }

    /**
     * @brief Accessor for sx
     * @return sx
     */
    AmiVectorArray const& getStateSensitivity() const { return sx_; }

    /**
     * @brief Accessor for nroots
     * @return nroots
     */
    std::vector<int> const& getNumberOfRoots() const { return nroots_; }

    /**
     * @brief Get information on the discontinuities encountered so far.
     * @return The vector of discontinuities.
     */
    std::vector<Discontinuity> const& getDiscontinuities() const {
        return discs_;
    }

    /**
     * @brief Accessor for dJzdx
     * @return dJzdx
     */
    std::vector<realtype> const& getDJzdx() const { return dJzdx_; }

    /**
     * @brief Accessor for it
     * @return it
     */
    int getCurrentTimeIteration() const { return it_; }

    /**
     * @brief Returns final time point for which simulations are available
     * @return time point
     */
    realtype getFinalTime() const { return final_state_.t; }

    /**
     * @brief Returns maximal event index for which simulations are available
     * @return index
     */
    int getEventCounter() const {
        return gsl::narrow<int>(event_states_.size()) - 1;
    }

    /**
     * @brief Returns maximal event index for which the timepoint is available
     * @return index
     */
    int getRootCounter() const { return gsl::narrow<int>(discs_.size()) - 1; }

    /**
     * @brief Retrieves the carbon copy of the simulation state variables at
     * the specified timepoint index
     * @param it timepoint index
     * @return state
     */
    SimulationState const& getSimulationStateTimepoint(int it) const {
        if (model->getTimepoint(it) == initial_state_.t)
            return getInitialSimulationState();
        auto map_iter = timepoint_states_.find(model->getTimepoint(it));
        Ensures(map_iter != timepoint_states_.end());
        return map_iter->second;
    };

    /**
     * @brief Retrieves the carbon copy of the simulation state variables at
     * the specified event index
     * @param iroot event index
     * @return SimulationState
     */
    SimulationState const& getSimulationStateEvent(int iroot) const {
        return event_states_.at(iroot);
    };

    /**
     * @brief Retrieves the carbon copy of the simulation state variables at the
     * initial timepoint
     * @return SimulationState
     */
    SimulationState const& getInitialSimulationState() const {
        return initial_state_;
    };

    /**
     * @brief Retrieves the carbon copy of the simulation state variables at the
     * final timepoint (or when simulation failed)
     * @return SimulationState
     */
    SimulationState const& getFinalSimulationState() const {
        return final_state_;
    };

    /**
     * @brief Return the preequilibration SteadystateProblem.
     * @return The preequilibration SteadystateProblem, if any.
     */
    SteadystateProblem* getPreequilibrationProblem() {
        if (preeq_problem_.has_value())
            return &*preeq_problem_;
        return nullptr;
    }

    /**
     * @brief Return the preequilibration SteadystateProblem.
     * @return The preequilibration SteadystateProblem, if any.
     */
    SteadystateProblem const* getPreequilibrationProblem() const {
        if (preeq_problem_.has_value())
            return &*preeq_problem_;
        return nullptr;
    }

    /**
     * @brief Return the postequilibration SteadystateProblem.
     * @return The postequilibration SteadystateProblem, if any.
     */
    SteadystateProblem* getPostequilibrationProblem() {
        if (posteq_problem_.has_value())
            return &*posteq_problem_;
        return nullptr;
    }

    /**
     * @brief Return the postequilibration SteadystateProblem.
     * @return The postequilibration SteadystateProblem, if any.
     */
    SteadystateProblem const* getPostequilibrationProblem() const {
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
     * @brief Initialize model and solver for presimulation or
     * the main simulation if there is no presimulation.
     */
    void initialize();

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

    /**
     * @brief Execute everything necessary for the handling of events
     *
     * @param tlastroot Reference to the timepoint of the last event
     * @param initial_event initial event flag
     */

    void handleEvent(realtype& tlastroot, bool initial_event);

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
    void storeEvent();

    /**
     * @brief Execute everything necessary for the handling of data points
     *
     * @param t measurement timepoint
     */
    void handleDataPoint(realtype t);

    /**
     * @brief checks whether there are any events to fill
     *
     * @param nmaxevent maximal number of events
     */
    bool checkEventsToFill(int nmaxevent) const {
        return std::any_of(
            nroots_.cbegin(), nroots_.cend(),
            [nmaxevent](int curNRoots) { return curNRoots < nmaxevent; }
        );
    };

    /**
     * @brief fills events at final timepoint if necessary
     *
     * @param nmaxevent maximal number of events
     */
    void fillEvents(int nmaxevent) {
        if (checkEventsToFill(nmaxevent)) {
            discs_.emplace_back(t_);
            storeEvent();
        }
    }

    /**
     * @brief Creates a carbon copy of the current simulation state variables
     * @return state
     */
    SimulationState getSimulationState();

    /** array of number of found roots for a certain event type
     * (dimension: ne) */
    std::vector<int> nroots_;

    /** array of values of the root function (dimension: ne) */
    std::vector<realtype> rootvals_;

    /** temporary rootval storage to check crossing in secondary event
     * (dimension: ne) */
    std::vector<realtype> rval_tmp_;

    /** Discontinuities encountered so far (dimension: dynamic) */
    std::vector<Discontinuity> discs_;

    /** Events that are waiting to be handled at the current timepoint. */
    EventQueue pending_events_;

    /** state derivative of event likelihood
     * (dimension nJ x nx x nMaxEvent, ordering =?) */
    std::vector<realtype> dJzdx_;

    /** current time */
    realtype t_;

    /**
     * @brief Array of flags indicating which root has been found.
     *
     * Array of length nr (ne) with the indices of the user functions gi found
     * to have a root. For i = 0, . . . ,nr 1 or -1 if gi has a root, and = 0
     * if not. See CVodeGetRootInfo for details.
     */
    std::vector<int> roots_found_;

    /** simulation states history at timepoints */
    std::map<realtype, SimulationState> timepoint_states_;

    /** simulation state history at events */
    std::vector<SimulationState> event_states_;

    /** simulation state after initialization*/
    SimulationState initial_state_;

    /** simulation state after simulation */
    SimulationState final_state_;

    /** state vector (dimension: nx_solver) */
    AmiVector x_;

    /** old state vector (dimension: nx_solver) */
    AmiVector x_old_;

    /** differential state vector (dimension: nx_solver) */
    AmiVector dx_;

    /** time derivative state vector (dimension: nx_solver) */
    AmiVector xdot_;

    /** old time derivative state vector (dimension: nx_solver) */
    AmiVector xdot_old_;

    /** sensitivity state vector array (dimension: nx_cl x nplist, row-major) */
    AmiVectorArray sx_;

    /** differential sensitivity state vector array
     * (dimension: nx_cl x nplist, row-major) */
    AmiVectorArray sdx_;

    /** sensitivity of the event timepoint (dimension: nplist) */
    std::vector<realtype> stau_;

    /** storage for last found root */
    realtype tlastroot_{0.0};

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
                // programm termination.
                fwd_->final_state_ = fwd_->getSimulationState();
                // if there is an associated output timepoint, also store it in
                // timepoint_states if it's not present there.
                // this may happen if there is an error just at
                // (or indistinguishably before) an output timepoint
                auto final_time = fwd_->getFinalTime();
                auto const timepoints = fwd_->model->getTimepoints();
                if (!fwd_->timepoint_states_.count(final_time)
                    && std::find(
                           timepoints.cbegin(), timepoints.cend(), final_time
                       ) != timepoints.cend()) {
                    fwd_->timepoint_states_[final_time] = fwd_->final_state_;
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

#endif // FORWARDPROBLEM_H
