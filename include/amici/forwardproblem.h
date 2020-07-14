#ifndef AMICI_FORWARDPROBLEM_H
#define AMICI_FORWARDPROBLEM_H

#include "amici/defines.h"
#include "amici/vector.h"
#include "amici/model.h"
#include "amici/misc.h"
#include "amici/sundials_matrix_wrapper.h"

#include <sundials/sundials_direct.h>
#include <vector>
#include <memory>

namespace amici {

class ExpData;
class Solver;
class SteadystateProblem;
class FinalStateStorer;

/**
 * @brief implements an exchange format to store and transfer the state of a simulation at a
 * specific timepoint.
 */
struct SimulationState{
    /** timepoint */
    realtype t;
    /** state variables */
    AmiVector x;
    /** state variables */
    AmiVector dx;
    /** state variable sensitivity */
    AmiVectorArray sx;
    /** state of the model that was used for simulation */
    ModelState state;
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
     * @param preeq preequilibration with which to initialize the forward problem,
     * pass nullptr for no initialization
     */
    ForwardProblem(const ExpData *edata, Model *model, Solver *solver,
                   const SteadystateProblem *preeq);

    ~ForwardProblem() = default;

    /** allow FinalStateStorer to access private members and functions */
    friend FinalStateStorer;

    /**
     * @brief Solve the forward problem.
     *
     * If forward sensitivities are enabled this will also compute sensitivies.
     */
    void workForwardProblem();

    /**
     * @brief computes adjoint updates dJydx according to provided model and expdata
     * @param model Model instance
     * @param edata experimental data
     */
    void getAdjointUpdates(Model &model, const ExpData &edata);

    /**
     * @brief Accessor for t
     * @return t
     */
    realtype getTime() const {
        return t;
    }

    /**
     * @brief Accessor for x
     * @return x
     */
    AmiVector const& getState() const {
        return x;
    }

    /**
     * @brief Accessor for dx
     * @return dx
     */
    AmiVector const& getStateDerivative() const {
        return dx;
    }

    /**
     * @brief Accessor for sx
     * @return sx
     */
    AmiVectorArray const& getStateSensitivity() const {
        return sx;
    }

    /**
     * @brief Accessor for x_disc
     * @return x_disc
     */
    std::vector<AmiVector> const& getStatesAtDiscontinuities() const {
        return x_disc;
    }

    /**
     * @brief Accessor for xdot_disc
     * @return xdot_disc
     */
    std::vector<AmiVector> const& getRHSAtDiscontinuities() const {
        return xdot_disc;
    }

    /**
     * @brief Accessor for xdot_old_disc
     * @return xdot_old_disc
     */
    std::vector<AmiVector> const& getRHSBeforeDiscontinuities() const {
        return xdot_old_disc;
    }

    /**
     * @brief Accessor for nroots
     * @return nroots
     */
    std::vector<int> const& getNumberOfRoots() const {
        return nroots;
    }

    /**
     * @brief Accessor for discs
     * @return discs
     */
    std::vector<realtype> const& getDiscontinuities() const {
        return discs;
    }

    /**
     * @brief Accessor for rootidx
     * @return rootidx
     */
    std::vector<std::vector<int>> const& getRootIndexes() const {
        return rootidx;
    }

    /**
     * @brief Accessor for dJydx
     * @return dJydx
     */
   std::vector<realtype> const& getDJydx() const {
        return dJydx;
    }

    /**
     * @brief Accessor for dJzdx
     * @return dJzdx
     */
    std::vector<realtype> const& getDJzdx() const {
        return dJzdx;
    }

    /**
     * @brief Accessor for pointer to x
     * @return &x
     */
    AmiVector *getStatePointer() {
        return &x;
    }

    /**
     * @brief Accessor for pointer to dx
     * @return &dx
     */
    AmiVector *getStateDerivativePointer() {
        return &dx;
    }

    /**
     * @brief accessor for pointer to sx
     * @return &sx
     */
    AmiVectorArray *getStateSensitivityPointer() {
        return &sx;
    }

    /**
     * @brief Accessor for pointer to sdx
     * @return &sdx
     */
    AmiVectorArray *getStateDerivativeSensitivityPointer() {
        return &sdx;
    }

    /**
     * @brief Accessor for it
     * @return it
     */
    int getCurrentTimeIteration() const {
        return it;
    }

    /**
     * @brief Returns maximal time point index for which simulations are available
     * @return index
     */
    int getTimepointCounter() const {
        return static_cast<int>(timepoints.size() - 1);
    }

    /**
     * @brief Returns maximal event index for which simulations are available
     * @return index
     */
    int getEventCounter() const {
        return static_cast<int>(event_states.size() - 1);
    }

    /**
     * @brief Returns maximal event index for which the timepoint is available
     * @return index
     */
    int getRootCounter() const {
        return static_cast<int>(discs.size() - 1);
    }

    /**
     * @brief Retrieves the carbon copy of the simulation state variables at
     * the specified timepoint index
     * @param it timpoint index
     * @return state
     */
    const SimulationState &getSimulationStateTimepoint(int it) const {
        return timepoint_states.find(timepoints.at(it))->second;
    };

    /**
     * @brief Retrieves the carbon copy of the simulation state variables at
     * the specified event index
     * @param iroot event index
     * @return SimulationState
     */
    const SimulationState &getSimulationStateEvent(int iroot) const {
        return event_states.at(iroot);
    };

    /**
     * @brief Retrieves the carbon copy of the simulation state variables at the
     * initial timepoint
     * @return SimulationState
     */
    const SimulationState &getInitialSimulationState() const {
        return initial_state;
    };

    /**
     * @brief Retrieves the carbon copy of the simulation state variables at the
     * final timepoint (or when simulation failed)
     * @return SimulationState
     */
    const SimulationState &getFinalSimulationState() const {
        return final_state;
    };

    /** pointer to model instance */
    Model *model;

    /** pointer to solver instance */
    Solver *solver;

    /** pointer to experimental data instance */
    const ExpData *edata;
  private:

    void handlePresimulation();

    /**
     * @brief Execute everything necessary for the handling of events
     *
     * @param tlastroot pointer to the timepoint of the last event
     */

    void handleEvent(realtype *tlastroot,bool seflag);

    /**
     * @brief Extract output information for events
     */
    void storeEvent();

    /**
     * @brief Execute everything necessary for the handling of data points
     *
     * @param it index of data point
     */
    void handleDataPoint(int it);

    /**
     * @brief Applies the event bolus to the current state
     *
     * @param model pointer to model specification object
     */
    void applyEventBolus();

    /**
     * @brief Applies the event bolus to the current sensitivities
     */
    void applyEventSensiBolusFSA();

    /**
     * @brief checks whether there are any events to fill
     *
     * @param nmaxevent maximal number of events
     */
    bool checkEventsToFill(int nmaxevent) const {
        return std::any_of(nroots.cbegin(), nroots.cend(),
                           [nmaxevent](int curNRoots) {
                return curNRoots < nmaxevent;
        });
    };

    /**
     * @brief fills events at final timepoint if necessary
     *
     * @param nmaxevent maximal number of events
     */
    void fillEvents(int nmaxevent) {
        if (checkEventsToFill(nmaxevent)) {
            discs.push_back(t);
            storeEvent();
        }
    }

    /**
     * @brief Creates a carbon copy of the current simulation state variables
     * @return state
     */
    SimulationState getSimulationState() const;

    /** array of index vectors (dimension ne) indicating whether the respective
     * root has been detected for all so far encountered discontinuities,
     * extended as needed (dimension: dynamic) */
    std::vector<std::vector<int>> rootidx;

    /** array of number of found roots for a certain event type
     * (dimension: ne) */
    std::vector<int> nroots;

    /** array of values of the root function (dimension: ne) */
    std::vector<realtype> rootvals;

    /** temporary rootval storage to check crossing in secondary event
     * (dimension: ne) */
    std::vector<realtype> rvaltmp;

    /** array containing the time-points of discontinuities
     * (dimension: nmaxevent x ne, ordering = ?) */
    std::vector<realtype> discs;

    /** array containing the index of discontinuities
     * (dimension: nmaxevent x ne, ordering = ?) */
    std::vector<realtype> irdiscs;

    /** array of state vectors (dimension nx) for all so far encountered
     * discontinuities, extended as needed (dimension dynamic) */
    std::vector<AmiVector> x_disc;

    /** array of differential state vectors (dimension nx) for all so far
     * encountered discontinuities, extended as needed (dimension dynamic) */
    std::vector<AmiVector> xdot_disc;

    /** array of old differential state vectors (dimension nx) for all so far
     * encountered discontinuities, extended as needed (dimension dynamic) */
    std::vector<AmiVector> xdot_old_disc;

    /** state derivative of data likelihood
     * (dimension nJ x nx x nt, ordering =?) */
    std::vector<realtype> dJydx;

    /** state derivative of event likelihood
     * (dimension nJ x nx x nMaxEvent, ordering =?) */
    std::vector<realtype> dJzdx;

    /** current time */
    realtype t;

    /** number of timepoints */
    std::vector<realtype> timepoints;

    /**
     * @brief Array of flags indicating which root has beend found.
     *
     * Array of length nr (ne) with the indices of the user functions gi found
     * to have a root. For i = 0, . . . ,nr 1 if gi has a root, and = 0 if not.
     */
    std::vector<int> rootsfound;

    /** simulation states history at timepoints  */
    std::map<realtype, SimulationState> timepoint_states;

    /** simulation state history at events*/
    std::vector<SimulationState> event_states;

    /** simulation state after initialization*/
    SimulationState initial_state;

    /** simulation state after simulation*/
    SimulationState final_state;

    /** state vector (dimension: nx_solver) */
    AmiVector x;

    /** old state vector (dimension: nx_solver) */
    AmiVector x_old;

    /** differential state vector (dimension: nx_solver) */
    AmiVector dx;

    /** old differential state vector (dimension: nx_solver) */
    AmiVector dx_old;

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

    /** storage for last found root */
    realtype tlastroot = 0.0;

    /** flag to indicate wheter solver was preeinitialized via preequilibration */
    bool preequilibrated = false;

    /** current iteration number for time index */
    int it;

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
    explicit FinalStateStorer(ForwardProblem *fwd) : fwd(fwd) {
    }

    FinalStateStorer &operator=(const FinalStateStorer &other) = delete;

    /**
     * @brief destructor, stores simulation state
     */
    ~FinalStateStorer() {
        if(fwd)
            fwd->final_state = fwd->getSimulationState();
    }
  private:
    ForwardProblem *fwd;
};

} // namespace amici

#endif // FORWARDPROBLEM_H
