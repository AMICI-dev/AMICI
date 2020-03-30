#ifndef AMICI_FORWARDPROBLEM_H
#define AMICI_FORWARDPROBLEM_H

#include "amici/defines.h"
#include "amici/vector.h"
#include "amici/sundials_matrix_wrapper.h"

#include <sundials/sundials_direct.h>
#include <vector>
#include <memory>

namespace amici {

class ReturnData;
class ExpData;
class Solver;
class Model;
class SteadystateProblem;

/**
 * @brief The ForwardProblem class groups all functions for solving the
 * forward problem.f
 */
class ForwardProblem {
  public:
    /**
     * @brief Constructor
     * @param rdata pointer to ReturnData instance
     * @param edata pointer to ExpData instance
     * @param model pointer to Model instance
     * @param solver pointer to Solver instance
     * @param preeq preequilibration with which to initialize the forward problem,
     * pass nullptr for no initialization
     */
    ForwardProblem(ReturnData *rdata, const ExpData *edata,
                   Model *model, Solver *solver,
                   const SteadystateProblem *preeq);

    ~ForwardProblem() = default;

    /**
     * @brief Solve the forward problem.
     *
     * If forward sensitivities are enabled this will also compute sensitivies.
     */
    void workForwardProblem();

    /**
     * @brief Accessor for t
     * @return t
     */
    realtype getTime() const {
        return t;
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
    AmiVectorArray const& getStatesAtDiscontinuities() const {
        return x_disc;
    }

    /**
     * @brief Accessor for xdot_disc
     * @return xdot_disc
     */
    AmiVectorArray const& getRHSAtDiscontinuities() const {
        return xdot_disc;
    }

    /**
     * @brief Accessor for xdot_old_disc
     * @return xdot_old_disc
     */
    AmiVectorArray const& getRHSBeforeDiscontinuities() const {
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
    std::vector<int> const& getRootIndexes() const {
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
     * @brief Accessor for iroot
     * @return iroot
     */
    int getRootCounter() const {
        return iroot;
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
    int getCurrentTimeIteration() {
        return it;
    }
    
    /**
     * @brief Accessor for x_timepoints
     * @return x_timepoints.at(it)
     */
    const AmiVector &getInitialState() const {
        return x0;
    }
    
    /**
     * @brief Accessor for sx_timepoints
     * @return sx_timepoints.at(it)
     */
    const AmiVectorArray &getInitialStateSensitivity() const {
        return sx0;
    }
    
    
    /**
     * @brief Accessor for x_timepoints
     * @param it timepoint index
     * @return x_timepoints.at(it)
     */
    const AmiVector &getStateTimePoint(int it) const {
        return x_timepoints.at(it);
    }
    
    /**
     * @brief Accessor for sx_timepoints
     * @param it timepoint index
     * @return sx_timepoints.at(it)
     */
    const AmiVectorArray &getStateSensitivityTimePoint(int it) const {
        return sx_timepoints.at(it);
    }
    
    /**
     * @brief Accessor for x_timepoints
     * @param it timepoint index
     * @return x_timepoints.at(it)
     */
    const AmiVector &getStateEvent(int ie) const {
        return x_events.at(ie);
    }
    
    /**
     * @brief Accessor for sx_timepoints
     * @param it timepoint index
     * @return sx_timepoints.at(it)
     */
    const AmiVectorArray &getStateSensitivityEvent(int ie) const {
        return sx_events.at(ie);
    }
    
    /**
     * @brief Compute updates for backwards (ajoint) problem
     * @return &sdx
     */
    void getAdjointUpdates();

    /** pointer to model instance */
    Model *model;

    /** pointer to return data instance */
    ReturnData *rdata;

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
     * @brief Evaluates the Jacobian and differential equation right hand side,
     * stores it in RetunData
     */
    void storeJacobianAndDerivativeInReturnData();

    /**
     * @brief Extract output information for events
     */
    void getEventOutput();

    /**
     * @brief Extracts event information for forward sensitivity analysis
     *
     * @param ie index of event type
     */
    void getEventSensisFSA(int ie);

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

    /** array of index which root has been found
     * (dimension: ne * ne * nmaxevent, ordering = ?) */
    std::vector<int> rootidx;

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

    /** current root index, will be increased during the forward solve and
     * decreased during backward solve */
    int iroot = 0;

    /** array of state vectors at discontinuities
     * (dimension nx x nMaxEvent * ne, ordering =?) */
    AmiVectorArray x_disc;

    /** array of differential state vectors at discontinuities
     * (dimension nx x nMaxEvent * ne, ordering =?) */
    AmiVectorArray xdot_disc;

    /** array of old differential state vectors at discontinuities
     * (dimension nx x nMaxEvent * ne, ordering =?) */
    AmiVectorArray xdot_old_disc;

    /** state derivative of data likelihood
     * (dimension nJ x nx x nt, ordering =?) */
    std::vector<realtype> dJydx;

    /** state derivative of event likelihood
     * (dimension nJ x nx x nMaxEvent, ordering =?) */
    std::vector<realtype> dJzdx;

    /** current time */
    realtype t;

    /**
     * @brief Array of flags indicating which root has beend found.
     *
     * Array of length nr (ne) with the indices of the user functions gi found
     * to have a root. For i = 0, . . . ,nr 1 if gi has a root, and = 0 if not.
     */
    std::vector<int> rootsfound;

    /** temporary storage of Jacobian, kept here to avoid reallocation
     * (dimension: nx x nx, col-major) */
    SUNMatrixWrapper Jtmp;

    /** state vector history at timepoints  */
    std::vector<AmiVector> x_timepoints;
    
    /** state sensitivity vector history at timepoints */
    std::vector<AmiVectorArray> sx_timepoints;
    
    /** state vector history at events*/
    std::vector<AmiVector> x_events;
    
    /** state sensitivity vector history at events*/
    std::vector<AmiVectorArray> sx_events;
    
    /** initial state */
    AmiVector x0;
    
    /** initial state sensitivity */
    AmiVectorArray sx0;
    
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


} // namespace amici

#endif // FORWARDPROBLEM_H
