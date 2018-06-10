#ifndef AMICI_FORWARDPROBLEM_H
#define AMICI_FORWARDPROBLEM_H

#include "amici/defines.h"
#include "amici/vector.h"

#include <sundials/sundials_direct.h>
#include <vector>

namespace amici {

class ReturnData;
class ExpData;
class Solver;
class Model;

/**
 * @brief The ForwardProblem class groups all functions for solving the
 * backwards problem.
 * Has only static members.
 */
class ForwardProblem {
  public:
    ForwardProblem(ReturnData *rdata, const ExpData *edata,
                   Model *model, Solver *solver);

    /** destructor */
    ~ForwardProblem() {
        DestroyMat(Jtmp);
    }
    
    void workForwardProblem();
    
    /** pointer to model instance */
    Model *model;
    /** pointer to return data instance */
    ReturnData *rdata;
    /** pointer to solver instance */
    Solver *solver;
    /** pointer to experimental data instance */
    const ExpData *edata;
    
    /** accessor for t
     * @return t
     */
    realtype getTime() const {
        return t;
    }
    
    /** accessor for sx
     * @return sx
     */
    AmiVectorArray getStateSensitivity() const {
        return sx;
    }
    
    /** accessor for x_disc
     * @return x_disc
     */
    AmiVectorArray getStatesAtDiscontinuities() const {
        return x_disc;
    }
    
    /** accessor for xdot_disc
     * @return xdot_disc
     */
    AmiVectorArray getRHSAtDiscontinuities() const {
        return xdot_disc;
    }
    
    /** accessor for xdot_old_disc
     * @return xdot_old_disc
     */
    AmiVectorArray getRHSBeforeDiscontinuities() const {
        return xdot_old_disc;
    }
    
    /** accessor for nroots
     * @return nroots
     */
    std::vector<int> getNumberOfRoots() const {
        return nroots;
    }
    
    /** accessor for discs
     * @return discs
     */
    std::vector<realtype> getDiscontinuities() const {
        return discs;
    }
    
    /** accessor for rootidx
     * @return rootidx
     */
    std::vector<int> getRootIndexes() const {
        return rootidx;
    }
    
    /** accessor for dJydx
     * @return dJydx
     */
   std::vector<realtype> getDJydx() const {
        return dJydx;
    }
    
    /** accessor for dJzdx
     * @return dJzdx
     */
    std::vector<realtype> getDJzdx() const {
        return dJzdx;
    }
    
    /** accessor for iroot
     * @return iroot
     */
    int getRootCounter() const {
        return iroot;
    }
    
    /** accessor for pointer to x
     * @return &x
     */
    AmiVector *getStatePointer() {
        return &x;
    }
    
    /** accessor for pointer to dx
     * @return &dx
     */
    AmiVector *getStateDerivativePointer() {
        return &dx;
    }
    
    /** accessor for pointer to sx
     * @return &sx
     */
    AmiVectorArray *getStateSensitivityPointer() {
        return &sx;
    }
    
    /** accessor for pointer to sdx
     * @return &sdx
     */
    AmiVectorArray *getStateDerivativeSensitivityPointer() {
        return &sdx;
    }
    
  private:

    void handleEvent(realtype *tlastroot,const bool seflag);

    void storeJacobianAndDerivativeInReturnData();

    void getEventOutput();

    void prepEventSensis(int ie);

    void getEventSensisFSA(int ie);

    void handleDataPoint(int it);

    void getDataOutput(int it);

    void prepDataSensis(int it);

    void getDataSensisFSA(int it);

    void applyEventBolus();

    void applyEventSensiBolusFSA();
    
    /** array of index which root has been found  (dimension: ne * ne * ne * nmaxevent, ordering = ?) */
    std::vector<int> rootidx;
    /** array of number of found roots for a certain event type (dimension: ne) */
    std::vector<int> nroots;
    /** array of values of the root function (dimension: ne) */
    std::vector<realtype> rootvals;
    /** temporary rootval storage to check crossing in secondary event (dimension: ne) */
    std::vector<realtype> rvaltmp;
    
    /** array containing the time-points of discontinuities (dimension: nmaxevent x ne, ordering = ?) */
    std::vector<realtype> discs;
    /** array containing the index of discontinuities (dimension: nmaxevent x ne, ordering = ?) */
    std::vector<realtype> irdiscs;
    
    /** current root index, will be increased during the forward solve and
     * decreased during backward solve */
    int iroot = 0;
    
    /** array of state vectors at discontinuities  (dimension nx x nMaxEvent * ne, ordering =?) */
    AmiVectorArray x_disc;
    /** array of differential state vectors at discontinuities (dimension nx x nMaxEvent * ne, ordering =?) */
    AmiVectorArray xdot_disc;
    /** array of old differential state vectors at discontinuities (dimension nx x nMaxEvent * ne, ordering =?) */
    AmiVectorArray xdot_old_disc;
    
    /** state derivative of data likelihood (dimension nJ x nx x nt, ordering =?) */
    std::vector<realtype> dJydx;
    /** state derivative of event likelihood (dimension nJ x nx x nMaxEvent, ordering =?) */
    std::vector<realtype> dJzdx;
    
    /** data likelihood */
    std::vector<realtype> Jy; // TODO remove?
    /** event likelihood */
    std::vector<realtype> Jz; // TODO remove?

    /** parameter derivative of data likelihood */ // TODO remove?
    std::vector<realtype> dJydp;
    /** parameter derivative of event likelihood */  // TODO remove?
    std::vector<realtype> dJzdp;
    
    /** current time */
    realtype t;
    
    /** array of flags indicating which root has beend found.
     *  array of length nr (ne) with the indices of the user functions gi found to
     * have a root. For i = 0, . . . ,nr 1 if gi has a root, and = 0 if not.
     */
    std::vector<int> rootsfound;

    /** temporary storage of Jacobian, kept here to avoid reallocation (dimension: nx x nx, col-major) */
    DlsMat Jtmp;
    
    /** state vector (dimension: nx) */
    AmiVector x;
    /** old state vector (dimension: nx) */
    AmiVector x_old;
    /** differential state vector (dimension: nx) */
    AmiVector dx;
    /** old differential state vector (dimension: nx) */
    AmiVector dx_old;
    /** time derivative state vector (dimension: nx) */
    AmiVector xdot;
    /** old time derivative state vector (dimension: nx) */
    AmiVector xdot_old;

    /** sensitivity state vector array (dimension: nx x nplist, row-major) */
    AmiVectorArray sx;
    /** differential sensitivity state vector array (dimension: nx x nplist, row-major) */
    AmiVectorArray sdx;
};


} // namespace amici

#endif // FORWARDPROBLEM_H
