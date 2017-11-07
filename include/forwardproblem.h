#ifndef AMICI_FORWARDPROBLEM_H
#define AMICI_FORWARDPROBLEM_H

#include "include/amici_defines.h"
#include "include/amici_vector.h"
#include <sundials/sundials_direct.h>
#include <vector>
namespace amici {

class UserData;
class TempData;
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
    ForwardProblem(const UserData *udata,
                   ReturnData *rdata, const ExpData *edata,
                   Model *model, Solver *solver);
    ~ForwardProblem() {
        DestroyMat(Jtmp);
    };
    
    
    void workForwardProblem();
    
    Model *model;
    ReturnData *rdata;
    Solver *solver;
    const UserData *udata;
    const ExpData *edata;
    
    /** array of index which root has been found */
    int *rootidx = nullptr;
    /** array of number of found roots for a certain event type */
    int *nroots = nullptr;
    /** array of values of the root function */
    realtype *rootvals = nullptr;
    /** temporary rootval storage to check crossing in secondary event */
    realtype *rvaltmp = nullptr;
    
    /** flag indicating whether a certain heaviside function should be active or
     not */
    realtype *h = nullptr;
    
    /** array containing the time-points of discontinuities*/
    realtype *discs = nullptr;
    /** array containing the index of discontinuities */
    realtype *irdiscs = nullptr;
    
    /** current root index, will be increased during the forward solve and
     * decreased during backward solve */
    int iroot = 0;
    
    /** array of state vectors at discontinuities*/
    AmiVector *x_disc;
    /** array of differential state vectors at discontinuities*/
    AmiVector *xdot_disc;
    /** array of old differential state vectors at discontinuities*/
    AmiVector *xdot_old_disc;
    
    /** state derivative of data likelihood */
    std::vector<double> dJydx;
    /** state derivative of event likelihood */
    std::vector<double> dJzdx;
    
    
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

    void updateHeaviside();
    
    /** data likelihood */
    std::vector<double> Jy;
    /** event likelihood */
    std::vector<double> Jz;

    /** parameter derivative of data likelihood */
    std::vector<double> dJydp;
    /** parameter derivative of event likelihood */
    std::vector<double> dJzdp;
    
    /** current time */
    realtype t;
    
    /** array of flags indicating which root has beend found.
     *  array of length nr with the indices of the user functions gi found to
     * have a
     *  root. For i = 0, . . . ,nr 1 if gi has a root, and = 0 if not.
     */
    int *rootsfound = nullptr;

    DlsMat Jtmp;
    
    /** state vector */
    AmiVector x;
    /** old state vector */
    AmiVector x_old;
    /** differential state vector */
    AmiVector dx;
    /** old differential state vector */
    AmiVector dx_old;
    /** time derivative state vector */
    AmiVector xdot;
    /** old time derivative state vector */
    AmiVector xdot_old;

    /** sensitivity state vector array */
    AmiVectorArray sx;
    /** differential sensitivity state vector array */
    AmiVectorArray sdx;
    
    ForwardProblem();
};


} // namespace amici

#endif // FORWARDPROBLEM_H
