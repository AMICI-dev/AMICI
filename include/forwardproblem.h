#ifndef AMICI_FORWARDPROBLEM_H
#define AMICI_FORWARDPROBLEM_H

#include "include/amici_defines.h"
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
                   Model *model, Solver *solver) :
    dJzdx(model->nJ() * model->nx() * udata->nmaxevent(), 0.0),
    dJydx(model->nJ() * model->nx() * udata->nt(), 0.0) {
        
    }
    
    static void workForwardProblem();
    
  private:

    static void handleEvent(realtype *tlastroot, const UserData *udata,
                           ReturnData *rdata, const ExpData *edata,
                           TempData *tdata, int seflag, Solver *solver,
                           Model *model);

    static void storeJacobianAndDerivativeInReturnData(TempData *tdata,
                                                      ReturnData *rdata,
                                                      Model *model);

    static void getEventOutput(const UserData *udata, ReturnData *rdata,
                              const ExpData *edata, TempData *tdata,
                              Model *model);

    static void prepEventSensis(int ie, ReturnData *rdata, const ExpData *edata,
                               TempData *tdata, Model *model);

    static void getEventSensisFSA(int ie, ReturnData *rdata,
                                 const ExpData *edata, TempData *tdata,
                                 Model *model);

    static void handleDataPoint(int it, const UserData *udata, ReturnData *rdata,
                               const ExpData *edata, TempData *tdata,
                               Solver *solver, Model *model);

    static void getDataOutput(int it, const UserData *udata, ReturnData *rdata,
                             const ExpData *edata, TempData *tdata,
                             Solver *solver, Model *model);

    static void prepDataSensis(int it, ReturnData *rdata, const ExpData *edata,
                              TempData *tdata, Model *model);

    static void getDataSensisFSA(int it, const UserData *udata,
                                ReturnData *rdata, const ExpData *edata,
                                TempData *tdata, Solver *solver, Model *model);

    static void applyEventBolus(TempData *tdata, Model *model);

    static void applyEventSensiBolusFSA(TempData *tdata, Model *model);

    static void updateHeaviside(TempData *tdata, const int ne);

  private:
    
    /** data likelihood */
    std::vector<double> Jy;
    /** event likelihood */
    std::vector<double> Jz;

    /** state derivative of data likelihood */
    std::vector<double> dJydx;
    /** parameter derivative of data likelihood */
    std::vector<double> dJydp;
    /** parameter derivative of event likelihood */
    std::vector<double> dJzdp;
    /** state derivative of event likelihood */
    std::vector<double> dJzdx;
    
    /** current time */
    realtype t;
    
    /** array of flags indicating which root has beend found.
     *  array of length nr with the indices of the user functions gi found to
     * have a
     *  root. For i = 0, . . . ,nr 1 if gi has a root, and = 0 if not.
     */
    int *rootsfound = nullptr;
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
    
    /** integer for indexing of backwards problems */
    int which = 0;
    
    /** array containing the time-points of discontinuities*/
    realtype *discs = nullptr;
    /** array containing the index of discontinuities */
    realtype *irdiscs = nullptr;
    
    /** current root index, will be increased during the forward solve and
     * decreased during backward solve */
    int iroot = 0;
    
    
    /** state vector */
    AmiVector x;
    /** old state vector */
    AmiVector x_old;
    /** array of state vectors at discontinuities*/
    AmiVector *x_disc;
    /** array of differential state vectors at discontinuities*/
    AmiVector *xdot_disc;
    /** array of old differential state vectors at discontinuities*/
    AmiVector *xdot_old_disc;
    /** differential state vector */
    AmiVector dx;
    /** old differential state vector */
    AmiVector dx_old;
    /** time derivative state vector */
    AmiVector xdot;
    /** old time derivative state vector */
    AmiVector xdot_old;
    /** adjoint state vector */
    AmiVector xB;
    /** old adjoint state vector */
    AmiVector xB_old;
    /** differential adjoint state vector */
    AmiVector dxB;
    /** quadrature state vector */
    AmiVector xQB;
    /** old quadrature state vector */
    AmiVector xQB_old;
    /** sensitivity state vector array */
    AmiVectorArray *sx;
    /** differential sensitivity state vector array */
    AmiVectorArray *sdx;
    
    ForwardProblem();
};


} // namespace amici

#endif // FORWARDPROBLEM_H
