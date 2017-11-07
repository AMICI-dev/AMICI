#ifndef AMICI_BACKWARDPROBLEM_H
#define AMICI_BACKWARDPROBLEM_H

#include "include/amici_defines.h"
#include "include/amici_vector.h"
#include <vector>

namespace amici {

class UserData;
class ReturnData;
class ExpData;
class Solver;
class Model;
class ForwardProblem;

//!  class to solve backwards problems.
/*!
  solves the backwards problem for adjoint sensitivity analysis and handles
  events and data-points
*/

class BackwardProblem {
  public:
    void workBackwardProblem();

    BackwardProblem(ForwardProblem fwd);

  private:
    
    Model *model;
    ReturnData *rdata;
    Solver *solver;
    const UserData *udata;
    
    void handleEventB(int iroot);
    
    void handleDataPointB(int it);
    
    void updateHeavisideB(int iroot);
    
    realtype getTnext(realtype *troot, int iroot, int it);
    
    /** current time */
    realtype t;
    /** parameter derivative of likelihood array */
    std::vector<double> llhS0;
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
    /** array of state vectors at discontinuities*/
    AmiVector *x_disc;
    /** array of differential state vectors at discontinuities*/
    AmiVector *xdot_disc;
    /** array of old differential state vectors at discontinuities*/
    AmiVector *xdot_old_disc;
    /** sensitivity state vector array */
    AmiVectorArray sx;
    /** array of number of found roots for a certain event type */
    int *nroots = nullptr;
    /** array containing the time-points of discontinuities*/
    realtype *discs = nullptr;
    /** array containing the index of discontinuities */
    realtype *irdiscs = nullptr;
    /** index of the backward problem */
    int which = 0;
    /** flag indicating whether a certain heaviside function should be active or
     not */
    realtype *h = nullptr;
    /** current root index, will be increased during the forward solve and
     * decreased during backward solve */
    int iroot = 0;
    /** array of index which root has been found */
    int *rootidx = nullptr;
    
    /** state derivative of data likelihood */
    std::vector<double> dJydx;
    /** state derivative of event likelihood */
    std::vector<double> dJzdx;

    
    
    
};

} // namespace amici

#endif // BACKWARDPROBLEM_H
