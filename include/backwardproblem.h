#ifndef AMICI_BACKWARDPROBLEM_H
#define AMICI_BACKWARDPROBLEM_H

#include "include/amici_defines.h"
#include "include/amici_vector.h"
#include <vector>

namespace amici {

class UserData;
class TempData;
class ReturnData;
class ExpData;
class Solver;
class Model;

//!  class to solve backwards problems.
/*!
  solves the backwards problem for adjoint sensitivity analysis and handles
  events and data-points
*/

class BackwardProblem {
  public:
    static void workBackwardProblem();

    BackwardProblem(const UserData *udata,
                    ReturnData *rdata, const ExpData *edata,
                    Model *model, Solver *solver) :
    xB(model->nx), xB_old(model->nx),
    dxB(model->nx), xQB(model->nx), xQB_old(model->nx)
    {
      // TODO: TBD
    };

  private:
    
    Model *model;
    ReturnData *rdata;
    Solver *solver;
    
    static void handleEventB(int iroot);
    
    static void handleDataPointB(int it);
    
    static void updateHeavisideB(int iroot, int ne);
    
    static realtype getTnext(realtype *troot, int iroot, realtype *tdata,
                             int it);
    
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
    /** array of number of found roots for a certain event type */
    int *nroots = nullptr;
    
};

} // namespace amici

#endif // BACKWARDPROBLEM_H
