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
    static void workBackwardProblem(const UserData *udata, TempData *tdata,
                                   ReturnData *rdata, Model *model);

    static void handleEventB(int iroot, TempData *tdata, Model *model);

    static void handleDataPointB(int it, ReturnData *rdata, TempData *tdata,
                                Solver *solver, Model *model);

    static void updateHeavisideB(int iroot, TempData *tdata, int ne);

    static realtype getTnext(realtype *troot, int iroot, realtype *tdata,
                             int it, Model *model);

  private:
    
    /** parameter derivative of likelihood array */
    std::vector<double> llhS0;
    
    BackwardProblem();
};

} // namespace amici

#endif // BACKWARDPROBLEM_H
