#ifndef AMICI_BACKWARDPROBLEM_H
#define AMICI_BACKWARDPROBLEM_H

#include "amici/defines.h"
#include "amici/vector.h"

#include <vector>

namespace amici {

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

    BackwardProblem(const ForwardProblem *fwd);
    
    /** accessor for t
     * @return t
     */
    realtype gett() const {
        return t;
    }
    
    /** accessor for which
     * @return which
     */
    int getwhich() const {
        return which;
    }
    
    /** accessor for pointer to which
     * @return which
     */
    int *getwhichptr() {
        return &which;
    }
    
    /** accessor for pointer to xB
     * @return &xB
     */
    AmiVector *getxBptr() {
        return &xB;
    }
    
    /** accessor for pointer to xQB
     * @return &xQB
     */
    AmiVector *getxQBptr() {
        return &xQB;
    }
    
    /** accessor for pointer to dxB
     * @return &dxB
     */
    AmiVector *getdxBptr() {
        return &dxB;
    }
    
    /** accessor for dJydx
     * @return dJydx
     */
    std::vector<realtype> const& getdJydx() const {
        return dJydx;
    }

  private:
        
    void handleEventB(int iroot);
    
    void handleDataPointB(int it);
    
    void updateHeavisideB(int iroot);
    
    realtype getTnext(std::vector<realtype> const& troot, const int iroot, const int it);
    
    Model *model;
    ReturnData *rdata;
    Solver *solver;

    /** current time */
    realtype t;
    /** parameter derivative of likelihood array */
    std::vector<realtype> llhS0;
    /** adjoint state vector */
    AmiVector xB;
    /** differential adjoint state vector */
    AmiVector dxB;
    /** quadrature state vector */
    AmiVector xQB;
    /** array of state vectors at discontinuities*/
    const AmiVectorArray x_disc;
    /** array of differential state vectors at discontinuities*/
    const AmiVectorArray xdot_disc;
    /** array of old differential state vectors at discontinuities*/
    const AmiVectorArray xdot_old_disc;
    /** sensitivity state vector array */
    AmiVectorArray sx;
    /** array of number of found roots for a certain event type */
    std::vector<int> nroots;
    /** array containing the time-points of discontinuities*/
    const std::vector<realtype> discs;
    /** array containing the index of discontinuities */
    const std::vector<realtype> irdiscs;
    /** index of the backward problem */
    int which = 0;
    /** current root index, will be increased during the forward solve and
     * decreased during backward solve */
    int iroot = 0;
    /** array of index which root has been found */
    const std::vector<int> rootidx;
    
    /** state derivative of data likelihood */
    const std::vector<realtype> dJydx;
    /** state derivative of event likelihood */
    const std::vector<realtype> dJzdx;
};

} // namespace amici

#endif // BACKWARDPROBLEM_H
