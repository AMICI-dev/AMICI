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
    /**
     * @brief Construct backward problem from forward problem
     * @param fwd pointer to corresponding forward problem
     * @return new BackwardProblem instance
     */
    explicit BackwardProblem(const ForwardProblem *fwd);

    /**
     * @brief Solve the backward problem.
     *
     * If adjoint sensitivities are enabled this will also compute
     * sensitivities. workForwardProblem must be called before this is
     * function is called.
     */
    void workBackwardProblem();

    /**
     * @brief Accessor for current time t
     * @return t
     */
    realtype gett() const {
        return t;
    }

    /**
     * @brief Accessor for which
     * @return which
     */
    int getwhich() const {
        return which;
    }

    /**
     * @brief Accessor for pointer to which
     * @return which
     */
    int *getwhichptr() {
        return &which;
    }

    /**
     * @brief Accessor for dJydx
     * @return dJydx
     */
    std::vector<realtype> const& getdJydx() const {
        return dJydx;
    }

  private:
    /**
     * @brief Execute everything necessary for the handling of events
     * for the backward problem
     *
     * @param iroot index of event @type int
     */
    void handleEventB(int iroot);

    /**
     * @brief Execute everything necessary for the handling of data
     * points for the backward problems
     *
     * @param it index of data point @type int
     */
    void handleDataPointB(int it);


    /**
     * @brief Compute the next timepoint to integrate to.
     *
     * This is the maximum of tdata and troot but also takes into account if
     * it<0 or iroot<0 where these expressions do not necessarily make sense.
     *
     * @param troot timepoint of next event @type realtype
     * @param iroot index of next event @type int
     * @param it index of next data point @type int
     * @param model pointer to model specification object @type Model
     * @return tnext next timepoint @type realtype
     */
    realtype getTnext(std::vector<realtype> const& troot, int iroot, int it);

    /**
     * @brief Compute likelihood sensitivities.
     */
    void computeLikelihoodSensitivities();

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
    AmiVectorArray sx0;
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
