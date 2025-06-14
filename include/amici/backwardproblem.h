#ifndef AMICI_BACKWARDPROBLEM_H
#define AMICI_BACKWARDPROBLEM_H

#include "amici/defines.h"
#include "amici/forwardproblem.h"
#include "amici/vector.h"

#include <vector>

namespace amici {

class ExpData;
class Solver;
class Model;
class ForwardProblem;
class SteadystateProblem;

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
     */
    explicit BackwardProblem(ForwardProblem& fwd);

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
    realtype gett() const { return t_; }

    /**
     * @brief Accessor for which
     * @return which
     */
    int getwhich() const { return which; }

    /**
     * @brief Accessor for pointer to which
     * @return which
     */
    int* getwhichptr() { return &which; }

    /**
     * @brief Accessor for dJydx
     * @return dJydx
     */
    std::vector<realtype> const& getdJydx() const { return dJydx_; }

    /**
     * @brief Accessor for xB
     * @return xB
     */
    AmiVector const& getAdjointState() const { return xB_; }

    /**
     * @brief Accessor for xQB
     * @return xQB
     */
    AmiVector const& getAdjointQuadrature() const { return xQB_; }

  private:
    void handlePostequilibration();

    /**
     * @brief Execute everything necessary for the handling of events
     * for the backward problem
     * @param disc The discontinuity to handle
     */
    void handleEventB(Discontinuity const& disc);

    /**
     * @brief Execute everything necessary for the handling of data
     * points for the backward problems
     *
     * @param it index of data point
     */
    void handleDataPointB(int it);

    /**
     * @brief Compute the next timepoint to integrate to.
     *
     * This is the maximum of tdata and troot but also takes into account if
     * it<0 or iroot<0 where these expressions do not necessarily make sense.
     *
     * @param it index of next data point
     * @return tnext next timepoint
     */
    realtype getTnext(int it);

    Model* model_;
    Solver* solver_;
    ExpData const* edata_;

    /** current time */
    realtype t_;
    /** adjoint state vector */
    AmiVector xB_;
    /** differential adjoint state vector */
    AmiVector dxB_;
    /** quadrature state vector */
    AmiVector xQB_;
    /** sensitivity state vector array */
    AmiVectorArray sx0_;
    /** array of number of found roots for a certain event type */
    std::vector<int> nroots_;
    /** array containing the time-points of discontinuities*/
    std::vector<Discontinuity> discs_;
    /** index of the backward problem */
    int which = 0;

    /** state derivative of data likelihood */
    std::vector<realtype> dJydx_;
    /** state derivative of event likelihood */
    std::vector<realtype> const dJzdx_;

    /** The preequilibration steadystate problem from the forward problem. */
    SteadystateProblem* preeq_problem_;

    /** The postequilibration steadystate problem from the forward problem. */
    SteadystateProblem* posteq_problem_;
};

} // namespace amici

#endif // BACKWARDPROBLEM_H
