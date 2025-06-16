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

/**
 * @brief The BwdSimWorkspace class is used to store temporary simulation
 * state during backward simulations.
 */
struct BwdSimWorkspace {
    /**
     * @brief Constructor
     * @param model The model for which to set up the workspace.
     * @param solver The solver for which to set up this workspace.
     */
    BwdSimWorkspace(
        gsl::not_null<Model*> model, gsl::not_null<Solver const*> solver
    );

    /** The model. */
    Model* model_;

    /** adjoint state vector */
    AmiVector xB_;
    /** differential adjoint state vector */
    AmiVector dxB_;
    /** quadrature state vector */
    AmiVector xQB_;

    /** array of number of found roots for a certain event type */
    std::vector<int> nroots_;
    /** array containing the time-points of discontinuities*/
    std::vector<Discontinuity> discs_;
    /** index of the backward problem */
    int which = 0;
};

/**
 * @brief The EventHandlingBwdSimulator class runs a backward simulation
 * and processes events and measurements general.
 */
class EventHandlingBwdSimulator {
  public:
    /**
     * @brief EventHandlingBwdSimulator constructor.
     * @param model The model to simulate.
     * @param solver The solver to use for the simulation.
     * @param ws The workspace to use for the simulation.
     */
    EventHandlingBwdSimulator(
        gsl::not_null<Model*> model, gsl::not_null<Solver*> solver,
        gsl::not_null<BwdSimWorkspace*> ws
    )
        : model_(model)
        , solver_(solver)
        , ws_(ws) {};

    /**
     * @brief Run the simulation.
     *
     * It will run the backward simulation from the initial time of this period
     * to the final timepoint of this period, handling events
     * and data points as they occur.
     *
     * Expects the model and the solver to be set up, and `ws` to be initialized
     * for this period.
     *
     * @param t_start The initial time of this period.
     * @param t_end The final time of this period.
     * @param it The index of the timepoint in `timepoints` to start with.
     * @param timepoints The output timepoints or measurement timepoints of
     * this period. This must contain at least the final timepoint of this
     * period.
     * @param dJydx State-derivative of data likelihood. Must be non-null if
     * there are any data points in this period.
     * @param dJzdx State-derivative of event likelihood. Must be non-null if
     * the model has any event-observables.
     */
    void
    run(realtype t_start, realtype t_end, realtype it,
        std::vector<realtype> const& timepoints,
        std::vector<realtype> const* dJydx, std::vector<realtype> const* dJzdx);

  private:
    /**
     * @brief Execute everything necessary for the handling of events
     * for the backward problem
     * @param disc The discontinuity to handle
     * @param dJzdx State-derivative of event likelihood
     */
    void
    handleEventB(Discontinuity const& disc, std::vector<realtype> const* dJzdx);

    /**
     * @brief Execute everything necessary for the handling of data
     * points for the backward problems
     *
     * @param it index of data point
     * @param dJydx State-derivative of data likelihood
     */
    void handleDataPointB(int it, std::vector<realtype> const* dJydx);

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

    /** The model to simulate. */
    Model* model_;

    /** The solver to use for the simulation. */
    Solver* solver_;

    /** The workspace to use for the simulation. */
    gsl::not_null<BwdSimWorkspace*> ws_;

    /** current time */
    realtype t_{0};
};

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
     * If adjoint sensitivities are enabled, this will also compute
     * sensitivities. workForwardProblem must be called before this function is
     * called.
     */
    void workBackwardProblem();

    /**
     * @brief Accessor for xB
     * @return xB
     */
    [[nodiscard]] AmiVector const& getAdjointState() const { return ws_.xB_; }

    /**
     * @brief Accessor for xQB
     * @return xQB
     */
    [[nodiscard]] AmiVector const& getAdjointQuadrature() const {
        return ws_.xQB_;
    }

  private:
    void handlePostequilibration();

    Model* model_;
    Solver* solver_;
    ExpData const* edata_;

    /** current time */
    realtype t_;
    /** sensitivity state vector array */
    AmiVectorArray sx0_;
    /** array of number of found roots for a certain event type */
    std::vector<int> nroots_;
    /** array containing the time-points of discontinuities*/
    std::vector<Discontinuity> discs_;

    /** state derivative of data likelihood */
    std::vector<realtype> dJydx_;
    /** state derivative of event likelihood */
    std::vector<realtype> const dJzdx_;

    /** The preequilibration steadystate problem from the forward problem. */
    SteadystateProblem* preeq_problem_;

    /** The postequilibration steadystate problem from the forward problem. */
    SteadystateProblem* posteq_problem_;

    BwdSimWorkspace ws_;

    EventHandlingBwdSimulator simulator_;
};

} // namespace amici

#endif // AMICI_BACKWARDPROBLEM_H
