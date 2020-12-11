#ifndef AMICI_SOLVER_H
#define AMICI_SOLVER_H

#include "amici/amici.h"
#include "amici/defines.h"
#include "amici/sundials_linsol_wrapper.h"
#include "amici/symbolic_functions.h"
#include "amici/vector.h"

#include <cmath>
#include <functional>
#include <memory>

namespace amici {

class ReturnData;
class ForwardProblem;
class BackwardProblem;
class Model;
class Solver;
class AmiciApplication;

extern AmiciApplication defaultContext;
} // namespace amici

// for serialization friend in Solver
namespace boost {
namespace serialization {
template <class Archive>
void serialize(Archive &ar, amici::Solver &s, unsigned int version);
}
} // namespace boost::serialization

namespace amici {

/**
 * The Solver class provides a generic interface to CVODES and IDAS solvers,
 * individual realizations are realized in the CVodeSolver and the IDASolver
 * class. All transient private/protected members (CVODES/IDAS memory, interface
 * variables and status flags) are specified as mutable and not included in
 * serialization or equality checks. No solver setting parameter should be
 * marked mutable.
 *
 * NOTE: Any changes in data members here must be propagated to copy ctor,
 * equality operator, serialization functions in serialization.h, and
 * amici::hdf5::readSolverSettingsFromHDF5 in hdf5.cpp.
 */
class Solver {
  public:
    Solver() = default;

    /**
     * @brief Constructor
     * @param app AMICI application context
     */
    Solver(AmiciApplication *app);

    /**
     * @brief Solver copy constructor
     * @param other
     */
    Solver(const Solver &other);

    virtual ~Solver() = default;

    /**
     * @brief Clone this instance
     * @return The clone
     */
    virtual Solver *clone() const = 0;

    /**
     * @brief runs a forward simulation until the specified timepoint
     *
     * @param tout next timepoint
     * @return status flag
     */
    int run(realtype tout) const;

    /**
     * @brief makes a single step in the simulation
     *
     * @param tout next timepoint
     * @return status flag
     */
    int step(realtype tout) const;

    /**
     * @brief runs a backward simulation until the specified timepoint
     *
     * @param tout next timepoint
     */
    void runB(realtype tout) const;

    /**
     * @brief Initializes the ami memory object and applies specified options
     * @param t0 initial timepoint
     * @param model pointer to the model instance
     * @param x0 initial states
     * @param dx0 initial derivative states
     * @param sx0 initial state sensitivities
     * @param sdx0 initial derivative state sensitivities
     */

    void setup(realtype t0, Model *model, const AmiVector &x0,
               const AmiVector &dx0, const AmiVectorArray &sx0,
               const AmiVectorArray &sdx0) const;

    /**
     * @brief Initializes the AMI memory object for the backwards problem
     * @param which index of the backward problem, will be set by this routine
     * @param tf final timepoint (initial timepoint for the bwd problem)
     * @param model pointer to the model instance
     * @param xB0 initial adjoint states
     * @param dxB0 initial adjoint derivative states
     * @param xQB0 initial adjoint quadratures
     */

    void setupB(int *which, realtype tf, Model *model, const AmiVector &xB0,
                const AmiVector &dxB0, const AmiVector &xQB0) const;

    /**
     * @brief Initializes the ami memory for quadrature computation
     * @param t0 initial timepoint
     * @param model pointer to the model instance
     * @param x0 initial states
     * @param dx0 initial derivative states
     * @param xB0 initial adjoint states
     * @param dxB0 initial derivative adjoint states
     * @param xQ0 initial quadrature vector
     */

    void setupSteadystate(const realtype t0, Model *model, const AmiVector &x0,
                          const AmiVector &dx0, const AmiVector &xB0,
                          const AmiVector &dxB0, const AmiVector &xQ0) const;

    /**
     * @brief Reinitializes state and respective sensitivities (if necessary) according
     * to changes in fixedParameters
     *
     * @param model pointer to the model instance
     */
    void updateAndReinitStatesAndSensitivities(Model *model);

    /**
     * getRootInfo extracts information which event occurred
     *
     * @param rootsfound array with flags indicating whether the respective
     * event occurred
     */
    virtual void getRootInfo(int *rootsfound) const = 0;

    /**
     * @brief Calculates consistent initial conditions, assumes initial
     * states to be correct (DAE only)
     *
     * @param tout1 next timepoint to be computed (sets timescale)
     */
    virtual void calcIC(realtype tout1) const = 0;

    /**
     * @brief Calculates consistent initial conditions for the backwards
     * problem, assumes initial states to be correct (DAE only)
     *
     * @param which identifier of the backwards problem
     * @param tout1 next timepoint to be computed (sets timescale)
     */
    virtual void calcICB(int which, realtype tout1) const = 0;

    /**
     * @brief Solves the backward problem until a predefined timepoint
     * (adjoint only)
     *
     * @param tBout timepoint until which simulation should be performed
     * @param itaskB task identifier, can be CV_NORMAL or CV_ONE_STEP
     */
    virtual void solveB(realtype tBout, int itaskB) const = 0;

    /**
     * @brief Disable rootfinding
     */
    virtual void turnOffRootFinding() const = 0;

    /**
     * @brief Return current sensitivity method
     * @return method enum
     */
    SensitivityMethod getSensitivityMethod() const;

    /**
     * @brief Set sensitivity method
     * @param sensi_meth
     */
    void setSensitivityMethod(SensitivityMethod sensi_meth);

    /**
     * @brief Return current sensitivity method during preequilibration
     * @return method enum
     */
    SensitivityMethod getSensitivityMethodPreequilibration() const;

    /**
     * @brief Set sensitivity method for preequilibration
     * @param sensi_meth_preeq
     */
    void setSensitivityMethodPreequilibration(SensitivityMethod sensi_meth_preeq);

    /**
     * @brief Disable forward sensitivity integration (used in steady state sim)
     */
    void switchForwardSensisOff() const;

    /**
     * @brief Get maximum number of allowed Newton steps for steady state
     * computation
     * @return
     */
    int getNewtonMaxSteps() const;

    /**
     * @brief Set maximum number of allowed Newton steps for steady state
     * computation
     * @param newton_maxsteps
     */
    void setNewtonMaxSteps(int newton_maxsteps);

    /**
     * @brief Get if model preequilibration is enabled
     * @return
     */
    bool getPreequilibration() const;

    /**
     * @brief Enable/disable model preequilibration
     * @param require_preequilibration
     */
    void setPreequilibration(bool require_preequilibration);

    /**
     * @brief Get maximum number of allowed linear steps per Newton step for
     * steady state computation
     * @return
     */
    int getNewtonMaxLinearSteps() const;

    /**
     * @brief Set maximum number of allowed linear steps per Newton step for
     * steady state computation
     * @param newton_maxlinsteps
     */
    void setNewtonMaxLinearSteps(int newton_maxlinsteps);

    /**
     * @brief Get a state of the damping factor used in the Newton solver
     * @return
     */
    NewtonDampingFactorMode getNewtonDampingFactorMode() const;

    /**
     * @brief Turn on/off a damping factor in the Newton method
     * @param dampingFactorMode
     */
    void setNewtonDampingFactorMode(NewtonDampingFactorMode dampingFactorMode);

    /**
     * @brief Get a lower bound of the damping factor used in the Newton solver
     * @return
     */
    double getNewtonDampingFactorLowerBound() const;

    /**
     * @brief Set a lower bound of the damping factor in the Newton solver
     * @param dampingFactorLowerBound
     */
    void setNewtonDampingFactorLowerBound(double dampingFactorLowerBound);

    /**
     * @brief Get sensitivity order
     * @return sensitivity order
     */
    SensitivityOrder getSensitivityOrder() const;

    /**
     * @brief Set the sensitivity order
     * @param sensi sensitivity order
     */
    void setSensitivityOrder(SensitivityOrder sensi);

    /**
     * @brief Get the relative tolerances for the forward problem
     *
     * Same tolerance is used for the backward problem if not specified
     * differently via setRelativeToleranceASA.
     *
     * @return relative tolerances
     */
    double getRelativeTolerance() const;

    /**
     * @brief Sets the relative tolerances for the forward problem
     *
     * Same tolerance is used for the backward problem if not specified
     * differently via setRelativeToleranceASA.
     *
     * @param rtol relative tolerance (non-negative number)
     */
    void setRelativeTolerance(double rtol);

    /**
     * @brief Get the absolute tolerances for the forward problem
     *
     * Same tolerance is used for the backward problem if not specified
     * differently via setAbsoluteToleranceASA.
     *
     * @return absolute tolerances
     */
    double getAbsoluteTolerance() const;

    /**
     * @brief Sets the absolute tolerances for the forward problem
     *
     * Same tolerance is used for the backward problem if not specified
     * differently via setAbsoluteToleranceASA.
     *
     * @param atol absolute tolerance (non-negative number)
     */
    void setAbsoluteTolerance(double atol);

    /**
     * @brief Returns the relative tolerances for the forward sensitivity
     * problem
     * @return relative tolerances
     */
    double getRelativeToleranceFSA() const;

    /**
     * @brief Sets the relative tolerances for the forward sensitivity problem
     * @param rtol relative tolerance (non-negative number)
     */
    void setRelativeToleranceFSA(double rtol);

    /**
     * @brief Returns the absolute tolerances for the forward sensitivity
     * problem
     * @return absolute tolerances
     */
    double getAbsoluteToleranceFSA() const;

    /**
     * @brief Sets the absolute tolerances for the forward sensitivity problem
     * @param atol absolute tolerance (non-negative number)
     */
    void setAbsoluteToleranceFSA(double atol);

    /**
     * @brief Returns the relative tolerances for the adjoint sensitivity
     * problem
     * @return relative tolerances
     */
    double getRelativeToleranceB() const;

    /**
     * @brief Sets the relative tolerances for the adjoint sensitivity problem
     * @param rtol relative tolerance (non-negative number)
     */
    void setRelativeToleranceB(double rtol);

    /**
     * @brief Returns the absolute tolerances for the backward problem for
     * adjoint sensitivity analysis
     * @return absolute tolerances
     */
    double getAbsoluteToleranceB() const;

    /**
     * @brief Sets the absolute tolerances for the backward problem for
     * adjoint sensitivity analysis
     * @param atol absolute tolerance (non-negative number)
     */
    void setAbsoluteToleranceB(double atol);

    /**
     * @brief Returns the relative tolerance for the quadrature problem
     * @return relative tolerance
     */
    double getRelativeToleranceQuadratures() const;

    /**
     * @brief sets the relative tolerance for the quadrature problem
     * @param rtol relative tolerance (non-negative number)
     */
    void setRelativeToleranceQuadratures(double rtol);

    /**
     * @brief returns the absolute tolerance for the quadrature problem
     * @return absolute tolerance
     */
    double getAbsoluteToleranceQuadratures() const;

    /**
     * @brief sets the absolute tolerance for the quadrature problem
     * @param atol absolute tolerance (non-negative number)
     */
    void setAbsoluteToleranceQuadratures(double atol);

    /**
     * @brief returns the relative tolerance for the steady state problem
     * @return relative tolerance
     */
    double getRelativeToleranceSteadyState() const;

    /**
     * @brief sets the relative tolerance for the steady state problem
     * @param rtol relative tolerance (non-negative number)
     */
    void setRelativeToleranceSteadyState(double rtol);

    /**
     * @brief returns the absolute tolerance for the steady state problem
     * @return absolute tolerance
     */
    double getAbsoluteToleranceSteadyState() const;

    /**
     * @brief sets the absolute tolerance for the steady state problem
     * @param atol absolute tolerance (non-negative number)
     */
    void setAbsoluteToleranceSteadyState(double atol);

    /**
     * @brief returns the relative tolerance for the sensitivities of the
     * steady state problem
     * @return relative tolerance
     */
    double getRelativeToleranceSteadyStateSensi() const;

    /**
     * @brief sets the relative tolerance for the sensitivities of the
     * steady state problem
     * @param rtol relative tolerance (non-negative number)
     */
    void setRelativeToleranceSteadyStateSensi(double rtol);

    /**
     * @brief returns the absolute tolerance for the sensitivities of the
     * steady state problem
     * @return absolute tolerance
     */
    double getAbsoluteToleranceSteadyStateSensi() const;

    /**
     * @brief sets the absolute tolerance for the sensitivities of the
     * steady state problem
     * @param atol absolute tolerance (non-negative number)
     */
    void setAbsoluteToleranceSteadyStateSensi(double atol);

    /**
     * @brief returns the maximum number of solver steps for the forward
     * problem
     * @return maximum number of solver steps
     */
    long int getMaxSteps() const;

    /**
     * @brief sets the maximum number of solver steps for the forward problem
     * @param maxsteps maximum number of solver steps (positive number)
     */
    void setMaxSteps(long int maxsteps);

    /**
     * @brief returns the maximum number of solver steps for the backward
     * problem
     * @return maximum number of solver steps
     */
    long int getMaxStepsBackwardProblem() const;

    /**
     * @brief sets the maximum number of solver steps for the backward problem
     *
     * @param maxsteps maximum number of solver steps (non-negative number)
     *
     * @note default behaviour (100 times the value for the forward problem) can
     * be restored by passing maxsteps=0
     */
    void setMaxStepsBackwardProblem(long int maxsteps);

    /**
     * @brief returns the linear system multistep method
     * @return linear system multistep method
     */
    LinearMultistepMethod getLinearMultistepMethod() const;

    /**
     * @brief sets the linear system multistep method
     * @param lmm linear system multistep method
     */
    void setLinearMultistepMethod(LinearMultistepMethod lmm);

    /**
     * @brief returns the nonlinear system solution method
     * @return
     */
    NonlinearSolverIteration getNonlinearSolverIteration() const;

    /**
     * @brief sets the nonlinear system solution method
     * @param iter nonlinear system solution method
     */
    void setNonlinearSolverIteration(NonlinearSolverIteration iter);

    /**
     * @brief getInterpolationType
     * @return
     */
    InterpolationType getInterpolationType() const;

    /**
     * @brief sets the interpolation of the forward solution that is used for
     * the backwards problem
     * @param interpType interpolation type
     */
    void setInterpolationType(InterpolationType interpType);

    /**
     * @brief Gets KLU / SuperLUMT state ordering mode
     *
     * @return State-ordering as integer according to
     * SUNLinSolKLU::StateOrdering or SUNLinSolSuperLUMT::StateOrdering (which
     * differ).
     */
    int getStateOrdering() const;

    /**
     * @brief Sets KLU / SuperLUMT state ordering mode
     *
     * This only applies when linsol is set to LinearSolver::KLU or
     * LinearSolver::SuperLUMT. Mind the difference between
     * SUNLinSolKLU::StateOrdering and SUNLinSolSuperLUMT::StateOrdering.
     * @param ordering state ordering
     */
    void setStateOrdering(int ordering);

    /**
     * @brief returns stability limit detection mode
     * @return stldet can be false (deactivated) or true (activated)
     */
    bool getStabilityLimitFlag() const;

    /**
     * @brief set stability limit detection mode
     * @param stldet can be false (deactivated) or true (activated)
     */
    void setStabilityLimitFlag(bool stldet);

    /**
     * @brief getLinearSolver
     * @return
     */
    LinearSolver getLinearSolver() const;

    /**
     * @brief setLinearSolver
     * @param linsol
     */
    void setLinearSolver(LinearSolver linsol);

    /**
     * @brief returns the internal sensitivity method
     * @return internal sensitivity method
     */
    InternalSensitivityMethod getInternalSensitivityMethod() const;

    /**
     * @brief sets the internal sensitivity method
     * @param ism internal sensitivity method
     */
    void setInternalSensitivityMethod(InternalSensitivityMethod ism);

    /**
     * @brief returns the ReturnData reporting mode
     * @return ReturnData reporting mode
     */
    RDataReporting getReturnDataReportingMode() const;

    /**
     * @brief sets the ReturnData reporting mode
     * @param rdrm ReturnData reporting mode
     */
    void setReturnDataReportingMode(RDataReporting rdrm);

    /**
     * @brief write solution from forward simulation
     * @param t time
     * @param x state
     * @param dx derivative state
     * @param sx state sensitivity
     * @param xQ quadrature
     */
    void writeSolution(realtype *t, AmiVector &x, AmiVector &dx,
                       AmiVectorArray &sx, AmiVector &xQ) const;

    /**
     * @brief write solution from forward simulation
     * @param t time
     * @param xB adjoint state
     * @param dxB adjoint derivative state
     * @param xQB adjoint quadrature
     * @param which index of adjoint problem
     */
    void writeSolutionB(realtype *t, AmiVector &xB, AmiVector &dxB,
                        AmiVector &xQB, int which) const;

    /**
     * @brief Access state solution at time t
     * @param t time
     * @return x or interpolated solution dky
     */
    const AmiVector &getState(realtype t) const;

    /**
     * @brief Access derivative state solution at time t
     * @param t time
     * @return dx or interpolated solution dky
     */
    const AmiVector &getDerivativeState(realtype t) const;

    /**
     * @brief Access state sensitivity solution at time t
     * @param t time
     * @return (interpolated) solution sx
     */
    const AmiVectorArray &getStateSensitivity(realtype t) const;

    /**
     * @brief Access adjoint solution at time t
     * @param which adjoint problem index
     * @param t time
     * @return (interpolated) solution xB
     */
    const AmiVector &getAdjointState(int which, realtype t) const;

    /**
     * @brief Access adjoint derivative solution at time t
     * @param which adjoint problem index
     * @param t time
     * @return (interpolated) solution dxB
     */
    const AmiVector &getAdjointDerivativeState(int which, realtype t) const;

    /**
     * @brief Access adjoint quadrature solution at time t
     * @param which adjoint problem index
     * @param t time
     * @return (interpolated) solution xQB
     */
    const AmiVector &getAdjointQuadrature(int which, realtype t) const;

    /**
     * @brief Access quadrature solution at time t
     * @param t time
     * @return (interpolated) solution xQ
     */
    const AmiVector &getQuadrature(realtype t) const;

    /**
     * @brief Reinitializes the states in the solver after an event occurrence
     *
     * @param t0 reinitialization timepoint
     * @param yy0 initial state variables
     * @param yp0 initial derivative state variables (DAE only)
     */
    virtual void reInit(realtype t0, const AmiVector &yy0,
                        const AmiVector &yp0) const = 0;

    /**
     * @brief Reinitializes the state sensitivities in the solver after an
     * event occurrence
     *
     * @param yyS0 new state sensitivity
     * @param ypS0 new derivative state sensitivities (DAE only)
     */
    virtual void sensReInit(const AmiVectorArray &yyS0,
                            const AmiVectorArray &ypS0) const = 0;

    /**
     * @brief Switches off computation of  state sensitivities without
     * deallocating the memory for sensitivities
     */
    virtual void sensToggleOff() const = 0;

    /**
     * @brief Reinitializes the adjoint states after an event occurrence
     *
     * @param which identifier of the backwards problem
     * @param tB0 reinitialization timepoint
     * @param yyB0 new adjoint state
     * @param ypB0 new adjoint derivative state
     */
    virtual void reInitB(int which, realtype tB0, const AmiVector &yyB0,
                         const AmiVector &ypB0) const = 0;

    /**
     * @brief Reinitialize the adjoint states after an event occurrence
     *
     * @param which identifier of the backwards problem
     * @param yQB0 new adjoint quadrature state
     */
    virtual void quadReInitB(int which, const AmiVector &yQB0) const = 0;

    /**
     * @brief current solver timepoint
     * @return t
     */
    realtype gett() const;

    /**
     * @brief Reads out the CPU time needed for forward solve
     * @return cpu_time
     */
    realtype getCpuTime() const;

    /**
     * @brief Reads out the CPU time needed for backward solve
     * @return cpu_timeB
     */
    realtype getCpuTimeB() const;

    /**
     * @brief number of states with which the solver was initialized
     * @return x.getLength()
     */
    int nx() const;

    /**
     * @brief number of parameters with which the solver was initialized
     * @return sx.getLength()
     */
    int nplist() const;

    /**
     * @brief number of quadratures with which the solver was initialized
     * @return xQB.getLength()
     */
    int nquad() const;

    /**
     * @brief check if FSA is being computed
     * @return flag
     */
    bool computingFSA() const {
        return getSensitivityOrder() >= SensitivityOrder::first &&
        getSensitivityMethod() == SensitivityMethod::forward && nplist() > 0;
    }

    /**
     * @brief check if ASA is being computed
     * @return flag
     */
    bool computingASA() const {
        return getSensitivityOrder() >= SensitivityOrder::first &&
        getSensitivityMethod() == SensitivityMethod::adjoint && nplist() > 0;
    }

    /**
     * @brief Resets vectors containing diagnosis information
     */
    void resetDiagnosis() const;

    /**
     * @brief Stores diagnosis information from solver memory block for forward problem
     */
    void storeDiagnosis() const;

    /**
     * @brief Stores diagnosis information from solver memory block for backward problem
     *
     * @param which identifier of the backwards problem
     */
    void storeDiagnosisB(int which) const;

    /**
     * @brief Accessor ns
     * @return ns
     */
    std::vector<int> const& getNumSteps() const {
        return ns_;
    }

    /**
     * @brief Accessor nsB
     * @return nsB
     */
    std::vector<int> const& getNumStepsB() const {
        return nsB_;
    }

    /**
     * @brief Accessor nrhs
     * @return nrhs
     */
    std::vector<int> const& getNumRhsEvals() const {
        return nrhs_;
    }

    /**
     * @brief Accessor nrhsB
     * @return nrhsB
     */
    std::vector<int> const& getNumRhsEvalsB() const {
        return nrhsB_;
    }

    /**
     * @brief Accessor netf
     * @return netf
     */
    std::vector<int> const& getNumErrTestFails() const {
        return netf_;
    }

    /**
     * @brief Accessor netfB
     * @return netfB
     */
    std::vector<int> const& getNumErrTestFailsB() const {
        return netfB_;
    }

    /**
     * @brief Accessor nnlscf
     * @return nnlscf
     */
    std::vector<int> const& getNumNonlinSolvConvFails() const {
        return nnlscf_;
    }

    /**
     * @brief Accessor nnlscfB
     * @return nnlscfB
     */
    std::vector<int> const& getNumNonlinSolvConvFailsB() const {
        return nnlscfB_;
    }

    /**
     * @brief Accessor order
     * @return order
     */
    std::vector<int> const& getLastOrder() const {
        return order_;
    }

    /**
     * @brief Serialize Solver (see boost::serialization::serialize)
     * @param ar Archive to serialize to
     * @param s Data to serialize
     * @param version Version number
     */
    template <class Archive>
    friend void boost::serialization::serialize(Archive &ar, Solver &s,
                                                unsigned int version);

    /**
     * @brief Check equality of data members excluding solver memory
     * @param a
     * @param b
     * @return
     */
    friend bool operator==(const Solver &a, const Solver &b);

    /** AMICI context */
    AmiciApplication *app = &defaultContext;

  protected:
    /**
     * @brief Sets a timepoint at which the simulation will be stopped
     *
     * @param tstop timepoint until which simulation should be performed
     */
    virtual void setStopTime(realtype tstop) const = 0;

    /**
     * @brief Solves the forward problem until a predefined timepoint
     *
     * @param tout timepoint until which simulation should be performed
     * @param itask task identifier, can be CV_NORMAL or CV_ONE_STEP
     * @return status flag indicating success of execution
     */
    virtual int solve(realtype tout, int itask) const = 0;

    /**
     * @brief Solves the forward problem until a predefined timepoint
     * (adjoint only)
     *
     * @param tout timepoint until which simulation should be performed
     * @param itask task identifier, can be CV_NORMAL or CV_ONE_STEP
     * @param ncheckPtr pointer to a number that counts the internal
     * checkpoints
     * @return status flag indicating success of execution
     */
    virtual int solveF(realtype tout, int itask, int *ncheckPtr) const = 0;

    /**
     * @brief reInitPostProcessF postprocessing of the solver memory after a
     * discontinuity in the forward problem
     * @param tnext next timepoint (defines integration direction)
     */
    virtual void reInitPostProcessF(realtype tnext) const = 0;

    /**
     * @brief reInitPostProcessB postprocessing of the solver memory after a
     * discontinuity in the backward problem
     * @param tnext next timepoint (defines integration direction)
     */
    virtual void reInitPostProcessB(realtype tnext) const = 0;

    /**
     * @brief extracts the state sensitivity at the current timepoint from
     * solver memory and writes it to the sx member variable
     */
    virtual void getSens() const = 0;

    /**
     * @brief extracts the adjoint state at the current timepoint from
     * solver memory and writes it to the xB member variable
     * @param which index of the backwards problem
     */
    virtual void getB(int which) const = 0;

    /**
     * @brief extracts the adjoint quadrature state at the current timepoint
     * from solver memory and writes it to the xQB member variable
     * @param which index of the backwards problem
     */
    virtual void getQuadB(int which) const = 0;

    /**
     * @brief extracts the quadrature at the current timepoint from solver
     * memory and writes it to the xQ member variable
     *
     * @param t timepoint for quadrature extraction
     */
    virtual void getQuad(realtype &t) const = 0;

    /**
     * @brief Initializes the states at the specified initial timepoint
     *
     * @param t0 initial timepoint
     * @param x0 initial states
     * @param dx0 initial derivative states
     */
    virtual void init(realtype t0, const AmiVector &x0,
                      const AmiVector &dx0) const = 0;

    /**
     * @brief Initializes the states at the specified initial timepoint
     *
     * @param t0 initial timepoint
     * @param x0 initial states
     * @param dx0 initial derivative states
     */
    virtual void initSteadystate(realtype t0, const AmiVector &x0,
                                 const AmiVector &dx0) const = 0;

    /**
     * @brief Initializes the forward sensitivities
     * @param sx0 initial states sensitivities
     * @param sdx0 initial derivative states sensitivities
     */
    virtual void sensInit1(const AmiVectorArray &sx0,
                           const AmiVectorArray &sdx0) const = 0;

    /**
     * @brief Initialize the adjoint states at the specified final timepoint
     *
     * @param which identifier of the backwards problem
     * @param tf final timepoint
     * @param xB0 initial adjoint state
     * @param dxB0 initial adjoint derivative state
     */
    virtual void binit(int which, realtype tf, const AmiVector &xB0,
                       const AmiVector &dxB0) const = 0;

    /**
     * @brief Initialize the quadrature states at the specified final timepoint
     *
     * @param which identifier of the backwards problem
     * @param xQB0 initial adjoint quadrature state
     */
    virtual void qbinit(int which, const AmiVector &xQB0) const = 0;

    /**
     * @brief Initializes the rootfinding for events
     *
     * @param ne number of different events
     */
    virtual void rootInit(int ne) const = 0;

    /**
     * @brief Initalize non-linear solver for sensitivities
     * @param model Model instance
     */
    void initializeNonLinearSolverSens(const Model *model) const;

    /**
     * @brief Set the dense Jacobian function
     */
    virtual void setDenseJacFn() const = 0;

    /**
     * @brief sets the sparse Jacobian function
     */
    virtual void setSparseJacFn() const = 0;

    /**
     * @brief sets the banded Jacobian function
     */
    virtual void setBandJacFn() const = 0;

    /**
     * @brief sets the Jacobian vector multiplication function
     */
    virtual void setJacTimesVecFn() const = 0;

    /**
     * @brief sets the dense Jacobian function
     *
     * @param which identifier of the backwards problem
     */
    virtual void setDenseJacFnB(int which) const = 0;

    /**
     * @brief sets the sparse Jacobian function
     *
     * @param which identifier of the backwards problem
     */
    virtual void setSparseJacFnB(int which) const = 0;

    /**
     * @brief sets the banded Jacobian function
     *
     * @param which identifier of the backwards problem
     */
    virtual void setBandJacFnB(int which) const = 0;

    /**
     * @brief sets the Jacobian vector multiplication function
     *
     * @param which identifier of the backwards problem
     */
    virtual void setJacTimesVecFnB(int which) const = 0;

    /**
     * @brief sets the sparse Jacobian function for backward steady state case
     */
    virtual void setSparseJacFn_ss() const = 0;

    /**
     * @brief Create specifies solver method and initializes solver memory for
     * the forward problem
     */
    virtual void allocateSolver() const = 0;

    /**
     * @brief sets scalar relative and absolute tolerances for the forward
     * problem
     *
     * @param rtol relative tolerances
     * @param atol absolute tolerances
     */
    virtual void setSStolerances(double rtol, double atol) const = 0;

    /**
     * @brief activates sets scalar relative and absolute tolerances for the
     * sensitivity variables
     *
     * @param rtol relative tolerances
     * @param atol array of absolute tolerances for every sensitivity variable
     */
    virtual void setSensSStolerances(double rtol, const double *atol) const = 0;

    /**
     * SetSensErrCon specifies whether error control is also enforced for
     * sensitivities for the forward problem
     *
     * @param error_corr activation flag
     */
    virtual void setSensErrCon(bool error_corr) const = 0;

    /**
     * @brief Specifies whether error control is also enforced for the
     * backward quadrature problem
     *
     * @param which identifier of the backwards problem
     * @param flag activation flag
     */
    virtual void setQuadErrConB(int which, bool flag) const = 0;

    /**
     * @brief Specifies whether error control is also enforced for the
     * forward quadrature problem
     *
     * @param flag activation flag
     */
    virtual void setQuadErrCon(bool flag) const = 0;

    /**
     * @brief Attaches the error handler function (errMsgIdAndTxt)
     * to the solver
     *
     */
    virtual void setErrHandlerFn() const = 0;

    /**
     * @brief Attaches the user data instance (here this is a Model) to the
     * forward problem
     *
     * @param model Model instance
     */
    virtual void setUserData(Model *model) const = 0;

    /**
     * @brief attaches the user data instance (here this is a Model) to the
     * backward problem
     *
     * @param which identifier of the backwards problem
     * @param model Model instance
     */
    virtual void setUserDataB(int which, Model *model) const = 0;

    /**
     * @brief specifies the maximum number of steps for the forward
     * problem
     *
     * @param mxsteps number of steps
     * @note in contrast to the SUNDIALS method, this sets the overall maximum, not the maximum between output times.
     */
    virtual void setMaxNumSteps(long int mxsteps) const = 0;

    /**
     * @brief specifies the maximum number of steps for the forward
     * problem
     *
     * @param which identifier of the backwards problem
     * @param mxstepsB number of steps
     * @note in contrast to the SUNDIALS method, this sets the overall maximum, not the maximum between output times.
     */
    virtual void setMaxNumStepsB(int which, long int mxstepsB) const = 0;

    /**
     * @brief activates stability limit detection for the forward
     * problem
     *
     * @param stldet flag for stability limit detection (TRUE or FALSE)
     *
     */
    virtual void setStabLimDet(int stldet) const = 0;

    /**
     * @brief activates stability limit detection for the backward
     * problem
     *
     * @param which identifier of the backwards problem
     * @param stldet flag for stability limit detection (TRUE or FALSE)
     *
     */
    virtual void setStabLimDetB(int which, int stldet) const = 0;

    /**
     * @brief specify algebraic/differential components (DAE only)
     *
     * @param model model specification
     */
    virtual void setId(const Model *model) const = 0;

    /**
     * @brief deactivates error control for algebraic components (DAE only)
     *
     * @param flag deactivation flag
     */
    virtual void setSuppressAlg(bool flag) const = 0;

    /**
     * @brief specifies the scaling and indexes for sensitivity
     * computation
     *
     * @param p parameters
     * @param pbar parameter scaling constants
     * @param plist parameter index list
     */
    virtual void setSensParams(const realtype *p, const realtype *pbar,
                               const int *plist) const = 0;

    /**
     * @brief interpolates the (derivative of the) solution at the requested
     * timepoint
     *
     * @param t timepoint
     * @param k derivative order
     */
    virtual void getDky(realtype t, int k) const = 0;

    /**
     * @brief interpolates the (derivative of the) solution at the requested
     * timepoint
     *
     * @param t timepoint
     * @param k derivative order
     * @param which index of backward problem
     */
    virtual void getDkyB(realtype t, int k, int which) const = 0;

    /**
     * @brief interpolates the (derivative of the) solution at the requested
     * timepoint
     *
     * @param t timepoint
     * @param k derivative order
     */
    virtual void getSensDky(realtype t, int k) const = 0;

    /**
     * @brief interpolates the (derivative of the) solution at the requested
     * timepoint
     *
     * @param t timepoint
     * @param k derivative order
     * @param which index of backward problem
     */
    virtual void getQuadDkyB(realtype t, int k, int which) const = 0;

    /**
     * @brief interpolates the (derivative of the) solution at the requested
     * timepoint
     *
     * @param t timepoint
     * @param k derivative order
     */
    virtual void getQuadDky(realtype t, int k) const = 0;

    /**
     * @brief initializes the adjoint problem
     */
    virtual void adjInit() const = 0;

    /**
     * @brief initializes the quadratures
     * @param xQ0 vector with initial values for xQ
     */
    virtual void quadInit(const AmiVector &xQ0) const = 0;

    /**
     * @brief Specifies solver method and initializes solver memory for the
     * backward problem
     *
     * @param which identifier of the backwards problem
     */
    virtual void allocateSolverB(int *which) const = 0;

    /**
     * @brief sets relative and absolute tolerances for the backward
     * problem
     *
     * @param which identifier of the backwards problem
     * @param relTolB relative tolerances
     * @param absTolB absolute tolerances
     */
    virtual void setSStolerancesB(int which, realtype relTolB,
                                  realtype absTolB) const = 0;

    /**
     * @brief sets relative and absolute tolerances for the quadrature
     * backward problem
     *
     * @param which identifier of the backwards problem
     * @param reltolQB relative tolerances
     * @param abstolQB absolute tolerances
     */
    virtual void quadSStolerancesB(int which, realtype reltolQB,
                                   realtype abstolQB) const = 0;

    /**
     * @brief sets relative and absolute tolerances for the quadrature problem
     *
     * @param reltolQB relative tolerances
     * @param abstolQB absolute tolerances
     */
    virtual void quadSStolerances(realtype reltolQB,
                                  realtype abstolQB) const = 0;

    /**
     * @brief reports the number of solver steps
     *
     * @param ami_mem pointer to the solver memory instance (can be from
     * forward or backward problem)
     * @param numsteps output array
     */
    virtual void getNumSteps(const void *ami_mem, long int *numsteps) const = 0;

    /**
     * @brief reports the number of right hand evaluations
     *
     * @param ami_mem pointer to the solver memory instance (can be from
     * forward or backward problem)
     * @param numrhsevals output array
     */
    virtual void getNumRhsEvals(const void *ami_mem,
                                long int *numrhsevals) const = 0;

    /**
     * @brief reports the number of local error test failures
     *
     * @param ami_mem pointer to the solver memory instance (can be from
     * forward or backward problem)
     * @param numerrtestfails output array
     */
    virtual void getNumErrTestFails(const void *ami_mem,
                                    long int *numerrtestfails) const = 0;

    /**
     * @brief reports the number of nonlinear convergence failures
     *
     * @param ami_mem pointer to the solver memory instance (can be from
     * forward or backward problem)
     * @param numnonlinsolvconvfails output array
     */
    virtual void
    getNumNonlinSolvConvFails(const void *ami_mem,
                              long int *numnonlinsolvconvfails) const = 0;

    /**
     * @brief Reports the order of the integration method during the
     * last internal step
     *
     * @param ami_mem pointer to the solver memory instance (can be from
     * forward or backward problem)
     * @param order output array
     */
    virtual void getLastOrder(const void *ami_mem, int *order) const = 0;

    /**
     * @brief Initializes and sets the linear solver for the forward problem
     *
     * @param model pointer to the model object
     */
    void initializeLinearSolver(const Model *model) const;

    /**
     * @brief Sets the non-linear solver
     */
    void initializeNonLinearSolver() const;

    /**
     * @brief Sets the linear solver for the forward problem
     */
    virtual void setLinearSolver() const = 0;

    /**
     * @brief Sets the linear solver for the backward problem
     * @param which index of the backward problem
     */
    virtual void setLinearSolverB(int which) const = 0;

    /**
     * @brief Set the non-linear solver for the forward problem
     */
    virtual void setNonLinearSolver() const = 0;

    /**
     * @brief Set the non-linear solver for the backward problem
     * @param which index of the backward problem
     */
    virtual void setNonLinearSolverB(int which) const = 0;

    /**
     * @brief Set the non-linear solver for sensitivities
     */
    virtual void setNonLinearSolverSens() const = 0;

    /**
     * @brief Initializes the linear solver for the backward problem
     *
     * @param model pointer to the model object
     * @param which index of the backward problem
     */

    void initializeLinearSolverB(const Model *model, int which) const;

    /**
     * @brief Initializes the non-linear solver for the backward problem
     * @param which index of the backward problem
     */
    void initializeNonLinearSolverB(int which) const;

    /**
     * Accessor function to the model stored in the user data
     *
     * @return user data model
     */
    virtual const Model *getModel() const = 0;

    /**
     * @brief checks whether memory for the forward problem has been allocated
     *
     * @return proxy for solverMemory->(cv|ida)_MallocDone
     */
    bool getInitDone() const;

    /**
     * @brief checks whether memory for forward sensitivities has been allocated
     *
     * @return proxy for solverMemory->(cv|ida)_SensMallocDone
     */
    bool getSensInitDone() const;

    /**
     * @brief checks whether memory for forward interpolation has been allocated
     *
     * @return proxy for solverMemory->(cv|ida)_adjMallocDone
     */
    bool getAdjInitDone() const;

    /**
     * @brief checks whether memory for the backward problem has been allocated
     * @param which adjoint problem index
     * @return proxy for solverMemoryB->(cv|ida)_MallocDone
     */
    bool getInitDoneB(int which) const;

    /**
     * @brief checks whether memory for backward quadratures has been allocated
     * @param which adjoint problem index
     * @return proxy for solverMemoryB->(cv|ida)_QuadMallocDone
     */
    bool getQuadInitDoneB(int which) const;

    /**
     * @brief checks whether memory for quadratures has been allocated
     * @return proxy for solverMemory->(cv|ida)_QuadMallocDone
     */
    bool getQuadInitDone() const;

    /**
     * @brief attaches a diagonal linear solver to the forward problem
     */
    virtual void diag() const = 0;

    /**
     * @brief attaches a diagonal linear solver to the backward problem
     *
     * @param which identifier of the backwards problem
     */
    virtual void diagB(int which) const = 0;

    /**
     * @brief resets solverMemory and solverMemoryB
     * @param nx new number of state variables
     * @param nplist new number of sensitivity parameters
     * @param nquad new number of quadratures (only differs from nplist for
     * higher order sensitivity computation)
     */
    void resetMutableMemory(int nx, int nplist, int nquad) const;

    /**
     * @brief Retrieves the solver memory instance for the backward problem
     *
     * @param which identifier of the backwards problem
     * @param ami_mem pointer to the forward solver memory instance
     * @return A (void *) pointer to the CVODES memory allocated for the
     * backward problem.
     */
    virtual void *getAdjBmem(void *ami_mem, int which) const = 0;

    /**
     * @brief updates solver tolerances according to the currently specified
     * member variables
     */
    void applyTolerances() const;

    /**
     * @brief updates FSA solver tolerances according to the currently
     * specified member variables
     */
    void applyTolerancesFSA() const;

    /**
     * @brief updates ASA solver tolerances according to the currently
     * specified member variables
     *
     * @param which identifier of the backwards problem
     */
    void applyTolerancesASA(int which) const;

    /**
     * @brief updates ASA quadrature solver tolerances according to the
     * currently specified member variables
     *
     * @param which identifier of the backwards problem
     */
    void applyQuadTolerancesASA(int which) const;

    /**
     * @brief updates quadrature solver tolerances according to the
     * currently specified member variables
     */
    void applyQuadTolerances() const;

    /**
     * @brief updates all sensitivity solver tolerances according to the
     * currently specified member variables
     */
    void applySensitivityTolerances() const;

    /** pointer to solver memory block */
    mutable std::unique_ptr<void, std::function<void(void *)>> solver_memory_;

    /** pointer to solver memory block */
    mutable std::vector<std::unique_ptr<void, std::function<void(void *)>>>
        solver_memory_B_;

    /** internal sensitivity method flag used to select the sensitivity solution
     * method. Only applies for Forward Sensitivities. */
    InternalSensitivityMethod ism_ {InternalSensitivityMethod::simultaneous};

    /** specifies the linear multistep method.
     */
    LinearMultistepMethod lmm_ {LinearMultistepMethod::BDF};

    /**
     * specifies the type of nonlinear solver iteration
     */
    NonlinearSolverIteration iter_ {NonlinearSolverIteration::newton};

    /** interpolation type for the forward problem solution which
     * is then used for the backwards problem.
     */
    InterpolationType interp_type_ {InterpolationType::hermite};

    /** maximum number of allowed integration steps */
    long int maxsteps_ {10000};

    /** linear solver for the forward problem */
    mutable std::unique_ptr<SUNLinSolWrapper> linear_solver_;

    /** linear solver for the backward problem */
    mutable std::unique_ptr<SUNLinSolWrapper> linear_solver_B_;

    /** non-linear solver for the forward problem */
    mutable std::unique_ptr<SUNNonLinSolWrapper> non_linear_solver_;

    /** non-linear solver for the backward problem */
    mutable std::unique_ptr<SUNNonLinSolWrapper> non_linear_solver_B_;

    /** non-linear solver for the sensitivities */
    mutable std::unique_ptr<SUNNonLinSolWrapper> non_linear_solver_sens_;

    /** flag indicating whether the forward solver has been called */
    mutable bool solver_was_called_F_ {false};

    /** flag indicating whether the backward solver has been called */
    mutable bool solver_was_called_B_ {false};

    /**
     * @brief sets that memory for the forward problem has been allocated
     */
    void setInitDone() const;

    /**
     * @brief sets that memory for forward sensitivities has been allocated
     */
    void setSensInitDone() const;

    /**
     * @brief sets that memory for forward sensitivities has not been allocated
     */
    void setSensInitOff() const;

    /**
     * @brief sets that memory for forward interpolation has been allocated
     */
    void setAdjInitDone() const;

    /**
     * @brief sets that memory for the backward problem has been allocated
     * @param which adjoint problem index
     */
    void setInitDoneB(int which) const;

    /**
     * @brief sets that memory for backward quadratures has been allocated
     * @param which adjoint problem index
     */
    void setQuadInitDoneB(int which) const;

    /**
     * @brief sets that memory for quadratures has been allocated
     */
    void setQuadInitDone() const;

    /**
     * @brief Sets sensitivity method (for simulation or preequilibration)
     * @param sensi_meth new value for sensi_meth[_preeq]
     * @param preequilibration flag indicating preequilibration or simulation
     */
    void checkSensitivityMethod(const SensitivityMethod sensi_meth,
                                bool preequilibration) const;

    /** state (dimension: nx_solver) */
    mutable AmiVector x_ {0};

    /** state interface variable (dimension: nx_solver) */
    mutable AmiVector dky_ {0};

    /** state derivative dummy (dimension: nx_solver) */
    mutable AmiVector dx_ {0};

    /** state sensitivities interface variable (dimension: nx_solver x nplist)
     */
    mutable AmiVectorArray sx_ {0, 0};
    /** state derivative sensitivities dummy (dimension: nx_solver x nplist)
     */
    mutable AmiVectorArray sdx_ {0, 0};

    /** adjoint state interface variable (dimension: nx_solver) */
    mutable AmiVector xB_ {0};

    /** adjoint derivative dummy variable (dimension: nx_solver) */
    mutable AmiVector dxB_ {0};

    /** adjoint quadrature interface variable (dimension: nJ x nplist) */
    mutable AmiVector xQB_ {0};

    /** forward quadrature interface variable (dimension: nx_solver) */
    mutable AmiVector xQ_ {0};

    /** integration time of the forward problem */
    mutable realtype t_ {std::nan("")};

    /** flag to force reInitPostProcessF before next call to solve */
    mutable bool force_reinit_postprocess_F_ {false};

    /** flag to force reInitPostProcessB before next call to solveB */
    mutable bool force_reinit_postprocess_B_ {false};

  private:

    /**
     * @brief applies total number of steps for next solver call
     */
    void apply_max_num_steps() const;

    /**
     * @brief applies total number of steps for next backwards solver call
     */
    void apply_max_num_steps_B() const;


    /** method for sensitivity computation */
    SensitivityMethod sensi_meth_ {SensitivityMethod::forward};

    /** method for sensitivity computation in preequilibration */
    SensitivityMethod sensi_meth_preeq_ {SensitivityMethod::forward};

    /** flag controlling stability limit detection */
    booleantype stldet_ {true};

    /** state ordering */
    int ordering_ {static_cast<int>(SUNLinSolKLU::StateOrdering::AMD)};

    /** maximum number of allowed Newton steps for steady state computation */
    long int newton_maxsteps_ {0L};

    /** maximum number of allowed linear steps per Newton step for steady state
     * computation */
    long int newton_maxlinsteps_ {0L};

    /** Damping factor state used int the Newton method */
    NewtonDampingFactorMode newton_damping_factor_mode_
        {NewtonDampingFactorMode::on};

    /** Lower bound of the damping factor. */
    realtype newton_damping_factor_lower_bound_ {1e-8};

    /** Enable model preequilibration */
    bool requires_preequilibration_ {false};

    /** linear solver specification */
    LinearSolver linsol_ {LinearSolver::KLU};

    /** absolute tolerances for integration */
    realtype atol_ {1e-16};

    /** relative tolerances for integration */
    realtype rtol_ {1e-8};

    /** absolute tolerances for forward sensitivity integration */
    realtype atol_fsa_ {NAN};

    /** relative tolerances for forward sensitivity integration */
    realtype rtol_fsa_ {NAN};

    /** absolute tolerances for adjoint sensitivity integration */
    realtype atolB_ {NAN};

    /** relative tolerances for adjoint sensitivity integration */
    realtype rtolB_ {NAN};

    /** absolute tolerances for backward quadratures */
    realtype quad_atol_ {1e-12};

    /** relative tolerances for backward quadratures */
    realtype quad_rtol_ {1e-8};

    /** absolute tolerances for steadystate computation */
    realtype ss_atol_ {NAN};

    /** relative tolerances for steadystate computation */
    realtype ss_rtol_ {NAN};

    /** absolute tolerances for steadystate computation */
    realtype ss_atol_sensi_ {NAN};

    /** relative tolerances for steadystate computation */
    realtype ss_rtol_sensi_ {NAN};

    RDataReporting rdata_mode_ {RDataReporting::full};

    /** CPU time, forward solve */
    mutable realtype cpu_time_ {0.0};

    /** CPU time, backward solve */
    mutable realtype cpu_timeB_ {0.0};

    /** maximum number of allowed integration steps for backward problem */
    long int maxstepsB_ {0L};

    /** flag indicating whether sensitivities are supposed to be computed */
    SensitivityOrder sensi_ {SensitivityOrder::none};

    /** flag indicating whether init was called */
    mutable bool initialized_ {false};

    /** flag indicating whether sensInit1 was called */
    mutable bool sens_initialized_ {false};

    /** flag indicating whether adjInit was called */
    mutable bool adj_initialized_ {false};

    /** flag indicating whether (forward) quadInit was called */
    mutable bool quad_initialized_ {false};

    /** vector of flags indicating whether binit was called for respective
     which */
    mutable std::vector<bool> initializedB_{false};

    /** vector of flags indicating whether qbinit was called for respective
     which */
    mutable std::vector<bool> initializedQB_{false};

    /** number of checkpoints in the forward problem */
    mutable int ncheckPtr_ {0};

    /** number of integration steps forward problem (dimension: nt) */
    mutable std::vector<int> ns_;

    /** number of integration steps backward problem (dimension: nt) */
    mutable std::vector<int> nsB_;

    /** number of right hand side evaluations forward problem (dimension: nt) */
    mutable std::vector<int> nrhs_;

    /** number of right hand side evaluations backward problem (dimension: nt) */
    mutable std::vector<int> nrhsB_;

    /** number of error test failures forward problem (dimension: nt) */
    mutable std::vector<int> netf_;

    /** number of error test failures backward problem (dimension: nt) */
    mutable std::vector<int> netfB_;

    /**
     * number of linear solver convergence failures forward problem (dimension:
     * nt) */
    mutable std::vector<int> nnlscf_;

    /**
     * number of linear solver convergence failures backward problem (dimension:
     * nt) */
    mutable std::vector<int> nnlscfB_;

    /** employed order forward problem (dimension: nt) */
    mutable std::vector<int> order_;
};

bool operator==(const Solver &a, const Solver &b);

/**
 * @brief Extracts diagnosis information from solver memory block and
 * passes them to the specified output function
 *
 * @param error_code error identifier
 * @param module name of the module in which the error occurred
 * @param function name of the function in which the error occurred
 * @param msg error message
 * @param eh_data amici::Solver as void*
 */
void wrapErrHandlerFn(int error_code, const char *module, const char *function,
                      char *msg, void *eh_data);

} // namespace amici

#endif // AMICISOLVER_H
