#ifndef AMICI_SOLVER_H
#define AMICI_SOLVER_H

#include "amici/defines.h"
#include "amici/sundials_linsol_wrapper.h"
#include "amici/symbolic_functions.h"
#include "amici/vector.h"

#include <functional>
#include <memory>

namespace amici {

class ReturnData;
class ForwardProblem;
class BackwardProblem;
class Model;
class Solver;
} // namespace amici

// for serialization friend in Solver
namespace boost {
namespace serialization {
template <class Archive>
void serialize(Archive &ar, amici::Solver &u, const unsigned int version);
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
     * @param tout next timepooint
     * @return status flag
     */
    int run(const realtype tout) const;
    
    /**
     * @brief makes a single step in the simulation
     *
     * @param tout next timepooint
     * @return status flag
     */
    int step(const realtype tout) const;
    
    /**
     * @brief runs a backward simulation until the specified timepoint
     *
     * @param tout next timepooint
     * @return status flag
     */
    void runB(const realtype tout) const;
    
    /**
     * @brief Initialises the ami memory object and applies specified options
     * @param model pointer to the model object
     * @param x0 initial states
     * @param dx0 initial derivative states
     * @param sx0 initial state sensitivities
     * @param sdx0 initial derivative state sensitivities
     */

    void setup(const realtype t0, Model *model, const AmiVector &x0,
               const AmiVector &dx0, const AmiVectorArray &sx0,
               const AmiVectorArray &sdx0) const;

    /**
     * @brief Initialises the AMI memory object for the backwards problem
     * @param which index of the backward problem, will be set by this routine
     * @param tf final timepoint (initial timepoint for the bwd problem)
     * @param model pointer to the model object
     * @param xB0 initial adjoint states
     * @param dxB0 initial adjoint derivative states
     * @param xQB0 initial adjoint quadratures
     */

    void setupB(int *which, const realtype tf, Model *model,
                const AmiVector &xB0, const AmiVector &dxB0,
                const AmiVector &xQB0) const;

    /**
     * @brief Extracts diagnosis information from solver memory block and
     * writes them into the return data object
     *
     * @param it time-point index
     * @param rdata pointer to the return data object
     */
    void getDiagnosis(const int it, ReturnData *rdata) const;

    /**
     * @brief Extracts diagnosis information from solver memory block and
     * writes them into the return data object for the backward problem
     *
     * @param it time-point index
     * @param rdata pointer to the return data object
     * @param which identifier of the backwards problem
     */
    void getDiagnosisB(const int it, ReturnData *rdata, const int which) const;

    /**
     * getRootInfo extracts information which event occured
     *
     * @param rootsfound array with flags indicating whether the respective
     * event occured
     */
    virtual void getRootInfo(int *rootsfound) const = 0;

    /**
     * @brief Calculates consistent initial conditions, assumes initial
     * states to be correct (DAE only)
     *
     * @param tout1 next timepoint to be computed (sets timescale)
     */
    virtual void calcIC(const realtype tout1) const = 0;

    /**
     * @brief Calculates consistent initial conditions for the backwards
     * problem, assumes initial states to be correct (DAE only)
     *
     * @param which identifier of the backwards problem
     * @param tout1 next timepoint to be computed (sets timescale)
     */
    virtual void calcICB(const int which, const realtype tout1) const = 0;

    /**
     * @brief Solves the backward problem until a predefined timepoint
     * (adjoint only)
     *
     * @param tBout timepoint until which simulation should be performed
     * @param itaskB task identifier, can be CV_NORMAL or CV_ONE_STEP
     */
    virtual void solveB(const realtype tBout, const int itaskB) const = 0;

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
    void setNewtonMaxSteps(const int newton_maxsteps);

    /**
     * @brief Get if preequilibration of model via Newton solver is enabled
     * @return
     */
    bool getNewtonPreequilibration() const;

    /**
     * @brief Enable/disable preequilibration of model via Newton solver
     * @param newton_preeq
     */
    void setNewtonPreequilibration(const bool newton_preeq);

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
    void setNewtonMaxLinearSteps(const int newton_maxlinsteps);

    /**
     * @brief Get sensitvity order
     * @return sensitivity order
     */
    SensitivityOrder getSensitivityOrder() const;

    /**
     * @brief Set the sensitvity order
     * @param sensi sensitivity order
     */
    void setSensitivityOrder(const SensitivityOrder sensi);

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
    void setRelativeTolerance(const double rtol);

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
    void setAbsoluteTolerance(const double atol);

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
    void setRelativeToleranceFSA(const double rtol);

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
    void setAbsoluteToleranceFSA(const double atol);

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
    void setRelativeToleranceB(const double rtol);

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
    void setAbsoluteToleranceB(const double atol);

    /**
     * @brief Returns the relative tolerance for the quadrature problem
     * @return relative tolerance
     */
    double getRelativeToleranceQuadratures() const;

    /**
     * @brief sets the relative tolerance for the quadrature problem
     * @param rtol relative tolerance (non-negative number)
     */
    void setRelativeToleranceQuadratures(const double rtol);

    /**
     * @brief returns the absolute tolerance for the quadrature problem
     * @return absolute tolerance
     */
    double getAbsoluteToleranceQuadratures() const;

    /**
     * @brief sets the absolute tolerance for the quadrature problem
     * @param atol absolute tolerance (non-negative number)
     */
    void setAbsoluteToleranceQuadratures(const double atol);

    /**
     * @brief returns the relative tolerance for the steady state problem
     * @return relative tolerance
     */
    double getRelativeToleranceSteadyState() const;

    /**
     * @brief sets the relative tolerance for the steady state problem
     * @param rtol relative tolerance (non-negative number)
     */
    void setRelativeToleranceSteadyState(const double rtol);

    /**
     * @brief returns the absolute tolerance for the steady state problem
     * @return absolute tolerance
     */
    double getAbsoluteToleranceSteadyState() const;

    /**
     * @brief sets the absolute tolerance for the steady state problem
     * @param atol absolute tolerance (non-negative number)
     */
    void setAbsoluteToleranceSteadyState(const double atol);

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
    void setRelativeToleranceSteadyStateSensi(const double rtol);

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
    void setAbsoluteToleranceSteadyStateSensi(const double atol);

    /**
     * @brief returns the maximum number of solver steps for the forward
     * problem
     * @return maximum number of solver steps
     */
    long int getMaxSteps() const;

    /**
     * @brief sets the maximum number of solver steps for the forward problem
     * @param maxsteps maximum number of solver steps (non-negative number)
     */
    void setMaxSteps(const long int maxsteps);

    /**
     * @brief sets the maximum number of solver steps for the forward problem
     * @param maxsteps maximum number of solver steps (non-negative number)
     */
    void setMaxSteps(const int maxsteps);

    /**
     * @brief returns the maximum number of solver steps for the backward
     * problem
     * @return maximum number of solver steps
     */
    long int getMaxStepsBackwardProblem() const;

    /**
     * @brief sets the maximum number of solver steps for the backward problem
     * @param maxsteps maximum number of solver steps (non-negative number)
     */
    void setMaxStepsBackwardProblem(const int maxsteps);

    /**
     * @brief sets the maximum number of solver steps for the backward problem
     * @param maxsteps maximum number of solver steps (non-negative number)
     */
    void setMaxStepsBackwardProblem(const long int maxsteps);

    /**
     * @brief returns the linear system multistep method
     * @return linear system multistep method
     */
    LinearMultistepMethod getLinearMultistepMethod() const;

    /**
     * @brief sets the linear system multistep method
     * @param lmm linear system multistep method
     */
    void setLinearMultistepMethod(const LinearMultistepMethod lmm);

    /**
     * @brief returns the nonlinear system solution method
     * @return
     */
    NonlinearSolverIteration getNonlinearSolverIteration() const;

    /**
     * @brief sets the nonlinear system solution method
     * @param iter nonlinear system solution method
     */
    void setNonlinearSolverIteration(const NonlinearSolverIteration iter);

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
    void setInterpolationType(const InterpolationType interpType);

    /**
     * @brief sets KLU state ordering mode
     * @return
     */
    StateOrdering getStateOrdering() const;

    /**
     * @brief sets KLU state ordering mode (only applies when linsol is set to
     * amici.AMICI_KLU)
     * @param ordering state ordering
     */
    void setStateOrdering(const StateOrdering ordering);

    /**
     * @brief returns stability limit detection mode
     * @return stldet can be amici.FALSE (deactivated) or amici.TRUE (activated)
     */
    booleantype getStabilityLimitFlag() const;

    /**
     * @brief set stability limit detection mode
     * @param stldet can be amici.FALSE (deactivated) or amici.TRUE (activated)
     */
    void setStabilityLimitFlag(const booleantype stldet);

    /**
     * @brief getLinearSolver
     * @return
     */
    LinearSolver getLinearSolver() const;

    /**
     * @brief setLinearSolver
     * @param linsol
     */
    void setLinearSolver(const LinearSolver linsol);

    /**
     * @brief returns the internal sensitivity method
     * @return internal sensitivity method
     */
    InternalSensitivityMethod getInternalSensitivityMethod() const;

    /**
     * @brief sets the internal sensitivity method
     * @param ism internal sensitivity method
     */
    void setInternalSensitivityMethod(const InternalSensitivityMethod ism);
    
    /**
     * @brief write solution from forward simulation
     * @param t time
     * @param x state
     * @param dx derivative state
     * @param sx state sensitivity
     */
    void writeSolution(realtype *t, AmiVector &x, AmiVector &dx,
                       AmiVectorArray &sx) const;
    
    /**
     * @brief write solution from forward simulation
     * @param t time
     * @param xB adjoint state
     * @param dxB adjoint derivative state
     * @param xQB adjoint quadrature
     * @param which index of adjoint problem
     */
    void writeSolutionB(realtype *t, AmiVector &xB, AmiVector &dxB,
                        AmiVector &xQB, const int which) const;
    
    /**
     * @brief Access state solution at time t
     * @param t time
     * @return x or interpolated solution dky
     */
    const AmiVector &getState(const realtype t) const;
    
    /**
     * @brief Access derivative state solution at time t
     * @param t time
     * @return dx or interpolated solution dky
     */
    const AmiVector &getDerivativeState(const realtype t) const;

    /**
     * @brief Access state sensitivity solution at time t
     * @param t time
     * @return (interpolated) solution sx
     */
    const AmiVectorArray &getStateSensitivity(const realtype t) const;

    /**
     * @brief Access adjoint solution at time t
     * @param t time
     * @return (interpolated) solution xB
     */
    const AmiVector &getAdjointState(const int which, const realtype t) const;
    
    /**
     * @brief Access adjoint derivative solution at time t
     * @param t time
     * @return (interpolated) solution dxB
     */
    const AmiVector &getAdjointDerivativeState(const int which,
                                               const realtype t) const;

    /**
     * @brief Access adjoint quadrature solution at time t
     * @param t time
     * @return (interpolated) solution xQB
     */
    const AmiVector &getAdjointQuadrature(const int which, const realtype t) const;
    
    /**
     * @brief Reinitializes the states in the solver after an event occurence
     *
     * @param t0 initial timepoint
     * @param yy0 inital state variables
     * @param yp0 initial derivative state variables (DAE only)
     */
    virtual void reInit(const realtype t0, const AmiVector &yy0,
                        const AmiVector &yp0) const = 0;
    
    /**
     * @brief Reinitializes the state sensitivites in the solver after an
     * event occurence
     *
     * @param yS0 new state sensitivity
     * @param ypS0 new derivative state sensitivities (DAE only)
     */
    virtual void sensReInit(const AmiVectorArray &yyS0,
                            const AmiVectorArray &ypS0) const = 0;
    
    /**
     * @brief Reinitializes the adjoint states after an event occurence
     *
     * @param which identifier of the backwards problem
     * @param yQB0 new adjoint state
     * @param yQB0 new adjoint derivative state
     */
    virtual void reInitB(const int which, const realtype tB0,
                         const AmiVector &yyB0, const AmiVector &ypB0) const = 0;
    
    /**
     * @brief Reinitialize the adjoint states after an event occurence
     *
     * @param which identifier of the backwards problem
     * @param yQB0 new adjoint quadrature state
     */
    virtual void quadReInitB(const int which, const AmiVector &yQB0) const = 0;

    const realtype gett() const;
    
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
     * @brief Serialize Solver (see boost::serialization::serialize)
     * @param ar Archive to serialize to
     * @param r Data to serialize
     * @param version Version number
     */
    template <class Archive>
    friend void boost::serialization::serialize(Archive &ar, Solver &r,
                                                const unsigned int version);

    /**
     * @brief Check equality of data members excluding solver memory
     * @param a
     * @param b
     * @return
     */
    friend bool operator==(const Solver &a, const Solver &b);

  protected:
    /**
     * @brief Sets a timepoint at which the simulation will be stopped
     *
     * @param tstop timepoint until which simulation should be performed
     */
    virtual void setStopTime(const realtype tstop) const = 0;
    
    /**
     * @brief Solves the forward problem until a predefined timepoint
     *
     * @param tout timepoint until which simulation should be performed
     * @param itask task identifier, can be CV_NORMAL or CV_ONE_STEP
     * @return status flag indicating success of execution
     */
    virtual int solve(const realtype tout, const int itask) const = 0;
    
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
    virtual int solveF(const realtype tout, const int itask,
                       int *ncheckPtr) const = 0;
    
    /**
     * @brief reInitPostProcessF
     * @param tnext next timepoint (defines integration direction)
     */
    virtual void reInitPostProcessF(const realtype tnext) const = 0;

    /**
     * @brief reInitPostProcessB
     * @param tnext next timepoint (defines integration direction)
     */
    virtual void reInitPostProcessB(const realtype tnext) const = 0;

    /**
     * @brief sets sx to the state sensitivity at the current timepoint
     */
    virtual void getSens() const = 0;

    /**
     * @brief sets xB to the adjoint state at the current timepoint
     * @param which index of the backwards problem
     */
    virtual void getB(const int which) const = 0;

    /**
     * @brief sets xQB to the adjoint quadrature state at the current timepoint
     * @param which index of the backwards problem
     */
    virtual void getQuadB(const int which) const = 0;

    /**
     * @brief Initialises the states at the specified initial timepoint
     *
     * @param t0 initial timepoint
     * @param x0 initial states
     * @param dx0 initial derivative states
     */
    virtual void init(const realtype t0, const AmiVector &x0,
                      const AmiVector &dx0) const = 0;

    /**
     * @brief initialises the forward sensitivities
     * @param sx0 initial states semsitivities
     * @param sdx0 initial derivative states sensitivities
     */
    virtual void sensInit1(const AmiVectorArray &sx0,
                           const AmiVectorArray &sdx0) const = 0;

    /**
     * @brief Initialise the adjoint states at the specified final timepoint
     *
     * @param which identifier of the backwards problem
     * @param tf final timepoint
     * @param xB0 initial adjoint state
     * @param dxB0 initial adjoint derivative state
     */
    virtual void binit(const int which, const realtype tf, const AmiVector &xB0,
                       const AmiVector &dxB0) const = 0;

    /**
     * @brief Initialise the quadrature states at the specified final timepoint
     *
     * @param which identifier of the backwards problem
     * @param qQB0 intial adjoint quadrature state
     */
    virtual void qbinit(const int which, const AmiVector &xQB0) const = 0;

    /**
     * @brief Initialises the rootfinding for events
     *
     * @param ne number of different events
     */
    virtual void rootInit(const int ne) const = 0;

    /**
     * @brief Initalize non-linear solver for sensitivities
     * @param x
     * @param model
     */
    void initalizeNonLinearSolverSens(const Model *model) const;

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
    virtual void setDenseJacFnB(const int which) const = 0;

    /**
     * @brief sets the sparse Jacobian function
     *
     * @param which identifier of the backwards problem
     */
    virtual void setSparseJacFnB(const int which) const = 0;

    /**
     * @brief sets the banded Jacobian function
     *
     * @param which identifier of the backwards problem
     */
    virtual void setBandJacFnB(const int which) const = 0;

    /**
     * @brief sets the Jacobian vector multiplication function
     *
     * @param which identifier of the backwards problem
     */
    virtual void setJacTimesVecFnB(const int which) const = 0;

    /**
     * @brief Extracts diagnosis information from solver memory block and
     * writes them into the return data object for the backward problem
     *
     * @param error_code error identifier
     * @param module name of the module in which the error occured
     * @param function name of the function in which the error occured @type
     * char
     * @param msg error message
     * @param eh_data unused input
     */
    static void wrapErrHandlerFn(int error_code, const char *module,
                                 const char *function, char *msg,
                                 void *eh_data);

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
    virtual void setSStolerances(const double rtol,
                                 const double atol) const = 0;

    /**
     * @brief activates sets scalar relative and absolute tolerances for the
     * sensitivity variables
     *
     * @param rtol relative tolerances
     * @param atol array of absolute tolerances for every sensitivy variable
     */
    virtual void setSensSStolerances(const double rtol,
                                     const double *atol) const = 0;

    /**
     * SetSensErrCon specifies whether error control is also enforced for
     * sensitivities for the forward problem
     *
     * @param error_corr activation flag
     */
    virtual void setSensErrCon(const bool error_corr) const = 0;

    /**
     * @brief Specifies whether error control is also enforced for the
     * backward quadrature problem
     *
     * @param which identifier of the backwards problem
     * @param flag activation flag
     */
    virtual void setQuadErrConB(const int which, const bool flag) const = 0;

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
     * @param model Model instance,
     */
    virtual void setUserData(Model *model) const = 0;

    /**
     * @brief attaches the user data instance (here this is a Model) to the
     * backward problem
     *
     * @param which identifier of the backwards problem
     * @param model Model instance,
     */
    virtual void setUserDataB(const int which, Model *model) const = 0;

    /**
     * @brief specifies the maximum number of steps for the forward
     * problem
     *
     * @param mxsteps number of steps
     */
    virtual void setMaxNumSteps(long int mxsteps) const = 0;

    /**
     * @brief specifies the maximum number of steps for the forward
     * problem
     *
     * @param which identifier of the backwards problem
     * @param mxstepsB number of steps
     */
    virtual void setMaxNumStepsB(const int which, long int mxstepsB) const = 0;

    /**
     * @brief activates stability limit detection for the forward
     * problem
     *
     * @param stldet flag for stability limit detection (TRUE or FALSE)
     *
     */
    virtual void setStabLimDet(const int stldet) const = 0;

    /**
     * @brief activates stability limit detection for the backward
     * problem
     *
     * @param which identifier of the backwards problem
     * @param stldet flag for stability limit detection (TRUE or FALSE)
     *
     */
    virtual void setStabLimDetB(const int which, const int stldet) const = 0;

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
     * @param p paramaters
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
     * @param dky interpolated solution
     */
    virtual void getDky(const realtype t, const int k) const = 0;

    /**
     * @brief interpolates the (derivative of the) solution at the requested
     * timepoint
     *
     * @param t timepoint
     * @param k derivative order
     * @param dky interpolated solution
     * @param which index of backward problem
     */
    virtual void getDkyB(const realtype t, const int k,
                         const int which) const = 0;

    /**
     * @brief interpolates the (derivative of the) solution at the requested
     * timepoint
     *
     * @param t timepoint
     * @param k derivative order
     * @param dky interpolated solution
     */
    virtual void getSensDky(const realtype t, const int k) const = 0;

    /**
     * @brief interpolates the (derivative of the) solution at the requested
     * timepoint
     *
     * @param t timepoint
     * @param k derivative order
     * @param dky interpolated solution
     * @param which index of backward problem
     */
    virtual void getQuadDkyB(const realtype t, const int k,
                             const int which) const = 0;

    /**
     * @brief initializes the adjoint problem
     *
     */
    virtual void adjInit() const = 0;

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
    virtual void setSStolerancesB(const int which, const realtype relTolB,
                                  const realtype absTolB) const = 0;

    /**
     * @brief sets relative and absolute tolerances for the quadrature
     * backward problem
     *
     * @param which identifier of the backwards problem
     * @param reltolQB relative tolerances
     * @param abstolQB absolute tolerances
     */
    virtual void quadSStolerancesB(const int which, const realtype reltolQB,
                                   const realtype abstolQB) const = 0;

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
     * @param x
     */
    void initializeLinearSolver(const Model *model) const;

    /**
     * @brief Sets the non-linear solver
     * @param x
     */
    void initializeNonLinearSolver() const;

    /**
     * @brief Sets the linear solver for the forward problem
     */
    virtual void setLinearSolver() const = 0;

    /**
     * @brief Sets the linear solver for the backward problem
     * @param which
     */
    virtual void setLinearSolverB(const int which) const = 0;

    /**
     * @brief Set the non-linear solver for the forward problem
     */
    virtual void setNonLinearSolver() const = 0;

    /**
     * @brief Set the non-linear solver for the backward problem
     * @param which
     */
    virtual void setNonLinearSolverB(const int which) const = 0;

    /**
     * @brief Set the non-linear solver for sensitivities
     */
    virtual void setNonLinearSolverSens() const = 0;

    /**
     * @brief Initializes the linear solver for the backward problem
     *
     * @param model pointer to the model object
     * @param xB
     * @param which index of the backward problem
     */

    void initializeLinearSolverB(const Model *model, const int which) const;

    /**
     * @brief Initializes the non-linear solver for the backward problem
     * @param xB
     * @param which
     */
    void initializeNonLinearSolverB(const int which) const;

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
    bool getInitDoneB(const int which) const;

    /**
     * @brief checks whether memory for backward quadratures has been allocated
     * @param which adjoint problem index
     * @return proxy for solverMemoryB->(cv|ida)_QuadMallocDone
     */
    bool getQuadInitDoneB(const int which) const;

    /**
     * @brief attaches a diagonal linear solver to the forward problem
     */
    virtual void diag() const = 0;

    /**
     * @brief attaches a diagonal linear solver to the backward problem
     *
     * @param which identifier of the backwards problem
     */
    virtual void diagB(const int which) const = 0;

    /**
     * @brief resets solverMemory and solverMemoryB
     */
    void resetMutableMemory(const int nx, const int nplist, const int nquad) const;

    /**
     * @brief retrieves the solver memory instance for the backward problem
     *
     * @param which identifier of the backwards problem
     * @param ami_mem pointer to the forward solver memory instance
     * @return ami_memB pointer to the backward solver memory instance
     */
    virtual void *getAdjBmem(void *ami_mem, const int which) const = 0;

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
    void applyTolerancesASA(const int which) const;

    /**
     * @brief updates ASA quadrature solver tolerances according to the
     * currently specified member variables
     *
     * @param which identifier of the backwards problem
     */
    void applyQuadTolerancesASA(const int which) const;

    /**
     * @brief updates all senstivivity solver tolerances according to the
     * currently specified member variables
     */
    void applySensitivityTolerances() const;

    /** pointer to solver memory block */
    mutable std::unique_ptr<void, std::function<void(void *)>> solverMemory;

    /** pointer to solver memory block */
    mutable std::vector<std::unique_ptr<void, std::function<void(void *)>>>
        solverMemoryB;

    /** internal sensitivity method flag used to select the sensitivity solution
     * method. Only applies for Forward Sensitivities. */
    InternalSensitivityMethod ism = InternalSensitivityMethod::simultaneous;

    /** specifies the linear multistep method.
     */
    LinearMultistepMethod lmm = LinearMultistepMethod::BDF;

    /**
     * specifies the type of nonlinear solver iteration
     */
    NonlinearSolverIteration iter = NonlinearSolverIteration::newton;

    /** interpolation type for the forward problem solution which
     * is then used for the backwards problem.
     */
    InterpolationType interpType = InterpolationType::hermite;

    /** maximum number of allowed integration steps */
    long int maxsteps = 10000;

    /** linear solver for the forward problem */
    mutable std::unique_ptr<SUNLinSolWrapper> linearSolver;
    /** linear solver for the backward problem */
    mutable std::unique_ptr<SUNLinSolWrapper> linearSolverB;

    /** non-linear solver for the forward problem */
    mutable std::unique_ptr<SUNNonLinSolWrapper> nonLinearSolver;
    /** non-linear solver for the backward problem */
    mutable std::unique_ptr<SUNNonLinSolWrapper> nonLinearSolverB;
    /** non-linear solver for the sensitivities */
    mutable std::unique_ptr<SUNNonLinSolWrapper> nonLinearSolverSens;

    /** flag indicating whether the forward solver has been called */
    mutable bool solverWasCalledF = false;
    /** flag indicating whether the backward solver has been called */
    mutable bool solverWasCalledB = false;

    /**
     * @brief sets that memory for the forward problem has been allocated
     */
    void setInitDone() const;

    /**
     * @brief sets that memory for forward sensitivities has been allocated
     */
    void setSensInitDone() const;

    /**
     * @brief sets that memory for forward interpolation has been allocated
     */
    void setAdjInitDone() const;

    /**
     * @brief sets that memory for the backward problem has been allocated
     * @param which adjoint problem index
     */
    void setInitDoneB(const int which) const;

    /**
     * @brief sets that memory for backward quadratures has been allocated
     * @param which adjoint problem index
     */
    void setQuadInitDoneB(const int which) const;

    /** state (dimension: nx_solver) */
    mutable AmiVector x = AmiVector(0);
    /** state interface variable (dimension: nx_solver) */
    mutable AmiVector dky = AmiVector(0);
    /** state derivative dummy (dimension: nx_solver) */
    mutable AmiVector dx = AmiVector(0);
    /** state sensititivities interface variable (dimension: nx_solver x nplist)
     */
    mutable AmiVectorArray sx = AmiVectorArray(0, 0);
    /** state derivative sensititivities dummy (dimension: nx_solver x nplist)
     */
    mutable AmiVectorArray sdx = AmiVectorArray(0, 0);

    /** adjoint state interface variable (dimension: nx_solver) */
    mutable AmiVector xB = AmiVector(0);
    /** adjoint derivative dummy variable (dimension: nx_solver) */
    mutable AmiVector dxB = AmiVector(0);
    /** adjoint quadrature interface variable (dimension: nJ x nplist) */
    mutable AmiVector xQB = AmiVector(0);

    /** integration time of the forward problem */
    mutable realtype t;

    /** flag to force reInitPostProcessF before next call to solve */
    mutable bool forceReInitPostProcessF = false;

    /** flag to force reInitPostProcessB before next call to solveB */
    mutable bool forceReInitPostProcessB = false;

  private:
    
    /** method for sensitivity computation */
    SensitivityMethod sensi_meth = SensitivityMethod::forward;

    /** flag controlling stability limit detection */
    booleantype stldet = true;

    /** state ordering */
    StateOrdering ordering = StateOrdering::AMD;

    /** maximum number of allowed Newton steps for steady state computation */
    long int newton_maxsteps = 0;

    /** maximum number of allowed linear steps per Newton step for steady state
     * computation */
    long int newton_maxlinsteps = 0;

    /** Preequilibration of model via Newton solver? */
    bool newton_preeq = false;

    /** linear solver specification */
    LinearSolver linsol = LinearSolver::KLU;

    /** absolute tolerances for integration */
    realtype atol = 1e-16;

    /** relative tolerances for integration */
    realtype rtol = 1e-8;

    /** absolute tolerances for forward sensitivity integration */
    realtype atol_fsa = NAN;

    /** relative tolerances for forward sensitivity integration */
    realtype rtol_fsa = NAN;

    /** absolute tolerances for adjoint sensitivity integration */
    realtype atolB = NAN;

    /** relative tolerances for adjoint sensitivity integration */
    realtype rtolB = NAN;

    /** absolute tolerances for backward quadratures */
    realtype quad_atol = 1e-12;

    /** relative tolerances for backward quadratures */
    realtype quad_rtol = 1e-8;

    /** absolute tolerances for steadystate computation */
    realtype ss_atol = NAN;

    /** relative tolerances for steadystate computation */
    realtype ss_rtol = NAN;

    /** absolute tolerances for steadystate computation */
    realtype ss_atol_sensi = NAN;

    /** relative tolerances for steadystate computation */
    realtype ss_rtol_sensi = NAN;

    /** maximum number of allowed integration steps for backward problem */
    long int maxstepsB = 0;

    /** flag indicating whether sensitivities are supposed to be computed */
    SensitivityOrder sensi = SensitivityOrder::none;

    /** flag indicating whether init was called */
    mutable bool initialized = false;
    /** flag indicating whether sensInit1 was called */
    mutable bool sensInitialized = false;
    /** flag indicating whether adjInit was called */
    mutable bool adjInitialized = false;
    /** vector of flags indicating whether binit was called for respective
     which */
    mutable std::vector<bool> initializedB{false};
    /** vector of flags indicating whether qbinit was called for respective
     which */
    mutable std::vector<bool> initializedQB{false};
    
    /* number of checkpoints in the forward problem */
    mutable int ncheckPtr;
};

bool operator==(const Solver &a, const Solver &b);

} // namespace amici

#endif // AMICISOLVER_H
