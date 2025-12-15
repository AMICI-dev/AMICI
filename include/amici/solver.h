#ifndef AMICI_SOLVER_H
#define AMICI_SOLVER_H

#include "amici/defines.h"
#include "amici/logging.h"
#include "amici/misc.h"
#include "amici/model_state.h"
#include "amici/sundials_linsol_wrapper.h"
#include "amici/vector.h"

#include <chrono>
#include <cmath>
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
namespace boost::serialization {
template <class Archive>
void serialize(Archive& ar, amici::Solver& s, unsigned int version);
} // namespace boost::serialization

namespace amici {

/**
 * The Solver class provides a generic interface to CVODES and IDAS solvers,
 * individual realizations are realized in the CVodeSolver and the IDASolver
 * class. All transient private/protected members (CVODES/IDAS memory, interface
 * variables and status flags) are specified as mutable and not included in
 * serialization or equality checks. No solver setting parameter should be
 * marked mutable.
 */
/*
 * NOTE: Any changes in data members here must be propagated to copy ctor,
 * equality operator, serialization functions in serialization.h, and
 * amici::hdf5::(read/write)SolverSettings(From/To)HDF5 in hdf5.cpp.
 */
class Solver {
  public:
    /** Type of what is passed to Sundials solvers as user_data */
    using user_data_type = std::pair<Model*, Solver const*>;
    /** Type of the function to free a raw sundials solver pointer */
    using free_solver_ptr = std::function<void(void*)>;
    /**
     * @brief Default constructor
     */
    Solver() = default;

    /**
     * @brief Solver copy constructor
     *
     * @param other
     */
    Solver(Solver const& other);

    virtual ~Solver() = default;

    /**
     * @brief Clone this instance
     *
     * @return The clone
     */
    virtual Solver* clone() const = 0;

    /**
     * @brief Get SUNDIALS context
     *
     * @return context
     */
    SUNContext get_sun_context() const;

    /**
     * @brief Get the name of this class.
     *
     * @return Class name.
     */
    virtual std::string get_class_name() const = 0;

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
    void run_b(realtype tout) const;

    /**
     * @brief Initializes the ami memory object and applies specified options
     *
     * @param t0 initial timepoint
     * @param model pointer to the model instance
     * @param x0 initial states
     * @param dx0 initial derivative states
     * @param sx0 initial state sensitivities
     * @param sdx0 initial derivative state sensitivities
     */

    void setup(
        realtype t0, Model* model, AmiVector const& x0, AmiVector const& dx0,
        AmiVectorArray const& sx0, AmiVectorArray const& sdx0
    ) const;

    /**
     * @brief Initializes the AMI memory object for the backwards problem
     *
     * @param which index of the backward problem, will be set by this routine
     * @param tf final timepoint (initial timepoint for the bwd problem)
     * @param model pointer to the model instance
     * @param xB0 initial adjoint states
     * @param dxB0 initial adjoint derivative states
     * @param xQB0 initial adjoint quadratures
     */

    void setup_b(
        int* which, realtype tf, Model* model, AmiVector const& xB0,
        AmiVector const& dxB0, AmiVector const& xQB0
    ) const;

    /**
     * @brief Initializes the ami memory for quadrature computation
     *
     * @param t0 initial timepoint
     * @param model pointer to the model instance
     * @param x0 initial states
     * @param dx0 initial derivative states
     * @param xB0 initial adjoint states
     * @param dxB0 initial derivative adjoint states
     * @param xQ0 initial quadrature vector
     */

    void setup_steady_state(
        realtype t0, Model* model, AmiVector const& x0, AmiVector const& dx0,
        AmiVector const& xB0, AmiVector const& dxB0, AmiVector const& xQ0
    ) const;

    /**
     * @brief Reinitializes state and respective sensitivities (if necessary)
     * according to changes in fixedParameters
     *
     * @param model pointer to the model instance
     */
    void update_and_reinit_states_and_sensitivities(Model* model) const;

    /**
     * @brief Extract information which event occurred.
     *
     * @param rootsfound array with flags indicating whether the respective
     * event occurred
     */
    virtual void get_root_info(int* rootsfound) const = 0;

    /**
     * @brief Calculates consistent initial conditions, assumes initial
     * states to be correct (DAE only)
     *
     * @param tout1 next timepoint to be computed (sets timescale)
     */
    virtual void calc_ic(realtype tout1) const = 0;

    /**
     * @brief Calculates consistent initial conditions for the backwards
     * problem, assumes initial states to be correct (DAE only)
     *
     * @param which identifier of the backwards problem
     * @param tout1 next timepoint to be computed (sets timescale)
     */
    virtual void calc_ic_b(int which, realtype tout1) const = 0;

    /**
     * @brief Solves the backward problem until a predefined timepoint
     * (adjoint only)
     *
     * @param tBout timepoint until which simulation should be performed
     * @param itaskB task identifier, can be CV_NORMAL or CV_ONE_STEP
     */
    virtual void solve_b(realtype tBout, int itaskB) const = 0;

    /**
     * @brief Disable rootfinding
     */
    virtual void turn_off_root_finding() const = 0;

    /**
     * @brief Return current sensitivity method
     *
     * @return method enum
     */
    SensitivityMethod get_sensitivity_method() const;

    /**
     * @brief Set sensitivity method
     *
     * @param sensi_meth
     */
    void set_sensitivity_method(SensitivityMethod sensi_meth);

    /**
     * @brief Return current sensitivity method during preequilibration
     *
     * @return method enum
     */
    SensitivityMethod get_sensitivity_method_pre_equilibration() const;

    /**
     * @brief Set sensitivity method for preequilibration
     *
     * @param sensi_meth_preeq
     */
    void set_sensitivity_method_pre_equilibration(
        SensitivityMethod sensi_meth_preeq
    );

    /**
     * @brief Disable forward sensitivity integration (used in steady state sim)
     */
    void switch_forward_sensis_off() const;

    /**
     * @brief Get maximum number of allowed Newton steps for steady state
     * computation
     *
     * @return
     */
    int get_newton_max_steps() const;

    /**
     * @brief Set maximum number of allowed Newton steps for steady state
     * computation
     *
     * @param newton_maxsteps
     */
    void set_newton_max_steps(int newton_maxsteps);

    /**
     * @brief Get a state of the damping factor used in the Newton solver
     *
     * @return
     */
    NewtonDampingFactorMode get_newton_damping_factor_mode() const;

    /**
     * @brief Turn on/off a damping factor in the Newton method
     *
     * @param dampingFactorMode
     */
    void
    set_newton_damping_factor_mode(NewtonDampingFactorMode dampingFactorMode);

    /**
     * @brief Get a lower bound of the damping factor used in the Newton solver
     *
     * @return
     */
    double get_newton_damping_factor_lower_bound() const;

    /**
     * @brief Set a lower bound of the damping factor in the Newton solver
     *
     * @param dampingFactorLowerBound
     */
    void set_newton_damping_factor_lower_bound(double dampingFactorLowerBound);

    /**
     * @brief Get sensitivity order
     *
     * @return sensitivity order
     */
    SensitivityOrder get_sensitivity_order() const;

    /**
     * @brief Set the sensitivity order
     *
     * @param sensi sensitivity order
     */
    void set_sensitivity_order(SensitivityOrder sensi);

    /**
     * @brief Get the relative tolerances for the forward problem
     *
     * Same tolerance is used for the backward problem if not specified
     * differently via \ref Solver::set_relative_tolerance_b.
     *
     * @return relative tolerances
     */
    double get_relative_tolerance() const;

    /**
     * @brief Sets the relative tolerances for the forward problem
     *
     * Same tolerance is used for the backward problem if not specified
     * differently via \ref Solver::set_relative_tolerance_b.
     *
     * @param rtol relative tolerance (non-negative number)
     */
    void set_relative_tolerance(double rtol);

    /**
     * @brief Get the absolute tolerances for the forward problem
     *
     * Same tolerance is used for the backward problem if not specified
     * differently via \ref Solver::set_absolute_tolerance_b.
     *
     * @return absolute tolerances
     */
    double get_absolute_tolerance() const;

    /**
     * @brief Sets the absolute tolerances for the forward problem
     *
     * Same tolerance is used for the backward problem if not specified
     * differently via \ref Solver::set_absolute_tolerance_b.
     *
     * @param atol absolute tolerance (non-negative number)
     */
    void set_absolute_tolerance(double atol);

    /**
     * @brief Returns the relative tolerances for the forward sensitivity
     * problem
     *
     * @return relative tolerances
     */
    double get_relative_tolerance_fsa() const;

    /**
     * @brief Sets the relative tolerances for the forward sensitivity problem
     *
     * @param rtol relative tolerance (non-negative number)
     */
    void set_relative_tolerance_fsa(double rtol);

    /**
     * @brief Returns the absolute tolerances for the forward sensitivity
     * problem
     *
     * @return absolute tolerances
     */
    double get_absolute_tolerance_fsa() const;

    /**
     * @brief Sets the absolute tolerances for the forward sensitivity problem
     *
     * @param atol absolute tolerance (non-negative number)
     */
    void set_absolute_tolerance_fsa(double atol);

    /**
     * @brief Returns the relative tolerances for the adjoint sensitivity
     * problem
     *
     * @return relative tolerances
     */
    double get_relative_tolerance_b() const;

    /**
     * @brief Sets the relative tolerances for the adjoint sensitivity problem
     *
     * @param rtol relative tolerance (non-negative number)
     */
    void set_relative_tolerance_b(double rtol);

    /**
     * @brief Returns the absolute tolerances for the backward problem for
     * adjoint sensitivity analysis
     *
     * @return absolute tolerances
     */
    double get_absolute_tolerance_b() const;

    /**
     * @brief Sets the absolute tolerances for the backward problem for
     * adjoint sensitivity analysis
     *
     * @param atol absolute tolerance (non-negative number)
     */
    void set_absolute_tolerance_b(double atol);

    /**
     * @brief Returns the relative tolerance for the quadrature problem
     *
     * @return relative tolerance
     */
    double get_relative_tolerance_quadratures() const;

    /**
     * @brief sets the relative tolerance for the quadrature problem
     *
     * @param rtol relative tolerance (non-negative number)
     */
    void set_relative_tolerance_quadratures(double rtol);

    /**
     * @brief returns the absolute tolerance for the quadrature problem
     *
     * @return absolute tolerance
     */
    double get_absolute_tolerance_quadratures() const;

    /**
     * @brief sets the absolute tolerance for the quadrature problem
     *
     * @param atol absolute tolerance (non-negative number)
     */
    void set_absolute_tolerance_quadratures(double atol);

    /**
     * @brief returns the steady state simulation tolerance factor.
     *
     * Steady state simulation tolerances are the product of the simulation
     * tolerances and this factor, unless manually set with
     * `set(Absolute/Relative)ToleranceSteadyState()`.
     *
     * @return steady state simulation tolerance factor
     */
    double get_steady_state_tolerance_factor() const;

    /**
     * @brief set the steady state simulation tolerance factor.
     *
     * Steady state simulation tolerances are the product of the simulation
     * tolerances and this factor, unless manually set with
     * `set(Absolute/Relative)ToleranceSteadyState()`.
     *
     * @param factor tolerance factor (non-negative number)
     */
    void set_steady_state_tolerance_factor(double factor);

    /**
     * @brief returns the relative tolerance for the steady state problem
     *
     * @return relative tolerance
     */
    double get_relative_tolerance_steady_state() const;

    /**
     * @brief sets the relative tolerance for the steady state problem
     *
     * @param rtol relative tolerance (non-negative number)
     */
    void set_relative_tolerance_steady_state(double rtol);

    /**
     * @brief returns the absolute tolerance for the steady state problem
     *
     * @return absolute tolerance
     */
    double get_absolute_tolerance_steady_state() const;

    /**
     * @brief sets the absolute tolerance for the steady state problem
     *
     * @param atol absolute tolerance (non-negative number)
     */
    void set_absolute_tolerance_steady_state(double atol);

    /**
     * @brief returns the steady state sensitivity simulation tolerance factor.
     *
     * Steady state sensitivity simulation tolerances are the product of the
     * sensitivity simulation tolerances and this factor, unless manually set
     * with `set_{absolute,relative}_tolerance_steady_state_sensi()`.
     *
     * @return steady state simulation tolerance factor
     */
    double get_steady_state_sensi_tolerance_factor() const;

    /**
     * @brief set the steady state sensitivity simulation tolerance factor.
     *
     * Steady state sensitivity simulation tolerances are the product of the
     * sensitivity simulation tolerances and this factor, unless manually set
     * with `set_{absolute,relative}_tolerance_steady_state_sensi()`.
     *
     * @param factor tolerance factor (non-negative number)
     */
    void set_steady_state_sensi_tolerance_factor(double factor);

    /**
     * @brief returns the relative tolerance for the sensitivities of the
     * steady state problem
     *
     * @return relative tolerance
     */
    double get_relative_tolerance_steady_state_sensi() const;

    /**
     * @brief sets the relative tolerance for the sensitivities of the
     * steady state problem
     * @param rtol relative tolerance (non-negative number)
     */
    void set_relative_tolerance_steady_state_sensi(double rtol);

    /**
     * @brief returns the absolute tolerance for the sensitivities of the
     * steady state problem
     *
     * @return absolute tolerance
     */
    double get_absolute_tolerance_steady_state_sensi() const;

    /**
     * @brief sets the absolute tolerance for the sensitivities of the
     * steady state problem
     *
     * @param atol absolute tolerance (non-negative number)
     */
    void set_absolute_tolerance_steady_state_sensi(double atol);

    /**
     * @brief returns the maximum number of solver steps for the forward
     * problem
     *
     * @return maximum number of solver steps
     */
    long int get_max_steps() const;

    /**
     * @brief sets the maximum number of solver steps for the forward problem
     *
     * @param maxsteps maximum number of solver steps (positive number)
     */
    void set_max_steps(long int maxsteps);

    /**
     * @brief Returns the maximum time allowed for integration
     *
     * @return Time in seconds
     */
    double get_max_time() const;

    /**
     * @brief Set the maximum CPU time allowed for integration
     *
     * @param maxtime Time in seconds. Zero means infinite time.
     */
    void set_max_time(double maxtime);

    /**
     * @brief Start timer for tracking integration time
     */
    void start_timer() const;

    /**
     * @brief Check whether maximum integration time was exceeded
     *
     * @param interval Only check the time every ``interval`` ths call to avoid
     * potentially relatively expensive syscalls

     * @return True if the maximum integration time was exceeded,
     * false otherwise.
     */
    bool time_exceeded(int interval = 1) const;

    /**
     * @brief returns the maximum number of solver steps for the backward
     * problem
     *
     * @return maximum number of solver steps
     */
    long int get_max_steps_backward_problem() const;

    /**
     * @brief sets the maximum number of solver steps for the backward problem
     *
     * @param maxsteps maximum number of solver steps (non-negative number)
     *
     * @note default behaviour (100 times the value for the forward problem) can
     * be restored by passing maxsteps=0
     */
    void set_max_steps_backward_problem(long int maxsteps);

    /**
     * @brief returns the linear system multistep method
     *
     * @return linear system multistep method
     */
    LinearMultistepMethod get_linear_multistep_method() const;

    /**
     * @brief sets the linear system multistep method
     *
     * @param lmm linear system multistep method
     */
    void set_linear_multistep_method(LinearMultistepMethod lmm);

    /**
     * @brief returns the nonlinear system solution method
     *
     * @return
     */
    NonlinearSolverIteration get_non_linear_solver_iteration() const;

    /**
     * @brief sets the nonlinear system solution method
     *
     * @param iter nonlinear system solution method
     */
    void set_non_linear_solver_iteration(NonlinearSolverIteration iter);

    /**
     * @brief get_interpolation_type
     *
     * @return
     */
    InterpolationType get_interpolation_type() const;

    /**
     * @brief sets the interpolation of the forward solution that is used for
     * the backwards problem
     *
     * @param interpType interpolation type
     */
    void set_interpolation_type(InterpolationType interpType);

    /**
     * @brief Gets KLU / SuperLUMT state ordering mode
     *
     * @return State-ordering as integer according to
     * SUNLinSolKLU::StateOrdering or SUNLinSolSuperLUMT::StateOrdering (which
     * differ).
     */
    int get_state_ordering() const;

    /**
     * @brief Sets KLU / SuperLUMT state ordering mode
     *
     * This only applies when linsol is set to LinearSolver::KLU or
     * LinearSolver::SuperLUMT. Mind the difference between
     * SUNLinSolKLU::StateOrdering and SUNLinSolSuperLUMT::StateOrdering.
     *
     * @param ordering state ordering
     */
    void set_state_ordering(int ordering);

    /**
     * @brief returns stability limit detection mode
     *
     * @return stldet can be false (deactivated) or true (activated)
     */
    bool get_stability_limit_flag() const;

    /**
     * @brief set stability limit detection mode
     *
     * @param stldet can be false (deactivated) or true (activated)
     */
    void set_stability_limit_flag(bool stldet);

    /**
     * @brief get_linear_solver
     *
     * @return
     */
    LinearSolver get_linear_solver() const;

    /**
     * @brief set_linear_solver
     *
     * @param linsol
     */
    void set_linear_solver(LinearSolver linsol);

    /**
     * @brief returns the internal sensitivity method
     *
     * @return internal sensitivity method
     */
    InternalSensitivityMethod get_internal_sensitivity_method() const;

    /**
     * @brief sets the internal sensitivity method
     *
     * @param ism internal sensitivity method
     */
    void set_internal_sensitivity_method(InternalSensitivityMethod ism);

    /**
     * @brief returns the ReturnData reporting mode
     *
     * @return ReturnData reporting mode
     */
    RDataReporting get_return_data_reporting_mode() const;

    /**
     * @brief sets the ReturnData reporting mode
     *
     * @param rdrm ReturnData reporting mode
     */
    void set_return_data_reporting_mode(RDataReporting rdrm);

    /**
     * @brief write solution from forward simulation
     *
     * @param t time
     * @param x state
     * @param dx derivative state
     * @param sx state sensitivity
     * @param xQ quadrature
     */
    void write_solution(
        realtype& t, AmiVector& x, AmiVector& dx, AmiVectorArray& sx,
        AmiVector& xQ
    ) const;

    /**
     * @brief write solution from forward simulation
     *
     * @param t time
     * @param x state
     * @param dx derivative state
     * @param sx state sensitivity
     */
    void write_solution(
        realtype& t, AmiVector& x, AmiVector& dx, AmiVectorArray& sx
    ) const;

    /**
     * @brief write solution from forward simulation
     *
     * @param sol solution state
     */
    void write_solution(SolutionState& sol) const;

    /**
     * @brief write solution from forward simulation
     *
     * @param t Time for which to retrieve the solution
     * (interpolated if necessary). Must be greater than or equal to
     * the initial timepoint and less than or equal to the current timepoint.
     * @param sol solution state
     */
    void write_solution(realtype t, SolutionState& sol) const;

    /**
     * @brief write solution from backward simulation
     *
     * @param t time
     * @param xB adjoint state
     * @param dxB adjoint derivative state
     * @param xQB adjoint quadrature
     * @param which index of adjoint problem
     */
    void write_solution_b(
        realtype& t, AmiVector& xB, AmiVector& dxB, AmiVector& xQB, int which
    ) const;

    /**
     * @brief Access state solution at time t
     *
     * @param t time
     * @return x or interpolated solution dky
     */
    AmiVector const& get_state(realtype t) const;

    /**
     * @brief Access derivative state solution at time t
     *
     * @param t time
     * @return dx or interpolated solution dky
     */
    AmiVector const& get_derivative_state(realtype t) const;

    /**
     * @brief Access state sensitivity solution at time t
     *
     * @param t time
     * @return (interpolated) solution sx
     */
    AmiVectorArray const& get_state_sensitivity(realtype t) const;

    /**
     * @brief Access adjoint solution at time t
     *
     * @param which adjoint problem index
     * @param t time
     * @return (interpolated) solution xB
     */
    AmiVector const& get_adjoint_state(int which, realtype t) const;

    /**
     * @brief Access adjoint derivative solution at time t
     *
     * @param which adjoint problem index
     * @param t time
     * @return (interpolated) solution dxB
     */
    AmiVector const& get_adjoint_derivative_state(int which, realtype t) const;

    /**
     * @brief Access adjoint quadrature solution at time t
     *
     * @param which adjoint problem index
     * @param t time
     * @return (interpolated) solution xQB
     */
    AmiVector const& get_adjoint_quadrature(int which, realtype t) const;

    /**
     * @brief Access quadrature solution at time t
     *
     * @param t time
     * @return (interpolated) solution xQ
     */
    AmiVector const& get_quadrature(realtype t) const;

    /**
     * @brief Reinitializes the states in the solver after an event occurrence
     *
     * @param t0 reinitialization timepoint
     * @param yy0 initial state variables
     * @param yp0 initial derivative state variables (DAE only)
     */
    virtual void
    reinit(realtype t0, AmiVector const& yy0, AmiVector const& yp0) const
        = 0;

    /**
     * @brief Reinitializes the state sensitivities in the solver after an
     * event occurrence
     *
     * @param yyS0 new state sensitivity
     * @param ypS0 new derivative state sensitivities (DAE only)
     */
    virtual void
    sens_reinit(AmiVectorArray const& yyS0, AmiVectorArray const& ypS0) const
        = 0;

    /**
     * @brief Switches off computation of  state sensitivities without
     * deallocating the memory for sensitivities
     */
    virtual void sens_toggle_off() const = 0;

    /**
     * @brief Reinitializes the adjoint states after an event occurrence
     *
     * @param which identifier of the backwards problem
     * @param tB0 reinitialization timepoint
     * @param yyB0 new adjoint state
     * @param ypB0 new adjoint derivative state
     */
    virtual void reinit_b(
        int which, realtype tB0, AmiVector const& yyB0, AmiVector const& ypB0
    ) const
        = 0;

    /**
     * @brief Reinitialize the adjoint states after an event occurrence
     *
     * @param which identifier of the backwards problem
     * @param yQB0 new adjoint quadrature state
     */
    virtual void reinit_quad_b(int which, AmiVector const& yQB0) const = 0;

    /**
     * @brief current solver timepoint
     *
     * @return t
     */
    realtype get_t() const;

    /**
     * @brief Reads out the CPU time needed for forward solve
     *
     * @return cpu_time
     */
    realtype get_cpu_time() const;

    /**
     * @brief Reads out the CPU time needed for backward solve
     *
     * @return cpu_timeB
     */
    realtype get_cpu_time_b() const;

    /**
     * @brief number of state variables with which the solver was initialized
     *
     * @return x.getLength()
     */
    int nx() const;

    /**
     * @brief number of parameters with which the solver was initialized
     *
     * @return sx.getLength()
     */
    int nplist() const;

    /**
     * @brief number of quadratures with which the solver was initialized
     *
     * @return xQB.getLength()
     */
    int nquad() const;

    /**
     * @brief check if FSA is being computed
     *
     * @return flag
     */
    bool computing_fsa() const {
        return get_sensitivity_order() >= SensitivityOrder::first
               && get_sensitivity_method() == SensitivityMethod::forward
               && nplist() > 0 && nx() > 0;
    }

    /**
     * @brief check if ASA is being computed
     *
     * @return flag
     */
    bool computing_asa() const {
        return get_sensitivity_order() >= SensitivityOrder::first
               && get_sensitivity_method() == SensitivityMethod::adjoint
               && nplist() > 0;
    }

    /**
     * @brief Resets vectors containing diagnosis information
     */
    void reset_diagnosis() const;

    /**
     * @brief Stores diagnosis information from solver memory block for forward
     * problem
     */
    void store_diagnosis() const;

    /**
     * @brief Stores diagnosis information from solver memory block for backward
     * problem
     *
     * @param which identifier of the backwards problem
     */
    void store_diagnosis_b(int which) const;

    /**
     * @brief Number of integration steps for the forward problem
     * (dimension: nt).
     *
     * @return ns
     */
    std::vector<int> const& get_num_steps() const { return ns_; }

    /**
     * @brief Number of integration steps for the backward problem
     * (dimension: nt).
     *
     * @return nsB
     */
    std::vector<int> const& get_num_steps_b() const { return nsB_; }

    /**
     * @brief Number of right hand side evaluations for the forward problem
     * (dimension: nt)
     *
     * @return nrhs
     */
    std::vector<int> const& get_num_rhs_evals() const { return nrhs_; }

    /**
     * @brief Number of right hand side evaluations for the backward problem
     * (dimension: nt)
     *
     * @return nrhsB
     */
    std::vector<int> const& get_num_rhs_evals_b() const { return nrhsB_; }

    /**
     * @brief Number of error test failures during the forward simulation
     * (dimension: nt)
     *
     * @return netf
     */
    std::vector<int> const& get_num_err_test_fails() const { return netf_; }

    /**
     * @brief Number of error test failures during the backward simulation
     * (dimension: nt)
     *
     * @return netfB
     */
    std::vector<int> const& get_num_err_test_fails_b() const { return netfB_; }

    /**
     * @brief Number of non-linear solver convergence failures during the
     * forward simulation (dimension: nt).
     *
     * @return nnlscf
     */
    std::vector<int> const& get_num_non_lin_solv_conv_fails() const {
        return nnlscf_;
    }

    /**
     * @brief Number of non-linear solver convergence failures during the
     * backward simulation (dimension: nt).
     *
     * @return nnlscfB
     */
    std::vector<int> const& get_num_non_lin_solv_conv_fails_b() const {
        return nnlscfB_;
    }

    /**
     * @brief Accessor order
     *
     * @return order
     */
    std::vector<int> const& get_last_order() const { return order_; }

    /**
     * @brief Returns how convergence checks for steadystate computation are
     * performed.
     *
     * If activated, convergence checks are limited to every 25 steps
     * in the simulation solver to limit performance impact.
     *
     * @return boolean flag indicating newton step (true) or the right hand side
     * (false)
     */
    bool get_newton_step_steady_state_check() const {
        return newton_step_steadystate_conv_;
    }

    /**
     * @brief Returns how convergence checks for steadystate computation are
     * performed.
     *
     * @return boolean flag indicating state and sensitivity equations (true) or
     * only state variables (false).
     */
    bool get_sensi_steady_state_check() const {
        return check_sensi_steadystate_conv_;
    }

    /**
     * @brief Sets how convergence checks for steadystate computation are
     * performed.
     *
     * @param flag boolean flag to pick newton step (true) or the right hand
     * side (false, default)
     */
    void set_newton_step_steady_state_check(bool const flag) {
        newton_step_steadystate_conv_ = flag;
    }

    /**
     * @brief Sets for which variables convergence checks for steadystate
     * computation are performed.
     *
     * @param flag boolean flag to pick state and sensitivity equations (true,
     * default) or only state variables (false).
     */
    void set_sensi_steady_state_check(bool const flag) {
        check_sensi_steadystate_conv_ = flag;
    }

    /**
     * @brief Set the maximum number of nonlinear solver iterations permitted
     * per step.
     *
     * @param max_nonlin_iters maximum number of nonlinear solver iterations
     */
    void set_max_nonlin_iters(int max_nonlin_iters);

    /**
     * @brief Get the maximum number of nonlinear solver iterations permitted
     * per step.
     *
     * @return maximum number of nonlinear solver iterations
     */
    int get_max_nonlin_iters() const;

    /**
     * @brief Set the maximum number of nonlinear solver convergence failures
     * permitted per step.
     *
     * @param max_conv_fails maximum number of nonlinear solver convergence
     */
    void set_max_conv_fails(int max_conv_fails);

    /**
     * @brief Get the maximum number of nonlinear solver convergence failures
     * permitted per step.
     *
     * @return maximum number of nonlinear solver convergence
     */
    int get_max_conv_fails() const;

    /**
     * @brief Set constraints on the model state.
     *
     * See
     * https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#c.CVodeSetConstraints.
     *
     * @param constraints
     */
    void set_constraints(std::vector<realtype> const& constraints);

    /**
     * @brief Get constraints on the model state.
     *
     * @return constraints
     */
    std::vector<realtype> get_constraints() const {
        return constraints_.get_vector();
    }

    /**
     * @brief Set the maximum step size
     *
     * @param max_step_size maximum step size. `0.0` means no limit.
     */
    void set_max_step_size(realtype max_step_size);

    /**
     * @brief Get the maximum step size
     *
     * @return maximum step size
     */
    realtype get_max_step_size() const;

    /**
     * @brief Serialize Solver (see boost::serialization::serialize)
     *
     * @param ar Archive to serialize to
     * @param s Data to serialize
     * @param version Version number
     */
    template <class Archive>
    friend void boost::serialization::serialize(
        Archive& ar, Solver& s, unsigned int version
    );

    /**
     * @brief Check equality of data members excluding solver memory
     *
     * @param a
     * @param b
     * @return
     */
    friend bool operator==(Solver const& a, Solver const& b);

    /**
     * @brief Get the logger.
     *
     * @return Logger pointer or nullptr.
     */
    Logger* get_logger() const { return logger_; }

    /**
     * @brief Set the logger.
     *
     * @param logger Logger pointer or nullptr.
     */
    void set_logger(Logger* logger) { logger_ = logger; }

  protected:
    /**
     * @brief Sets a timepoint at which the simulation will be stopped
     *
     * @param tstop timepoint until which simulation should be performed
     */
    virtual void set_stop_time(realtype tstop) const = 0;

    /**
     * @brief Solves the forward problem until a predefined timepoint
     *
     * @param tout timepoint until which simulation should be performed
     * @param itask task identifier, can be \c CV_NORMAL or \c CV_ONE_STEP
     * @return status flag indicating success of execution
     */
    virtual int solve(realtype tout, int itask) const = 0;

    /**
     * @brief Solves the forward problem until a predefined timepoint
     * (adjoint only)
     *
     * @param tout timepoint until which simulation should be performed
     * @param itask task identifier, can be \c CV_NORMAL or \c CV_ONE_STEP
     * @param ncheckPtr pointer to a number that counts the internal
     * checkpoints
     * @return status flag indicating success of execution
     */
    virtual int solve_f(realtype tout, int itask, int* ncheckPtr) const = 0;

    /**
     * @brief Postprocess the solver memory after a discontinuity in the
     * forward problem
     *
     * @param tnext next timepoint (defines integration direction)
     */
    virtual void reinit_post_process_f(realtype tnext) const = 0;

    /**
     * @brief Postprocess the solver memory after a discontinuity in the
     * backward problem
     *
     * @param tnext next timepoint (defines integration direction)
     */
    virtual void reinit_post_process_b(realtype tnext) const = 0;

    /**
     * @brief extracts the state sensitivity at the current timepoint from
     * solver memory and writes it to the sx member variable
     */
    virtual void get_sens() const = 0;

    /**
     * @brief extracts the adjoint state at the current timepoint from
     * solver memory and writes it to the xB member variable
     *
     * @param which index of the backwards problem
     */
    virtual void get_b(int which) const = 0;

    /**
     * @brief extracts the adjoint quadrature state at the current timepoint
     * from solver memory and writes it to the xQB member variable
     *
     * @param which index of the backwards problem
     */
    virtual void get_quad_b(int which) const = 0;

    /**
     * @brief extracts the quadrature at the current timepoint from solver
     * memory and writes it to the xQ member variable
     *
     * @param t timepoint for quadrature extraction
     */
    virtual void get_quad(realtype& t) const = 0;

    /**
     * @brief Initializes the states at the specified initial timepoint
     *
     * @param t0 initial timepoint
     * @param x0 initial states
     * @param dx0 initial derivative states
     */
    virtual void
    init(realtype t0, AmiVector const& x0, AmiVector const& dx0) const
        = 0;

    /**
     * @brief Initializes the states at the specified initial timepoint
     *
     * @param t0 initial timepoint
     * @param x0 initial states
     * @param dx0 initial derivative states
     */
    virtual void init_steady_state(
        realtype t0, AmiVector const& x0, AmiVector const& dx0
    ) const
        = 0;

    /**
     * @brief Initializes the forward sensitivities
     *
     * @param sx0 initial states sensitivities
     * @param sdx0 initial derivative states sensitivities
     */
    virtual void
    sens_init_1(AmiVectorArray const& sx0, AmiVectorArray const& sdx0) const
        = 0;

    /**
     * @brief Initialize the adjoint states at the specified final timepoint
     *
     * @param which identifier of the backwards problem
     * @param tf final timepoint
     * @param xB0 initial adjoint state
     * @param dxB0 initial adjoint derivative state
     */
    virtual void b_init(
        int which, realtype tf, AmiVector const& xB0, AmiVector const& dxB0
    ) const
        = 0;

    /**
     * @brief Initialize the quadrature states at the specified final timepoint
     *
     * @param which identifier of the backwards problem
     * @param xQB0 initial adjoint quadrature state
     */
    virtual void qb_init(int which, AmiVector const& xQB0) const = 0;

    /**
     * @brief Initializes the root-finding for events
     *
     * @param ne number of different events
     */
    virtual void root_init(int ne) const = 0;

    /**
     * @brief Initialize non-linear solver for sensitivities
     *
     * @param model Model instance
     */
    void initialize_non_linear_solver_sens(Model const* model) const;

    /**
     * @brief Set the dense Jacobian function
     */
    virtual void set_dense_jac_fn() const = 0;

    /**
     * @brief sets the sparse Jacobian function
     */
    virtual void set_sparse_jac_fn() const = 0;

    /**
     * @brief sets the banded Jacobian function
     */
    virtual void set_band_jac_fn() const = 0;

    /**
     * @brief sets the Jacobian vector multiplication function
     */
    virtual void set_jac_times_vec_fn() const = 0;

    /**
     * @brief sets the dense Jacobian function
     *
     * @param which identifier of the backwards problem
     */
    virtual void set_dense_jac_fn_b(int which) const = 0;

    /**
     * @brief sets the sparse Jacobian function
     *
     * @param which identifier of the backwards problem
     */
    virtual void set_sparse_jac_fn_b(int which) const = 0;

    /**
     * @brief sets the banded Jacobian function
     *
     * @param which identifier of the backwards problem
     */
    virtual void set_band_jac_fn_b(int which) const = 0;

    /**
     * @brief sets the Jacobian vector multiplication function
     *
     * @param which identifier of the backwards problem
     */
    virtual void set_jac_times_vec_fn_b(int which) const = 0;

    /**
     * @brief sets the sparse Jacobian function for backward steady state case
     */
    virtual void set_sparse_jac_fn_ss() const = 0;

    /**
     * @brief Create specifies solver method and initializes solver memory for
     * the forward problem
     */
    virtual void allocate_solver() const = 0;

    /**
     * @brief sets scalar relative and absolute tolerances for the forward
     * problem
     *
     * @param rtol relative tolerances
     * @param atol absolute tolerances
     */
    virtual void set_ss_tolerances(double rtol, double atol) const = 0;

    /**
     * @brief activates sets scalar relative and absolute tolerances for the
     * sensitivity variables
     *
     * @param rtol relative tolerances
     * @param atol array of absolute tolerances for every sensitivity variable
     */
    virtual void set_sens_ss_tolerances(double rtol, double const* atol) const
        = 0;

    /**
     * @brief Set whether error control is also enforced for
     * sensitivities for the forward problem
     *
     * @param error_corr activation flag
     */
    virtual void set_sens_err_con(bool error_corr) const = 0;

    /**
     * @brief Specifies whether error control is also enforced for the
     * backward quadrature problem
     *
     * @param which identifier of the backwards problem
     * @param flag activation flag
     */
    virtual void set_quad_err_con_b(int which, bool flag) const = 0;

    /**
     * @brief Specifies whether error control is also enforced for the
     * forward quadrature problem
     *
     * @param flag activation flag
     */
    virtual void set_quad_err_con(bool flag) const = 0;

    /**
     * @brief Attaches the error handler function to the solver
     */
    virtual void set_err_handler_fn() const;

    /**
     * @brief Attaches the user data to the forward problem
     */
    virtual void set_user_data() const = 0;

    /**
     * @brief attaches the user data to the backward problem
     *
     * @param which identifier of the backwards problem
     */
    virtual void set_user_data_b(int which) const = 0;

    /**
     * @brief specifies the maximum number of steps for the forward
     * problem
     *
     * @param mxsteps number of steps
     * @note in contrast to the SUNDIALS method, this sets the overall maximum,
     * not the maximum between output times.
     */
    virtual void set_max_num_steps(long int mxsteps) const = 0;

    /**
     * @brief specifies the maximum number of steps for the forward
     * problem
     *
     * @param which identifier of the backwards problem
     * @param mxstepsB number of steps
     * @note in contrast to the SUNDIALS method, this sets the overall maximum,
     * not the maximum between output times.
     */
    virtual void set_max_num_steps_b(int which, long int mxstepsB) const = 0;

    /**
     * @brief activates stability limit detection for the forward
     * problem
     *
     * @param stldet flag for stability limit detection (TRUE or FALSE)
     *
     */
    virtual void set_stab_lim_det(int stldet) const = 0;

    /**
     * @brief activates stability limit detection for the backward
     * problem
     *
     * @param which identifier of the backwards problem
     * @param stldet flag for stability limit detection (TRUE or FALSE)
     *
     */
    virtual void set_stab_lim_det_b(int which, int stldet) const = 0;

    /**
     * @brief specify algebraic/differential components (DAE only)
     *
     * @param model model specification
     */
    virtual void set_id(Model const* model) const = 0;

    /**
     * @brief deactivates error control for algebraic components (DAE only)
     *
     * @param flag deactivation flag
     */
    virtual void set_suppress_alg(bool flag) const = 0;

    /**
     * @brief specifies the scaling and indexes for sensitivity
     * computation
     *
     * @param p parameters
     * @param pbar parameter scaling constants
     * @param plist parameter index list
     */
    virtual void set_sens_params(
        realtype const* p, realtype const* pbar, int const* plist
    ) const
        = 0;

    /**
     * @brief interpolates the (derivative of the) solution at the requested
     * timepoint
     *
     * @param t timepoint
     * @param k derivative order
     */
    virtual void get_dky(realtype t, int k) const = 0;

    /**
     * @brief interpolates the (derivative of the) solution at the requested
     * timepoint
     *
     * @param t timepoint
     * @param k derivative order
     * @param which index of backward problem
     */
    virtual void get_dky_b(realtype t, int k, int which) const = 0;

    /**
     * @brief interpolates the (derivative of the) solution at the requested
     * timepoint
     *
     * @param t timepoint
     * @param k derivative order
     */
    virtual void get_sens_dky(realtype t, int k) const = 0;

    /**
     * @brief interpolates the (derivative of the) solution at the requested
     * timepoint
     *
     * @param t timepoint
     * @param k derivative order
     * @param which index of backward problem
     */
    virtual void get_quad_dky_b(realtype t, int k, int which) const = 0;

    /**
     * @brief interpolates the (derivative of the) solution at the requested
     * timepoint
     *
     * @param t timepoint
     * @param k derivative order
     */
    virtual void get_quad_dky(realtype t, int k) const = 0;

    /**
     * @brief initializes the adjoint problem
     */
    virtual void adj_init() const = 0;

    /**
     * @brief initializes the quadratures#
     * @param xQ0 vector with initial values for xQ
     */
    virtual void quad_init(AmiVector const& xQ0) const = 0;

    /**
     * @brief Specifies solver method and initializes solver memory for the
     * backward problem
     *
     * @param which identifier of the backwards problem
     */
    virtual void allocate_solver_b(int* which) const = 0;

    /**
     * @brief sets relative and absolute tolerances for the backward
     * problem
     *
     * @param which identifier of the backwards problem
     * @param relTolB relative tolerances
     * @param absTolB absolute tolerances
     */
    virtual void
    set_ss_tolerances_b(int which, realtype relTolB, realtype absTolB) const
        = 0;

    /**
     * @brief sets relative and absolute tolerances for the quadrature
     * backward problem
     *
     * @param which identifier of the backwards problem
     * @param reltolQB relative tolerances
     * @param abstolQB absolute tolerances
     */
    virtual void
    quad_ss_tolerances_b(int which, realtype reltolQB, realtype abstolQB) const
        = 0;

    /**
     * @brief sets relative and absolute tolerances for the quadrature problem
     *
     * @param reltolQB relative tolerances
     * @param abstolQB absolute tolerances
     */
    virtual void quad_ss_tolerances(realtype reltolQB, realtype abstolQB) const
        = 0;

    /**
     * @brief reports the number of solver steps
     *
     * @param ami_mem pointer to the solver memory instance (can be from
     * forward or backward problem)
     * @param numsteps output array
     */
    virtual void get_num_steps(void const* ami_mem, long int* numsteps) const
        = 0;

    /**
     * @brief reports the number of right hand evaluations
     *
     * @param ami_mem pointer to the solver memory instance (can be from
     * forward or backward problem)
     * @param numrhsevals output array
     */
    virtual void
    get_num_rhs_evals(void const* ami_mem, long int* numrhsevals) const
        = 0;

    /**
     * @brief reports the number of local error test failures
     *
     * @param ami_mem pointer to the solver memory instance (can be from
     * forward or backward problem)
     * @param numerrtestfails output array
     */
    virtual void
    get_num_err_test_fails(void const* ami_mem, long int* numerrtestfails) const
        = 0;

    /**
     * @brief reports the number of nonlinear convergence failures
     *
     * @param ami_mem pointer to the solver memory instance (can be from
     * forward or backward problem)
     * @param numnonlinsolvconvfails output array
     */
    virtual void get_num_non_lin_solv_conv_fails(
        void const* ami_mem, long int* numnonlinsolvconvfails
    ) const
        = 0;

    /**
     * @brief Reports the order of the integration method during the
     * last internal step
     *
     * @param ami_mem pointer to the solver memory instance (can be from
     * forward or backward problem)
     * @param order output array
     */
    virtual void get_last_order(void const* ami_mem, int* order) const = 0;

    /**
     * @brief Initializes and sets the linear solver for the forward problem
     *
     * @param model pointer to the model object
     */
    void initialize_linear_solver(Model const* model) const;

    /**
     * @brief Sets the non-linear solver
     */
    void initialize_non_linear_solver() const;

    /**
     * @brief Sets the linear solver for the forward problem
     */
    virtual void set_linear_solver() const = 0;

    /**
     * @brief Sets the linear solver for the backward problem
     *
     * @param which index of the backward problem
     */
    virtual void set_linear_solver_b(int which) const = 0;

    /**
     * @brief Set the non-linear solver for the forward problem
     */
    virtual void set_non_linear_solver() const = 0;

    /**
     * @brief Set the non-linear solver for the backward problem
     *
     * @param which index of the backward problem
     */
    virtual void set_non_linear_solver_b(int which) const = 0;

    /**
     * @brief Set the non-linear solver for sensitivities
     */
    virtual void set_non_linear_solver_sens() const = 0;

    /**
     * @brief Initializes the linear solver for the backward problem
     *
     * @param model pointer to the model object
     * @param which index of the backward problem
     */

    void initialize_linear_solver_b(Model const* model, int which) const;

    /**
     * @brief Initializes the non-linear solver for the backward problem
     *
     * @param which index of the backward problem
     */
    void initialize_non_linear_solver_b(int which) const;

    /**
     * @brief Accessor function to the model stored in the user data
     *
     * @return user data model
     */
    virtual Model const* get_model() const = 0;

    /**
     * @brief checks whether memory for the forward problem has been allocated
     *
     * @return proxy for solverMemory->(cv|ida)_MallocDone
     */
    bool get_init_done() const;

    /**
     * @brief checks whether memory for forward sensitivities has been allocated
     *
     * @return proxy for solverMemory->(cv|ida)_SensMallocDone
     */
    bool get_sens_init_done() const;

    /**
     * @brief checks whether memory for forward interpolation has been allocated
     *
     * @return proxy for solverMemory->(cv|ida)_adjMallocDone
     */
    bool get_adj_init_done() const;

    /**
     * @brief checks whether memory for the backward problem has been allocated
     *
     * @param which adjoint problem index
     * @return proxy for solverMemoryB->(cv|ida)_MallocDone
     */
    bool get_init_done_b(int which) const;

    /**
     * @brief checks whether memory for backward quadratures has been allocated
     *
     * @param which adjoint problem index
     * @return proxy for solverMemoryB->(cv|ida)_QuadMallocDone
     */
    bool get_quad_init_done_b(int which) const;

    /**
     * @brief checks whether memory for quadratures has been allocated
     *
     * @return proxy for solverMemory->(cv|ida)_QuadMallocDone
     */
    bool get_quad_init_done() const;

    /**
     * @brief attaches a diagonal linear solver to the forward problem
     */
    virtual void diag() const = 0;

    /**
     * @brief attaches a diagonal linear solver to the backward problem
     *
     * @param which identifier of the backwards problem
     */
    virtual void diag_b(int which) const = 0;

    /**
     * @brief resets solverMemory and solverMemoryB
     *
     * @param nx new number of state variables
     * @param nplist new number of sensitivity parameters
     * @param nquad new number of quadratures (only differs from nplist for
     * higher order sensitivity computation)
     */
    void reset_mutable_memory(int nx, int nplist, int nquad) const;

    /**
     * @brief Retrieves the solver memory instance for the backward problem
     *
     * @param which identifier of the backwards problem
     * @param ami_mem pointer to the forward solver memory instance
     * @return A (void *) pointer to the CVODES memory allocated for the
     * backward problem.
     */
    virtual void* get_adj_b_mem(void* ami_mem, int which) const = 0;

    /**
     * @brief updates solver tolerances according to the currently specified
     * member variables
     */
    void apply_tolerances() const;

    /**
     * @brief updates FSA solver tolerances according to the currently
     * specified member variables
     */
    void apply_tolerances_fsa() const;

    /**
     * @brief updates ASA solver tolerances according to the currently
     * specified member variables
     *
     * @param which identifier of the backwards problem
     */
    void apply_tolerances_asa(int which) const;

    /**
     * @brief updates ASA quadrature solver tolerances according to the
     * currently specified member variables
     *
     * @param which identifier of the backwards problem
     */
    void apply_quad_tolerances_asa(int which) const;

    /**
     * @brief updates quadrature solver tolerances according to the
     * currently specified member variables
     */
    void apply_quad_tolerances() const;

    /**
     * @brief updates all sensitivity solver tolerances according to the
     * currently specified member variables
     */
    void apply_sensitivity_tolerances() const;

    /**
     * @brief Apply the constraints to the solver.
     */
    virtual void apply_constraints() const;

    /** SUNDIALS context */
    sundials::Context sunctx_;

    /** pointer to solver memory block */
    mutable std::unique_ptr<void, free_solver_ptr> solver_memory_;

    /** pointer to solver memory block */
    mutable std::vector<std::unique_ptr<void, free_solver_ptr>>
        solver_memory_B_;

    /** Sundials user_data_ */
    mutable user_data_type user_data_;

    /** internal sensitivity method flag used to select the sensitivity solution
     * method. Only applies for Forward Sensitivities. */
    InternalSensitivityMethod ism_{InternalSensitivityMethod::simultaneous};

    /** specifies the linear multistep method.
     */
    LinearMultistepMethod lmm_{LinearMultistepMethod::BDF};

    /**
     * specifies the type of nonlinear solver iteration
     */
    NonlinearSolverIteration iter_{NonlinearSolverIteration::newton};

    /** interpolation type for the forward problem solution which
     * is then used for the backwards problem.
     */
    InterpolationType interp_type_{InterpolationType::polynomial};

    /** maximum number of allowed integration steps */
    long int maxsteps_{10000};

    /** Maximum CPU-time for integration in seconds */
    std::chrono::duration<double, std::ratio<1>> maxtime_{0};

    /** Time at which solver timer was started */
    mutable CpuTimer simulation_timer_;

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
    mutable bool solver_was_called_F_{false};

    /** flag indicating whether the backward solver has been called */
    mutable bool solver_was_called_B_{false};

    /**
     * @brief sets that memory for the forward problem has been allocated
     */
    void set_init_done() const;

    /**
     * @brief sets that memory for forward sensitivities has been allocated
     */
    void set_sens_init_done() const;

    /**
     * @brief sets that memory for forward interpolation has been allocated
     */
    void set_adj_init_done() const;

    /**
     * @brief sets that memory for the backward problem has been allocated
     *
     * @param which adjoint problem index
     */
    void set_init_done_b(int which) const;

    /**
     * @brief sets that memory for backward quadratures has been allocated
     *
     * @param which adjoint problem index
     */
    void set_quad_init_done_b(int which) const;

    /**
     * @brief sets that memory for quadratures has been allocated
     */
    void set_quad_init_done() const;

    /**
     * @brief Sets sensitivity method (for simulation or preequilibration)
     *
     * @param sensi_meth new value for sensi_meth[_preeq]
     * @param preequilibration flag indicating preequilibration or simulation
     */
    void check_sensitivity_method(
        SensitivityMethod sensi_meth, bool preequilibration
    ) const;

    /**
     * @brief Apply the maximum number of nonlinear solver iterations permitted
     * per step.
     */
    virtual void apply_max_nonlin_iters() const = 0;

    /**
     * @brief Apply the maximum number of nonlinear solver convergence failures
     * permitted per step.
     */
    virtual void apply_max_conv_fails() const = 0;

    /**
     * @brief Apply the allowed maximum stepsize to the solver.
     */
    virtual void apply_max_step_size() const = 0;

    /** state (dimension: nx_solver) */
    mutable AmiVector x_{0, sunctx_};

    /** state interface variable (dimension: nx_solver) */
    mutable AmiVector dky_{0, sunctx_};

    /** state derivative dummy (dimension: nx_solver) */
    mutable AmiVector dx_{0, sunctx_};

    /** state sensitivities interface variable (dimension: nx_solver x nplist)
     */
    mutable AmiVectorArray sx_{0, 0, sunctx_};
    /** state derivative sensitivities dummy (dimension: nx_solver x nplist)
     */
    mutable AmiVectorArray sdx_{0, 0, sunctx_};

    /** adjoint state interface variable (dimension: nx_solver) */
    mutable AmiVector xB_{0, sunctx_};

    /** adjoint derivative dummy variable (dimension: nx_solver) */
    mutable AmiVector dxB_{0, sunctx_};

    /** adjoint quadrature interface variable (dimension: nJ x nplist) */
    mutable AmiVector xQB_{0, sunctx_};

    /** forward quadrature interface variable (dimension: nx_solver) */
    mutable AmiVector xQ_{0, sunctx_};

    /** integration time of the forward problem */
    mutable realtype t_{std::nan("")};

    /** flag to force reInitPostProcessF before next call to solve */
    mutable bool force_reinit_postprocess_F_{false};

    /** flag to force reInitPostProcessB before next call to solveB */
    mutable bool force_reinit_postprocess_B_{false};

    /** flag indicating whether sensInit1 was called */
    mutable bool sens_initialized_{false};

    /** Vector of constraints on the solution */
    mutable AmiVector constraints_;

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
    SensitivityMethod sensi_meth_{SensitivityMethod::forward};

    /** method for sensitivity computation in preequilibration */
    SensitivityMethod sensi_meth_preeq_{SensitivityMethod::forward};

    /** flag controlling stability limit detection */
    sunbooleantype stldet_{SUNTRUE};

    /** state ordering */
    int ordering_{static_cast<int>(SUNLinSolKLU::StateOrdering::AMD)};

    /** maximum number of allowed Newton steps for steady state computation */
    long int newton_maxsteps_{0L};

    /** Damping factor state used int the Newton method */
    NewtonDampingFactorMode newton_damping_factor_mode_{
        NewtonDampingFactorMode::on
    };

    /** Lower bound of the damping factor. */
    realtype newton_damping_factor_lower_bound_{1e-8};

    /** linear solver specification */
    LinearSolver linsol_{LinearSolver::KLU};

    /** absolute tolerances for integration */
    realtype atol_{1e-16};

    /** relative tolerances for integration */
    realtype rtol_{1e-8};

    /** absolute tolerances for forward sensitivity integration */
    realtype atol_fsa_{NAN};

    /** relative tolerances for forward sensitivity integration */
    realtype rtol_fsa_{NAN};

    /** absolute tolerances for adjoint sensitivity integration */
    realtype atolB_{NAN};

    /** relative tolerances for adjoint sensitivity integration */
    realtype rtolB_{NAN};

    /** absolute tolerances for backward quadratures */
    realtype quad_atol_{1e-12};

    /** relative tolerances for backward quadratures */
    realtype quad_rtol_{1e-8};

    /** steady state simulation tolerance factor */
    realtype ss_tol_factor_{1e2};

    /** absolute tolerances for steadystate computation */
    realtype ss_atol_{NAN};

    /** relative tolerances for steadystate computation */
    realtype ss_rtol_{NAN};

    /** steady state sensitivity simulation tolerance factor */
    realtype ss_tol_sensi_factor_{1e2};

    /** absolute tolerances for steadystate sensitivity computation */
    realtype ss_atol_sensi_{NAN};

    /** relative tolerances for steadystate sensitivity computation */
    realtype ss_rtol_sensi_{NAN};

    RDataReporting rdata_mode_{RDataReporting::full};

    /** whether newton step should be used for convergence steps */
    bool newton_step_steadystate_conv_{false};

    /** whether sensitivities should be checked for convergence to steadystate
     */
    bool check_sensi_steadystate_conv_{true};

    /** Maximum number of nonlinear solver iterations permitted per step */
    int max_nonlin_iters_{3};

    /** Maximum number of nonlinear solver convergence failures permitted per
     *  step */
    int max_conv_fails_{10};

    /** Maximum allowed step size */
    realtype max_step_size_{0.0};

    /** CPU time, forward solve */
    mutable realtype cpu_time_{0.0};

    /** CPU time, backward solve */
    mutable realtype cpu_time_b_{0.0};

    /** maximum number of allowed integration steps for backward problem */
    long int maxstepsB_{0L};

    /** flag indicating whether sensitivities are supposed to be computed */
    SensitivityOrder sensi_{SensitivityOrder::none};

    /** flag indicating whether init was called */
    mutable bool initialized_{false};

    /** flag indicating whether adjInit was called */
    mutable bool adj_initialized_{false};

    /** flag indicating whether (forward) quadInit was called */
    mutable bool quad_initialized_{false};

    /** vector of flags indicating whether binit was called for respective
     which */
    mutable std::vector<bool> initializedB_{false};

    /** vector of flags indicating whether qbinit was called for respective
     which */
    mutable std::vector<bool> initializedQB_{false};

    /** number of checkpoints in the forward problem */
    mutable int ncheckPtr_{0};

    /** number of integration steps forward problem (dimension: nt) */
    mutable std::vector<int> ns_;

    /** number of integration steps backward problem (dimension: nt) */
    mutable std::vector<int> nsB_;

    /** number of right hand side evaluations forward problem (dimension: nt) */
    mutable std::vector<int> nrhs_;

    /** number of right hand side evaluations backward problem (dimension: nt)
     */
    mutable std::vector<int> nrhsB_;

    /** number of error test failures forward problem (dimension: nt) */
    mutable std::vector<int> netf_;

    /** number of error test failures backward problem (dimension: nt) */
    mutable std::vector<int> netfB_;

    /**
     * number of non-linear solver convergence failures forward problem
     * (dimension: nt)
     */
    mutable std::vector<int> nnlscf_;

    /**
     * number of non-linear solver convergence failures backward problem
     * (dimension: nt)
     */
    mutable std::vector<int> nnlscfB_;

    /** employed order forward problem (dimension: nt) */
    mutable std::vector<int> order_;

    /** logger */
    Logger* logger_ = nullptr;
};

bool operator==(Solver const& a, Solver const& b);

} // namespace amici

#endif // AMICISOLVER_H
