%module solver

// Add necessary symbols to generated header
%{
#include "amici/solver.h"
using namespace amici;
%}

%rename(equals) operator==;

// remove functions that use AmiVector(Array) since that class anyways cannot
// be exposed in swig
%ignore get_adjoit_derivative_state;
%ignore get_adjoint_quadrature;
%ignore get_quadrature;
%ignore get_adjointState;
%ignore get_derivative_state;
%ignore get_state;
%ignore get_state_sensitivity;
%ignore quad_reinit_b;
%ignore reinit;
%ignore reinit_b;
%ignore sens_reinit;
%ignore setup;
%ignore setup_b;
%ignore setup_steady_state;
%ignore write_solution;
%ignore write_solution_b;
%ignore calc_ic;
%ignore calc_ic_b;
%ignore sens_toggle_off;
%ignore solve_b;
%ignore step;
%ignore run;
%ignore run_b;
%ignore reset_diagnosis;
%ignore store_diagnosis;
%ignore store_diagnosis_b;
%ignore turn_off_root_finding;
%ignore get_root_info;
%ignore update_and_reinit_states_and_sensitivities;
%ignore get_cpu_time;
%ignore get_cpu_time_b;
%ignore get_last_order;
%ignore get_num_err_test_fails;
%ignore get_num_err_test_fails_b;
%ignore get_num_non_lin_solv_conv_fails;
%ignore get_num_non_lin_solv_conv_fails_b;
%ignore get_num_rhs_evals;
%ignore get_num_rhs_evals_b;
%ignore get_num_steps;
%ignore get_num_steps_b;
%ignore gett;
%ignore startTimer;
%ignore switchForwardSensisOff;
%ignore timeExceeded;
%ignore getSunContext;
%ignore get_adjoint_state;
%ignore get_adjoint_derivative_state;
%ignore reinit_quad_b;
%ignore get_logger;
%ignore set_logger;
%ignore get_sun_context;

// Solver.__repr__
%pythoncode %{
def _solver_repr(self: "Solver"):
    return "\n".join([
        self.this.__repr__()[:-1],
        "  reporting_mode: "
        f"{RDataReporting(self.get_return_data_reporting_mode())!r}",
        f"  sens_meth: {SensitivityMethod(self.get_sensitivity_method())!r}",
        f"  sens_order: {SensitivityOrder(self.get_sensitivity_order())!r}",
        f"  sens_meth_preeq: {SensitivityMethod(self.get_sensitivity_method_pre_equilibration())!r}",
        f"  maxsteps: {self.get_max_steps()}",
        f"  maxtime: {self.get_max_time()}s",
        f"  abs_tol: {self.get_absolute_tolerance()}",
        f"  rel_tol: {self.get_relative_tolerance()}",
        f"  abs_tol_b: {self.get_absolute_tolerance_b()}",
        f"  rel_tol_b: {self.get_relative_tolerance_b()}",
        f"  abs_tol_fsa: {self.get_absolute_tolerance_fsa()}",
        f"  rel_tol_fsa: {self.get_relative_tolerance_fsa()}",
        f"  abs_tol_quad: {self.get_absolute_tolerance_quadratures()}",
        f"  rel_tol_quad: {self.get_relative_tolerance_quadratures()}",
        f"  abs_tol_ss: {self.get_absolute_tolerance_steady_state()}",
        f"  rel_tol_ss: {self.get_relative_tolerance_steady_state()}",
        f"  abs_tol_sss: {self.get_absolute_tolerance_steady_state_sensi()}",
        f"  rel_tol_sss: {self.get_relative_tolerance_steady_state_sensi()}",
        f"  int_sens_meth: {InternalSensitivityMethod(self.get_internal_sensitivity_method())!r}",
        f"  int_type: {InterpolationType(self.get_interpolation_type())!r}",
        f"  linsol: {LinearSolver(self.get_linear_solver())!r}",
        f"  lmm: {LinearMultistepMethod(self.get_linear_multistep_method())!r}",
        f"  newton_damp_mode: {NewtonDampingFactorMode(self.get_newton_damping_factor_mode())!r}",
        f"  newton_damp_lb: {self.get_newton_damping_factor_lower_bound()}",
        f"  newton_maxsteps: {self.get_newton_max_steps()}",
        f"  newton_ss_check: {self.get_newton_step_steady_state_check()}",
        f"  sens_ss_check: {self.get_sensi_steady_state_check()}",
        f"  interpolation_type: {InterpolationType(self.get_interpolation_type())!r}",
        f"  nlsol_iter: {NonlinearSolverIteration(self.get_non_linear_solver_iteration())!r}",
        f"  stability_limit: {self.get_stability_limit_flag()}",
        f"  state_ordering: {self.get_state_ordering()}",
        ">"
    ])

def _solver_reduce(self: "Solver"):
    """
    For now, we just store solver settings in a temporary HDF5 file.
    This is sufficient for multiprocessing use cases, but will not survive
    reboots and will not work in distributed (MPI) settings.
    This requires that amici was compiled with HDF5 support.
    """
    from amici.sim.sundials._swig_wrappers import restore_solver, write_solver_settings_to_hdf5
    from tempfile import NamedTemporaryFile
    import os
    with NamedTemporaryFile(suffix=".h5", delete=False) as tmpfile:
        tmpfilename = tmpfile.name
        write_solver_settings_to_hdf5(self, tmpfilename)

    return (
        restore_solver,
        (self.__class__, self.get_class_name(), tmpfilename,),
    )

%}
%extend amici::CVodeSolver {
%pythoncode %{
def __repr__(self):
    return _solver_repr(self)

def __reduce__(self):
    return _solver_reduce(self)

%}
};
%extend amici::IDASolver {
%pythoncode %{
def __repr__(self):
    return _solver_repr(self)

def __reduce__(self):
    return _solver_reduce(self)

%}
};

%extend std::unique_ptr<amici::Solver> {
%pythoncode %{
def __repr__(self):
    return _solver_repr(self)

def __deepcopy__(self, memo):
    return self.clone()

def __reduce__(self):
    return _solver_reduce(self)
%}
};

%extend amici::Solver {
%pythoncode %{
def __deepcopy__(self, memo):
    return self.clone()

def __reduce__(self):
    return _solver_reduce(self)

%}
};

%newobject amici::Solver::clone;
// Process symbols in header
%include "amici/solver.h"
