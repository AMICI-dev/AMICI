%module model

// Add necessary symbols to generated header
%{
#include "amici/model.h"
using namespace amici;
%}

// remove functions that use AmiVector(Array) since that class anyways cannot
// be exposed in swig
%ignore add_adjoint_quadrature_eventUpdate;
%ignore add_adjoint_state_event_update;
%ignore add_event_objective;
%ignore add_event_objective_regularization;
%ignore add_event_objective_sensitivity;
%ignore add_observable_objective;
%ignore add_observable_objective_sensitivity;
%ignore add_partial_event_objective_sensitivity;
%ignore add_partial_observable_objective_sensitivity;
%ignore add_state_event_update;
%ignore add_state_sensitivity_event_update;
%ignore fsx_rdata;
%ignore fx_rdata;
%ignore get_adjoint_state_event_update;
%ignore get_event_time_sensitivity;
%ignore get_adjoint_state_observable_update;
%ignore get_event;
%ignore get_events;
%ignore get_event_regularization;
%ignore get_event_regularization_sensitivity;
%ignore get_event_sensitivity;
%ignore get_event_time_sensitivity;
%ignore get_explicit_roots;
%ignore get_observable;
%ignore get_observable_sensitivity;
%ignore get_expression;
%ignore init_events;
%ignore reinit_events;
%ignore initialize;
%ignore initialize_b;
%ignore initialize_state_sensitivities;
%ignore initialize_state;
%ignore reinitialize;
%ignore ModelState;
%ignore get_model_state;
%ignore set_model_state;
%ignore fx0;
%ignore fx0_fixedParameters;
%ignore fsx0;
%ignore fsx0_fixedParameters;
%ignore get_dxdotdp;
%ignore get_dxdotdp_full;
%ignore check_inite;
%ignore fJrz;
%ignore fJy;
%ignore fJz;
%ignore fdJrzdsigma;
%ignore fdJrzdz;
%ignore fdJzdsigma;
%ignore fdJzdz;
%ignore fdJydsigma;
%ignore fdeltaqB;
%ignore fdeltasx;
%ignore fdeltax;
%ignore fdeltaxB;
%ignore fdrzdp;
%ignore fdrzdx;
%ignore fdsigmaydp;
%ignore fdsigmazdp;
%ignore fdydp;
%ignore fdydx;
%ignore fdzdp;
%ignore fdzdx;
%ignore frz;
%ignore fsdx0;
%ignore fsigmay;
%ignore fsigmaz;
%ignore fsrz;
%ignore fstau;
%ignore fsz;
%ignore fw;
%ignore fy;
%ignore fz;
%ignore update_heaviside;
%ignore update_heaviside_b;
%ignore get_event_sigma;
%ignore get_event_sigma_sensitivity;
%ignore get_observable_sigma;
%ignore get_observable_sigma_sensitivity;
%ignore get_unobserved_event_sensitivity;
%ignore fdsigmaydy;
%ignore fdspline_slopesdp;
%ignore fdspline_valuesdp;
%ignore fdtotal_cldp;
%ignore fdtotal_cldx_rdata;
%ignore fdx_rdatadp;
%ignore fdx_rdatadtcl;
%ignore fdx_rdatadx_solver;
%ignore fdsigmaydy;
%ignore get_steadystate_mask_av;
%ignore initialize_splines;
%ignore initialize_spline_sensitivities;
%ignore initialize_events;
%ignore reinit_explicit_roots;
%ignore add_adjoint_quadrature_event_update;
%ignore check_finite;

%newobject amici::Model::clone;

%rename(create_solver) amici::Model::get_solver;
%rename(_cpp_model_clone) amici::Model::clone;

%extend amici::Model {
%pythoncode %{
def clone(self):
    """Clone the model instance."""
    clone = self._cpp_model_clone()
    try:
        # copy module reference if present
        clone.module = self.module
    except Exception:
        pass

    return clone

def __deepcopy__(self, memo):
    return self.clone()

def __reduce__(self):
    from amici.sim.sundials._swig_wrappers import restore_model, get_model_settings, file_checksum

    return (
        restore_model,
        (
            self.get_name(),
            Path(self.module.__spec__.origin).parent,
            get_model_settings(self),
            file_checksum(self.module.extension_path),
        ),
        {}
    )

@overload
def simulate(
    self: AmiciModel,
    *,
    solver: Solver | None = None,
    edata: AmiciExpData | None = None,
    sensi_method: SensitivityMethod | str = None,
    sensi_order: SensitivityOrder | str = None,
) -> ReturnDataView: ...


@overload
def simulate(
    self: AmiciModel,
    *,
    solver: Solver | None = None,
    edata: AmiciExpDataVector | None = None,
    failfast: bool = True,
    num_threads: int = 1,
    sensi_method: SensitivityMethod | str = None,
    sensi_order: SensitivityOrder | str = None,
) -> list[ReturnDataView]: ...


def simulate(
    self: AmiciModel,
    *,
    solver: Solver | None = None,
    edata: AmiciExpData | AmiciExpDataVector | None = None,
    failfast: bool = True,
    num_threads: int = 1,
    sensi_method: SensitivityMethod | str = None,
    sensi_order: SensitivityOrder | str = None,
) -> ReturnDataView | list[ReturnDataView]:
    """Simulate model with given solver and experimental data.

    :param solver:
        Solver to use for simulation. Defaults to :meth:`Model.create_solver`.
    :param edata:
        Experimental data to use for simulation.
        A single :class:`ExpData` instance or a sequence of such instances.
        If `None`, no experimental data is used and the model is simulated
        as is.
    :param sensi_method:
        Sensitivity method to use for simulation.
        If `None`, the solver's current sensitivity method is used.
    :param sensi_order:
        Sensitivity order to use for simulation.
        If `None`, the solvers's current sensitivity order is used.
    :param failfast:
        Whether to stop simulations on first failure.
        Only relevant if `edata` is a sequence of :class:`ExpData` instances.
    :param num_threads:
        Number of threads to use for simulation.
        Only relevant if AMICI was compiled with OpenMP support and if `edata`
        is a sequence of :class:`ExpData` instances.
    :return:
        A single :class:`ReturnDataView` instance containing the simulation
        results if `edata` is a single :class:`ExpData` instance or `None`.
        If `edata` is a sequence of :class:`ExpData` instances, a list of
        :class:`ReturnDataView` instances is returned.
    """
    from amici.sim.sundials._swig_wrappers import _Model__simulate

    return _Model__simulate(
        self,
        solver=solver,
        edata=edata,
        failfast=failfast,
        num_threads=num_threads,
        sensi_method=sensi_method,
        sensi_order=sensi_order,
    )
%}
};

%extend std::unique_ptr<amici::Model> {
%pythoncode %{
def clone(self):
    """Clone the model instance."""
    clone = self._cpp_model_clone()
    try:
        # copy module reference if present
        clone.module = self.module
    except Exception:
        pass

    return clone

def __deepcopy__(self, memo):
    return self.clone()

def __reduce__(self):
    from amici.sim.sundials._swig_wrappers import restore_model, get_model_settings, file_checksum

    return (
        restore_model,
        (
            self.get_name(),
            Path(self.module.__spec__.origin).parent,
            get_model_settings(self),
            file_checksum(self.module.extension_path),
        ),
        {}
    )


@overload
def simulate(
    self: AmiciModel,
    *,
    solver: Solver | None = None,
    edata: AmiciExpData | None = None,
    sensi_method: SensitivityMethod | str = None,
    sensi_order: SensitivityOrder | str = None,
) -> ReturnDataView: ...


@overload
def simulate(
    self: AmiciModel,
    *,
    solver: Solver | None = None,
    edata: AmiciExpDataVector | None = None,
    failfast: bool = True,
    num_threads: int = 1,
    sensi_method: SensitivityMethod | str = None,
    sensi_order: SensitivityOrder | str = None,
) -> list[ReturnDataView]: ...


def simulate(
    self: AmiciModel,
    *,
    solver: Solver | None = None,
    edata: AmiciExpData | AmiciExpDataVector | None = None,
    failfast: bool = True,
    num_threads: int = 1,
    sensi_method: SensitivityMethod | str = None,
    sensi_order: SensitivityOrder | str = None,
) -> ReturnDataView | list[ReturnDataView]:
    """Simulate model with given solver and experimental data.

    :param solver:
        Solver to use for simulation. Defaults to :meth:`Model.get_solver`.
    :param edata:
        Experimental data to use for simulation.
        A single :class:`ExpData` instance or a sequence of such instances.
        If `None`, no experimental data is used and the model is simulated
        as is.
    :param sensi_method:
        Sensitivity method to use for simulation.
        If `None`, the solver's current sensitivity method is used.
    :param sensi_order:
        Sensitivity order to use for simulation.
        If `None`, the solvers's current sensitivity order is used.
    :param failfast:
        Whether to stop simulations on first failure.
        Only relevant if `edata` is a sequence of :class:`ExpData` instances.
    :param num_threads:
        Number of threads to use for simulation.
        Only relevant if AMICI was compiled with OpenMP support and if `edata`
        is a sequence of :class:`ExpData` instances.
    :return:
        A single :class:`ReturnDataView` instance containing the simulation
        results if `edata` is a single :class:`ExpData` instance or `None`.
        If `edata` is a sequence of :class:`ExpData` instances, a list of
        :class:`ReturnDataView` instances is returned.
    """
    from amici.sim.sundials._swig_wrappers import _Model__simulate

    return _Model__simulate(
        self,
        solver=solver,
        edata=edata,
        failfast=failfast,
        num_threads=num_threads,
        sensi_method=sensi_method,
        sensi_order=sensi_order,
    )
%}
};

// Process symbols in header
%include "amici/model.h"
