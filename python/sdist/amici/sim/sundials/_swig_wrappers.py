"""
Convenience wrappers for the swig interface.

While this functionality could be implemented in the swig interface
directly, it is easier to maintain and extend in pure Python.
"""

from __future__ import annotations

import contextlib
import logging
import warnings
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import amici
import amici._installation.amici as amici_swig
from amici._installation.amici import (
    AmiciExpData,
    AmiciExpDataVector,
    AmiciModel,
    AmiciSolver,
    RDataReporting,
    SensitivityMethod,
    SensitivityOrder,
    Solver,
    _get_ptr,
)
from amici.logging import get_logger

from . import ReturnDataView

logger = get_logger(__name__, log_level=logging.DEBUG)


__all__ = [
    "run_simulation",
    "run_simulations",
    "read_solver_settings_from_hdf5",
    "write_solver_settings_to_hdf5",
    "set_model_settings",
    "get_model_settings",
]


def run_simulation(
    model: AmiciModel,
    solver: AmiciSolver,
    edata: AmiciExpData | None = None,
) -> ReturnDataView:
    """
    Simulate a model with given solver and experimental data.

    :param model:
        Model instance

    :param solver:
        Solver instance, must be generated from
        :py:meth:`amici.amici.Model.create_solver`

    :param edata:
        ExpData instance (optional)

    :returns:
        ReturnData object with simulation results
    """
    if (
        model.ne > 0
        and solver.get_sensitivity_method()
        == amici_swig.SensitivityMethod.adjoint
        and solver.get_sensitivity_order() == amici_swig.SensitivityOrder.first
    ):
        warnings.warn(
            "Adjoint sensitivity analysis for models with discontinuous right hand sides (events/piecewise functions) has not been thoroughly tested. "
            "Sensitivities might be wrong. Tracked at https://github.com/AMICI-dev/AMICI/issues/18. "
            "Adjoint sensitivity analysis may work if the location of the discontinuity is not parameter-dependent, but we still recommend testing accuracy of gradients.",
            stacklevel=1,
        )

    rdata = amici_swig.run_simulation(
        _get_ptr(solver), _get_ptr(edata), _get_ptr(model)
    )
    _log_simulation(rdata)
    if solver.get_return_data_reporting_mode() == RDataReporting.full:
        _ids_and_names_to_rdata(rdata, model)
    return ReturnDataView(rdata)


def run_simulations(
    model: AmiciModel,
    solver: AmiciSolver,
    edata_list: AmiciExpDataVector,
    failfast: bool = True,
    num_threads: int = 1,
) -> list[ReturnDataView]:
    """
    Convenience wrapper for loops of amici.runAmiciSimulation

    :param model: Model instance
    :param solver: Solver instance, must be generated from Model.getSolver()
    :param edata_list: list of ExpData instances
    :param failfast: returns as soon as an integration failure is encountered
    :param num_threads: number of threads to use (only used if compiled
        with openmp)

    :returns: list of simulation results
    """
    if (
        model.ne > 0
        and solver.get_sensitivity_method()
        == amici_swig.SensitivityMethod.adjoint
        and solver.get_sensitivity_order() == amici_swig.SensitivityOrder.first
    ):
        warnings.warn(
            "Adjoint sensitivity analysis for models with discontinuous right hand sides (events/piecewise functions) has not been thoroughly tested. "
            "Sensitivities might be wrong. Tracked at https://github.com/AMICI-dev/AMICI/issues/18. "
            "Adjoint sensitivity analysis may work if the location of the discontinuity is not parameter-dependent, but we still recommend testing accuracy of gradients.",
            stacklevel=1,
        )

    edata_ptr_vector = amici_swig.ExpDataPtrVector(edata_list)
    rdata_ptr_list = amici_swig.run_simulations(
        _get_ptr(solver),
        edata_ptr_vector,
        _get_ptr(model),
        failfast,
        num_threads,
    )
    for rdata in rdata_ptr_list:
        _log_simulation(rdata)
        if solver.get_return_data_reporting_mode() == RDataReporting.full:
            _ids_and_names_to_rdata(rdata, model)

    return [ReturnDataView(r) for r in rdata_ptr_list]


def read_solver_settings_from_hdf5(
    file: str, solver: AmiciSolver, location: str | None = "solverSettings"
) -> None:
    """
    Apply solver settings from an HDF5 file to a Solver instance.

    :param file: hdf5 filename
    :param solver: Solver instance to which settings will be transferred
    :param location: location of solver settings in hdf5 file
    """
    amici_swig.read_solver_settings_from_hdf5(file, _get_ptr(solver), location)


def write_solver_settings_to_hdf5(
    solver: AmiciSolver,
    file: str | object,
    location: str | None = "solverSettings",
) -> None:
    """
    Write solver settings from a Solver instance to an HDF5 file.

    :param file: hdf5 filename
    :param solver: Solver instance from which settings will be stored
    :param location: location of solver settings in hdf5 file
    """
    amici_swig.write_solver_settings_to_hdf5(_get_ptr(solver), file, location)


# Values are suffixes of `get[...]` and `set[...]` `amici.Model` methods.
# If either the getter or setter is not named with this pattern, then the value
# is a tuple where the first and second elements are the getter and setter
# methods, respectively.
model_instance_settings = [
    # `set_parameter_{list,scale}` will clear initial state sensitivities, so
    #  `set_parameter_{list,scale}` has to be called first.
    "parameter_list",
    "parameter_scale",  # getter returns a SWIG object
    "add_sigma_residuals",
    "always_check_finite",
    "fixed_parameters",
    "initial_state",
    (
        "get_initial_state_sensitivities",
        "set_unscaled_initial_state_sensitivities",
    ),
    "minimum_sigma_residuals",
    ("n_max_event", "set_n_max_event"),
    "free_parameters",
    "reinitialization_state_idxs",
    "reinitialize_fixed_parameter_initial_states",
    "state_is_non_negative",
    "steady_state_computation_mode",
    "steady_state_sensitivity_mode",
    ("t0", "set_t0"),
    ("t0_preeq", "set_t0_preeq"),
    "timepoints",
    "steadystate_mask",
]


def get_model_settings(
    model: AmiciModel,
) -> dict[str, Any]:
    """Get model settings that are set independently of the compiled model.

    This function is used for serialization of model settings.
    It must only be used in combination with :func:`set_model_settings`.

    :param model: The AMICI model instance.

    :returns:
        Value to be passed to :func:`set_model_settings` to restore the model
        settings.
        Keys are AMICI model attributes, values are attribute values.
    """
    settings = {}
    for setting in model_instance_settings:
        getter = setting[0] if isinstance(setting, tuple) else f"get_{setting}"

        if (
            getter == "get_initial_state"
            and not model.has_custom_initial_state()
        ):
            settings[setting] = []
            continue
        if (
            getter == "get_initial_state_sensitivities"
            and not model.has_custom_initial_state_sensitivities()
        ):
            settings[setting] = []
            continue

        settings[setting] = getattr(model, getter)()
        # TODO `amici.Model.get_parameter_scale` returns a SWIG object instead
        # of a Python list/tuple.
        if setting == "parameter_scale":
            settings[setting] = tuple(settings[setting])
    return settings


def set_model_settings(
    model: AmiciModel,
    settings: dict[str, Any],
) -> None:
    """Set model settings.

    This function is used for deserialization of model settings.
    It must only be used in combination with :func:`get_model_settings`.

    :param model: The AMICI model instance.
    :param settings: The return value of :func:`get_model_settings`.
        Keys are callable attributes (setters) of an AMICI model,
        values are provided to the setters.
    """
    for setting, value in settings.items():
        setter = setting[1] if isinstance(setting, tuple) else f"set_{setting}"
        getattr(model, setter)(value)


def _log_simulation(rdata: amici_swig.ReturnData):
    """Extension warnings to Python logging."""
    amici_severity_to_logging = {
        amici_swig.LogSeverity_debug: logging.DEBUG,
        amici_swig.LogSeverity_warning: logging.WARNING,
        amici_swig.LogSeverity_error: logging.ERROR,
    }
    for msg in rdata.messages:
        condition = f"[{rdata.id}]" if rdata.id else ""
        logger.log(
            amici_severity_to_logging[msg.severity],
            f"{condition}[{msg.identifier}] {msg.message}",
        )


def _ids_and_names_to_rdata(
    rdata: amici_swig.ReturnData, model: amici_swig.Model
):
    """Copy entity IDs and names from a Model to ReturnData."""
    for entity_type in (
        "state",
        "observable",
        "expression",
        "free_parameter",
        "fixed_parameter",
    ):
        for name_or_id in ("ids", "names"):
            names_or_ids = getattr(model, f"get_{entity_type}_{name_or_id}")()
            setattr(
                rdata,
                f"{entity_type.lower()}_{name_or_id.lower()}",
                names_or_ids,
            )

    rdata.state_ids_solver = model.get_state_ids_solver()
    rdata.state_names_solver = model.get_state_names_solver()


@contextlib.contextmanager
def _solver_settings(solver, sensi_method=None, sensi_order=None):
    """Context manager to temporarily apply solver settings."""
    old_method = old_order = None

    if sensi_method is not None:
        old_method = solver.get_sensitivity_method()
        if isinstance(sensi_method, str):
            sensi_method = SensitivityMethod[sensi_method]
        solver.set_sensitivity_method(sensi_method)

    if sensi_order is not None:
        old_order = solver.get_sensitivity_order()
        if isinstance(sensi_order, str):
            sensi_order = SensitivityOrder[sensi_order]
        solver.set_sensitivity_order(sensi_order)

    try:
        yield solver
    finally:
        if old_method is not None:
            solver.set_sensitivity_method(old_method)
        if old_order is not None:
            solver.set_sensitivity_order(old_order)


def _Model__simulate(
    self: AmiciModel,
    *,
    solver: Solver | None = None,
    edata: AmiciExpData | AmiciExpDataVector | None = None,
    failfast: bool = True,
    num_threads: int = 1,
    sensi_method: SensitivityMethod | str = None,
    sensi_order: SensitivityOrder | str = None,
) -> ReturnDataView | list[ReturnDataView]:
    """
    For use in `swig/model.i` to avoid code duplication in subclasses.

    Keep in sync with `Model.simulate` and `ModelPtr.simulate`.

    """
    if solver is None:
        solver = self.create_solver()

    with _solver_settings(
        solver=solver, sensi_method=sensi_method, sensi_order=sensi_order
    ):
        if isinstance(edata, Sequence):
            return run_simulations(
                model=_get_ptr(self),
                solver=_get_ptr(solver),
                edata_list=edata,
                failfast=failfast,
                num_threads=num_threads,
            )

        return run_simulation(
            model=_get_ptr(self),
            solver=_get_ptr(solver),
            edata=_get_ptr(edata),
        )


def restore_model(
    module_name: str, module_path: Path, settings: dict, checksum: str = None
) -> amici.Model:
    """
    Recreate a model instance with given settings.

    For use in ModelPtr.__reduce__.

    :param module_name:
        Name of the model module.
    :param module_path:
        Path to the model module.
    :param settings:
        Model settings to be applied.
        See `set_model_settings` / `get_model_settings`.
    :param checksum:
        Checksum of the model extension to verify integrity.
    """
    from amici import import_model_module

    model_module = import_model_module(module_name, module_path)
    model = model_module.get_model()
    model.module = model_module._self
    set_model_settings(model, settings)

    if checksum is not None and checksum != file_checksum(
        model.module.extension_path
    ):
        raise RuntimeError(
            f"Model file checksum does not match the expected checksum "
            f"({checksum}). The model file may have been modified "
            f"after the model was pickled."
        )

    return model


def file_checksum(
    path: str | Path, algorithm: str = "sha256", chunk_size: int = 8192
) -> str:
    """
    Compute checksum for `path` using `algorithm` (e.g. 'md5', 'sha1', 'sha256').
    Returns the hexadecimal digest string.
    """
    import hashlib

    path = Path(path)
    h = hashlib.new(algorithm)
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def restore_edata(
    init_args: Sequence,
    simulation_parameter_dict: dict[str, Any],
) -> amici_swig.ExpData:
    """
    Recreate an ExpData instance.

    For use in ExpData.__reduce__.
    """
    edata = amici_swig.ExpData(*init_args)

    edata.pscale = amici_swig.parameter_scaling_from_int_vector(
        simulation_parameter_dict.pop("pscale")
    )
    for key, value in simulation_parameter_dict.items():
        if key == "timepoints":
            # timepoints are set during ExpData construction
            continue
        assert hasattr(edata, key)
        setattr(edata, key, value)
    return edata


def restore_solver(cls: type, cls_name: str, hdf5_file: str) -> Solver:
    """
    Recreate a Solver or SolverPtr instance from an HDF5 file.

    For use in Solver.__reduce__.

    :param cls:
        Class of the original object ({CVode,IDA}Solver or SolverPtr).
    :param cls_name:
        Name of the (pointed to) solver class ("CVodeSolver" or "IDASolver").
    :param hdf5_file:
        HDF5 file from which to read the solver settings.
    """
    from . import CVodeSolver, IDASolver

    if cls_name == "CVodeSolver":
        solver = CVodeSolver()
    elif cls_name == "IDASolver":
        solver = IDASolver()
    else:
        raise ValueError(f"Unknown solver class name: {cls_name}")

    if not issubclass(cls, Solver):
        solver = cls(solver)
    read_solver_settings_from_hdf5(hdf5_file, solver)
    return solver
