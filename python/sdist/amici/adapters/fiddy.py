"""
Adapters for using AMICI with the `fiddy <https://github.com/ICB-DCM/fiddy/>`__
package for finite difference checks.


.. note::

    Like fiddy, this module is experimental and subject to change.
"""

from __future__ import annotations

from collections.abc import Callable
from functools import partial
from inspect import signature
from typing import TYPE_CHECKING, Any

import numpy as np
import petab.v1 as petab
from fiddy import CachedFunction, Type
from petab.v1.C import LIN, LOG, LOG10

from amici.sim.sundials import (
    AmiciExpData,
    AmiciModel,
    AmiciSolver,
    Model,
    ReturnData,
    SensitivityOrder,
    run_simulation,
)
from amici.sim.sundials.petab.v1 import LLH, SLLH, create_edatas
from amici.sim.sundials.petab.v1._parameter_mapping import (
    create_parameter_mapping,
)

if TYPE_CHECKING:
    from amici.sim.sundials.petab import PetabSimulationResult, PetabSimulator

__all__ = [
    "run_simulation_to_function_and_derivative",
    "simulate_petab_to_function_and_derivative",
    "simulate_petab_v2_to_function_and_derivative",
    "output_labels_for_derivatives",
]

LOG_E_10 = np.log(10)


def _transform_gradient_lin_to_lin(gradient_value, _):
    return gradient_value


def _transform_gradient_lin_to_log(gradient_value, parameter_value):
    return gradient_value * parameter_value


def _transform_gradient_lin_to_log10(gradient_value, parameter_value):
    return gradient_value * (parameter_value * LOG_E_10)


transforms = {
    LIN: _transform_gradient_lin_to_lin,
    LOG: _transform_gradient_lin_to_log,
    LOG10: _transform_gradient_lin_to_log10,
}


all_rdata_derivatives = {
    "x": "sx",
    "x0": "sx0",
    "x_ss": "sx_ss",
    "y": "sy",
    "sigmay": "ssigmay",
    "z": "sz",
    "rz": "srz",
    "sigmaz": "ssigmaz",
    "llh": "sllh",
    "sllh": "s2llh",
    "res": "sres",
}

# The dimension of the AMICI ReturnData that contains parameters.
# Should be shifted to the last dimension to be compatible with fiddy.
derivative_parameter_dimension = {
    "sx": 1,
    "sx0": 0,
    "sx_ss": 0,
    "sy": 1,
    "ssigmay": 1,
    # 'sz'      : ???,
    "srz": 2,
    # 'ssigmaz' : ???,
    "sllh": 0,
    "s2llh": 1,
    "sres": 1,
}


def _rdata_array_transpose(array: np.ndarray, variable: str) -> tuple[int]:
    if array.size == 0:
        return array
    original_parameter_dimension = derivative_parameter_dimension[variable]
    return np.moveaxis(array, original_parameter_dimension, -1)


default_derivatives = {
    k: v
    for k, v in all_rdata_derivatives.items()
    if v not in ["sz", "srz", "ssigmaz", "s2llh"]
}

# Entities to id type mapping
_entity_ids_by_variable = {
    "x": "state",
    "x0": "state",
    "x_ss": "state",
    "y": "observable",
    "sigmay": "observable",
    "res": "observable",
}
# Entities that have a time
_has_timepoint_axis = {"x", "y", "sigmay", "res"}


def output_labels_for_derivatives(
    amici_model: AmiciModel,
    derivative_variables: list[str] = None,
    timepoints: list[float] = None,
) -> list[str]:
    """Per-flat-row labels for fiddy's `function`/`derivative`'s bundled output.

    :param amici_model: The AMICI model (for state/observable IDs).
    :param derivative_variables: Same meaning/default as
        :func:`run_simulation_to_function_and_derivative`.
    :param timepoints: Output timepoints, for variables with a timepoint
        axis. Defaults to `amici_model.get_timepoints()`.
    :return: One label per flat output row, in bundling order.
    :raises NotImplementedError: For a variable with no label source
        (``z``, ``rz``, ``sigmaz``, or second-order ``sllh``).
    """
    variables = list(
        default_derivatives
        if derivative_variables is None
        else derivative_variables
    )
    unsupported = [v for v in variables if v not in default_derivatives]
    if unsupported:
        raise NotImplementedError(
            f"No output labels available for {unsupported} -- only "
            f"{list(default_derivatives)} are supported."
        )
    if timepoints is None:
        timepoints = list(amici_model.get_timepoints())
    ids_by_kind = {
        "state": list(amici_model.get_state_ids()),
        "observable": list(amici_model.get_observable_ids()),
    }

    labels = []
    for variable in variables:
        if variable == "llh":
            labels.append("llh")
            continue
        entity_ids = ids_by_kind[_entity_ids_by_variable[variable]]
        if variable in _has_timepoint_axis:
            labels.extend(
                f"{variable}[t={t:g}, id={entity_id}]"
                for t in timepoints
                for entity_id in entity_ids
            )
        else:
            labels.extend(
                f"{variable}[id={entity_id}]" for entity_id in entity_ids
            )
    return labels


def run_simulation_to_function_and_derivative(
    amici_model: AmiciModel,
    *,
    cache: bool = False,
    free_parameter_ids: list[str] = None,
    amici_solver: AmiciSolver = None,
    amici_edata: AmiciExpData = None,
    derivative_variables: list[str] = None,
):
    """Convert `run_simulation` to a fiddy-checkable ``(function,
    derivative)`` pair, e.g. for :func:`fiddy.check_jacobian`.

    Both `function` and `derivative` return a dict keyed by
    `derivative_variables` (or `default_derivatives`' keys, if not given)
    -- one simulation output per key for `function` (`x`, `y`, `llh`, ...),
    its forward-sensitivity counterpart for `derivative` (`sx`, `sy`,
    `sllh`, ..., with the parameter axis moved last, and sliced/reordered to
    `free_parameter_ids` from each simulation's `rdata.plist`.

    :param amici_model:
        The AMICI model to simulate.
    :param amici_solver:
        The AMICI solver to use. If `None`, a new solver will be created from
        the model.
    :param amici_edata:
        The AMICI ExpData to use. If `None`, no data will be used.
    :param derivative_variables:
        The variables that derivatives will be computed or approximated for.
        See the keys of `all_rdata_derivatives` for options.
    :param free_parameter_ids:
        IDs for the values in the simulated free parameter vector. Each
        must be in the resolved `plist` (`amici_model` or `amici_edata`),
        or `derivative` raises `ValueError`.
    :param cache:
        Whether to cache the function calls.
    :returns: A tuple of `(function, derivative)`.
    """
    if amici_solver is None:
        amici_solver = amici_model.create_solver()
    if free_parameter_ids is None:
        free_parameter_ids = amici_model.get_free_parameter_ids()
    if amici_edata is not None and amici_edata.free_parameters:
        raise NotImplementedError(
            "Customization of parameter values inside AMICI ExpData."
        )
    chosen_derivatives = default_derivatives
    if derivative_variables is not None:
        chosen_derivatives = {
            k: all_rdata_derivatives[k] for k in derivative_variables
        }
    amici_free_parameter_ids = amici_model.get_free_parameter_ids()

    def run_amici_simulation(
        point: Type.POINT, order: SensitivityOrder
    ) -> ReturnData:
        problem_parameters = dict(zip(free_parameter_ids, point, strict=True))
        amici_model.set_free_parameter_by_id(problem_parameters)
        amici_solver.set_sensitivity_order(order)
        rdata = run_simulation(
            model=amici_model, solver=amici_solver, edata=amici_edata
        )
        return rdata

    def function(point: Type.POINT) -> dict[str, np.ndarray]:
        rdata = run_amici_simulation(point=point, order=SensitivityOrder.none)
        outputs = {}
        for variable in chosen_derivatives:
            value = getattr(rdata, variable)
            # AMICI represents a structurally empty field (e.g. `x` for a
            # model with zero states) as `None`, not an empty array --
            # `np.asarray(None, dtype=float)` would silently produce a 0-d
            # NaN scalar instead, which is both the wrong shape and would
            # spuriously fail fiddy's non-finite-value check.
            if value is not None:
                outputs[variable] = np.asarray(value, dtype=float)
        return outputs

    def derivative(point: Type.POINT) -> dict[str, np.ndarray]:
        rdata = run_amici_simulation(point=point, order=SensitivityOrder.first)
        rdata_free_parameter_ids = [
            amici_free_parameter_ids[i] for i in rdata.plist
        ]
        try:
            parameter_indices = [
                rdata_free_parameter_ids.index(parameter_id)
                for parameter_id in free_parameter_ids
            ]
        except ValueError as error:
            raise ValueError(
                f"{error}. `free_parameter_ids` requested a parameter "
                "whose sensitivity was not computed by this simulation "
                "-- check `amici_model.get_parameter_list()` and "
                "`amici_edata.plist` (if `amici_edata` is given, its own "
                "`plist` takes priority over the model's whenever it is "
                "non-empty)."
            ) from error
        outputs = {}
        for variable, derivative_variable in chosen_derivatives.items():
            value = getattr(rdata, derivative_variable)
            if value is not None:  # see `function`'s comment above
                outputs[variable] = _rdata_array_transpose(
                    array=np.asarray(value, dtype=float),
                    variable=derivative_variable,
                )[..., parameter_indices]
        return outputs

    if cache:
        function = CachedFunction(function)

    return function, derivative


def simulate_petab_to_function_and_derivative(
    petab_problem: petab.Problem,
    *,
    amici_model: Model,
    free_parameter_ids: list[str] = None,
    cache: bool = False,
    precreate_edatas: bool = True,
    precreate_parameter_mapping: bool = True,
    simulate_petab: Callable[[Any], str] = None,
    **kwargs,
) -> tuple[Type.FUNCTION, Type.FUNCTION]:
    """
    Convert :func:`amici.sim.sundials.petab.v1.simulate_petab`
    (PEtab v1 simulations) to a fiddy-checkable ``(function, derivative)``
    pair, e.g. for :func:`fiddy.check_gradient`.

    Note that all gradients are provided on linear scale. The correction from
    `'log10'` scale is automatically done.

    :param amici_model:
        The AMICI model to simulate.
    :param simulate_petab:
        A method to simulate PEtab problems with AMICI, e.g.
        `amici.petab_objective.simulate_petab`.
    :param free_parameter_ids:
        The IDs of the parameters, in the order that parameter values will
        be supplied. Defaults to `petab_problem.parameter_df.index`.
    :param petab_problem:
        The PEtab problem.
    :param cache:
        Whether to cache the function call.
    :param precreate_edatas:
        Whether to create the AMICI measurements object in advance, to save
        time.
    :param precreate_parameter_mapping:
        Whether to create the AMICI parameter mapping object in advance, to
        save time.
    :param kwargs:
        Passed to `simulate_petab`.
    :returns:
        A tuple of:

        * 1: A method to compute the function at a point.
        * 2: A method to compute the gradient at a point.
    """
    if free_parameter_ids is None:
        free_parameter_ids = list(petab_problem.parameter_df.index)

    if simulate_petab is None:
        from amici.sim.sundials.petab.v1._simulations import simulate_petab

    edatas = None
    if precreate_edatas:
        edatas = create_edatas(
            amici_model=amici_model,
            petab_problem=petab_problem,
            simulation_conditions=petab_problem.get_simulation_conditions_from_measurement_df(),
        )

    parameter_mapping = None
    if precreate_parameter_mapping:
        parameter_mapping = create_parameter_mapping(
            petab_problem=petab_problem,
            simulation_conditions=petab_problem.get_simulation_conditions_from_measurement_df(),
            scaled_parameters=kwargs.get(
                "scaled_parameters",
                (
                    signature(simulate_petab)
                    .parameters["scaled_parameters"]
                    .default
                ),
            ),
            amici_model=amici_model,
        )

    precreated_kwargs = {
        "edatas": edatas,
        "parameter_mapping": parameter_mapping,
        "petab_problem": petab_problem,
    }
    precreated_kwargs = {
        k: v for k, v in precreated_kwargs.items() if v is not None
    }

    amici_solver = kwargs.pop("solver", amici_model.create_solver())

    simulate_petab_partial = partial(
        simulate_petab,
        amici_model=amici_model,
        **precreated_kwargs,
        **kwargs,
    )

    def simulate_petab_full(point: Type.POINT, order: SensitivityOrder):
        problem_parameters = dict(zip(free_parameter_ids, point, strict=True))
        amici_solver.set_sensitivity_order(order)
        result = simulate_petab_partial(
            problem_parameters=problem_parameters,
            solver=amici_solver,
        )
        return result

    def function(point: Type.POINT):
        output = simulate_petab_full(point, order=SensitivityOrder.none)
        result = output[LLH]
        return np.array(result)

    def derivative(point: Type.POINT) -> Type.POINT:
        result = simulate_petab_full(point, order=SensitivityOrder.first)

        if result[SLLH] is None:
            raise RuntimeError("Simulation failed.")

        sllh = np.array(
            [result[SLLH][parameter_id] for parameter_id in free_parameter_ids]
        )
        return sllh

    if cache:
        function = CachedFunction(function)

    return function, derivative


def simulate_petab_v2_to_function_and_derivative(
    petab_simulator: PetabSimulator,
    *,
    free_parameter_ids: list[str] = None,
    cache: bool = False,
) -> tuple[Type.FUNCTION, Type.FUNCTION]:
    r"""Create a fiddy-checkable ``(function, derivative)`` pair for a
    `PetabSimulator`, e.g. for :func:`fiddy.check_gradient`.

    :param petab_simulator:
        The PEtab simulator to use.
    :param free_parameter_ids:
        The IDs of the parameters, in the order that parameter values will
        be supplied. Defaults to the estimated parameters of the PEtab problem.
    :param cache:
        Whether to cache the function call.
    :returns:
        tuple of:

        * 1: A method to compute the function at a point.
        * 2: A method to compute the gradient at a point.
    """
    if free_parameter_ids is None:
        free_parameter_ids = list(petab_simulator._petab_problem.x_free_ids)

    def simulate(
        point: Type.POINT, order: SensitivityOrder
    ) -> PetabSimulationResult:
        problem_parameters = dict(zip(free_parameter_ids, point, strict=True))
        petab_simulator.solver.set_sensitivity_order(order)

        result = petab_simulator.simulate(
            problem_parameters=problem_parameters,
        )
        return result

    def function(point: Type.POINT) -> np.ndarray:
        output = simulate(point, order=SensitivityOrder.none)
        result = output.llh
        return np.array(result)

    def derivative(point: Type.POINT) -> Type.POINT:
        result = simulate(point, order=SensitivityOrder.first)

        if result.sllh is None:
            raise RuntimeError("Simulation failed.")

        sllh = np.array(
            [result.sllh[parameter_id] for parameter_id in free_parameter_ids]
        )
        return sllh

    if cache:
        function = CachedFunction(function)

    return function, derivative
