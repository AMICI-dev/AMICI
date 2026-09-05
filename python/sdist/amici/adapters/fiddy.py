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


def run_simulation_to_function_and_derivative(
    amici_model: AmiciModel,
    *,
    cache: bool = True,
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
    `sllh`, ..., with the parameter axis moved last via
    :func:`_rdata_array_transpose`, and already sliced down to just
    `free_parameter_ids`, in that order -- AMICI's own sensitivity arrays
    are w.r.t. `amici_model.get_free_parameter_ids()`, which need not be
    the same set/order as `free_parameter_ids`, so this slicing happens
    once here rather than requiring every caller to redo it). fiddy's own
    :class:`fiddy.Function`/:func:`fiddy.check_jacobian` handle flattening
    and unbundling a dict-returning function internally -- no manual
    concatenation or index bookkeeping needed here.

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
        The IDs that correspond to the values in the free parameter vector that is
        simulated.
    :param cache:
        Whether to cache the function calls.
    :returns: A tuple of `(function, derivative)`.
    """
    if amici_solver is None:
        amici_solver = amici_model.create_solver()
    if free_parameter_ids is None:
        free_parameter_ids = amici_model.get_free_parameter_ids()
    if amici_edata is not None and amici_edata.free_parameters is not None:
        raise NotImplementedError(
            "Customization of parameter values inside AMICI ExpData."
        )
    chosen_derivatives = default_derivatives
    if derivative_variables is not None:
        chosen_derivatives = {
            k: all_rdata_derivatives[k] for k in derivative_variables
        }
    # AMICI's own sensitivity arrays are w.r.t. `amici_model`'s full free
    # parameter vector, which need not match `free_parameter_ids` (subset
    # and/or order) -- slice/reorder to `free_parameter_ids` once here.
    amici_free_parameter_ids = amici_model.get_free_parameter_ids()
    parameter_indices = [
        amici_free_parameter_ids.index(parameter_id)
        for parameter_id in free_parameter_ids
    ]

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
        return {
            variable: np.asarray(getattr(rdata, variable), dtype=float)
            for variable in chosen_derivatives
        }

    def derivative(point: Type.POINT) -> dict[str, np.ndarray]:
        rdata = run_amici_simulation(point=point, order=SensitivityOrder.first)
        return {
            variable: _rdata_array_transpose(
                array=np.asarray(
                    getattr(rdata, derivative_variable), dtype=float
                ),
                variable=derivative_variable,
            )[..., parameter_indices]
            for variable, derivative_variable in chosen_derivatives.items()
        }

    if cache:
        # Only `function` -- the one fiddy's own FD engine calls, and
        # calls repeatedly at the same point via its own caching-aware
        # batch dispatch -- benefits from this. `derivative` is called at
        # most a handful of times, each at a different (jittered) point,
        # so caching it has no practical benefit; worse, `CachedFunction`
        # is a `fiddy.Function` subclass, which always flattens a dict
        # return into a flat array -- silently breaking `derivative`'s
        # dict-shaped return for any caller expecting it back untouched
        # (e.g. `fiddy.check_jacobian`'s `expected` argument).
        function = CachedFunction(function)

    return function, derivative


def simulate_petab_to_function_and_derivative(
    petab_problem: petab.Problem,
    *,
    amici_model: Model,
    free_parameter_ids: list[str] = None,
    cache: bool = True,
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
        # Only `function` -- the one fiddy's own FD engine calls
        # repeatedly -- benefits from caching. `derivative` is called at
        # most a handful of times, each at a different (jittered) point,
        # so caching it has no practical benefit; also avoids relying on
        # `CachedFunction` (a `fiddy.Function` subclass, which always
        # flattens a dict return into a flat array) for a function whose
        # return shape a caller expects back untouched.
        function = CachedFunction(function)

    return function, derivative


def simulate_petab_v2_to_function_and_derivative(
    petab_simulator: PetabSimulator,
    *,
    free_parameter_ids: list[str] = None,
    cache: bool = True,
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
        # Only `function` -- the one fiddy's own FD engine calls
        # repeatedly -- benefits from caching. `derivative` is called at
        # most a handful of times, each at a different (jittered) point,
        # so caching it has no practical benefit; also avoids relying on
        # `CachedFunction` (a `fiddy.Function` subclass, which always
        # flattens a dict return into a flat array) for a function whose
        # return shape a caller expects back untouched.
        function = CachedFunction(function)

    return function, derivative
