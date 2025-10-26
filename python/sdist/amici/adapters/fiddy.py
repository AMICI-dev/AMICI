"""
Adapters for using AMICI with the `fiddy <https://github.com/ICB-DCM/fiddy/>`__
package for finite difference checks.

NOTE: Like fiddy, this module is experimental and subject to change.
"""

from collections.abc import Callable
from functools import partial
from inspect import signature
from typing import Any

import numpy as np
import petab.v1 as petab
from fiddy import CachedFunction, Type, fiddy_array
from petab.v1.C import LIN, LOG, LOG10

import amici
import amici.petab.simulations
from amici.petab.conditions import create_edatas
from amici.petab.parameter_mapping import create_parameter_mapping
from amici.petab.simulations import LLH, SLLH

from .. import ReturnData, SensitivityOrder

__all__ = [
    "run_simulation_to_cached_functions",
    "simulate_petab_to_cached_functions",
]

LOG_E_10 = np.log(10)


def transform_gradient_lin_to_lin(gradient_value, _):
    return gradient_value


def transform_gradient_lin_to_log(gradient_value, parameter_value):
    return gradient_value * parameter_value


def transform_gradient_lin_to_log10(gradient_value, parameter_value):
    return gradient_value * (parameter_value * LOG_E_10)


transforms = {
    LIN: transform_gradient_lin_to_lin,
    LOG: transform_gradient_lin_to_log,
    LOG10: transform_gradient_lin_to_log10,
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


def rdata_array_transpose(array: np.ndarray, variable: str) -> tuple[int]:
    if array.size == 0:
        return array
    original_parameter_dimension = derivative_parameter_dimension[variable]
    return np.moveaxis(array, original_parameter_dimension, -1)


def fiddy_array_transpose(array: np.ndarray, variable: str) -> tuple[int]:
    if array.size == 0:
        return array
    original_parameter_dimension = derivative_parameter_dimension[variable]
    return np.moveaxis(array, -1, original_parameter_dimension)


default_derivatives = {
    k: v
    for k, v in all_rdata_derivatives.items()
    if v not in ["sz", "srz", "ssigmaz", "s2llh"]
}


def run_simulation_to_cached_functions(
    amici_model: amici.AmiciModel,
    *,
    cache: bool = True,
    parameter_ids: list[str] = None,
    amici_solver: amici.AmiciSolver = None,
    amici_edata: amici.AmiciExpData = None,
    derivative_variables: list[str] = None,
):
    """Convert `amici.run_simulation` to fiddy functions.

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
    :param parameter_ids:
        The IDs that correspond to the values in the parameter vector that is
        simulated.
    :param cache:
        Whether to cache the function calls.
    :returns: function, derivatives and structure
    """
    if amici_solver is None:
        amici_solver = amici_model.create_solver()
    if parameter_ids is None:
        parameter_ids = amici_model.get_parameter_ids()
    if amici_edata is not None and amici_edata.parameters is not None:
        raise NotImplementedError(
            "Customization of parameter values inside AMICI ExpData."
        )
    chosen_derivatives = default_derivatives
    if derivative_variables is not None:
        chosen_derivatives = {
            k: all_rdata_derivatives[k] for k in derivative_variables
        }

    def run_amici_simulation(
        point: Type.POINT, order: SensitivityOrder
    ) -> ReturnData:
        problem_parameters = dict(zip(parameter_ids, point, strict=True))
        amici_model.set_parameter_by_id(problem_parameters)
        amici_solver.set_sensitivity_order(order)
        rdata = amici.run_simulation(
            model=amici_model, solver=amici_solver, edata=amici_edata
        )
        return rdata

    def function(point: Type.POINT):
        rdata = run_amici_simulation(point=point, order=SensitivityOrder.none)
        outputs = {
            variable: fiddy_array(getattr(rdata, variable))
            for variable in chosen_derivatives
        }
        rdata_flat = np.concatenate(
            [output.flat for output in outputs.values()]
        )
        return rdata_flat

    def derivative(point: Type.POINT, return_dict: bool = False):
        rdata = run_amici_simulation(point=point, order=SensitivityOrder.first)
        outputs = {
            variable: rdata_array_transpose(
                array=fiddy_array(getattr(rdata, derivative_variable)),
                variable=derivative_variable,
            )
            for variable, derivative_variable in chosen_derivatives.items()
        }
        rdata_flat = np.concatenate(
            [
                output_array.reshape(-1, output_array.shape[-1])
                for output_array in outputs.values()
            ],
            axis=0,
        )
        if return_dict:
            return outputs
        return rdata_flat

    if cache:
        function = CachedFunction(function)
        derivative = CachedFunction(derivative)

    # Get structure
    dummy_point = fiddy_array(
        [amici_model.get_parameter_by_id(par_id) for par_id in parameter_ids]
    )
    dummy_rdata = run_amici_simulation(
        point=dummy_point, order=SensitivityOrder.first
    )

    structures = {
        "function": {variable: None for variable in chosen_derivatives},
        "derivative": {variable: None for variable in chosen_derivatives},
    }
    function_position = 0
    derivative_position = 0
    for variable, derivative_variable in chosen_derivatives.items():
        function_array = fiddy_array(getattr(dummy_rdata, variable))
        derivative_array = fiddy_array(
            getattr(dummy_rdata, derivative_variable)
        )
        structures["function"][variable] = (
            function_position,
            function_position + function_array.size,
            function_array.shape,
        )
        structures["derivative"][variable] = (
            derivative_position,
            derivative_position + derivative_array.size,
            derivative_array.shape,
        )
        function_position += function_array.size
        derivative_position += derivative_array.size

    return function, derivative, structures


# (start, stop, shape)
TYPE_STRUCTURE = tuple[int, int, tuple[int, ...]]


def flatten(arrays: dict[str, Type.ARRAY]) -> Type.ARRAY:
    flattened_value = np.concatenate([array.flat for array in arrays.values()])
    return flattened_value


def reshape(
    array: Type.ARRAY,
    structure: TYPE_STRUCTURE,
    sensitivities: bool = False,
) -> dict[str, Type.ARRAY]:
    reshaped = {}
    for variable, (start, stop, shape) in structure.items():
        # array is currently "flattened" w.r.t. fiddy dimensions
        # hence, if sensis, reshape w.r.t. fiddy dimensions
        if sensitivities and (
            dimension0 := derivative_parameter_dimension.get(
                "s" + variable, False
            )
        ):
            shape = [
                size
                for dimension, size in enumerate(shape)
                if dimension != dimension0
            ] + [shape[dimension0]]

        array = array[start:stop]
        if array.size != 0:
            array = array.reshape(shape)

        # now reshape to AMICI dimensions
        if sensitivities and (
            derivative_parameter_dimension.get(f"s{variable}", False)
        ):
            array = fiddy_array_transpose(
                array=array,
                variable=f"s{variable}",
            )
        reshaped[variable] = array

    return reshaped


def simulate_petab_to_cached_functions(
    petab_problem: petab.Problem,
    amici_model: amici.Model,
    parameter_ids: list[str] = None,
    cache: bool = True,
    precreate_edatas: bool = True,
    precreate_parameter_mapping: bool = True,
    simulate_petab: Callable[[Any], str] = None,
    **kwargs,
) -> tuple[Type.FUNCTION, Type.FUNCTION]:
    """
    Convert :func:`amici.petab.simulations.simulate_petab`
    (PEtab v1 simulations) to fiddy functions.

    Note that all gradients are provided on linear scale. The correction from
    `'log10'` scale is automatically done.

    :param amici_model:
        The AMICI model to simulate.
    :param simulate_petab:
        A method to simulate PEtab problems with AMICI, e.g.
        `amici.petab_objective.simulate_petab`.
    :param parameter_ids:
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
    if parameter_ids is None:
        parameter_ids = list(petab_problem.parameter_df.index)

    if simulate_petab is None:
        simulate_petab = amici.petab.simulations.simulate_petab

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
        problem_parameters = dict(zip(parameter_ids, point, strict=True))
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
            [result[SLLH][parameter_id] for parameter_id in parameter_ids]
        )
        return sllh

    if cache:
        function = CachedFunction(function)
        derivative = CachedFunction(derivative)

    return function, derivative
