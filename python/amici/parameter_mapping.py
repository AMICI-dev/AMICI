"""
Parameter mapping
-----------------

When performing parameter inference, often parameters need to be mapped from
simulation to estimation parameters, and parameters can differ between
conditions. This can be handled using the `ParameterMapping`.

Note
~~~~

While the parameter mapping can be used direcly with AMICI, it was developed
for usage together with PEtab, for which the whole worklow of generating
the mapping is automatized.
"""


import numbers
from typing import Any, Dict, List, Union
from collections.abc import Sequence

import amici
import numpy as np
from petab.C import *  # noqa: F403


SingleParameterMapping = Dict[str, Union[numbers.Number, str]]
SingleScaleMapping = Dict[str, str]
AmiciModel = Union[amici.Model, amici.ModelPtr]


class ParameterMappingForCondition:
    """Parameter mapping for condition.

    Contains mappings for free parameters, fixed parameters, and fixed
    preequilibration parameters, both for parameters and scales.

    In the scale mappings, for each simulation parameter the scale
    on which the value is passed (and potentially gradients are to be
    returned) is given. In the parameter mappings, for each simulation
    parameter a corresponding optimization parameter (or a numeric value)
    is given.

    If a mapping is not passed, the parameter mappings are assumed to be empty,
    and if a scale mapping is not passed, all scales are set to linear.

    :param map_sim_var:
        Mapping for free simulation parameters.
    :param scale_map_sim_var:
        Scales for free simulation parameters.
    :param map_preeq_fix:
        Mapping for fixed preequilibration parameters.
    :param scale_map_preeq_fix:
        Scales for fixed preequilibration parameters.
    :param map_sim_fix:
        Mapping for fixed simulation parameters.
    :param scale_map_sim_fix:
        Scales for fixed simulation parameters.
    """

    def __init__(
            self,
            map_sim_var: SingleParameterMapping = None,
            scale_map_sim_var: SingleScaleMapping = None,
            map_preeq_fix: SingleParameterMapping = None,
            scale_map_preeq_fix: SingleScaleMapping = None,
            map_sim_fix: SingleParameterMapping = None,
            scale_map_sim_fix: SingleScaleMapping = None,
    ):
        if map_sim_var is None:
            map_sim_var = {}
        self.map_sim_var = map_sim_var

        if scale_map_sim_var is None:
            scale_map_sim_var = {key: LIN for key in map_sim_var}
        self.scale_map_sim_var = scale_map_sim_var

        if map_preeq_fix is None:
            map_preeq_fix = {}
        self.map_preeq_fix = map_preeq_fix

        if scale_map_preeq_fix is None:
            scale_map_preeq_fix = {key: LIN for key in map_preeq_fix}
        self.scale_map_preeq_fix = scale_map_preeq_fix

        if map_sim_fix is None:
            map_sim_fix = {}
        self.map_sim_fix = map_sim_fix

        if scale_map_sim_fix is None:
            scale_map_sim_fix = {key: LIN for key in map_sim_fix}
        self.scale_map_sim_fix = scale_map_sim_fix


class ParameterMapping(Sequence):
    """Parameter mapping for multiple conditions.

    This can be used like a list of `ParameterMappingForCondition`s.

    :param parameter_mappings:
        List of parameter mappings for specific conditions.
    """

    def __init__(
            self,
            parameter_mappings: List[ParameterMappingForCondition] = None
    ):
        super().__init__()
        if parameter_mappings is None:
            parameter_mappings = []
        self.parameter_mappings = parameter_mappings

    def __iter__(self):
        for mapping in self.parameter_mappings:
            yield mapping

    def __getitem__(self, item):
        return self.parameter_mappings[item]

    def __len__(self):
        return len(self.parameter_mappings)

    def append(
            self,
            parameter_mapping_for_condition: ParameterMappingForCondition):
        """Append a condition specific parameter mapping."""
        self.parameter_mappings.append(parameter_mapping_for_condition)


def fill_in_parameters(
        edatas: List[amici.ExpData],
        problem_parameters: Dict[str, numbers.Number],
        scaled_parameters: bool,
        parameter_mapping: ParameterMapping,
        amici_model: AmiciModel) -> None:
    """Fill fixed and dynamic parameters into the edatas (in-place).

    :param edatas:
        List of experimental datas :class:`amici.amici.ExpData` with
        everything except parameters filled.
    :param problem_parameters:
        Problem parameters as parameterId=>value dict. Only
        parameters included here will be set. Remaining parameters will
        be used as currently set in `amici_model`.
    :param scaled_parameters:
        If True, problem_parameters are assumed to be on the scale provided
        in the parameter mapping. If False, they are assumed
        to be in linear scale.
    :param parameter_mapping:
        Parameter mapping for all conditions.
    :param amici_model:
        AMICI model.
    """
    for edata, mapping_for_condition in zip(edatas, parameter_mapping):
        fill_in_parameters_for_condition(
            edata, problem_parameters, scaled_parameters,
            mapping_for_condition, amici_model)


def fill_in_parameters_for_condition(
        edata: amici.ExpData,
        problem_parameters: Dict[str, numbers.Number],
        scaled_parameters: bool,
        parameter_mapping: ParameterMappingForCondition,
        amici_model: AmiciModel) -> None:
    """Fill fixed and dynamic parameters into the edata for condition
    (in-place).

    :param edata:
        Experimental data object to fill parameters into.
    :param problem_parameters:
        Problem parameters as parameterId=>value dict. Only
        parameters included here will be set. Remaining parameters will
        be used as already set in `amici_model` and `edata`.
    :param scaled_parameters:
        If True, problem_parameters are assumed to be on the scale provided
        in the parameter mapping. If False, they
        are assumed to be in linear scale.
    :param parameter_mapping:
        Parameter mapping for current condition.
    :param amici_model:
        AMICI model
    """
    map_sim_var = parameter_mapping.map_sim_var
    scale_map_sim_var = parameter_mapping.scale_map_sim_var
    map_preeq_fix = parameter_mapping.map_preeq_fix
    scale_map_preeq_fix = parameter_mapping.scale_map_preeq_fix
    map_sim_fix = parameter_mapping.map_sim_fix
    scale_map_sim_fix = parameter_mapping.scale_map_sim_fix

    # Parameter mapping may contain parameter_ids as values, these *must*
    # be replaced

    def _get_par(model_par, value):
        """Replace parameter IDs in mapping dicts by values from
        problem_parameters where necessary"""
        if isinstance(value, str):
            # estimated parameter
            # (condition table overrides must have been handled already,
            # e.g. by the PEtab parameter mapping)
            return problem_parameters[value]
        if model_par in problem_parameters:
            # user-provided
            return problem_parameters[model_par]
        # constant value
        return value

    map_preeq_fix = {key: _get_par(key, val)
                     for key, val in map_preeq_fix.items()}
    map_sim_fix = {key: _get_par(key, val)
                   for key, val in map_sim_fix.items()}
    map_sim_var = {key: _get_par(key, val)
                   for key, val in map_sim_var.items()}

    # If necessary, (un)scale parameters
    if scaled_parameters:
        unscale_parameters_dict(map_preeq_fix, scale_map_preeq_fix)
        unscale_parameters_dict(map_sim_fix, scale_map_sim_fix)
    if not scaled_parameters:
        # We scale all parameters to the scale they are estimated on, and pass
        # that information to amici via edata.{parameters,pscale}.
        # The scaling is necessary to obtain correct derivatives.
        scale_parameters_dict(map_sim_var, scale_map_sim_var)
        # We can skip preequilibration parameters, because they are identical
        # with simulation parameters, and only the latter are used from here
        # on.

    ##########################################################################
    # variable parameters and parameter scale

    # parameter list from mapping dict
    parameters = [map_sim_var[par_id]
                  for par_id in amici_model.getParameterIds()]

    # scales list from mapping dict
    scales = [petab_to_amici_scale(scale_map_sim_var[par_id])
              for par_id in amici_model.getParameterIds()]

    if parameters:
        edata.parameters = parameters

    if scales:
        edata.pscale = amici.parameterScalingFromIntVector(scales)

    ##########################################################################
    # fixed parameters preequilibration
    if map_preeq_fix:
        fixed_pars_preeq = [map_preeq_fix[par_id]
                            for par_id in amici_model.getFixedParameterIds()]
        edata.fixedParametersPreequilibration = fixed_pars_preeq

    ##########################################################################
    # fixed parameters simulation
    if map_sim_fix:
        fixed_pars_sim = [map_sim_fix[par_id]
                          for par_id in amici_model.getFixedParameterIds()]
        edata.fixedParameters = fixed_pars_sim


def petab_to_amici_scale(petab_scale: str) -> int:
    """Convert petab scale id to amici scale id."""
    if petab_scale == LIN:
        return amici.ParameterScaling_none
    if petab_scale == LOG10:
        return amici.ParameterScaling_log10
    if petab_scale == LOG:
        return amici.ParameterScaling_ln
    raise ValueError(f"PEtab scale not recognized: {petab_scale}")


def amici_to_petab_scale(amici_scale: int) -> str:
    """Convert amici scale id to petab scale id."""
    if amici_scale == amici.ParameterScaling_none:
        return LIN
    if amici_scale == amici.ParameterScaling_log10:
        return LOG10
    if amici_scale == amici.ParameterScaling_ln:
        return LOG
    raise ValueError(f"AMICI scale not recognized: {amici_scale}")


def scale_parameter(value: numbers.Number,
                    petab_scale: str) -> numbers.Number:
    """Bring parameter from linear scale to target scale.

    :param value:
        Value to scale
    :param petab_scale:
        Target scale of ``value``

    :return:
        ``value`` on target scale
    """
    if petab_scale == LIN:
        return value
    if petab_scale == LOG10:
        return np.log10(value)
    if petab_scale == LOG:
        return np.log(value)
    raise ValueError(f"Unknown parameter scale {petab_scale}. "
                     f"Must be from {(LIN, LOG, LOG10)}")


def unscale_parameter(value: numbers.Number,
                      petab_scale: str) -> numbers.Number:
    """Bring parameter from scale to linear scale.

    :param value:
        Value to scale
    :param petab_scale:
        Target scale of ``value``

    :return:
        ``value`` on linear scale
    """
    if petab_scale == LIN:
        return value
    if petab_scale == LOG10:
        return np.power(10, value)
    if petab_scale == LOG:
        return np.exp(value)
    raise ValueError(f"Unknown parameter scale {petab_scale}. "
                     f"Must be from {(LIN, LOG, LOG10)}")


def scale_parameters_dict(
        value_dict: Dict[Any, numbers.Number],
        petab_scale_dict: Dict[Any, str]) -> None:
    """
    Bring parameters from linear scale to target scale.

    Bring values in ``value_dict`` from linear scale to the scale
    provided in ``petab_scale_dict`` (in-place).
    Both arguments are expected to have the same length and matching keys.

    :param value_dict:
        Values to scale

    :param petab_scale_dict:
        Target scales of ``values``
    """
    if not value_dict.keys() == petab_scale_dict.keys():
        raise AssertionError("Keys don't match.")

    for key, value in value_dict.items():
        value_dict[key] = scale_parameter(value, petab_scale_dict[key])


def unscale_parameters_dict(
        value_dict: Dict[Any, numbers.Number],
        petab_scale_dict: Dict[Any, str]) -> None:
    """
    Bring parameters from target scale to linear scale.

    Bring values in ``value_dict`` from linear scale to the scale
    provided in ``petab_scale_dict`` (in-place).
    Both arguments are expected to have the same length and matching keys.

    :param value_dict:
        Values to scale

    :param petab_scale_dict:
        Target scales of ``values``
    """
    if not value_dict.keys() == petab_scale_dict.keys():
        raise AssertionError("Keys don't match.")

    for key, value in value_dict.items():
        value_dict[key] = unscale_parameter(value, petab_scale_dict[key])
