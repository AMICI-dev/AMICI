"""
Parameter mapping
-----------------

When performing parameter inference, often parameters need to be mapped from
simulation to estimation parameters, and parameters can differ between
conditions. This can be handled using the `ParameterMapping`.

Note
~~~~

While the parameter mapping can be used directly with AMICI, it was developed
for usage together with PEtab, for which the whole workflow of generating
the mapping is automatized.
"""
from __future__ import annotations

import numbers
import warnings
from collections.abc import Sequence
from itertools import chain
from typing import Any, Dict, List, Set, Union

import amici
import numpy as np
from petab.C import *  # noqa: F403

SingleParameterMapping = Dict[str, Union[numbers.Number, str]]
SingleScaleMapping = Dict[str, str]

# some extra imports for backward-compatibility
from .petab.conditions import (  # noqa # pylint: disable=unused-import
    fill_in_parameters,
    fill_in_parameters_for_condition,
)


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

    def __repr__(self):
        return (
            f"{self.__class__.__name__}("
            f"map_sim_var={repr(self.map_sim_var)},"
            f"scale_map_sim_var={repr(self.scale_map_sim_var)},"
            f"map_preeq_fix={repr(self.map_preeq_fix)},"
            f"scale_map_preeq_fix={repr(self.scale_map_preeq_fix)},"
            f"map_sim_fix={repr(self.map_sim_fix)},"
            f"scale_map_sim_fix={repr(self.scale_map_sim_fix)})"
        )

    @property
    def free_symbols(self) -> Set[str]:
        """Get IDs of all (symbolic) parameters present in this mapping"""
        return {
            p
            for p in chain(
                self.map_sim_var.values(),
                self.map_preeq_fix.values(),
                self.map_sim_fix.values(),
            )
            if isinstance(p, str)
        }


class ParameterMapping(Sequence):
    r"""Parameter mapping for multiple conditions.

    This can be used like a list of :class:`ParameterMappingForCondition`\ s.

    :param parameter_mappings:
        List of parameter mappings for specific conditions.
    """

    def __init__(
        self, parameter_mappings: List[ParameterMappingForCondition] = None
    ):
        super().__init__()
        if parameter_mappings is None:
            parameter_mappings = []
        self.parameter_mappings = parameter_mappings

    def __iter__(self):
        yield from self.parameter_mappings

    def __getitem__(
        self, item
    ) -> Union[ParameterMapping, ParameterMappingForCondition]:
        result = self.parameter_mappings[item]
        if isinstance(result, ParameterMappingForCondition):
            return result
        return ParameterMapping(result)

    def __len__(self):
        return len(self.parameter_mappings)

    def append(
        self, parameter_mapping_for_condition: ParameterMappingForCondition
    ):
        """Append a condition specific parameter mapping."""
        self.parameter_mappings.append(parameter_mapping_for_condition)

    def __repr__(self):
        return f"{self.__class__.__name__}({repr(self.parameter_mappings)})"

    @property
    def free_symbols(self) -> Set[str]:
        """Get IDs of all (symbolic) parameters present in this mapping"""
        return set.union(*(mapping.free_symbols for mapping in self))


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


def scale_parameter(value: numbers.Number, petab_scale: str) -> numbers.Number:
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
    raise ValueError(
        f"Unknown parameter scale {petab_scale}. "
        f"Must be from {(LIN, LOG, LOG10)}"
    )


def unscale_parameter(
    value: numbers.Number, petab_scale: str
) -> numbers.Number:
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
    raise ValueError(
        f"Unknown parameter scale {petab_scale}. "
        f"Must be from {(LIN, LOG, LOG10)}"
    )


def scale_parameters_dict(
    value_dict: Dict[Any, numbers.Number], petab_scale_dict: Dict[Any, str]
) -> None:
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
    if value_dict.keys() != petab_scale_dict.keys():
        raise AssertionError("Keys don't match.")

    for key, value in value_dict.items():
        value_dict[key] = scale_parameter(value, petab_scale_dict[key])


def unscale_parameters_dict(
    value_dict: Dict[Any, numbers.Number], petab_scale_dict: Dict[Any, str]
) -> None:
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
    if value_dict.keys() != petab_scale_dict.keys():
        raise AssertionError("Keys don't match.")

    for key, value in value_dict.items():
        value_dict[key] = unscale_parameter(value, petab_scale_dict[key])
