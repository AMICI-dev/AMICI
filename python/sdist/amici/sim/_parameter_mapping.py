"""Parameter mapping for conditions.

Only code that is independent of JAX or SUNDIALS objects.
"""

from __future__ import annotations

import numbers
from collections.abc import Sequence
from itertools import chain

from petab.v1 import LIN

SingleParameterMapping = dict[str, numbers.Number | str]
SingleScaleMapping = dict[str, str]


class ParameterMappingForCondition:
    """Parameter mapping for condition.

    Contains mappings for free parameters, fixed parameters, and fixed
    pre-equilibration parameters, both for parameters and scales.

    In the scale mappings, for each simulation parameter, the scale
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
        Mapping for fixed pre-equilibration parameters.
    :param scale_map_preeq_fix:
        Scales for fixed pre-equilibration parameters.
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
    def free_symbols(self) -> set[str]:
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
        self, parameter_mappings: list[ParameterMappingForCondition] = None
    ):
        super().__init__()
        if parameter_mappings is None:
            parameter_mappings = []
        self.parameter_mappings = parameter_mappings

    def __iter__(self):
        yield from self.parameter_mappings

    def __getitem__(
        self, item
    ) -> ParameterMapping | ParameterMappingForCondition:
        result = self.parameter_mappings[item]
        if isinstance(result, ParameterMappingForCondition):
            return result
        return ParameterMapping(result)

    def __len__(self):
        return len(self.parameter_mappings)

    def append(
        self, parameter_mapping_for_condition: ParameterMappingForCondition
    ):
        """Append a condition-specific parameter mapping."""
        self.parameter_mappings.append(parameter_mapping_for_condition)

    def __repr__(self):
        return f"{self.__class__.__name__}({repr(self.parameter_mappings)})"

    @property
    def free_symbols(self) -> set[str]:
        """Get IDs of all (symbolic) parameters present in this mapping"""
        return set.union(*(mapping.free_symbols for mapping in self))
