"""Scales of PEtab v2 problem parameters.

Only code that is independent of JAX or SUNDIALS objects.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING, Literal, get_args

if TYPE_CHECKING:
    import petab.v2 as petabv2

#: Scale of a PEtab problem parameter: linear, natural logarithm, or
#: decadic logarithm.
ParameterScale = Literal["lin", "log", "log10"]


def get_parameter_scales(
    petab_problem: petabv2.Problem,
    parameter_scales: Mapping[str, ParameterScale] | None,
) -> dict[str, ParameterScale]:
    """Get the scales of all problem parameters.

    PEtab v2 does not specify on which scale parameters are estimated, so
    the scales are chosen by the user.

    :param petab_problem:
        The PEtab v2 problem.
    :param parameter_scales:
        The user-provided scales of (a subset of) the problem parameters.
        Parameters not included are on linear scale.
    :return:
        The scales of all problem parameters, in the order of
        ``Problem.x_ids``.
    :raises ValueError:
        If ``parameter_scales`` contains IDs that are not problem parameters,
        or invalid scales.
    """
    parameter_scales = dict(parameter_scales or {})
    if unknown := set(parameter_scales) - set(petab_problem.x_ids):
        raise ValueError(
            "Parameter scales were provided for parameters that are not "
            f"PEtab problem parameters: {sorted(unknown)}"
        )
    valid_scales = get_args(ParameterScale)
    for par_id, scale in parameter_scales.items():
        if scale not in valid_scales:
            raise ValueError(
                f"Invalid scale {scale!r} for parameter {par_id!r}. "
                f"Must be one of {', '.join(map(repr, valid_scales))}."
            )
    return {
        par_id: parameter_scales.get(par_id, "lin")
        for par_id in petab_problem.x_ids
    }
