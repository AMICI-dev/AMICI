from __future__ import annotations

import numbers
from typing import Any

import numpy as np
from petab.v1 import LIN, LOG, LOG10

from amici._installation.amici import ParameterScaling


def petab_to_amici_scale(petab_scale: str) -> int:
    """Convert petab scale id to amici scale id."""
    if petab_scale == LIN:
        return ParameterScaling.none
    if petab_scale == LOG10:
        return ParameterScaling.log10
    if petab_scale == LOG:
        return ParameterScaling.ln
    raise ValueError(f"PEtab scale not recognized: {petab_scale}")


def amici_to_petab_scale(amici_scale: int) -> str:
    """Convert amici scale id to petab scale id."""
    if amici_scale == ParameterScaling.none:
        return LIN
    if amici_scale == ParameterScaling.log10:
        return LOG10
    if amici_scale == ParameterScaling.ln:
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
    value_dict: dict[Any, numbers.Number], petab_scale_dict: dict[Any, str]
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
    value_dict: dict[Any, numbers.Number], petab_scale_dict: dict[Any, str]
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
