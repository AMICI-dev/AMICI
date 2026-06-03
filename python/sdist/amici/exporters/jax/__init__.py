"""
Code generation for JAX models for simulation with diffrax solvers.

This module provides an interface to generate and use AMICI models with JAX.
Please note that this module is experimental, the API may substantially change
in the future. Use at your own risk and do not expect backward compatibility.
"""

from .nn import (
    BatchNorm,
    Flatten,
    InstanceNorm,
    cat,
    generate_equinox,
    tanhshrink,
)

__all__ = [
    "BatchNorm",
    "Flatten",
    "InstanceNorm",
    "generate_equinox",
    "tanhshrink",
    "cat",
]
