"""
JAX
---

This module provides an interface to generate and use AMICI models with JAX. Please note that this module is
experimental, the API may substantially change in the future. Use at your own risk and do not expect backward
compatibility.
"""

from warnings import warn

from amici.jax.model import JAXModel
from amici.jax.petab import (
    JAXProblem,
    ReturnValue,
    petab_simulate,
    run_simulations,
)

warn(
    "The JAX module is experimental and the API may change in the future.",
    ImportWarning,
    stacklevel=2,
)

__all__ = [
    "JAXModel",
    "JAXProblem",
    "run_simulations",
    "petab_simulate",
    "ReturnValue",
]
