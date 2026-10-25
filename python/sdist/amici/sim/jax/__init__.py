"""
Functionality for simulating JAX-based AMICI models.

This module provides an interface to generate and use AMICI models with JAX.
Please note that this module is experimental, the API may substantially change
in the future. Use at your own risk and do not expect backward compatibility.
"""

from .model import JAXModel
from .petab import (
    JAXProblem,
    ReturnValue,
    petab_simulate,
    run_simulations,
)

__all__ = [
    "JAXModel",
    "JAXProblem",
    "run_simulations",
    "petab_simulate",
    "ReturnValue",
]
