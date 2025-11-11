"""Unified model import interface for AMICI.

This subpackage aggregates functionality for importing models from various
formats into AMICI.

Re-exported symbols:
- SBML:    `SbmlImporter`
- PySB:    `pysb2amici`, `pysb2jax`, `ode_model_from_pysb_importer`
- Antimony:`antimony2sbml`, `antimony2amici`
- Utils:   `MeasurementChannel`, noise helpers, `RESERVED_SYMBOLS`,
           `ObservableTransformation`

All imports are unconditional and will raise ImportError if a required
backend is missing.
"""

from __future__ import annotations

from .antimony import antimony2amici, antimony2sbml
from .pysb import ode_model_from_pysb_importer, pysb2amici, pysb2jax

# Unconditional imports of submodules / symbols
from .sbml import SbmlImporter
from .utils import (
    RESERVED_SYMBOLS,
    MeasurementChannel,
    ObservableTransformation,
    noise_distribution_to_cost_function,
    noise_distribution_to_observable_transformation,
)

__all__ = [
    # SBML
    "SbmlImporter",
    # PySB
    "pysb2amici",
    "pysb2jax",
    "ode_model_from_pysb_importer",
    # Antimony
    "antimony2sbml",
    "antimony2amici",
    # Utils
    "MeasurementChannel",
    "noise_distribution_to_cost_function",
    "noise_distribution_to_observable_transformation",
    "RESERVED_SYMBOLS",
    "ObservableTransformation",
]
