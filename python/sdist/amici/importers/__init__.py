"""Unified model import interface for AMICI.

This subpackage aggregates functionality for importing models from various
formats into AMICI.

Re-exported symbols (imported conditionally where dependencies are optional):
- SBML:    `SbmlImporter` (requires `libsbml` and related deps)
- PySB:    `pysb2amici`, `pysb2jax`, `ode_model_from_pysb_importer` (requires `pysb`)
- Antimony:`antimony2sbml`, `antimony2amici` (requires `antimony` Python package)
- Utils:   `MeasurementChannel`, noise helpers, `RESERVED_SYMBOLS`,
           `ObservableTransformation`
"""

from __future__ import annotations

# Utils are internal and safe to import unconditionally
from .utils import (
    RESERVED_SYMBOLS,
    MeasurementChannel,
    ObservableTransformation,
    noise_distribution_to_cost_function,
    noise_distribution_to_observable_transformation,
)

__all__: list[str] = [
    # Utils
    "MeasurementChannel",
    "noise_distribution_to_cost_function",
    "noise_distribution_to_observable_transformation",
    "RESERVED_SYMBOLS",
    "ObservableTransformation",
]

# SBML importer (optional)
try:
    from .sbml import (
        SbmlImporter,  # noqa: E402,F401
        sbml2jax,  # noqa: E402,F401
        smbl2amici,  # noqa: E402,F401
    )

    __all__ += [
        "sbml2amici",
        "sbml2jax",
        "SbmlImporter",
    ]
except Exception:
    # Leave unavailable if optional dependency is missing
    pass

# PySB importer (optional)
try:
    from .pysb import (
        ode_model_from_pysb_importer,  # noqa: E402,F401
        pysb2amici,  # noqa: E402,F401
        pysb2jax,  # noqa: E402,F401
    )

    __all__ += [
        "pysb2amici",
        "pysb2jax",
        "ode_model_from_pysb_importer",
    ]
except Exception:
    pass

# Antimony importer (optional)
try:
    from .antimony import antimony2amici, antimony2sbml  # noqa: E402,F401

    __all__ += ["antimony2sbml", "antimony2amici"]
except Exception:
    pass

# BNGL importer (optional)
try:
    from .bngl import bngl2amici  # noqa: E402,F401

    __all__.append("bngl2amici")
except Exception:
    pass
