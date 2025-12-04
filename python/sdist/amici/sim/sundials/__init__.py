"""
Functionality for simulating AMICI models using SUNDIALS solvers.

This module provides the relevant objects for simulating AMICI models
using :term:`SUNDIALS` solvers such as :term:`CVODES` and :term:`IDAS` and
for analyzing the simulation results.
"""

import os
import warnings
from types import ModuleType
from typing import Protocol, runtime_checkable

import amici
from amici import amici_path, import_model_module

#: boolean indicating if this is the full package with swig interface or
#  the raw package without extension
has_clibs: bool = any(
    os.path.isfile(os.path.join(amici_path, "_installation", wrapper))
    for wrapper in ["amici.py", "amici_without_hdf5.py"]
)
#: boolean indicating if amici was compiled with hdf5 support
hdf5_enabled: bool = False

__all__ = [
    "import_model_module",
    "ModelModule",
]

if has_clibs:
    # Import SWIG module and swig-dependent submodules
    #  if required and available

    # prevent segfaults under pytest
    #  see also:
    #  https://github.com/swig/swig/issues/2881
    #  https://github.com/AMICI-dev/AMICI/issues/2565
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            category=DeprecationWarning,
            message="builtin type .* has no __module__ attribute",
        )
        # The swig-generated Python module
        from amici._installation import amici as amici_swig_py

    # TODO: selective import to avoid loading unneeded symbols
    from amici._installation.amici import *
    from amici._installation.amici import _SWIG_VERSION as _SWIG_VERSION

    # has to be done before importing read_solver_settings_from_hdf5
    #  from .swig_wrappers
    hdf5_enabled = "read_solver_settings_from_hdf5" in dir()
    # These modules require the swig interface and other dependencies
    from ._numpy import ExpDataView as ExpDataView
    from ._numpy import ReturnDataView as ReturnDataView
    from ._numpy import evaluate as evaluate

    # Import after import of the swig module as this is supposed to
    #  shadow some swig symbols
    from ._swig_wrappers import *

    @runtime_checkable
    class ModelModule(Protocol):  # noqa: F811
        """Type of AMICI-generated model modules.

        To enable static type checking."""

        def get_model(self) -> amici_swig_py.Model:
            """Create a model instance."""
            ...

    AmiciModel = amici_swig_py.Model | amici_swig_py.ModelPtr

    from . import _swig_wrappers

    __all__ += [
        *amici._installation.amici.__all__,
        "ExpDataView",
        "ReturnDataView",
        "evaluate",
        "ModelModule",
        *_swig_wrappers.__all__,
        "AmiciModel",
    ]

    # expose the swig module itself
    amici = amici_swig_py
else:
    ModelModule = ModuleType
