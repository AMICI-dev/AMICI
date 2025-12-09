"""AMICI-generated module for model model_steadystate_py"""

import sys
import warnings
from pathlib import Path

import amici

# Ensure we are binary-compatible, see #556
if "1.0.0.dev" != amici.__version__:
    raise amici.AmiciVersionError(
        f"Cannot use model `model_steadystate_py` in {Path(__file__).parent}, "
        "generated with amici==1.0.0.dev, "
        f"together with amici=={amici.__version__} "
        "which is currently installed. To use this model, install "
        "amici==1.0.0.dev or re-import the model with the amici "
        "version currently installed."
    )

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

    model_steadystate_py = amici._module_from_path(
        "model_steadystate_py.model_steadystate_py",
        Path(__file__).parent / "model_steadystate_py.py",
    )

for var in dir(model_steadystate_py):
    if not var.startswith("_"):
        globals()[var] = getattr(model_steadystate_py, var)

try:
    # _self: this module; will be set during import
    #  via amici.import_model_module
    model_steadystate_py._model_module = _self  # noqa: F821
except NameError:
    # when the model package is imported via `import`
    model_steadystate_py._model_module = sys.modules[__name__]

__version__ = "0.1.0"
