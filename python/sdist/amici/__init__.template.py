"""AMICI-generated module for model TPL_MODELNAME"""

import sys
from pathlib import Path
import amici
import warnings

# Ensure we are binary-compatible, see #556
if "TPL_AMICI_VERSION" != amici.__version__:
    raise amici.AmiciVersionError(
        f"Cannot use model `TPL_MODELNAME` in {Path(__file__).parent}, "
        "generated with amici==TPL_AMICI_VERSION, "
        f"together with amici=={amici.__version__} "
        "which is currently installed. To use this model, install "
        "amici==TPL_AMICI_VERSION or re-import the model with the amici "
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

    TPL_MODELNAME = amici._module_from_path(
        "TPL_MODELNAME.TPL_MODELNAME",
        Path(__file__).parent / "TPL_MODELNAME.py",
    )

for var in dir(TPL_MODELNAME):
    if not var.startswith("_"):
        globals()[var] = getattr(TPL_MODELNAME, var)
get_model = TPL_MODELNAME.getModel

try:
    # _self: this module; will be set during import
    #  via amici.import_model_module
    TPL_MODELNAME._model_module = _self  # noqa: F821
except NameError:
    # when the model package is imported via `import`
    TPL_MODELNAME._model_module = sys.modules[__name__]

__version__ = "TPL_PACKAGE_VERSION"
