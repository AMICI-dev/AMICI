"""AMICI-generated module for model TPL_MODELNAME"""

import datetime
import os
import sys
from pathlib import Path
from typing import TYPE_CHECKING
import amici


if TYPE_CHECKING:
    from amici.jax import JAXModel

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

TPL_MODELNAME = amici._module_from_path(
    "TPL_MODELNAME.TPL_MODELNAME", Path(__file__).parent / "TPL_MODELNAME.py"
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


def get_jax_model() -> "JAXModel":
    # If the model directory was meanwhile overwritten, this would load the
    #  new version, which would not match the previously imported extension.
    #  This is not allowed, as it would lead to inconsistencies.
    jax_py_file = Path(__file__).parent / "jax.py"
    jax_py_file = jax_py_file.resolve()
    t_imported = TPL_MODELNAME._get_import_time()  # noqa: protected-access
    t_modified = os.path.getmtime(jax_py_file)
    if t_imported < t_modified:
        t_imp_str = datetime.datetime.fromtimestamp(t_imported).isoformat()
        t_mod_str = datetime.datetime.fromtimestamp(t_modified).isoformat()
        raise RuntimeError(
            f"Refusing to import {jax_py_file} which was changed since "
            f"TPL_MODELNAME was imported. This is to avoid inconsistencies "
            "between the different model implementations.\n"
            f"Imported at {t_imp_str}\nModified at {t_mod_str}.\n"
            "Import the module with a different name or restart the "
            "Python kernel."
        )
    jax = amici._module_from_path("jax", jax_py_file)
    return jax.JAXModel_TPL_MODELNAME()


__version__ = "TPL_PACKAGE_VERSION"
