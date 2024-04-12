"""AMICI-generated module for model TPL_MODELNAME"""

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

from .TPL_MODELNAME import *  # noqa: F403, F401
from .TPL_MODELNAME import getModel as get_model  # noqa: F401


def get_jax_model() -> "JAXModel":
    from .jax import JAXModel_TPL_MODELNAME

    return JAXModel_TPL_MODELNAME()


__version__ = "TPL_PACKAGE_VERSION"
