"""AMICI-generated module for model TPL_MODELNAME"""

import amici
from pathlib import Path

# Ensure we are binary-compatible, see #556
if 'TPL_AMICI_VERSION' != amici.__version__:
    raise amici.AmiciVersionError(
        f'Cannot use model `TPL_MODELNAME` in {Path(__file__).parent}, '
        'generated with amici==TPL_AMICI_VERSION, '
        f'together with amici=={amici.__version__} '
        'which is currently installed. To use this model, install '
        'amici==TPL_AMICI_VERSION or re-import the model with the amici '
        'version currently installed.'
    )

from TPL_MODELNAME._TPL_MODELNAME import *

__version__ = 'TPL_PACKAGE_VERSION'
