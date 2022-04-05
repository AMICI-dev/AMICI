"""AMICI-generated module for model TPL_MODELNAME"""

import amici

# Ensure we are binary-compatible, see #556
if 'TPL_AMICI_VERSION' != amici.__version__:
    raise amici.AmiciVersionError(
        'Cannot use model `TPL_MODELNAME`, generated with amici=='
        f'TPL_AMICI_VERSION, together with amici=={amici.__version__} '
        'which is currently installed. To use this model, install '
        'amici==TPL_AMICI_VERSION or re-import the model with the amici '
        'version currently installed.'
    )

from TPL_MODELNAME._TPL_MODELNAME import *

__version__ = 'TPL_PACKAGE_VERSION'
