"""AMICI-generated module for model TPL_MODELNAME"""

import amici

# Ensure we are binary-compatible, see #556
if 'TPL_AMICI_VERSION' != amici.__version__:
    raise RuntimeError('Cannot use model TPL_MODELNAME, generated with AMICI '
                       'version TPL_AMICI_VERSION, together with AMICI version'
                       f' {amici.__version__} which is present in your '
                       'PYTHONPATH. Install the AMICI package matching the '
                       'model version or regenerate the model with the AMICI '
                       'currently in your path.')

from TPL_MODELNAME._TPL_MODELNAME import *

__version__ = 'TPL_PACKAGE_VERSION'
