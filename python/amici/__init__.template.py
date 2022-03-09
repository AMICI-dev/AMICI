"""AMICI-generated module for model TPL_MODELNAME"""

# Ensure we are binary-compatible, see #556
#  (to be removed once the amici extension is removed from the amici base
#  package)
try:
    import amici
    if amici.__version__ != 'TPL_AMICI_VERSION':
        from warnings import warn
        warn('Model TPL_MODELNAME was generated with AMICI '
             'version TPL_AMICI_VERSION. Currently installed AMICI version is '
             f' {amici.__version__}. This model must not be used in '
             'combination with any C++ extension objects from the amici base '
             'package.')
except ModuleNotFoundError:
    # AMICI is not installed, so no problem here
    pass

from TPL_MODELNAME._TPL_MODELNAME import *
from .core import *

__version__ = 'TPL_PACKAGE_VERSION'
