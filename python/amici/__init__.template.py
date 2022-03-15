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

import os


def _imported_from_setup() -> bool:
    """Check whether this module is imported from `setup.py`"""

    from inspect import getouterframes, currentframe

    # in case we are imported from setup.py, this will be the AMICI package
    # root directory (otherwise it is most likely the Python library directory,
    # we are not interested in)
    package_root = os.path.realpath(os.path.dirname(os.path.dirname(__file__)))

    for frame in getouterframes(currentframe(), context=0):
        # Need to compare the full path, in case a user tries to import AMICI
        # from a module `*setup.py`. Will still cause trouble if some package
        # requires the AMICI extension during its installation, but seems
        # unlikely...
        frame_path = os.path.realpath(os.path.expanduser(frame.filename))
        if frame_path == os.path.join(package_root, 'setup.py'):
            return True

    return False


if not _imported_from_setup():
    from TPL_MODELNAME._TPL_MODELNAME import *

__version__ = 'TPL_PACKAGE_VERSION'
