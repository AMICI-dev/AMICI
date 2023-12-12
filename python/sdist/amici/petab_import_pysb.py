import warnings

from .petab.pysb_import import *  # noqa: F401, F403

# DEPRECATED - DON'T ADD ANYTHING NEW HERE

warnings.warn(
    "Importing amici.petab_import_pysb is deprecated. Use `amici.petab.pysb_import` instead.",
    DeprecationWarning,
)
