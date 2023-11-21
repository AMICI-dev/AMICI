"""
PEtab Import
------------
Import a model in the :mod:`petab` (https://github.com/PEtab-dev/PEtab) format
into AMICI.
"""
import warnings

warnings.warn(
    "Importing amici.petab_import is deprecated. Use `amici.petab` instead.",
    DeprecationWarning,
)

from .petab.import_helpers import (
    get_fixed_parameters,
    get_observation_model,
    petab_noise_distributions_to_amici,
    petab_scale_to_amici_scale,
)

# DEPRECATED - DON'T ADD ANYTHING NEW HERE
from .petab.petab_import import (
    check_model,
    import_model,
    import_model_sbml,
    import_petab_problem,
)
from .petab.sbml_import import species_to_parameters
