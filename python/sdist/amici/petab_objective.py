# THIS FILE IS TO BE REMOVED - DON'T ADD ANYTHING HERE!

import warnings

warnings.warn(
    f"Importing {__name__} is deprecated. Use `amici.petab.simulations` instead.",
    DeprecationWarning,
)

from .petab.parameter_mapping import create_parameter_mapping  # noqa: F401
from .petab.simulations import create_edatas  # noqa: F401
from .petab.simulations import (  # noqa: F401
    aggregate_sllh,
    rdatas_to_measurement_df,
    rdatas_to_simulation_df,
    rescale_sensitivity,
    simulate_petab,
)
