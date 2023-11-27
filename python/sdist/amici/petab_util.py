"""Various helper functions for working with PEtab problems."""
import warnings

from .petab.util import PREEQ_INDICATOR_ID, get_states_in_condition_table

warnings.warn(
    f"Importing {__name__} is deprecated. Use `amici.petab.util` instead.",
    DeprecationWarning,
)
