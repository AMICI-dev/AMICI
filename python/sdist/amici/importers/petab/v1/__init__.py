"""PEtab v1 import related code."""

# ID of model parameter that is to be added to SBML model to indicate
#  preequilibration
PREEQ_INDICATOR_ID = "preequilibration_indicator"

from ._petab_import import import_petab_problem

__all__ = [
    "import_petab_problem",
]
