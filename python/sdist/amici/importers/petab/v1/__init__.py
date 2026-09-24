"""PEtab v1 import related code.

.. note::
    PEtab v1 never formally specified a grammar for math expressions
    (observable/noise formulas, condition table entries, ...). For
    parsing such expressions, this module therefore uses PEtab v2's
    math grammar (:func:`petab.v2.math.sympify_petab`) instead, since
    that one is well-defined.
"""

# ID of model parameter that is to be added to SBML model to indicate
#  preequilibration
PREEQ_INDICATOR_ID = "preequilibration_indicator"

from ._petab_import import import_petab_problem

__all__ = [
    "import_petab_problem",
]
