"""
Functionality for simulating PEtab problems.

This submodule contains all functionality related to simulating PEtab problems.

See also :mod:`amici.importers.petab` for importing PEtab problems.

For simulating
`PEtab v2 <https://petab.readthedocs.io/en/latest/v2/documentation_data_format.html#>`__
problems, the relevant classes are:

* :class:`PetabSimulator`: Simulate PEtab problems with AMICI.
* :class:`ExperimentManager`: Create :class:`amici.ExpData` objects for PEtab
  experiments.

See :doc:`/examples/example_petab/petab_v2` for example usage.
Note that the PEtab v2 API is still under development and may change in future
releases.

The PEtab v1 legacy functionality is still available under
:mod:`amici.simulations.sim.sundials.petab.v1`.
Note that this functionality will be deprecated once the PEtab v2 import is
stable.
"""

from ._v2 import ExperimentManager, PetabSimulator
from .v1 import (
    EDATAS,
    LLH,
    RDATAS,
    RES,
    S2LLH,
    SLLH,
    SRES,
)

__all__ = [
    "PetabSimulator",
    "ExperimentManager",
    "EDATAS",
    "LLH",
    "RDATAS",
    "RES",
    "S2LLH",
    "SLLH",
    "SRES",
]
