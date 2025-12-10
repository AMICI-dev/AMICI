"""
PEtab import related code.

This submodule contains all functionality related to importing PEtab problems.

For importing and simulating
`PEtab v2 <https://petab.readthedocs.io/en/latest/v2/documentation_data_format.html#>`__
problems, the relevant classes are:

* :class:`PetabImporter`: Import a PEtab problem as an AMICI model.

See :doc:`/examples/example_petab/petab_v2` for example usage.
Note that the PEtab v2 API is still under development and may change in future
releases.

The PEtab legacy v1 import functionality is still available under
:mod:`amici.importers.petab.v1`.
Note that this functionality will be deprecated once the PEtab v2 import is
stable.
Most PEtab v1 problems can be also imported using the PEtab v2 import by
passing a :class:`petab.v1.Problem` instance to the PEtab v2 import functions.
"""

# FIXME: for some tests (petab-sciml, maybe petab-v1-pysb) we still rely on an
#   old PEtab version on which the petab v2 import does not work.
#   Once those tests are updated, we can remove this try-except block.
try:
    from ._petab_importer import *  # noqa: F403, F401
except ImportError:
    pass
