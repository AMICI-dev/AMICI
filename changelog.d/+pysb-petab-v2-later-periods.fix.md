Fixed PEtab v2 import of PySB models with condition changes in any period
but the first one, including the main period after pre-equilibration, which
failed with a `TypeError`. Species targets (via the mapping table) are now
supported in all periods, with target values that may depend on the current
model state in later periods. First-period species targets no longer require
a `pysb.Initial` with an identical pattern, as long as the species is part of
the generated reaction network, and the same species can be changed in
multiple experiments. Changes of `pysb.Expression`s and `pysb.Observable`s,
which are not valid condition targets according to the PEtab v2
specification and were silently ignored, now raise an error, as do changes of
parameters that compartment sizes depend on. The JAX backend raises a
`NotImplementedError` for PySB problems with changes after the first period.
