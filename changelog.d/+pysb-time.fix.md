PySB import now supports `pysb.time` (e.g., from BNGL `time()`).
Previously, models using it failed to compile,
and expressions depending on it were evaluated only once.
As in SBML import, PySB observables named `time` are rejected.
In PEtab v2 problems with PySB models,
`time` in observable and noise formulas now refers to the simulation time
instead of being turned into a parameter with value 1 (#3592).
Validating such problems requires a libpetab release newer than 0.10.0.
