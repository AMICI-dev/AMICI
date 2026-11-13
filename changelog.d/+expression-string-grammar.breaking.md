Math expression strings that are not part of an SBML model or PEtab problem
(`MeasurementChannel` formulas and sigmas, custom noise distributions,
spline arguments, `amici.sim.sundials.evaluate`) are now all parsed by
the new `amici.importers.utils.sympify_str`, which documents the grammar.
Names that sympy would resolve to constants or other non-function objects
are now plain symbols. In particular, `pi`, `E`, `I`, `oo`, `nan`, `true`,
and `false` no longer denote constants or booleans; use, e.g., `exp(1)` or
`acos(-1)` instead.
SBML import now raises a `ValueError` if a `MeasurementChannel` formula or
sigma contains symbols that are neither model entities, observables, nor
`observableParameter{n}`/`noiseParameter{n}` placeholders, instead of
generating model code that fails to compile.
Relates to #3281.
