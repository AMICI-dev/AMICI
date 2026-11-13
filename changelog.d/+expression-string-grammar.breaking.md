Math expression strings that are not part of an SBML model or PEtab problem
(`MeasurementChannel` formulas and sigmas, custom noise distributions,
spline arguments, `amici.sim.sundials.evaluate`) are now all parsed by
the new `amici.importers.utils.sympify_str`, which documents the grammar.
Names that are not followed by `(` are now always plain symbols, even if
sympy defines an object or function with that name. In particular, `pi`,
`E`, `I`, `oo`, `nan`, `true`, and `false` no longer denote constants or
booleans; use, e.g., `exp(1)` or `acos(-1)` instead. Names such as `re`,
`ff`, or `factor` can now be used as symbols. `gamma`, `beta`, and `zeta`
are symbols unless called as functions.
SBML import now raises a `ValueError` if a `MeasurementChannel` formula,
sigma, or custom noise distribution contains symbols that are neither model
entities, observables, nor `observableParameter{n}`/`noiseParameter{n}`
placeholders. Previously, such symbols remained undefined, and model import
or compilation failed later. Use `time` for the model time; `t` is no
longer accepted.
`observableParameter{n}` placeholders are now also recognized in
`MeasurementChannel`s without a sigma.
`amici.importers.utils.noise_distribution_to_cost_function` now returns a
function that takes the sympy symbols for the observable, measurement, and
sigma and returns a sympy expression, instead of a formula string.
Relates to #3281.
