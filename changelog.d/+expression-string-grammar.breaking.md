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
Use `time` for the model time in `MeasurementChannel` formulas and sigmas.
SBML import now checks measurement channels for symbols that the generated
model code can't access, and raises a `ValueError` early instead of
producing model code that fails to compile or computes wrong results.
Detected cases: symbols in `MeasurementChannel` formulas or sigmas that are
neither model entities, the model time, the channel's own observable (in
the sigma), nor placeholders declared for that channel; symbols in
custom noise distributions other than model parameters, the channel's own
observable, measurement, and sigma, and, for time-resolved channels, other
time-resolved observables; and sigmas of event observables that depend on
time-resolved observables, or vice versa.
`amici.importers.utils.noise_distribution_to_cost_function` now returns a
function that takes the sympy symbols for the observable, measurement, and
sigma and returns a sympy expression, instead of a formula string.
Relates to #3281.
