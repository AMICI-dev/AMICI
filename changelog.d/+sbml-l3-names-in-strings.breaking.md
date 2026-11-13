For SBML import, the SBML L3 formula names `INF`, `NaN`, `rem`, `avogadro`,
`exponentiale`, and `log10` no longer have a special meaning in
`MeasurementChannel` formulas and sigmas, neither in strings nor as sympy
symbols with these names. Strings are parsed by
`amici.importers.utils.sympify_str`. Unless the model has an entity with
that ID, `INF`, `NaN`, `avogadro`, and `exponentiale` now result in an
unknown-symbol error, and `log10(...)` in an unsupported-function error.
`rem(...)` is sympy's polynomial remainder now.
