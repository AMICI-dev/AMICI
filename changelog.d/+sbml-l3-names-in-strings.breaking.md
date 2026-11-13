For SBML import, the SBML L3 formula names `INF`, `NaN`, `rem`, `avogadro`,
`exponentiale`, and `log10` no longer have a special meaning in
`MeasurementChannel` formulas and sigmas, neither in strings nor as sympy
symbols with these names. Strings are parsed with sympy's grammar, where,
e.g., `rem` is the polynomial remainder. To use SBML L3 math, convert it to
a sympy expression first, e.g., using `sbmlmath`.
