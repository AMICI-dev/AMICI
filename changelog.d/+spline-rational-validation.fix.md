Fixed spline node/value validation (monotonicity, positivity) being
silently skipped for non-integer rational literals round-tripped through
SBML, due to AMICI's SBML math writer representing them as an unevaluated
division rather than a proper rational number. It now uses the
`sbmlmath` package, like the SBML math reading path already did (#2146).
