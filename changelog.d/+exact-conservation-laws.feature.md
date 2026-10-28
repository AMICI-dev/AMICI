Experimental conservation laws for non-constant species (enabled via
`AMICI_EXPERIMENTAL_SBML_NONCONST_CLS`) are now computed in exact rational
arithmetic. Previously, the floating-point computation of the left nullspace
of the stoichiometric matrix could miss conservation laws, fail with an
`AssertionError`, or silently return wrong ones, even for small integer
stoichiometries. Non-integer stoichiometries are interpreted as their
shortest decimal representation (e.g., `0.1` as `1/10`). This changes the
generated code for affected models (e.g., `1.0*x` becomes `x`, `0.5*x`
becomes `(1.0/2.0)*x`, and floating-point noise such as `10.000000000000012`
is gone). This is also considerably faster for large networks (#3366).
