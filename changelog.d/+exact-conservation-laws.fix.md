Conservation laws for non-constant species are now computed in exact rational
arithmetic. Previously, the floating-point computation of the left nullspace
of the stoichiometric matrix could miss conservation laws, fail with an
`AssertionError`, or silently return wrong ones, even for small integer
stoichiometries. Non-integer stoichiometries are interpreted as their
shortest decimal representation (e.g., `0.1` as `1/10`). This is also
considerably faster for large networks (#3366).
