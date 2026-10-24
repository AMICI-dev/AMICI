Fixed adjoint sensitivities for pre-equilibration or presimulation if
the Heaviside variables, static expressions, or conservation law totals at
their end differ from those during the main simulation, e.g., due to a
`piecewise` condition on a time-dependent or fixed-parameter-dependent value
or due to differing fixed parameters. This affected the adjoint right-hand
side and event assignments depending on such switches (#3443), and, for
static expressions, also outputs and likelihoods after adjoint
pre-equilibration (#3448).
