Fixed the algebraic states at the start of the main simulation or
presimulation after a previous period for DAE models. They were made
consistent with the Heaviside variables of the previous period, which, e.g.,
for a `piecewise` condition that switches only due to the changed fixed
parameters, could give the algebraic states of the wrong branch in `x0` and
at the initial timepoint (#3459).
