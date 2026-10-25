Fixed uncompilable model code for noise formulas (sigmas) that reference
expressions, such as targets of SBML assignment rules, for both observables
and event observables. These expressions are now replaced by their
definitions (#3377, #3384). Noise formulas that depend on quantities that are
not available for their evaluation are now rejected at import with an
informative error message instead of failing during compilation. These are
piecewise functions, splines, conservation law totals, states that cannot be
expressed in terms of the respective observable, and, for event observables,
states and the event observable itself (#3389).
