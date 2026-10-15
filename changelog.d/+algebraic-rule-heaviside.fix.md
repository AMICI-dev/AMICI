Fixed missing event handling for `piecewise` functions in algebraic
rules (#3358). The solver was not stopped at the switching points, and for
state-dependent conditions, forward sensitivities were wrong.
