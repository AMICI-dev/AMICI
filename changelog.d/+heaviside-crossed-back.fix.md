Fixed wrong simulation results or failures if an event assignment reverts a
trigger or `piecewise` condition that changed at the same time, e.g.,
`at x < thr: x = 1.5 * x` (#3425).
