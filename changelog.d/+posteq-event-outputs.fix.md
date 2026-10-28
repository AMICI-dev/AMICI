Fixed simulations of models with event outputs failing with
`vector::_M_range_check` if post-equilibration is used and any event occurred
fewer than `n_max_event` times until the last finite output timepoint, and
crashing with adjoint sensitivities and pre-simulation (#3413). The outputs of
events that did not occur were filled in at the end of each simulation period
instead of only at the last output timepoint of the main simulation. If the
last output timepoint was the initial time, they were thus also filled in
during pre-simulation, and their contributions to adjoint sensitivities were
counted twice. With post-equilibration, the outputs of events that did not
occur until the last finite output timepoint are not filled in (#3460).
