Fixed wrong gradients when the main simulation used adjoint
sensitivities, pre-equilibration used forward sensitivities, and
fixed-parameter-dependent state reinitialization was enabled.
`ReturnData::process_backward_problem()` misclassified this
combination as "pre-equilibration with adjoint sensitivities via
backward simulation" whenever pre-equilibration happened to require
actual integration (e.g. due to the default
`SteadyStateSensitivityMode::integrationOnly`), and then used the
pre-equilibrium steady state's sensitivity instead of the
reinitialized main-simulation initial-state sensitivity as the
gradient's boundary term (#3278).
