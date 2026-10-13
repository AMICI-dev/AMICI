Fixed `ForwardProblem::handle_presimulation()` rejecting presimulation
with state reinitialization for any sensitivity method, even though
only adjoint sensitivity analysis is actually affected. This was
masked by the `num_state_reinits()` bug above and incorrectly blocked
forward-sensitivity presimulation with reinitialization.
