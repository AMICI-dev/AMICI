Fixed `DEModel.num_state_reinits()` always returning 0, regardless of
whether the model actually has solver states with fixed-parameter-
dependent initial values. This fed into `Model::nx_reinit()`, which
guards two currently-unsupported combinations (adjoint preequilibration
with state reinitialization, see #1156; and the analogous case for
adjoint presimulation) — since the count was always 0, those guards
never triggered, silently allowing unsupported (and potentially
incorrect) simulations to run instead of raising an error. Models that
hit either combination now correctly fail with a clear error message.
The guards now only apply if reinitialization is actually enabled for
the respective condition (#3277).
