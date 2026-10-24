Fixed `DEModel.num_state_reinits()` always returning 0, regardless of
whether the model actually has solver states with fixed-parameter-
dependent initial values. This fed into `Model::nx_reinit()`, which
guards adjoint presimulation with state reinitialization, which is not
supported. Since the count was always 0, this guard never triggered,
silently allowing unsupported (and potentially incorrect) simulations
to run instead of raising an error. Such simulations now correctly fail
with a clear error message. The guard now only applies if
reinitialization is actually enabled for the respective condition
(#3277).
