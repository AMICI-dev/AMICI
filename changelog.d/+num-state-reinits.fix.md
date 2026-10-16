Fixed `DEModel.num_state_reinits()` always returning 0, regardless of
whether the model actually has solver states with fixed-parameter-
dependent initial values. This fed into `Model::nx_reinit()`, which
guards adjoint presimulation with state reinitialization (still
unsupported) and, before it was implemented (#1156), also adjoint
preequilibration with state reinitialization — since the count was
always 0, these guards never triggered, silently allowing the
unsupported combination to run instead of raising an error. The guards
now only apply if reinitialization is actually enabled for the
respective condition (#3277).
