Fixed incorrect adjoint sensitivities for models with multiple events
that can trigger at the same time point. Previously, all such
simultaneously-triggered events incorrectly shared a single pre-/post-event
state snapshot for the adjoint update (#2805).
