Fixed adjoint sensitivities for simultaneous events, where an earlier event
changes a state that a later one depends on (#3414). The event-time terms of
all events at the same time are now applied at the trigger state. Generated
model code changed (`fdeltaxB`, `fdeltaqB`, and new functions), so previously
imported models must be re-imported; using them with the new version is not
supported.
