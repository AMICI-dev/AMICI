Fixed wrong forward sensitivities after state-dependent events in DAE
models. The event sensitivity updates used the DAE residual instead of the
time derivative of the states, and the post-event time derivative did not
account for the change of the algebraic states (#3341). Also, without
sensitivities, algebraic states at output times coinciding with an event
are now consistent with the post-event differential states.
