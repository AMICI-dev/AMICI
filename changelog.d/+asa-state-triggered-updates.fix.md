Fixed incorrect adjoint sensitivities for models with SBML events that
are both state-triggered and state-updating. The trigger-time-sensitivity
terms were evaluated using the post-event state instead of the pre-event
state. The state trajectory was unaffected (#3257).
