Fixed adjoint sensitivities for models with event outputs. If an event
occurred fewer than `n_max_event` times, the adjoint simulation failed with
"Mismatch in Heaviside size" (#3360). Otherwise, the contributions of event
measurements to the gradient were taken from the next event occurrence, so
the first occurrence did not contribute, and if an event occurred
`n_max_event` times, data was read out of bounds.
