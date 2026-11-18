Fixed the sensitivities of event observables (`sz`) and their contribution
to `sllh` (forward and adjoint), which did not account for the dependence of
the event time on the parameters, e.g., for `z = x` at `x < c`, or for
`z = k * time` at `time > k` (#3361, #3551). This also applies to the event
observables of events triggered by other events' assignments, for which the
`SECONDARY_EVENT` warning is no longer issued for models without algebraic
states.
Also fixed adjoint sensitivity analysis for models with event observables
but without (non-event) observables, which failed with
`Incorrect buffer size`.
