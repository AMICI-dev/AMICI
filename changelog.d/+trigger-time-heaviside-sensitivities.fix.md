Fixed forward and adjoint sensitivities for events with
`useValuesFromTriggerTime=true` whose assignments depend on a `piecewise`
switch that flips between the trigger and the execution of the event
(#3432). The sensitivities used the switch values at execution time.
