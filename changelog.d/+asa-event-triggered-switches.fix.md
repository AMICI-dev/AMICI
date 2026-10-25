Fixed adjoint sensitivities for events and `piecewise` switches whose
trigger only depends on variables changed by event assignments (#3298).
The contribution of the triggering event's time sensitivity was missing.
