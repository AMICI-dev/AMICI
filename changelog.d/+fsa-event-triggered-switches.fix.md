Fixed forward sensitivities for events and `piecewise` switches whose
trigger only depends on variables changed by event assignments, or on time
only via such switches (#3298).
The contribution of the triggering event's time sensitivity was missing.
Models need to be re-imported.
