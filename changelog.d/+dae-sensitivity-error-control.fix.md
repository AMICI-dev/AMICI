Fixed forward sensitivities of DAE models being integrated without
sensitivity error control after the first event or discontinuity (e.g., a
switch in a piecewise function) (#3464). If this happened during
presimulation, it affected the entire main simulation. This affected the
staggered corrector (`InternalSensitivityMethod.staggered`), which may now
take more steps after events. For the default simultaneous corrector, IDAS
does not include sensitivities in the local error test at all; this is now
documented, together with the recommendation to use the staggered corrector
for error-controlled DAE forward sensitivities.
