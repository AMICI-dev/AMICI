Setting `InternalSensitivityMethod.staggered1` for DAE models now raises an
error. IDAS does not support this method, and the simulation failed with an
unspecific error in `IDASensInit`.
