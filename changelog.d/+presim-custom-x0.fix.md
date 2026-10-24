Fixed a custom initial state being overwritten by the state reinitialization
at the start of presimulation without pre-equilibration. States are now
reinitialized only at the start of the main simulation in that case, as
documented.
