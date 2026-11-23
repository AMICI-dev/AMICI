Fixed `ExpData(model)` ignoring the model's state reinitialization setting
(`Model.set_reinitialize_fixed_parameter_initial_states`) at the start of
presimulation; states were only reinitialized at the start of the main
simulation (#3458). Relatedly, `ExpData.reinitialization_state_idxs_presim`
and `ExpData.reinitialization_state_idxs_sim` are now also applied if no
presimulation or main simulation fixed parameters, respectively, are
provided in the `ExpData`, as is already the case for
`ExpData.reinitialize_fixed_parameter_initial_states`.
This changes results only for simulations with pre-equilibration and
presimulation and state reinitialization enabled, or with these index lists
set but the respective fixed parameters not provided.
Adjoint sensitivity analysis with presimulation now raises an error if
`ExpData.reinitialization_state_idxs_presim` is set without
`ExpData.fixed_parameters_presimulation` and these states are
reinitialized, since reinitialization at the start of presimulation is not
yet supported for adjoint sensitivity analysis.
