Adjoint pre-equilibration is now supported with fixed-parameter-dependent
reinitialization of non-constant states
(`reinitialize_fixed_parameter_initial_states`) for ODE models without
conservation laws. Previously, this failed with `AMICI_NOT_IMPLEMENTED`.
For models with conservation laws or DAE models, it fails with a clear
error message (#1156).
