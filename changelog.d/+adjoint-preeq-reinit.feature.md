Adjoint preequilibration is now supported for models with
fixed-parameter-dependent solver-state reinitialization
(`reinitialize_fixed_parameter_initial_states`); this previously raised
`AMICI_NOT_IMPLEMENTED`. Models with conservation laws, remain unsupported
for this combination and now raise a clear error instead (#1156).
