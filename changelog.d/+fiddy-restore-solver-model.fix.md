The `fiddy` adapter's `*_to_function_and_derivative` helpers no longer
leave the caller's solver at `SensitivityOrder.none` and, for
`run_simulation_to_function_and_derivative`, the model at the parameter
values of the last evaluated point. Both are now restored after each
evaluation (#3299).
