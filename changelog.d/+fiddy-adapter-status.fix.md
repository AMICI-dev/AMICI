`amici.adapters.fiddy.run_simulation_to_function_and_derivative` now raises
`fiddy.FunctionEvaluationError` if the simulation fails. Previously, the
simulation status was ignored, so a failed simulation could pass a finite
difference check if it still returned finite values.
