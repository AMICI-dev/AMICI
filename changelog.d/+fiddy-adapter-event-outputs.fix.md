`amici.adapters.fiddy.run_simulation_to_function_and_derivative` now
supports the derivatives of the event outputs `z`, `rz` and `sigmaz`.
Previously, the derivatives of `z` and `sigmaz` failed with a `KeyError`,
and the derivative of `rz` was taken along the event output axis instead of
the parameter axis, which gave wrong values or an `IndexError` (#3571).
