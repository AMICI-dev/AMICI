New functions `amici.sim.sundials.get_solver_settings` and
`amici.sim.sundials.set_solver_settings` get and apply all settings of a
`Solver` as a dictionary, like `get_model_settings` and `set_model_settings`
for models. Pickling a `Solver` now uses these settings instead of a
temporary HDF5 file, so it also works when AMICI is built without HDF5, and
it no longer leaves temporary files behind.

`IDASolver.set_linear_solver` now accepts a linear solver, as
`CVodeSolver.set_linear_solver` does. Previously, an internal overload
without arguments hid it on `IDASolver` instances.
