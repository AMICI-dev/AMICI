Fixed the JAX backend for PEtab v2 problems converted by libpetab's
`ExperimentsToSbmlConverter` after PEtab-dev/libpetab-python#579, where
experiments without model-encoded changes have no experiment-specific
condition. Simulating them failed with `ValueError: Experiment '...' has no
dynamic period with a condition id`, and PEtab SciML array inputs could be
taken from the wrong experiment, giving wrong results. Converted periods are
now matched to the original ones by experiment id instead of condition ids.
The signature of `JAXProblem.load_reinitialisation` changed to
`(experiment, is_preeq)`. To identify experiments in the results of
`amici.sim.jax.petab.run_simulations`, use `experiment_ids` instead of
`dynamic_conditions`.
