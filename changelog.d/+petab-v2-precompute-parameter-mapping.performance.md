`amici.sim.sundials.petab.ExperimentManager` now derives the
experiment-specific part of the parameter mapping (placeholder overrides,
sensitivity parameter list, parameter scales) only once per experiment
instead of on every call of `apply_parameters`, and `PetabSimulator`
aggregates the log-likelihood gradient with precomputed index mappings.
This reduces the overhead of a PEtab v2 objective function evaluation
compared to PEtab v1, especially for problems with many measurements.
