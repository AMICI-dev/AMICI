Fixed PySB import with `pysb_model_has_obs_and_noise=True`: sigmas were
assigned to observables in the order of the sigma expressions instead of the
order of the observables, and sigma expressions that were also used in the
model dynamics caused an error.
