PySB import:
* Observables are now ordered as in the `observation_model` list, instead of
  in the order of the corresponding PySB observables and expressions.
* Observable names are now taken from `MeasurementChannel.name`, if set,
  instead of the name of the PySB observable or expression. For PEtab
  problems with PySB models, these are the observable names from the
  PEtab problem, if set.
* A `MeasurementChannel` without formula whose ID is not the name of a PySB
  observable or expression now raises a `ValueError`, instead of being
  silently ignored. Duplicate channel IDs and event-resolved channels now
  raise an error, instead of being silently ignored or treated as
  time-resolved, respectively.
* `pysb_model_has_obs_and_noise` is now ignored. PySB expressions used as
  observable or sigma always remain model expressions, i.e., they are
  included in the model's expression IDs. This also applies to PEtab
  problems with PySB models.
