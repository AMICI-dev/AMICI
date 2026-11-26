PEtab v2: `ExperimentManager.apply_parameters` now raises `NotImplementedError`
for timepoint-specific placeholder overrides
instead of silently applying the last override (#3476, #3602).
This only affected an `ExperimentManager` used with a PEtab problem
other than the imported one.
