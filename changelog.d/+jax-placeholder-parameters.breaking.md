`SbmlImporter.sbml2jax` no longer detects placeholder parameters by name.
PEtab-based imports are not affected.
Placeholders `noiseParameter{n}` and `observableParameter{n}` in the
`observation_model` formulas must now be declared on the respective
`MeasurementChannel` via the new `noise_placeholders` and
`observable_placeholders` arguments.
Otherwise, they are not treated as placeholders.
The PEtab importers declare them automatically.
Model entities with IDs of this form are no longer treated as placeholders
by any importer; declaring a placeholder with the ID of a model entity is
an error.
Declared placeholders are rejected by `sbml2amici`, the PySB importer, and
for event-resolved measurement channels, which don't support them yet
([#3553](https://github.com/AMICI-dev/AMICI/pull/3553)).
