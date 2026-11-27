PySB import: `MeasurementChannel.formula` and `MeasurementChannel.sigma`
can now be arbitrary formulas referring to `time` and the names of model
parameters, expressions, and observables. For measurement channels with a
formula, IDs that are `time` or such a name are currently not supported.
Custom noise distributions referring to names that aren't available in the
log-likelihood now raise a `ValueError` at import, instead of producing
model code that fails to compile.
