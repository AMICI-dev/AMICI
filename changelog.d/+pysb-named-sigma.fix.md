Fixed PySB import with `MeasurementChannel.sigma` referring to a
`pysb.Expression`, which failed since v0.31.0. Sigma expressions now also
remain usable in the rest of the model.
See [#3538](https://github.com/AMICI-dev/AMICI/issues/3538).
