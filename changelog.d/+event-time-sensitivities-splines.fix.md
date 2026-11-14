Fixed the import of models with a time-varying spline and an event trigger
that depends on a model expression, which failed with a
`NotImplementedError` when generating sensitivity code, even if no trigger
depends on the spline.
