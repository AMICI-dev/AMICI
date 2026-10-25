Model extensions are now built with the `swig` executable of the
current Python environment, if available, instead of the first `swig`
on `PATH`, which could differ from the SWIG version used to build
AMICI. An explicitly set `SWIG` environment variable is still
respected. A SWIG version mismatch is now also reported before the
model is compiled, not only when the model module is imported.
