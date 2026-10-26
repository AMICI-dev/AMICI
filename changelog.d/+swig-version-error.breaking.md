A SWIG version mismatch between AMICI and a model is now an error instead of
a warning. Compiling a model with a SWIG executable of a different version
than the one used for building AMICI raises `amici.AmiciVersionError` before
the build starts, and importing such a model raises it as well. Model and
AMICI extensions share SWIG runtime type information, which differs between
SWIG versions. Set the `SWIG` environment variable to a matching SWIG
executable, or rebuild AMICI.
