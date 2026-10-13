Fixed splines being treated as unconditionally time-dependent via
fragile string matching, rather than through AMICI's normal symbolic
dependency tracking, when determining which expressions need to be
recomputed at every simulation step vs. only once. Models combining a
spline with an event trigger that depends on it now raise a clear error
at import time instead of silently generating incorrect C++, since
simulating such models is not yet supported (#3245).
