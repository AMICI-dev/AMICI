Creating a `ReturnDataView` or `ExpDataView` is much faster. The array shape
of a field is now read from the C++ object only when that field is first
accessed, instead of reading all shapes when the view is created. This speeds
up `run_simulations` for many experimental conditions, e.g. in PEtab
problems with thousands of experiments.
