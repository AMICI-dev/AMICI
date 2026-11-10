`DoubleVector` and `IntVector` now implement `__array__`, so
`numpy.array`/`numpy.asarray` copy their data in C++ instead of building the
array element by element via the Python iterator. This makes accessing the
fields of a `ReturnDataView` (e.g. `x`, `w`, `sx`, `FIM`) about three orders
of magnitude faster, which mattered for models with many expressions or
sensitivity parameters, and for PEtab problems with many experiments
([#3510](https://github.com/AMICI-dev/AMICI/issues/3510)).
