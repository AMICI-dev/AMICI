Fixed `TypeError` when indexing `DoubleVector` or `IntVector` (e.g.,
`rdata.ptr.x[0]`) if another SWIG module such as `libsbml` was imported
before the amici extension, which is always the case after `import amici`.
The amici extension and model extensions now use their own SWIG runtime type
table. Models need to be re-imported
([#3511](https://github.com/AMICI-dev/AMICI/issues/3511)).
