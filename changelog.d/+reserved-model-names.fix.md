Model names that coincide with the name of a generated source file (e.g.,
`x0` or `main`, case-insensitive) are now rejected at import. Previously,
model files were overwritten, and the import failed during compilation or
when loading the model module.
