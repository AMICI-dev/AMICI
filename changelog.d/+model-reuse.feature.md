Model import now stores a hash of the import inputs (model, import options,
relevant package versions and environment variables) in the model directory,
and only reuses a previously imported model if the inputs are unchanged
(#1834).
Previously, the PEtab importers reused existing models in the output
directory even if, for example, the PEtab problem had changed.
`SbmlImporter.sbml2amici`, `pysb2amici`, and the respective JAX functions
can reuse previously imported models via the new `skip_if_unchanged`
argument.
Setting the environment variable `AMICI_NO_MODEL_REUSE=1` disables reusing
models.
See the documentation section "Reusing imported models" for details and
limitations.
