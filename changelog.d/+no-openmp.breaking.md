AMICI no longer uses OpenMP. `run_simulations(..., num_threads=...)` and
`Model.simulate(..., num_threads=...)` now always simulate multiple
conditions in parallel using native threads, also on platforms where AMICI
previously was built without OpenMP, such as macOS without `libomp`.
OpenMP-specific compiler and linker flags are no longer needed for
installing AMICI or compiling models. `compiled_with_openmp` is deprecated
and always returns `False`.
