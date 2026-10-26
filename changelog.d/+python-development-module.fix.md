Building AMICI and AMICI models only requires the Python headers, no longer
the Python library (CMake's `Development.Module` instead of `Development`),
so they can be built with Python installations that do not ship
`libpython`, such as the manylinux images used for building wheels.
