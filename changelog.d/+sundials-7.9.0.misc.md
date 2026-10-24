Updated the bundled SUNDIALS from 7.1.1 to
[7.9.0](https://sundials.readthedocs.io/en/latest/Changelog_link.html).
Relevant upstream changes include fixes to the CVODES/IDAS handling of
inequality constraints (`Solver.set_constraints`), which now fail a step
after 10 constraint failures within a single step attempt, and a fix of the
IDAS quadrature predictor after computing consistent initial conditions,
which also applies to the backward quadratures after discontinuities in DAE
models. Simulation results may therefore differ slightly from previous
versions.

For editable installations (`pip install -e ...`), clean out the temporary
build directories (at least `ThirdParty/sundials/build/` and
`python/sdist/build/`) before rebuilding.
