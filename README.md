# About AMICI


AMICI provides a multilanguage (Python, C++, Matlab) interface for the SUNDIALS solvers CVODES (for ordinary differential equations) and IDAS (for algebraic differential equations). AMICI allows the user to read differential equation models specified as SBML and automatically compiles such models as .mex simulation files, C++ executables or python modules. In contrast to the SUNDIALSTB interface, all necessary functions are transformed into native C++ code, which allows for a significantly faster simulation. Beyond forward integration, the compiled simulation file also allows for forward sensitivity analysis, steady state sensitivity analysis and adjoint sensitivity analysis for likelihood based output functions.

The interface was designed to provide routines for efficient gradient computation in parameter estimation of biochemical reaction models but is also applicable to a wider range of differential equation constrained optimization problems.

## Publications

[Fröhlich, F., Kaltenbacher, B., Theis, F. J., & Hasenauer, J. (2017). Scalable Parameter Estimation for Genome-Scale Biochemical Reaction Networks. Plos Computational Biology, 13(1), e1005331. doi: 10.1371/journal.pcbi.1005331](https://doi.org/10.1371/journal.pcbi.1005331)

[Fröhlich, F., Theis, F. J., Rädler, J. O., & Hasenauer, J. (2017). Parameter estimation for dynamical systems with discrete events and logical operations. Bioinformatics, 33(7), 1049-1056. doi: 10.1093/bioinformatics/btw764](https://doi.org/10.1093/bioinformatics/btw764)

## Current build status

[![Build Status](https://travis-ci.org/ICB-DCM/AMICI.svg?branch=master)](https://travis-ci.org/ICB-DCM/AMICI)
[![CodeCov](https://codecov.io/gh/ICB-DCM/AMICI/branch/master/graph/badge.svg)](https://codecov.io/gh/ICB-DCM/AMICI)
