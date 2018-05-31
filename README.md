# About AMICI


AMICI provides a multilanguage (Python, C++, Matlab) interface for the SUNDIALS solvers CVODES (for ordinary differential equations) and IDAS (for algebraic differential equations). AMICI allows the user to read differential equation models specified as SBML and automatically compiles such models as .mex simulation files, C++ executables or python modules. In contrast to the SUNDIALSTB interface, all necessary functions are transformed into native C++ code, which allows for a significantly faster simulation. Beyond forward integration, the compiled simulation file also allows for forward sensitivity analysis, steady state sensitivity analysis and adjoint sensitivity analysis for likelihood based output functions.

The interface was designed to provide routines for efficient gradient computation in parameter estimation of biochemical reaction models but is also applicable to a wider range of differential equation constrained optimization problems.

## Current build status

[![Build Status](https://travis-ci.org/ICB-DCM/AMICI.svg?branch=master)](https://travis-ci.org/ICB-DCM/AMICI)
[![CodeCov](https://codecov.io/gh/ICB-DCM/AMICI/branch/master/graph/badge.svg)](https://codecov.io/gh/ICB-DCM/AMICI)


