About AMICI
==========

Summary: Maltab Interface for CVODES and IDAS which does all the model processing for you.

AMICI provides a MATLAB interface for the SUNDIALS solvers CVODES (for ordinary differential equations) and IDAS (for algebraic differential equations). AMICI allows the user to specify differential equation models in terms of symbolic variables in MATLAB and automatically compiles such models as .mex simulation files. In contrast to the SUNDIALSTB interface, all necessary functions are transformed into native C code, which allows for a significantly faster compilation. Beyond forward integration, the compiled simulation file also allows for forward sensitivity analysis, steady state sensitivity analysis and adjoint sensitivity analysis for likelihood based output functions.

The interface was designed to provide routines for efficient gradient computation in parameter estimation of biochemical reaction models but is also applicable to a wider range of differential equation constrained optimization problems.

Current build status
====================

Linux/Mac: [![TravisCI](https://travis-ci.org/ICB-DCM/AMICI.svg?branch=master)](https://travis-ci.org/ICB-DCM/AMICI)
Windows: [![AppVeyor](https://ci.appveyor.com/api/projects/status/ob315laj1i6i3om3?svg=true)](https://ci.appveyor.com/project/FFroehlich/amici)

FAQ
===

Q: My model fails to build.
A: Remove the corresponding model directory located in AMICI/models/*yourmodelname* and compile again.

Q: It still does not compile.
A: Make an [https://github.com/ICB-DCM/AMICI/issues](issue) and we will have a look.

Q: The simulation/sensitivities I get are incorrect.
A: There are some known issues, especially with adjoint sensitivities and events. If your particular problem is not featured in the Make an [https://github.com/ICB-DCM/AMICI/issues](issues), please add it!

 
