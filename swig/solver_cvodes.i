%module solver_cvodes

// Add necessary symbols to generated header
%{
#include "amici/solver_cvodes.h"
using namespace amici;
%}

%newobject amici::CVodeSolver::clone;
%feature("notabstract") amici::CVodeSolver;

// Process symbols in header
%include "amici/solver_cvodes.h"
