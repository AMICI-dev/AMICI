%module solver_idas

// Add necessary symbols to generated header
%{
#include "amici/solver_idas.h"
using namespace amici;
%}

%newobject amici::IDASolver::clone;

// Process symbols in header
%include "amici/solver_idas.h"
