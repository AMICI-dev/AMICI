%module solver_cvodes

// Add necessary symbols to generated header
%{
#include "amici/solver_cvodes.h"
using namespace amici;
%}

// Process symbols in header
%include "amici/solver_cvodes.h"
