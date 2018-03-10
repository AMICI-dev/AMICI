%module solver

// Add necessary symbols to generated header
%{
#include "amici_solver.h"
using namespace amici;
%}

// Process symbols in header
%include "amici_solver.h"
