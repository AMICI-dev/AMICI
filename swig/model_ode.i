%module model_ode

// Add necessary symbols to generated header
%{
#include "amici_model_ode.h"
using namespace amici;
%}

%include std_unique_ptr.i

wrap_unique_ptr(SolverPtr, amici::Solver)

// Process symbols in header

%include "amici_model_ode.h"
