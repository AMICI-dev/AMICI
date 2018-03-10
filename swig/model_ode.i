%module model_ode

// Add necessary symbols to generated header
%{
#include "amici_model_ode.h"
using namespace amici;
%}

// Process symbols in header
%ignore getSolver;
%include "amici_model_ode.h"
