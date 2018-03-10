%module model_dae

// Add necessary symbols to generated header
%{
#include "amici_model_dae.h"
using namespace amici;
%}

// Process symbols in header
%ignore getSolver;
%include "amici_model_dae.h"
