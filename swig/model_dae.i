%module model_dae

// Add necessary symbols to generated header
%{
#include "amici/model_dae.h"
using namespace amici;
%}

%ignore fM;

// Process symbols in header

%include "amici/model_dae.h"
