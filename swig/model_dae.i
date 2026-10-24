%module model_dae

// Add necessary symbols to generated header
%{
#include "amici/model_dae.h"
using namespace amici;
%}

%ignore fM;
%ignore fdxdotdt;

// Process symbols in header

%include "amici/model_dae.h"
