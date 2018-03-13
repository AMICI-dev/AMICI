%module model

// Add necessary symbols to generated header
%{
#include "amici_model.h"
using namespace amici;
%}

%ignore getSolver;

// Process symbols in header

%include "amici_model.h"
