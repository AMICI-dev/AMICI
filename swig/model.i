%module model

// Add necessary symbols to generated header
%{
#include "amici_model.h"
using namespace amici;
%}

// Process symbols in header
%ignore getSolver;
%include "amici_model.h"
