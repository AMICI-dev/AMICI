%module model

// Add necessary symbols to generated header
%{
#include "amici/model.h"
using namespace amici;
%}

// Process symbols in header

%include "amici/model.h"
