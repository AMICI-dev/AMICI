%module model

// Add necessary symbols to generated header
%{
#include "amici/model.h"
using namespace amici;
%}

%import "amici/abstract_model.h"
// Process symbols in header
%include "amici/model.h"
