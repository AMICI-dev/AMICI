%module edata

// Add necessary symbols to generated header
%{
#include "amici/edata.h"
using namespace amici;
%}

// Process symbols in header
%include "amici/edata.h"
