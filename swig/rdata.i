%module rdata

// Add necessary symbols to generated header
%{
#include "amici/rdata.h"
using namespace amici;
%}

%ignore processSimulationObjects;
// Process symbols in header
%include "amici/rdata.h"
