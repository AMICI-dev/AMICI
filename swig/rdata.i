%module rdata

// Add necessary symbols to generated header
%{
#include "amici/rdata.h"
using namespace amici;
%}

%ignore processSimulationObjects;
%ignore ModelContext;

// Process symbols in header
%include "amici/rdata.h"
