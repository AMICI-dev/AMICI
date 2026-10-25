%module rdata

// Add necessary symbols to generated header
%{
#include "amici/rdata.h"
using namespace amici;
%}

%ignore process_simulation_objects;
// ModelContext is already ignored globally in amici.i

// Process symbols in header
%include "amici/rdata.h"
