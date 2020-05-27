%module misc

%ignore amici::printfToString;

// Add necessary symbols to generated header
%{
#include "amici/misc.h"
using namespace amici;
%}

// Process symbols in header
%include "amici/misc.h"
