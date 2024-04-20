%module misc

%ignore amici::printfToString;
%ignore amici::regexErrorToString;
%ignore amici::writeSlice;
%ignore ContextManager;
%ignore amici::scaleParameters;
%ignore amici::unscaleParameters;

// Add necessary symbols to generated header
%{
#include "amici/misc.h"
using namespace amici;
%}

// Process symbols in header
%include "amici/misc.h"
