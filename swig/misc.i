%module misc

%ignore amici::printf_to_string;
%ignore amici::regex_error_to_string;
%ignore amici::write_slice;
%ignore ContextManager;
%ignore amici::scale_parameters;
%ignore amici::unscale_parameters;

// Add necessary symbols to generated header
%{
#include "amici/misc.h"
using namespace amici;
%}

// Process symbols in header
%include "amici/misc.h"
