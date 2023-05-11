%module model_ode

// Add necessary symbols to generated header
%{
#include "amici/model_ode.h"
using namespace amici;
%}

%ignore fJvB;
%ignore fxBdot;
%ignore fqBdot;
%ignore fqBdot_ss;


// Process symbols in header

%include "amici/model_ode.h"
