%module amici
%include edata.i
%include rdata.i
%include model.i
%include solver.i
%include model_ode.i
%include model_dae.i
%include solver_idas.i
%include solver_cvodes.i

// Add necessary symbols to generated header
%{
#include "amici.h"
using namespace amici;
%}

%include std_unique_ptr.i

wrap_unique_ptr(ReturnDataPtr, amici::ReturnData)

// Process symbols in header
%include "amici.h"
