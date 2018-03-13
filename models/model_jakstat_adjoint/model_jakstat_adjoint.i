%module model_jakstat_adjoint

// Add necessary symbols to generated header
%{
#include "wrapfunctions.h"
%}

%include ../swig/std_unique_ptr.i
%include ../swig/rdata.i

wrap_unique_ptr(ModelPtr, amici::Model)
wrap_unique_ptr(ReturnDataPtr, amici::ReturnData)

// Process symbols in header
%include "wrapfunctions.h"
