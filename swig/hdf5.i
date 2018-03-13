%module hdf5

// Add necessary symbols to generated header
%{
#include "amici_hdf5.h"
using namespace amici;
%}

%include std_unique_ptr.i

wrap_unique_ptr(ExpDataPtr, amici::ExpData)


// Process symbols in header

%include "amici_hdf5.h"
