%module hdf5

// Add necessary symbols to generated header
%{
#include "amici/hdf5.h"
using namespace amici;
%}

// Process symbols in header
%include "amici/hdf5.h"
