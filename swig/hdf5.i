%module hdf5

// Add necessary symbols to generated header
%{
#ifndef AMICI_SWIG_WITHOUT_HDF5
#include "amici/hdf5.h"
using namespace amici;
#endif
%}

%wrapper %{
#ifndef AMICI_SWIG_WITHOUT_HDF5
%}

// Process symbols in header
%include "amici/hdf5.h"

%wrapper %{
#endif // AMICI_SWIG_WITHOUT_HDF5
%}
