%module hdf5

%rename("$ignore", regextarget=1, fullname=1) "amici::hdf5::.*$";
%rename("%s") amici::hdf5::read_model_data_from_hdf5;
%rename("%s") amici::hdf5::read_exp_data_from_hdf5;
%rename("%s") amici::hdf5::read_solver_settings_from_hdf5;
%rename("%s") amici::hdf5::write_return_data_to_hdf5;
%rename("%s") amici::hdf5::write_exp_data_to_hdf5;
%rename("%s") amici::hdf5::write_solver_settings_to_hdf5;

// Add necessary symbols to generated header
%{
#ifndef AMICI_SWIG_WITHOUT_HDF5
#include "amici/hdf5.h"
using namespace amici;
#endif
%}

// Process symbols in header
%include "amici/hdf5.h"
