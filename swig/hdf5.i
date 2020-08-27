%module hdf5

%rename("$ignore", regextarget=1, fullname=1) "amici::hdf5::.*$";
%rename("%s") amici::hdf5::readModelDataFromHDF5;
%rename("%s") amici::hdf5::readSimulationExpData;
%rename("%s") amici::hdf5::readSolverSettingsFromHDF5;
%rename("%s") amici::hdf5::writeReturnData;
%rename("%s") amici::hdf5::writeSimulationExpData;
%rename("%s") amici::hdf5::writeSolverSettingsToHDF5;

// Add necessary symbols to generated header
%{
#ifndef AMICI_SWIG_WITHOUT_HDF5
#include "amici/hdf5.h"
using namespace amici;
#endif
%}

// Process symbols in header
%include "amici/hdf5.h"
