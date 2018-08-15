%module amici

%{
#define SWIG_FILE_WITH_INIT
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <vector>
#include <stdexcept>
%}
%init %{
import_array();
%}

%{
static_assert (sizeof(double) == sizeof (npy_double), "Numpy double size mismatch");
static_assert (sizeof(int) == sizeof (npy_int), "Numpy integer size mismatch");

#include "../swig/stdvec2numpy.h"
using namespace amici;
%}
%include "stdvec2numpy.h"

typedef double realtype;

%include <stl.i>

%include std_unique_ptr.i
wrap_unique_ptr(SolverPtr, amici::Solver)
wrap_unique_ptr(ReturnDataPtr, amici::ReturnData)
wrap_unique_ptr(ModelPtr, amici::Model)
wrap_unique_ptr(ExpDataPtr, amici::ExpData)

// Include before any other header which uses enums defined there
%include "amici/defines.h"


%include edata.i
%include rdata.i

%include solver.i
%include solver_idas.i
%include solver_cvodes.i
%include model.i
%include model_ode.i
%include model_dae.i

#ifndef AMICI_SWIG_WITHOUT_HDF5
%include hdf5.i
#endif

// Add necessary symbols to generated header
%{
#include "amici/amici.h"
using namespace amici;
%}

// Process symbols in header
%include "amici/amici.h"

// Expose vectors
%template(DoubleVector) std::vector<realtype>;
%template(IntVector) std::vector<int>;
%template(ParameterScalingVector) std::vector<amici::ParameterScaling>;
%template(StringVector) std::vector<std::string>;
