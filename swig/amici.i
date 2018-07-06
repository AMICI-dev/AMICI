%module amici

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
%template(ParameterScalingVector) std::vector<amici::AMICI_parameter_scaling>;
%template(StringVector) std::vector<std::string>;
