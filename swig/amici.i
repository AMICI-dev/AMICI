%module amici

typedef double realtype;

%include <std_string.i>
%include <std_vector.i>

%include std_unique_ptr.i
wrap_unique_ptr(SolverPtr, amici::Solver)
wrap_unique_ptr(ReturnDataPtr, amici::ReturnData)
wrap_unique_ptr(ModelPtr, amici::Model)
wrap_unique_ptr(ExpDataPtr, amici::ExpData)

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
%include "amici/defines.h"

namespace std
{
    %template(DoubleVector) vector<realtype>;
    %template(IntVector) vector<int>;
}
