%module amici

%include <std_string.i>
%include <std_vector.i>
%include std_unique_ptr.i

namespace std
{
    %template(DoubleVector) vector<double>;
    %template(DoubleVector) vector<realtype>;
    %template(IntVector) vector<int>;
}

%include edata.i
%include rdata.i

%include solver.i
%include solver_idas.i
%include solver_cvodes.i
%include model.i
%include model_ode.i
%include model_dae.i
%include hdf5.i


// Add necessary symbols to generated header
%{
#include "amici.h"
using namespace amici;
%}



wrap_unique_ptr(ReturnDataPtr, amici::ReturnData)

// Process symbols in header
%include "amici.h"
