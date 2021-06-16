%module amici

%include <exception.i>
%exception {
    try {
        $action
    } catch(std::exception const& ex) {
        SWIG_exception_fail(SWIG_RuntimeError, ex.what());
    } catch(...) {
        SWIG_exception_fail(SWIG_RuntimeError, "Unknown exception occurred");
    }
}

%{
#define SWIG_FILE_WITH_INIT
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <vector>
#include <stdexcept>
%}

// numpy
%init %{
import_array();
%}

// Expose vectors
%include <stl.i>
%template(DoubleVector) std::vector<double>;
%template(IntVector) std::vector<int>;
%template(BoolVector) std::vector<bool>;
%template(StringVector) std::vector<std::string>;
%template(StringDoubleMap) std::map<std::string, double>;

// Let numpy access std::vector
%{
static_assert (sizeof(double) == sizeof (npy_double), "Numpy double size mismatch");
static_assert (sizeof(int) == sizeof (npy_int), "Numpy integer size mismatch");

#include "../swig/stdvec2numpy.h"
using namespace amici;
%}
%include "stdvec2numpy.h"

typedef double realtype;

%include std_unique_ptr.i
wrap_unique_ptr(SolverPtr, amici::Solver)
wrap_unique_ptr(ReturnDataPtr, amici::ReturnData)
wrap_unique_ptr(ModelPtr, amici::Model)
wrap_unique_ptr(ExpDataPtr, amici::ExpData)

%naturalvar amici::SimulationParameters::x0;
%naturalvar amici::SimulationParameters::sx0;
%naturalvar amici::SimulationParameters::parameters;
%naturalvar amici::SimulationParameters::pscale;
%naturalvar amici::SimulationParameters::plist;
%naturalvar amici::SimulationParameters::fixedParameters;
%naturalvar amici::SimulationParameters::fixedParametersPreequilibration;
%naturalvar amici::SimulationParameters::fixedParametersPresimulation;
%naturalvar amici::SimulationParameters::reinitialization_state_idxs_sim;
%naturalvar amici::SimulationParameters::reinitialization_state_idxs_presim;

%ignore amici::ModelContext;
%ignore amici::ContextManager;
%ignore amici::ModelState;
%ignore amici::ModelStateDerived;

// Include before any other header which uses enums defined there
%include "amici/defines.h"

%include "amici/model_dimensions.h"
%include "amici/model_state.h"
%include "amici/simulation_parameters.h"

%include abstract_model.i
%include misc.i
%include edata.i
%include solver.i
%include solver_idas.i
%include solver_cvodes.i
%include model.i
%include model_ode.i
%include model_dae.i
%include rdata.i

#ifndef AMICI_SWIG_WITHOUT_HDF5
%include hdf5.i
#endif

// Return Python list of raw pointers instead of std::vector<std::unique_ptr> which is a huge pain
%typemap(out) std::vector<std::unique_ptr<amici::ReturnData>> %{
    $result = PyList_New($1.size());
    for (int i = 0; i < (int) $1.size(); i++) {
        PyObject *o = SWIG_NewPointerObj($1.at(i).release(), $descriptor(amici::ReturnData*), SWIG_POINTER_OWN);
        PyList_SetItem($result, i, o);
    }
%}

// Add necessary symbols to generated header
// Ignore due to https://github.com/swig/swig/issues/1643
%ignore amici::AmiciApplication::warningF;
%ignore amici::AmiciApplication::errorF;
%{
#include "amici/amici.h"
using namespace amici;
%}

// Prevent using ValueWrapper, but don't expose unique_ptr vector
%ignore std::vector<std::unique_ptr<amici::ReturnData>>;
%template(ReturnDataPtrVector) std::vector<std::unique_ptr<amici::ReturnData>>;

// Process symbols in header
%include "amici/amici.h"

// Expose vectors
%template(ExpDataPtrVector) std::vector<amici::ExpData*>;


// Convert integer values to enum class
// defeats the purpose of enum class, but didn't find a better way to allow for
// vectors of enum class types in python
%{
namespace amici {
std::vector<amici::ParameterScaling> parameterScalingFromIntVector(std::vector<int> const& intVec) {
    std::vector<amici::ParameterScaling> result(intVec.size());
    for (int i = 0; i < (int) result.size(); ++i) {
        result[i] = static_cast<amici::ParameterScaling>(intVec[i]);
    }
    return result;
}
}; // namespace amici
%}
namespace amici {
    std::vector<amici::ParameterScaling> parameterScalingFromIntVector(std::vector<int> const& intVec);
}
%template(ParameterScalingVector) std::vector<amici::ParameterScaling>;


// Add function to check if amici was compiled with OpenMP
%{
namespace amici {
/** AMICI extension was compiled with OpenMP? */
bool compiledWithOpenMP() {
#if defined(_OPENMP)
    return true;
#else
    return false;
#endif
}
};
%}
namespace amici {
bool compiledWithOpenMP();
}

%pythoncode %{
from enum import IntEnum
def enum(prefix):
    values = {k:v for k,v in globals().items() if k.startswith(prefix + '_')}
    values = {k[len(prefix)+1:]:v for k,v in values.items()}
    return IntEnum(prefix, values)
ParameterScaling = enum('ParameterScaling')
SecondOrderMode = enum('SecondOrderMode')
SensitivityOrder = enum('SensitivityOrder')
SensitivityMethod = enum('SensitivityMethod')
LinearSolver = enum('LinearSolver')
InternalSensitivityMethod = enum('InternalSensitivityMethod')
InterpolationType = enum('InterpolationType')
LinearMultistepMethod = enum('LinearMultistepMethod')
NonlinearSolverIteration = enum('NonlinearSolverIteration')
SteadyStateSensitivityMode = enum('SteadyStateSensitivityMode')
SteadyStateStatus = enum('SteadyStateStatus')
NewtonDampingFactorMode = enum('NewtonDampingFactorMode')
FixedParameterContext = enum('FixedParameterContext')
RDataReporting = enum('RDataReporting')
%}

%template(SteadyStateStatusVector) std::vector<amici::SteadyStateStatus>;

// add module docstring and import additional types for typehints
%pythonbegin %{
"""
Core C++ bindings
-----------------
This module encompasses the complete public C++ API of AMICI, which was
exposed via swig. All functions listed here are directly accessible in the
main amici package, i.e., :py:class:`amici.amici.ExpData` is available as
``amici.ExpData``.
Usage of functions and classes from the base ``amici`` package is
generally recommended as they often include convenience wrappers that avoid
common pitfalls when accessing C++ types from python and implement some
nonstandard type conversions.
"""

from typing import Iterable
%}
