%module model_ode

// Add necessary symbols to generated header
%{
#include "amici_model_ode.h"
using namespace amici;
%}

%ignore getSolver;

// Process symbols in header

%include "amici_model_ode.h"
%rename(getSolver) getSolverSwig;
%extend amici::Model_ODE {
    amici::Solver *getSolverSwig() {
        return new amici::CVodeSolver();
    }
}

