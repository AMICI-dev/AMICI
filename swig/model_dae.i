%module model_dae

// Add necessary symbols to generated header
%{
#include "amici_model_dae.h"
using namespace amici;
%}


%ignore getSolver;

// Process symbols in header

%include "amici_model_dae.h"
%rename(getSolver) getSolverSwig;
%extend amici::Model_DAE {
    amici::Solver *getSolverSwig() {
        return new amici::IDASolver();
    }
}
