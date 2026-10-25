%module solver_idas

// Add necessary symbols to generated header
%{
#include "amici/solver_idas.h"
using namespace amici;
%}

%ignore reinit_post_process_f;
%ignore reinit_post_process_b;
%ignore reinit_quad_b;
%ignore quad_ss_tolerances_b;
%ignore quad_ss_tolerances;
%ignore solve;
%ignore solve_f;
%ignore get_dky;
%ignore get_sens;
%ignore get_sens_dky;
%ignore get_b;
%ignore get_dky_b;
%ignore get_quad_b;
%ignore get_quad_dky_b;
%ignore get_quad;
%ignore get_quad_dky;

%newobject amici::IDASolver::clone;
%feature("notabstract") amici::IDASolver;

// Required for SWIG 4.2.0 https://github.com/AMICI-dev/AMICI/issues/2275
%extend amici::IDASolver {
    IDASolver() {
        return new IDASolver();
    }

    IDASolver(Solver const& solver) {
        return new IDASolver(dynamic_cast<IDASolver const&>(solver));
    }
}

// Process symbols in header
%include "amici/solver_idas.h"
