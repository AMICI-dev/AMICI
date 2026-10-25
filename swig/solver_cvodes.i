%module solver_cvodes

// Add necessary symbols to generated header
%{
#include "amici/solver_cvodes.h"
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

%newobject amici::CVodeSolver::clone;
%feature("notabstract") amici::CVodeSolver;

// Required with SWIG 4.2.0 https://github.com/AMICI-dev/AMICI/issues/2275
%extend amici::CVodeSolver {
    CVodeSolver() {
        return new CVodeSolver();
    }

    CVodeSolver(Solver const& solver) {
        return new CVodeSolver(dynamic_cast<CVodeSolver const&>(solver));
    }
}

// Process symbols in header
%include "amici/solver_cvodes.h"
