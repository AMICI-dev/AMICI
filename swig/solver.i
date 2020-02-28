%module solver

// Add necessary symbols to generated header
%{
#include "amici/solver.h"
using namespace amici;
%}

%rename(equals) operator==;

%ignore getAdjointDerivativeState;
%ignore getAdjointQuadrature;
%ignore getAdjointState;
%ignore getDerivativeState;
%ignore getState;
%ignore getStateSensitivity;
%ignore quadReInitB;
%ignore reInit;
%ignore reInitB;
%ignore sensReInit;
%ignore setup;
%ignore setupB;
%ignore writeSolution;
%ignore writeSolutionB;

// Process symbols in header
%include "amici/solver.h"
