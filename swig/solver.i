%module solver

// Add necessary symbols to generated header
%{
#include "amici/solver.h"
using namespace amici;
%}

%rename(equals) operator==;

// remove functions that use AmiVector(Array) since that class anyways cannot
// be exposed in swig
%ignore getAdjointDerivativeState;
%ignore getAdjointQuadrature;
%ignore getQuadrature;
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
%ignore setupSteadystate;
%ignore writeSolution;
%ignore writeSolutionB;
%ignore calcIC;
%ignore calcICB;
%ignore sensToggleOff;
%ignore solveB;
%ignore step;
%ignore run;
%ignore runB;
%ignore resetDiagnosis;
%ignore storeDiagnosis;
%ignore storeDiagnosisB;
%ignore turnOffRootFinding;
%ignore getRootInfo;
%ignore updateAndReinitStatesAndSensitivities;


%newobject amici::Solver::clone;
// Process symbols in header
%include "amici/solver.h"
