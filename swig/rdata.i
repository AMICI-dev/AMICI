%module rdata

// Add necessary symbols to generated header
%{
#include "amici/rdata.h"
using namespace amici;
%}

%ignore invalidate;
%ignore invalidateLLH;
%ignore invalidateSLLH;
%ignore applyChainRuleFactorToSimulationResults;
%ignore storeJacobianAndDerivativeInReturnData;
%ignore initializeObjectiveFunction;
%ignore processPreEquilibration;
%ignore processPostEquilibration;
%ignore processForwardProblem;
%ignore processBackwardProblem;
%ignore processSolver;

// Process symbols in header
%include "amici/rdata.h"
