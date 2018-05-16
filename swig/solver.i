%module solver

// Add necessary symbols to generated header
%{
#include "amici/solver.h"
using namespace amici;
%}

%rename(equals) operator==;

// Process symbols in header
%include "amici/solver.h"
