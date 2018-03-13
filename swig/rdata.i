%module rdata

// Add necessary symbols to generated header
%{
#include "amici/rdata.h"
using namespace amici;
%}

// Process symbols in header
%include "amici/rdata.h"
