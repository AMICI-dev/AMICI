%module rdata

// Add necessary symbols to generated header
%{
#include "rdata.h"
using namespace amici;
%}

// Process symbols in header
%include "rdata.h"
