%module logging

%nodefaultctor amici::LogItem;

// Add necessary symbols to generated header
%{
#include "amici/logging.h"
%}

// Process symbols in header
%include "amici/logging.h"
