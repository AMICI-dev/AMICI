%module model_jakstat_adjoint_o2

// Add necessary symbols to generated header

%{
#include "wrapfunctions.h"
%}

%rename(getModel) getModelSwig;
%inline %{
 amici::Model_ODE *getModelSwig() {
    return new Model_model_jakstat_adjoint_o2();
}
%}
%ignore getModel;



// Process symbols in header
%include "wrapfunctions.h"
