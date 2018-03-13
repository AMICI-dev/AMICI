%module model_dirac

// Add necessary symbols to generated header

%{
#include "wrapfunctions.h"
%}

%rename(getModel) getModelSwig;
%inline %{
 amici::Model_ODE *getModelSwig() {
    return new Model_model_dirac();
}
%}
%ignore getModel;



// Process symbols in header
%include "wrapfunctions.h"
