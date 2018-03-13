%module model_robertson

// Add necessary symbols to generated header

%{
#include "wrapfunctions.h"
%}

%rename(getModel) getModelSwig;
%inline %{
 amici::Model_DAE *getModelSwig() {
    return new Model_model_robertson();
}
%}
%ignore getModel;



// Process symbols in header
%include "wrapfunctions.h"
