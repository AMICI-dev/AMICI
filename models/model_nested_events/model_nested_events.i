%module model_nested_events

// Add necessary symbols to generated header

%{
#include "wrapfunctions.h"
%}

%rename(getModel) getModelSwig;
%inline %{
 amici::Model_ODE *getModelSwig() {
    return new Model_model_nested_events();
}
%}
%ignore getModel;



// Process symbols in header
%include "wrapfunctions.h"
