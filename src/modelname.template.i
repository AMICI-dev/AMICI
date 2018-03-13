%module TPL_MODELNAME

// Add necessary symbols to generated header

%{
#include "wrapfunctions.h"
%}

%rename(getModel) getModelSwig;
%inline %{
 amici::Model_TPL_MODELTYPE *getModelSwig() {
    return new Model_TPL_MODELNAME();
}
%}
%ignore getModel;



// Process symbols in header
%include "wrapfunctions.h"
