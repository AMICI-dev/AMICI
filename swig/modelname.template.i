%module TPL_MODELNAME
%import amici.i
// Add necessary symbols to generated header

%{
#include "wrapfunctions.h"
#include "amici/model_ode.h"
#include "amici/model_dae.h"
using namespace amici;
%}

// store the time a module was imported
%{
#include <chrono>
static std::chrono::time_point<std::chrono::system_clock> _module_import_time;

static double _get_import_time() {
    auto epoch = _module_import_time.time_since_epoch();
    return std::chrono::duration<double>(epoch).count();
}
%}

static double _get_import_time();

%init %{
    _module_import_time = std::chrono::system_clock::now();
%}


// Make model module accessible from the model
%feature("pythonappend") amici::generic_model::getModel %{
    if '.' in __name__:
        import sys
        val.module = sys.modules['.'.join(__name__.split('.')[:-1])]
%}


// Process symbols in header
%include "wrapfunctions.h"
