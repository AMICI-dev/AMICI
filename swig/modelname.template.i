%define MODULEIMPORT
"
import amici
import datetime
import importlib.util
import os
import sysconfig
from pathlib import Path

ext_suffix = sysconfig.get_config_var('EXT_SUFFIX')
_TPL_MODELNAME = amici._module_from_path(
    'TPL_MODELNAME._TPL_MODELNAME' if __package__ or '.' in __name__
    else '_TPL_MODELNAME',
    Path(__file__).parent / f'_TPL_MODELNAME{ext_suffix}',
)

def _get_import_time():
    return _TPL_MODELNAME._get_import_time()

t_imported = _get_import_time()
t_modified = os.path.getmtime(__file__)
if t_imported < t_modified:
    t_imp_str = datetime.datetime.fromtimestamp(t_imported).isoformat()
    t_mod_str = datetime.datetime.fromtimestamp(t_modified).isoformat()
    module_path = Path(__file__).resolve()
    raise RuntimeError(
        f'Cannot import extension for TPL_MODELNAME from '
        f'{module_path}, because an extension in the same location '
        f'has already been imported, but the file was modified on '
        f'disk. \\nImported at {t_imp_str}\\nModified at {t_mod_str}.\\n'
        'Import the module with a different name or restart the '
        'Python kernel.'
    )
"
%enddef

%module(package="TPL_MODELNAME",moduleimport=MODULEIMPORT) TPL_MODELNAME

%pythoncode %{
# the model-package __init__.py module (will be set during import)
_model_module = None


%}

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
        val.module = _model_module
%}


// Process symbols in header
%include "wrapfunctions.h"
