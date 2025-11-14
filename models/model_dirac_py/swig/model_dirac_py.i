%define MODULEIMPORT
"
import amici
import datetime
import importlib.util
import os
import sysconfig
from pathlib import Path

ext_suffix = sysconfig.get_config_var('EXT_SUFFIX')
extension_path = Path(__file__).parent / f'_model_dirac_py{ext_suffix}'
_model_dirac_py = amici._module_from_path(
    'model_dirac_py._model_dirac_py' if __package__ or '.' in __name__
    else '_model_dirac_py',
    extension_path,
)

def _get_import_time():
    return _model_dirac_py._get_import_time()

t_imported = _get_import_time()
t_modified = os.path.getmtime(__file__)
if t_imported < t_modified:
    t_imp_str = datetime.datetime.fromtimestamp(t_imported).isoformat()
    t_mod_str = datetime.datetime.fromtimestamp(t_modified).isoformat()
    module_path = Path(__file__).resolve()
    raise RuntimeError(
        f'Cannot import extension for model_dirac_py from '
        f'{module_path}, because an extension in the same location '
        f'has already been imported, but the file was modified on '
        f'disk. \\nImported at {t_imp_str}\\nModified at {t_mod_str}.\\n'
        'Import the module with a different name or restart the '
        'Python kernel.'
    )
"
%enddef

%module(package="model_dirac_py",moduleimport=MODULEIMPORT) model_dirac_py

// store swig version
%constant int SWIG_VERSION_MAJOR = (SWIG_VERSION >> 16);
%constant int SWIG_VERSION_MINOR = ((SWIG_VERSION >> 8) & 0xff);
%constant int SWIG_VERSION_PATCH = (SWIG_VERSION & 0xff);

%pythoncode %{
# SWIG version used to build the model extension as `(major, minor, patch)`
_SWIG_VERSION = (SWIG_VERSION_MAJOR, SWIG_VERSION_MINOR, SWIG_VERSION_PATCH)

if (amici_swig := amici.amici._SWIG_VERSION) != (model_swig := _SWIG_VERSION):
    import warnings
    warnings.warn(
        f"SWIG version mismatch between amici ({amici_swig}) and model "
        f"({model_swig}). This may lead to unexpected behavior. "
        "In that case, please recompile the model with swig=="
        f"{amici_swig[0]}.{amici_swig[1]}.{amici_swig[2]} or rebuild amici "
        f"with swig=={model_swig[0]}.{model_swig[1]}.{model_swig[2]}.",
        RuntimeWarning,
        stacklevel=2,
    )
%}

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
static std::chrono::time_point<std::chrono::system_clock> _module_import_time = std::chrono::system_clock::now();

static double _get_import_time() {
    auto epoch = _module_import_time.time_since_epoch();
    return std::chrono::duration<double>(epoch).count();
}
%}

static double _get_import_time();

%init %{
    // NOTE: from SWIG 4.4.0 onwards, %init code is executed every time the
    // module is executed - not only on first import
    // This code ends up in `SWIG_mod_exec`.
%}


// Make model module accessible from the model
%feature("pythonappend") amici::generic_model::get_model %{
    if '.' in __name__:
        val.module = _model_module
%}


// Process symbols in header
%include "wrapfunctions.h"
