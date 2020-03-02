"""
AMICI
-----
The AMICI Python module provides functionality for importing SBML models and
turning them into C++ Python extensions.

Getting started:
```
# creating a extension module for an SBML model:
import amici
amiSbml = amici.SbmlImporter('mymodel.sbml')
amiSbml.sbml2amici('modelName', 'outputDirectory')

# using the created module (set python path)
import modelName
help(modelName)
```

:var amici_path:
    absolute root path of the amici repository
:var amiciSwigPath:
    absolute path of the amici swig directory
:var amiciSrcPath:
    absolute path of the amici source directory
:var amiciModulePath:
    absolute root path of the amici module
:var hdf5_enabled:
    boolean indicating if amici was compiled with hdf5 support
:var has_clibs:
    boolean indicating if this is the full package with swig interface or
    the raw package without
:var capture_cstdout:
    context to redirect C/C++ stdout to python stdout if python stdout was
    redirected (doing nothing if not redirected).
"""

import os
import re
import sys
from contextlib import suppress


def _get_amici_path():
    """
    Determine package installation path, or, if used directly from git
    repository, get repository root
    """
    basedir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    if os.path.exists(os.path.join(basedir, '.git')):
        return os.path.abspath(basedir)
    return os.path.dirname(__file__)


def _get_commit_hash():
    """Get commit hash from file"""
    basedir = os.path.dirname(os.path.dirname(os.path.dirname(amici_path)))
    commitfile = next(
        (
            file for file in [
                os.path.join(basedir, '.git', 'FETCH_HEAD'),
                os.path.join(basedir, '.git', 'ORIG_HEAD'),
            ]
            if os.path.isfile(file)
        ),
        None
    )

    if commitfile:
        with open(commitfile) as f:
            return str(re.search(r'^([\w]*)', f.read().strip()).group())
    return 'unknown'


# redirect C/C++ stdout to python stdout if python stdout is redirected,
# e.g. in ipython notebook
capture_cstdout = suppress
if sys.stdout != sys.__stdout__:
    try:
        from wurlitzer import sys_pipes as capture_cstdout
    except ModuleNotFoundError:
        pass

hdf5_enabled = False
has_clibs = False

try:
    # Try importing AMICI SWIG-interface with HDF5 enabled
    from . import amici
    from .amici import *
    hdf5_enabled = True
    has_clibs = True
except ModuleNotFoundError:
    # import from setuptools or installation with `--no-clibs`
    pass
except (ImportError, AttributeError) as e:
    # No such module exists or there are some dynamic linking problems
    if isinstance(e, AttributeError) or "cannot import name" in str(e):
        # No such module exists (ImportError),
        #  or python tries to import a HDF5 function from the non-hdf5
        #  swig interface (AttributeError):
        #  try importing AMICI SWIG-interface without HDF5
        try:
            from . import amici_without_hdf5 as amici
            from .amici_without_hdf5 import *
            has_clibs = True
        except ModuleNotFoundError:
            # import from setuptools or installation with `--no-clibs`
            pass
        except ImportError as e:
            if "cannot import name" in str(e):
                # No such module exists
                # this probably means, the model was imported during setuptools
                # `setup` or after an installation with `--no-clibs`.
                pass
            else:
                # Probably some linking problem that we don't want to hide
                raise e
    else:
        # Probably some linking problem that we don't want to hide
        raise e

# Initialize AMICI paths
amici_path = _get_amici_path()
amiciSwigPath = os.path.join(amici_path, 'swig')
amiciSrcPath = os.path.join(amici_path, 'src')
amiciModulePath = os.path.dirname(__file__)

# Get version number from file
with open(os.path.join(amici_path, 'version.txt')) as f:
    __version__ = f.read().strip()

__commit__ = _get_commit_hash()

try:
    # These module require the swig interface and other dependencies which will
    # be installed if the the AMICI package was properly installed. If not,
    # AMICI was probably imported from setup.py and we don't need those.
    from .sbml_import import SbmlImporter, assignmentRules2observables
    from .numpy import ReturnDataView, ExpDataView
    from .pandas import getEdataFromDataFrame, \
        getDataObservablesAsDataFrame, getSimulationObservablesAsDataFrame, \
        getSimulationStatesAsDataFrame, getResidualsAsDataFrame
    from .ode_export import ODEModel, ODEExporter
except (ModuleNotFoundError, ImportError) as e:
    # import from setuptools or installation with `--no-clibs`
    if has_clibs:
        # cannot raise as we may also end up here when installing from an
        # already installed in-source installation without all requirements
        # installed (e.g. fresh virtualenv)
        print(f'Suppressing error {str(e)}')


def runAmiciSimulation(model, solver, edata=None):
    """ Convenience wrapper around amici.runAmiciSimulation (generated by swig)

    Arguments:
        model: Model instance
        solver: Solver instance, must be generated from Model.getSolver()
        edata: ExpData instance (optional)

    Returns:
        ReturnData object with simulation results

    Raises:

    """
    if edata and edata.__class__.__name__ == 'ExpDataPtr':
        edata = edata.get()

    with capture_cstdout():
        rdata = amici.runAmiciSimulation(solver.get(), edata, model.get())
    return numpy.ReturnDataView(rdata)


def ExpData(*args):
    """ Convenience wrapper for ExpData constructors

    Arguments:
        args: arguments

    Returns:
        ExpData Instance

    Raises:

    """
    if isinstance(args[0], ReturnDataView):
        return amici.ExpData(args[0]['ptr'].get(), *args[1:])
    elif isinstance(args[0], ExpDataPtr):
        # the *args[:1] should be empty, but by the time you read this,
        # the constructor signature may have changed and you are glad this
        # wrapper did not break.
        return amici.ExpData(args[0].get(), *args[:1])
    elif isinstance(args[0], ModelPtr):
        return amici.ExpData(args[0].get())
    else:
        return amici.ExpData(*args)


def runAmiciSimulations(model, solver, edata_list, failfast=True,
                        num_threads=1):
    """ Convenience wrapper for loops of amici.runAmiciSimulation

    Arguments:
        model: Model instance
        solver: Solver instance, must be generated from Model.getSolver()
        edata_list: list of ExpData instances
        failfast: returns as soon as an integration failure is encountered
        num_threads: number of threads to use
                     (only used if compiled with openmp)

    Returns:
        list of ReturnData objects with simulation results

    Raises:

    """
    with capture_cstdout():
        edata_ptr_vector = amici.ExpDataPtrVector(edata_list)
        rdata_ptr_list = amici.runAmiciSimulations(solver.get(),
                                                   edata_ptr_vector,
                                                   model.get(),
                                                   failfast,
                                                   num_threads)
    return [numpy.ReturnDataView(r) for r in rdata_ptr_list]
