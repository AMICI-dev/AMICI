"""
AMICI
-----

The AMICI Python module provides functionality for importing SBML or PySB
models and turning them into C++ Python extensions.

:var amici_path:
    absolute root path of the amici repository or Python package
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
"""

import importlib
import os
import re
import sys
from contextlib import suppress, contextmanager
from types import ModuleType as ModelModule
from typing import Optional, Union, Sequence, List


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
                os.path.join(basedir, '.git', 'ORIG_HEAD'), ]
            if os.path.isfile(file)
        ),
        None
    )

    if commitfile:
        with open(commitfile) as f:
            return str(re.search(r'^([\w]*)', f.read().strip()).group())
    return 'unknown'


def _imported_from_setup() -> bool:
    """Check whether this module is imported from `setup.py`"""

    from inspect import getouterframes, currentframe

    # in case we are imported from setup.py, this will be the AMICI package
    # root directory (otherwise it is most likely the Python library directory,
    # we are not interested in)
    package_root = os.path.realpath(os.path.dirname(os.path.dirname(__file__)))

    for frame in getouterframes(currentframe(), context=0):
        # Need to compare the full path, in case a user tries to import AMICI
        # from a module `*setup.py`. Will still cause trouble if some package
        # requires the AMICI extension during its installation, but seems
        # unlikely...
        frame_path = os.path.realpath(os.path.expanduser(frame.filename))
        if frame_path == os.path.join(package_root, 'setup.py'):
            return True

    return False


try:
    from wurlitzer import sys_pipes
except ModuleNotFoundError:
    sys_pipes = suppress


@contextmanager
def _capture_cstdout():
    """Redirect C/C++ stdout to python stdout if python stdout is redirected,
    e.g. in ipython notebook"""
    if sys.stdout == sys.__stdout__:
        yield
    else:
        with sys_pipes():
            yield


# Initialize AMICI paths
amici_path = _get_amici_path()
amiciSwigPath = os.path.join(amici_path, 'swig')
amiciSrcPath = os.path.join(amici_path, 'src')
amiciModulePath = os.path.dirname(__file__)

has_clibs = any([os.path.isfile(os.path.join(amici_path, wrapper))
                 for wrapper in ['amici.py', 'amici_without_hdf5.py']])

AmiciModel = Union['amici.Model', 'amici.ModelPtr']
AmiciSolver = Union['amici.Solver', 'amici.SolverPtr']
AmiciExpData = Union['amici.ExpData', 'amici.ExpDataPtr']
AmiciReturnData = Union['amici.ReturnData', 'amici.ReturnDataPtr']
AmiciExpDataVector = Union['amici.ExpDataPtrVector', Sequence[AmiciExpData]]

# Get version number from file
with open(os.path.join(amici_path, 'version.txt')) as f:
    __version__ = f.read().strip()

__commit__ = _get_commit_hash()

# Import SWIG module and swig-dependent submodules if required and available
if not _imported_from_setup():
    if has_clibs:
        from . import amici
        from .amici import *

        # These module require the swig interface and other dependencies
        from .numpy import ReturnDataView, ExpDataView
        from .pandas import (
            getEdataFromDataFrame,
            getDataObservablesAsDataFrame,
            getSimulationObservablesAsDataFrame,
            getSimulationStatesAsDataFrame,
            getResidualsAsDataFrame
        )

    # These modules don't require the swig interface
    from .sbml_import import SbmlImporter, assignmentRules2observables
    from .ode_export import ODEModel, ODEExporter

    try:
        # Requires Python>=3.8
        from typing import Protocol

        class ModelModule(Protocol):
            """Enable Python static type checking for AMICI-generated model modules"""
            def getModel(self) -> amici.Model:
                pass
    except ImportError:
        pass

hdf5_enabled = 'readSolverSettingsFromHDF5' in dir()


def _get_ptr(obj: Union[AmiciModel, AmiciExpData, AmiciSolver, AmiciReturnData]
             ) -> Union['amici.Model', 'amici.ExpData', 'amici.Solver',
                        'amici.ReturnData']:
    """
    Convenience wrapper that returns the smart pointer pointee, if applicable

    :param obj:
        Potential smart pointer

    :returns:
        Non-smart pointer
    """
    if isinstance(obj, (amici.ModelPtr, amici.ExpDataPtr, amici.SolverPtr,
                        amici.ReturnDataPtr)):
        return obj.get()
    return obj


def runAmiciSimulation(
        model: AmiciModel,
        solver: AmiciSolver,
        edata: Optional[AmiciExpData] = None
) -> 'numpy.ReturnDataView':
    """
    Convenience wrapper around :py:func:`amici.amici.runAmiciSimulation`
    (generated by swig)

    :param model:
        Model instance

`   :param solver:
        Solver instance, must be generated from
        :py:meth:`amici.amici.Model.getSolver`

    :param edata:
        ExpData instance (optional)

    :returns:
        ReturnData object with simulation results
    """
    with _capture_cstdout():
        rdata = amici.runAmiciSimulation(_get_ptr(solver), _get_ptr(edata),
                                         _get_ptr(model))
    return numpy.ReturnDataView(rdata)


def ExpData(*args) -> 'amici.ExpData':
    """
    Convenience wrapper for :py:class:`amici.amici.ExpData` constructors

    :param args: arguments

    :returns: ExpData Instance
    """
    if isinstance(args[0], ReturnDataView):
        return amici.ExpData(_get_ptr(args[0]['ptr']), *args[1:])
    elif isinstance(args[0], (amici.ExpData, amici.ExpDataPtr)):
        # the *args[:1] should be empty, but by the time you read this,
        # the constructor signature may have changed and you are glad this
        # wrapper did not break.
        return amici.ExpData(_get_ptr(args[0]), *args[1:])
    elif isinstance(args[0], (amici.Model, amici.ModelPtr)):
        return amici.ExpData(_get_ptr(args[0]))
    else:
        return amici.ExpData(*args)


def runAmiciSimulations(
        model: AmiciModel,
        solver: AmiciSolver,
        edata_list: AmiciExpDataVector,
        failfast: bool = True,
        num_threads: int = 1,
) -> List['numpy.ReturnDataView']:
    """
    Convenience wrapper for loops of amici.runAmiciSimulation

    :param model: Model instance
    :param solver: Solver instance, must be generated from Model.getSolver()
    :param edata_list: list of ExpData instances
    :param failfast: returns as soon as an integration failure is encountered
    :param num_threads: number of threads to use (only used if compiled
        with openmp)

    :returns: list of simulation results
    """
    with _capture_cstdout():
        edata_ptr_vector = amici.ExpDataPtrVector(edata_list)
        rdata_ptr_list = amici.runAmiciSimulations(_get_ptr(solver),
                                                   edata_ptr_vector,
                                                   _get_ptr(model),
                                                   failfast,
                                                   num_threads)
    return [numpy.ReturnDataView(r) for r in rdata_ptr_list]


def readSolverSettingsFromHDF5(
        file: str,
        solver: AmiciSolver,
        location: Optional[str] = 'solverSettings'
) -> None:
    """
    Convenience wrapper for :py:func:`amici.readSolverSettingsFromHDF5`

    :param file: hdf5 filename
    :param solver: Solver instance to which settings will be transferred
    :param location: location of solver settings in hdf5 file
    """
    amici.readSolverSettingsFromHDF5(file, _get_ptr(solver), location)


def writeSolverSettingsToHDF5(
        solver: AmiciSolver,
        file: Union[str, object],
        location: Optional[str] = 'solverSettings'
) -> None:
    """
    Convenience wrapper for :py:func:`amici.amici.writeSolverSettingsToHDF5`

    :param file: hdf5 filename, can also be object created by
        :py:func:`amici.amici.createOrOpenForWriting`
    :param solver: Solver instance from which settings will stored
    :param location: location of solver settings in hdf5 file
    """
    amici.writeSolverSettingsToHDF5(_get_ptr(solver), file, location)


class add_path:
    """Context manager for temporarily changing PYTHONPATH"""

    def __init__(self, path: str):
        self.path: str = path

    def __enter__(self):
        if self.path:
            sys.path.insert(0, self.path)

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            sys.path.remove(self.path)
        except ValueError:
            pass


def import_model_module(module_name: str,
                        module_path: Optional[str] = None) -> ModelModule:
    """
    Import Python module of an AMICI model

    :param module_name:
        Name of the python package of the model
    :param module_path:
        Absolute or relative path of the package directory
    :return:
        The model module
    """

    # ensure we will find the newly created module
    importlib.invalidate_caches()

    if not os.path.isdir(module_path):
        raise ValueError(f"module_path '{module_path}' is not a directory.")

    module_path = os.path.abspath(module_path)

    # module already loaded?
    if module_name in sys.modules:
        # if a module with that name is already in sys.modules, we remove it,
        # along with all other modules from that package. otherwise, there
        # will be trouble if two different models with the same name are to
        # be imported.
        del sys.modules[module_name]
        # collect first, don't delete while iterating
        to_unload = {loaded_module_name for loaded_module_name in
                     sys.modules.keys() if
                     loaded_module_name.startswith(f"{module_name}.")}
        for m in to_unload:
            del sys.modules[m]

    with add_path(module_path):
        return importlib.import_module(module_name)
