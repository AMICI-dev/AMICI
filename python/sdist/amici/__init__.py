"""
AMICI
-----

The AMICI Python module provides functionality for importing SBML or PySB
models and turning them into C++ Python extensions.
"""

import contextlib
import importlib.util
import importlib
import os
import re
import sys
from pathlib import Path
from types import ModuleType
from typing import Any
from collections.abc import Callable
import warnings


def _get_amici_path():
    """
    Determine package installation path, or, if used directly from git
    repository, get repository root
    """
    basedir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    if os.path.exists(os.path.join(basedir, ".git")):
        return os.path.abspath(basedir)
    return os.path.dirname(__file__)


def _get_commit_hash():
    """Get commit hash from file"""
    basedir = os.path.dirname(os.path.dirname(os.path.dirname(amici_path)))
    commitfile = next(
        (
            file
            for file in [
                os.path.join(basedir, ".git", "FETCH_HEAD"),
                os.path.join(basedir, ".git", "ORIG_HEAD"),
            ]
            if os.path.isfile(file)
        ),
        None,
    )

    if commitfile:
        with open(commitfile) as f:
            return str(re.search(r"^([\w]*)", f.read().strip()).group())
    return "unknown"


def _imported_from_setup() -> bool:
    """Check whether this module is imported from `setup.py`"""

    from inspect import currentframe, getouterframes
    from os import sep

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
        if frame_path == os.path.join(
            package_root, "setup.py"
        ) or frame_path.endswith(f"{sep}setuptools{sep}build_meta.py"):
            return True

    return False


# Initialize AMICI paths
#: absolute root path of the amici repository or Python package
amici_path = _get_amici_path()
#: absolute path of the amici swig directory
amiciSwigPath = os.path.join(amici_path, "swig")
#: absolute path of the amici source directory
amiciSrcPath = os.path.join(amici_path, "src")
#: absolute root path of the amici module
amiciModulePath = os.path.dirname(__file__)
#: boolean indicating if this is the full package with swig interface or
#  the raw package without extension
has_clibs: bool = any(
    os.path.isfile(os.path.join(amici_path, wrapper))
    for wrapper in ["amici.py", "amici_without_hdf5.py"]
)
#: boolean indicating if amici was compiled with hdf5 support
hdf5_enabled: bool = False

# Get version number from file
with open(os.path.join(amici_path, "version.txt")) as f:
    __version__ = f.read().strip()

__commit__ = _get_commit_hash()

# Import SWIG module and swig-dependent submodules if required and available
if not _imported_from_setup():
    if has_clibs:
        # prevent segfaults under pytest
        #  see also:
        #  https://github.com/swig/swig/issues/2881
        #  https://github.com/AMICI-dev/AMICI/issues/2565
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                category=DeprecationWarning,
                message="builtin type .* has no __module__ attribute",
            )
            from . import amici

        from .amici import *

        # has to be done before importing readSolverSettingsFromHDF5
        #  from .swig_wrappers
        hdf5_enabled = "readSolverSettingsFromHDF5" in dir()
        # These modules require the swig interface and other dependencies
        from .numpy import ExpDataView, ReturnDataView  # noqa: F401
        from .pandas import *
        from .swig_wrappers import *

    # These modules don't require the swig interface
    from typing import Protocol, runtime_checkable

    from .de_export import DEExporter  # noqa: F401
    from .sbml_import import (  # noqa: F401
        SbmlImporter,
        assignmentRules2observables,
    )

    try:
        from .jax import JAXModel
    except (ImportError, ModuleNotFoundError):
        JAXModel = object

    @runtime_checkable
    class ModelModule(Protocol):  # noqa: F811
        """Type of AMICI-generated model modules.

        To enable static type checking."""

        def getModel(self) -> amici.Model:
            """Create a model instance."""
            ...

        def get_model(self) -> amici.Model:
            """Create a model instance."""
            ...

    AmiciModel = Union[amici.Model, amici.ModelPtr]
else:
    ModelModule = ModuleType


class add_path:
    """Context manager for temporarily changing PYTHONPATH.

    Add a path to the PYTHONPATH for the duration of the context manager.
    """

    def __init__(self, path: str | Path):
        self.path: str = str(path)

    def __enter__(self):
        if self.path:
            sys.path.insert(0, self.path)

    def __exit__(self, exc_type, exc_value, traceback):
        with contextlib.suppress(ValueError):
            sys.path.remove(self.path)


class set_path:
    """Context manager for temporarily changing PYTHONPATH.

    Set the PYTHONPATH to a given path for the duration of the context manager.
    """

    def __init__(self, path: str | Path):
        self.path: str = str(path)

    def __enter__(self):
        self.orginal_path = sys.path.copy()
        sys.path = [self.path]

    def __exit__(self, exc_type, exc_value, traceback):
        sys.path = self.orginal_path


def _module_from_path(module_name: str, module_path: Path | str) -> ModuleType:
    """Import a module from a given path.

    Import a module from a given path. The module is not added to
    `sys.modules`. The `_self` attribute of the module is set to the module
    itself.

    :param module_name:
        Name of the module.
    :param module_path:
        Path to the module file. Absolute or relative to the current working
        directory.
    """
    module_path = Path(module_path).resolve()
    if not module_path.is_file():
        raise ModuleNotFoundError(f"Module file not found: {module_path}")
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    module = importlib.util.module_from_spec(spec)
    module._self = module
    spec.loader.exec_module(module)
    return module


def import_model_module(
    module_name: str, module_path: Path | str
) -> ModelModule:
    """
    Import Python module of an AMICI model

    :param module_name:
        Name of the python package of the model
    :param module_path:
        Absolute or relative path of the package directory
    :return:
        The model module
    """
    model_root = str(module_path)

    # ensure we will find the newly created module
    importlib.invalidate_caches()

    if not os.path.isdir(module_path):
        raise ValueError(f"module_path '{model_root}' is not a directory.")

    module_path = Path(model_root, module_name, "__init__.py")

    # We may want to import a matlab-generated model where the extension
    #  is in a different directory. This is not a regular use case. It's only
    #  used in the amici tests and can be removed at any time.
    #  The models (currently) use the default swig-import and require
    #  modifying sys.path.
    module_path_matlab = Path(model_root, f"{module_name}.py")
    if not module_path.is_file() and module_path_matlab.is_file():
        with set_path(model_root):
            # prevent segfaults under pytest
            #  see also:
            #  https://github.com/swig/swig/issues/2881
            #  https://github.com/AMICI-dev/AMICI/issues/2565
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    category=DeprecationWarning,
                    message="builtin type .* has no __module__ attribute",
                )
                return _module_from_path(module_name, module_path_matlab)

    module = _module_from_path(module_name, module_path)
    module._self = module
    return module


class AmiciVersionError(RuntimeError):
    """Error thrown if an AMICI model is loaded that is incompatible with
    the installed AMICI base package"""

    pass


def _get_default_argument(func: Callable, arg: str) -> Any:
    """Get the default value of the given argument in the given function."""
    import inspect

    signature = inspect.signature(func)
    if (
        default := signature.parameters[arg].default
    ) is not inspect.Parameter.empty:
        return default
    raise ValueError(f"No default value for argument {arg} of {func}.")
