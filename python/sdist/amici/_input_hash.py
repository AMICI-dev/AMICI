"""Hashes of the inputs of model imports.

A model import is determined by the model, the import options, and the
environment, i.e., the versions of AMICI and of the packages used for code
generation, some environment variables, and, for compiled models, the build
configuration. The hash of these inputs is stored in the marker of the model
directory (see :mod:`amici._model_dir`), so that a later import with the same
inputs can reuse the existing model instead of importing it again.

Inputs that can't be represented reliably, e.g., arbitrary callables, make
the hash ``None``, so that the model is never reused.
Setting the environment variable :data:`NO_MODEL_REUSE_ENV_VAR` to
``1``, ``on``, or ``true`` prevents reusing any model.

Limitations: Changes to the inputs that are not covered here are not
detected, e.g., changes to Python modules imported by a PySB model,
to external tools like BioNetGen or the compiler, or to system libraries.
"""

from __future__ import annotations

import enum
import functools
import hashlib
import importlib.metadata
import json
import logging
import math
import os
import platform
import sys
from pathlib import Path
from typing import Any

import sympy as sp

from .logging import get_logger

logger = get_logger(__name__, logging.WARNING)

#: Environment variable that disables reusing previously imported models
NO_MODEL_REUSE_ENV_VAR = "AMICI_NO_MODEL_REUSE"

#: Placeholder for inputs that can't be hashed
UNHASHABLE = "<unhashable>"

#: Arguments of the import functions that don't affect the generated model
NON_HASHED_ARGS = frozenset(
    {"self", "output_dir", "verbose", "compile", "skip_if_unchanged"}
)

# import options that are hashed as sets, since the importers only check
#  membership; remove an option here if its order starts to matter
_UNORDERED_OPTIONS = frozenset({"fixed_parameters", "hardcode_symbols"})

# environment variables that affect building model extensions
_BUILD_ENV_VARS = (
    "AMICI_BLAS_USE_SCIPY_OPENBLAS",
    "AMICI_CXXFLAGS",
    "AMICI_LDFLAGS",
    "BLAS_CFLAGS",
    "BLAS_LIBS",
    "CC",
    "CMAKE_PREFIX_PATH",
    "CXX",
    "CXXFLAGS",
    "ENABLE_AMICI_DEBUGGING",
    "ENABLE_GCOV_COVERAGE",
    "LDFLAGS",
    "MKL_INCDIR",
    "MKL_LIB",
    "MKLROOT",
    "SWIG",
)

# parts of the amici package that affect code generation
_CODEGEN_SOURCES = (
    "_input_hash.py",
    "_model_dir.py",
    "_symbolic",
    "constants.py",
    "exporters",
    "importers",
)


class _UnhashableError(Exception):
    """An input can't be represented reliably."""


def compute_input_hash(
    model: str | None,
    options: dict[str, Any],
    *,
    packages: tuple[str, ...],
    build: bool,
) -> tuple[str | None, dict[str, Any]]:
    """Compute the hash of the inputs of a model import.

    :param model:
        Deterministic serialization of the model, e.g., the SBML document,
        or ``None`` if the model can't be serialized.
    :param options:
        Import options, i.e., anything else that affects the generated code.
    :param packages:
        Distribution names of the packages whose versions affect the
        generated code, in addition to sympy.
    :param build:
        Whether the build configuration of model extensions is relevant,
        i.e., for SUNDIALS models, even if not compiled right away.
    :returns:
        The hash, and the hashed inputs, with the model replaced by its hash.
        The hash is ``None`` if any input can't be hashed or if reusing
        models is disabled.
    """
    hashable = True
    normalized_options = {}
    for key, value in sorted(options.items()):
        if key in _UNORDERED_OPTIONS and value is not None:
            value = set(value)
        try:
            normalized_options[key] = _normalize(value)
        except (_UnhashableError, OSError) as e:
            logger.debug(f"Cannot hash import option {key!r}: {e}")
            normalized_options[key] = UNHASHABLE
            hashable = False

    inputs = {
        "model": UNHASHABLE
        if model is None
        else hashlib.sha256(model.encode()).hexdigest(),
        "options": normalized_options,
        "environment": _get_environment(("sympy", *packages), build=build),
    }
    if not hashable or model is None or _model_reuse_disabled():
        return None, inputs

    return _hash(inputs), inputs


def _model_reuse_disabled() -> bool:
    """Check whether reusing models is disabled via the environment."""
    return os.environ.get(NO_MODEL_REUSE_ENV_VAR, "").lower() in (
        "1",
        "on",
        "true",
    )


def _hash(value: Any) -> str:
    """Hash a JSON-serializable value."""
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, allow_nan=False).encode()
    ).hexdigest()


def _normalize(value: Any) -> Any:
    """Convert a value to a deterministic JSON-serializable representation.

    :raises _UnhashableError: If the value can't be represented reliably.
    """
    from .importers.utils import MeasurementChannel, _default_simplify

    if value is None or isinstance(value, bool | int | str):
        return value
    if isinstance(value, float):
        # distinguish from equal int values and handle non-finite values
        return {"float": value.hex() if math.isfinite(value) else repr(value)}
    if isinstance(value, enum.Enum):
        return {
            "enum": type(value).__qualname__,
            "value": _normalize(value.value),
        }
    if isinstance(value, sp.Basic):
        return {"sympy": sp.srepr(_replace_dummies(value))}
    if isinstance(value, list | tuple):
        return [_normalize(x) for x in value]
    if isinstance(value, set | frozenset):
        return sorted(
            (_normalize(x) for x in value), key=lambda x: json.dumps(x)
        )
    if isinstance(value, dict):
        # lists instead of tuples, to compare equal after a JSON round trip
        items = [[_normalize(k), _normalize(v)] for k, v in value.items()]
        return {
            "dict": sorted(items, key=lambda kv: json.dumps(kv[0])),
        }
    if isinstance(value, Path):
        # the generated code may depend on the file contents
        if not value.is_file():
            raise _UnhashableError(f"{value} is not a file")
        return {
            "path": str(value),
            "sha256": hashlib.sha256(value.read_bytes()).hexdigest(),
        }
    if value is _default_simplify:
        return {"callable": "default_simplify"}
    if isinstance(value, MeasurementChannel) or _is_de_model_component(value):
        return {
            "class": type(value).__qualname__,
            "attributes": _normalize(vars(value)),
        }
    if hasattr(value, "model_dump") and hasattr(type(value), "model_fields"):
        # pydantic models
        return {
            "class": type(value).__qualname__,
            "attributes": _normalize(value.model_dump(mode="json")),
        }
    raise _UnhashableError(f"unsupported type {type(value).__qualname__}")


def _replace_dummies(expr: sp.Basic) -> sp.Basic:
    """Replace dummy symbols by regular symbols.

    The representation of dummy symbols, which include, e.g., PySB model
    components, contains a process-dependent index. They are replaced by
    regular symbols named after their class and name.
    """
    return expr.xreplace(
        {
            dummy: sp.Symbol(
                f"{type(dummy).__qualname__}:{dummy.name}",
                **dummy.assumptions0,
            )
            for dummy in expr.atoms(sp.Dummy)
        }
    )


def _is_de_model_component(value: Any) -> bool:
    """Check whether a value is a model component, e.g., an event."""
    from ._symbolic.de_model_components import ModelQuantity

    return isinstance(value, ModelQuantity)


def _get_environment(packages: tuple[str, ...], build: bool) -> dict:
    """Get the parts of the environment that affect an import."""
    from . import __version__
    from .exporters.sundials.cxxcodeprinter import _extract_cse_from_env

    environment = {
        "amici_version": __version__,
        "amici_sources": _get_amici_source_digest(),
        "python": f"{platform.python_implementation()}-"
        f"{sys.version_info.major}.{sys.version_info.minor}",
        "packages": {
            package: _get_package_version(package)
            for package in sorted(packages)
        },
        "extract_cse": _extract_cse_from_env(),
    }
    if build:
        environment["build"] = {
            "env": {
                var: os.environ.get(var) for var in sorted(_BUILD_ENV_VARS)
            },
            "libamici": _get_amici_library_stats(),
        }
    return environment


def _get_package_version(package: str) -> str | None:
    """Get the installed version of a package, or ``None``."""
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return None


@functools.cache
def _get_amici_source_digest() -> str:
    """Hash the amici sources that affect code generation.

    This detects changes to editable or modified installations which do not
    change the AMICI version.
    """
    package_dir = Path(__file__).parent
    digest = hashlib.sha256()
    for entry in _CODEGEN_SOURCES:
        path = package_dir / entry
        files = [path] if path.is_file() else path.rglob("*")
        for file in sorted(files):
            if (
                not file.is_file()
                or "__pycache__" in file.parts
                or file.suffix == ".pyc"
            ):
                continue
            digest.update(file.relative_to(package_dir).as_posix().encode())
            digest.update(b"\0")
            digest.update(hashlib.sha256(file.read_bytes()).digest())
    return digest.hexdigest()


def _get_amici_library_stats() -> list[list]:
    """Get size and modification time of the AMICI C++ libraries.

    Model extensions are built against those libraries. This detects
    rebuilds of editable installations.
    """
    lib_dir = Path(__file__).parent / "_installation"
    stats = []
    for file in sorted(lib_dir.glob("lib*/*amici*")):
        if file.is_file():
            stat = file.stat()
            stats.append(
                [
                    file.relative_to(lib_dir).as_posix(),
                    stat.st_size,
                    stat.st_mtime_ns,
                ]
            )
    return stats


def get_changed_inputs(old: dict | None, new: dict) -> list[str]:
    """Get the names of the inputs that differ between two imports.

    :param old: The inputs of the previous import, if known.
    :param new: The inputs of the current import.
    :returns: Sorted names of changed inputs, e.g., ``options.simplify``.
    """
    if not isinstance(old, dict):
        return ["<unknown>"]
    changed = []
    for key in sorted(old.keys() | new.keys()):
        old_value, new_value = old.get(key), new.get(key)
        if isinstance(old_value, dict) and isinstance(new_value, dict):
            changed.extend(
                f"{key}.{sub_key}"
                for sub_key in sorted(old_value.keys() | new_value.keys())
                if old_value.get(sub_key) != new_value.get(sub_key)
                or sub_key not in old_value
                or sub_key not in new_value
            )
        elif old_value != new_value:
            changed.append(key)
    return changed
