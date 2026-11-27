"""Locate the build tools of the current Python environment."""

import os
import shutil
import sysconfig
from collections.abc import MutableMapping


def get_env_swig() -> str | None:
    """Get the path to the SWIG executable of the current Python environment.

    :returns: The path to the ``swig`` executable in the scripts directory of
        the running interpreter, or ``None`` if there is none.
    """
    # Not the binary in `swig.BIN_DIR`: it requires `SWIG_LIB`, which the
    # wrapper script sets.
    return shutil.which("swig", path=sysconfig.get_path("scripts"))


def _get_env_tool_dir(name: str) -> str | None:
    """Get the directory of a build tool of the current Python environment.

    :param name: ``cmake`` or ``ninja``.
    :returns: The directory containing the executable provided by the
        respective Python package, or by the scripts directory of the
        running interpreter, or ``None`` if there is none.
    """
    candidates = []
    try:
        if name == "cmake":
            import cmake

            candidates.append(cmake.CMAKE_BIN_DIR)
        elif name == "ninja":
            import ninja

            candidates.append(ninja.BIN_DIR)
    except Exception:  # noqa: BLE001, S110
        # e.g., a broken installation: fall back to the scripts directory
        pass
    candidates.append(sysconfig.get_path("scripts"))

    for candidate in candidates:
        if candidate and shutil.which(name, path=candidate):
            return candidate
    return None


def use_env_build_tools(
    env: MutableMapping[str, str] | None = None,
) -> None:
    """Use the build tools of the current Python environment.

    Prepends the directories of the CMake and Ninja executables of the
    current Python environment to ``PATH``, unless ``PATH`` already resolves
    to them, so that they are used even if the environment is not activated.
    The Ninja directory is usually the environment's scripts directory, whose
    other executables then also take precedence.
    Sets ``SWIG`` to the SWIG executable of the current environment, unless
    already set.

    :param env: The environment variables to update. Defaults to
        :data:`os.environ`.
    """
    if env is None:
        env = os.environ

    path = env.get("PATH", os.defpath)
    prepend = []
    for name in ("cmake", "ninja"):
        if (tool_dir := _get_env_tool_dir(name)) is None:
            continue
        env_tool = shutil.which(name, path=tool_dir)
        current_tool = shutil.which(name, path=path)
        if current_tool is not None and os.path.samefile(
            current_tool, env_tool
        ):
            continue
        prepend.append(tool_dir)
    # ninja's directory may be the scripts directory, which also contains a
    #  cmake wrapper, so it must come after cmake's directory
    # an empty PATH entry would mean the current directory
    env["PATH"] = os.pathsep.join([*prepend, path] if path else prepend)

    if "SWIG" not in env and (swig := get_env_swig()):
        env["SWIG"] = swig
