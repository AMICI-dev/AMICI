"""
Functionality for building the C++ extensions of an amici-created model
package.
"""

import os
import re
import shutil
import subprocess
import sys
import sysconfig
import warnings
from pathlib import Path


def build_model_extension(
    package_dir: str | Path,
    verbose: bool | int | None = False,
    compiler: str | None = None,
    extra_msg: str | None = None,
) -> None:
    """
    Compile the model extension of an amici-created model package.

    The SWIG executable is taken from the ``SWIG`` environment variable,
    if set, otherwise from the scripts directory of the current Python
    environment, if available, otherwise from ``PATH``.

    :param package_dir:
        Directory of the model package to be compiled. I.e., the directory
        containing the `setup.py` file.

    :param verbose:
        Make model compilation verbose.

    :param compiler:
        Absolute path to the compiler executable to be used to build the Python
        extension, e.g. ``/usr/bin/clang``.

    :param extra_msg:
        Additional message to be printed in case of a failed build.
    """
    # setup.py assumes it is run from within the model directory
    package_dir = Path(package_dir)
    script_args = [sys.executable, package_dir / "setup.py"]

    if verbose:
        script_args.append("--verbose")
    else:
        script_args.append("--quiet")

    script_args.extend(
        [
            "build_ext",
            f"--build-lib={package_dir}",
            # This is generally not required, but helps to reduce the path
            # length of intermediate build files, that may easily become
            # problematic on Windows, due to its ridiculous 255-character path
            # length limit.
            f"--build-temp={package_dir / 'build'}",
        ]
    )

    env = os.environ.copy()
    if compiler is not None:
        # CMake will use the compiler specified in the CXX environment variable
        env["CXX"] = compiler

    if "SWIG" not in env and (swig := _get_env_swig()):
        # Use the SWIG of the current Python environment, which is most
        # likely the one used for building amici, instead of the first one
        # on PATH
        env["SWIG"] = swig
    if "SWIG" in env:
        _check_swig_version(env["SWIG"])

    # distutils.core.run_setup looks nicer, but does not let us check the
    # result easily
    try:
        result = subprocess.run(
            script_args,
            cwd=str(package_dir),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=True,
            env=env,
        )
    except subprocess.CalledProcessError as e:
        print(e.output.decode("utf-8"))
        print("Failed building the model extension.")
        if extra_msg:
            print(f"Note: {extra_msg}")
        raise

    if verbose:
        print(result.stdout.decode("utf-8"))


def _get_env_swig() -> str | None:
    """Get the path to the SWIG executable of the current Python environment.

    :returns: The path to the ``swig`` executable in the scripts directory of
        the running interpreter, or ``None`` if there is none.
    """
    return shutil.which("swig", path=sysconfig.get_path("scripts"))


def _check_swig_version(swig: str) -> None:
    """Warn if the given SWIG version differs from the one used for amici.

    Mismatching SWIG versions may lead to subtle errors at runtime.
    """
    try:
        from amici.sim.sundials import _SWIG_VERSION as amici_swig
    except ImportError:
        return

    try:
        output = subprocess.run(
            [swig, "-version"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError) as e:
        warnings.warn(
            f"Failed to determine the SWIG version of `{swig}`: {e}",
            RuntimeWarning,
            stacklevel=3,
        )
        return

    if not (match := re.search(r"SWIG Version (\d+)\.(\d+)\.(\d+)", output)):
        warnings.warn(
            f"Failed to determine the SWIG version of `{swig}`.",
            RuntimeWarning,
            stacklevel=3,
        )
        return

    model_swig = tuple(map(int, match.groups()))
    if model_swig != amici_swig:
        warnings.warn(
            f"SWIG version mismatch between amici "
            f"({'.'.join(map(str, amici_swig))}) and the SWIG executable "
            f"`{swig}` ({'.'.join(map(str, model_swig))}) used for "
            "building the model extension. This may lead to unexpected "
            "behavior. Set the `SWIG` environment variable to a matching "
            "SWIG executable.",
            RuntimeWarning,
            stacklevel=3,
        )
