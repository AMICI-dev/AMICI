"""Setuptools file for creating AMICI module

This file is based on setuptools alone and does not require CMake.
All sources are compiled anew.

This file expects to be run from within its directory.

Non-python-package requirements:
- swig>=3.0
- Optional: hdf5 libraries and headers
"""

import os
import sys
from pathlib import Path

from cmake_build_extension import CMakeExtension
from setuptools import setup
import shutil

# Add containing directory to path, as we need some modules from the AMICI
# package already for installation
sys.path.insert(0, os.path.dirname(__file__))

from amici.custom_commands import (
    AmiciBuildCMakeExtension,
    AmiciBuildPy,
    AmiciDevelop,
    AmiciInstall,
    AmiciInstallLib,
    AmiciSDist,
)


def get_extensions():
    """Get required C(++) extensions for build_ext"""
    # CMake prefix path for finding FindXXX.cmake to find SuiteSparse
    #  components
    install_dir = Path(__file__).parent / "amici"
    prefix_path = install_dir.absolute()
    AmiciBuildCMakeExtension.extend_cmake_prefix_path(str(prefix_path))

    # Used by all extensions
    global_cmake_configure_options = [
        "-DCMAKE_VERBOSE_MAKEFILE=ON",
        f"-DCMAKE_MODULE_PATH={prefix_path.as_posix()}",
    ]

    debug_build = os.getenv("ENABLE_AMICI_DEBUGGING", "").lower() in [
        "1",
        "true",
    ] or os.getenv("ENABLE_GCOV_COVERAGE", "").lower() in ["1", "true"]
    build_type = "Debug" if debug_build else "Release"

    # SuiteSparse
    suitesparse = CMakeExtension(
        name="SuiteSparse",
        install_prefix="amici",
        source_dir="amici/ThirdParty/SuiteSparse",
        cmake_build_type=build_type,
        cmake_configure_options=[
            *global_cmake_configure_options,
            "-DCMAKE_POSITION_INDEPENDENT_CODE=ON",
            "-DBUILD_SHARED_LIBS=OFF",
            "-DBUILD_TESTING=OFF",
            # Building SuiteSparse_config does not require a BLAS
            #  we just set BLAS_LIBRARIES to skip the search,
            #  the value is not used
            # "-DBLA_VENDOR=All",
            "-DBLAS_LIBRARIES=dummy",
            "-DSUITESPARSE_USE_64BIT_BLAS=ON",
            "-DSUITESPARSE_ENABLE_PROJECTS=amd;btf;colamd;klu",
            "-DSUITESPARSE_USE_CUDA=OFF",
            "-DSUITESPARSE_USE_FORTRAN=OFF",
            "-DSUITESPARSE_USE_PYTHON=OFF",
            "-DSUITESPARSE_USE_OPENMP=OFF",
            "-DSUITESPARSE_CONFIG_USE_OPENMP=OFF",
            "-DCHOLMOD_CAMD=OFF",
            "-DKLU_USE_CHOLMOD=OFF",
        ],
    )
    cmake_prefix_path = os.getenv("CMAKE_PREFIX_PATH", "").split(os.pathsep)
    # We need the potentially temporary and unpredictable build path
    #  to use artifacts from other extensions here. `${build_dir}` will
    #  be replaced by the actual path by `AmiciBuildCMakeExtension`
    #  before being passed to CMake.
    cmake_prefix_path.append("${build_dir}/amici")

    # SUNDIALS
    sundials = CMakeExtension(
        name="sundials",
        install_prefix="amici",
        source_dir="amici/ThirdParty/sundials",
        cmake_build_type=build_type,
        cmake_configure_options=[
            *global_cmake_configure_options,
            "-DBUILD_ARKODE=OFF",
            "-DBUILD_CVODE=OFF",
            "-DBUILD_IDA=OFF",
            "-DBUILD_KINSOL=OFF",
            "-DBUILD_SHARED_LIBS=OFF",
            "-DBUILD_STATIC_LIBS=ON",
            "-DBUILD_NVECTOR_MANYVECTOR=OFF",
            "-DEXAMPLES_ENABLE_C=OFF",
            "-DEXAMPLES_INSTALL=OFF",
            "-DENABLE_KLU=ON",
            f"-DCMAKE_PREFIX_PATH='{';'.join(cmake_prefix_path)}'",
        ],
    )

    # If scipy_openblas64 is installed, we make its cmake configuration
    # available
    try:
        import scipy_openblas64  # noqa: F401

        cmake_prefix_path.append("${build_dir}/amici/lib/cmake/openblas")
    except ImportError:
        pass

    # AMICI
    amici_ext = CMakeExtension(
        name="amici",
        install_prefix="amici",
        source_dir="amici",
        cmake_build_type=build_type,
        cmake_configure_options=[
            *global_cmake_configure_options,
            "-Werror=dev"
            # Turn warnings in to errors on GitHub Actions,
            #  match original repo and forks with default name
            if os.environ.get("GITHUB_REPOSITORY", "").endswith("/AMICI")
            else "-Wno-error=dev",
            "-DAMICI_PYTHON_BUILD_EXT_ONLY=ON",
            f"-DPython3_EXECUTABLE={Path(sys.executable).as_posix()}",
            f"-DCMAKE_PREFIX_PATH='{';'.join(cmake_prefix_path)}'",
        ],
    )

    # monkey patch a pre-cmake hook into the extension
    # see `amici.custom_commands.AmiciBuildCMakeExtension.build_extension`
    def amici_pre_cmake(ext, binary_dir):
        """Pre-cmake hook for AMICI extension"""
        # If scipy_openblas64 is installed, we make its cmake configuration
        # available
        try:
            import scipy_openblas64

            # we are not supposed to have a run-time dependency on the
            # scipy_openblas64 package, so we copy the relevant files
            # to the installation directory
            # (see warnings at https://github.com/MacPython/openblas-libs)
            binary_dir = Path(binary_dir)
            shutil.copytree(
                scipy_openblas64.get_lib_dir(),
                binary_dir / "amici" / "lib",
                dirs_exist_ok=True,
            )
            shutil.copytree(
                scipy_openblas64.get_include_dir(),
                binary_dir / "amici" / "include",
                dirs_exist_ok=True,
            )
        except ImportError:
            pass

    amici_ext.amici_pre_cmake = amici_pre_cmake

    # Order matters!
    return [suitesparse, sundials, amici_ext]


def main():
    ext_modules = get_extensions()

    # handle parallel building
    # Note: can be empty to use all hardware threads
    if (parallel_jobs := os.environ.get("AMICI_PARALLEL_COMPILE")) is not None:
        os.environ["CMAKE_BUILD_PARALLEL_LEVEL"] = parallel_jobs
    else:
        os.environ["CMAKE_BUILD_PARALLEL_LEVEL"] = "1"

    # Install
    setup(
        cmdclass={
            "install": AmiciInstall,
            "sdist": AmiciSDist,
            "build_ext": AmiciBuildCMakeExtension,
            "install_lib": AmiciInstallLib,
            "develop": AmiciDevelop,
            "build_py": AmiciBuildPy,
        },
        ext_modules=ext_modules,
    )


if __name__ == "__main__":
    main()
