"""AMICI model package setup"""

import os
import sys
from pathlib import Path

from amici import _get_amici_path
from amici.custom_commands import AmiciBuildCMakeExtension
from cmake_build_extension import CMakeExtension
from setuptools import find_namespace_packages, setup
import importlib.metadata


def get_extension() -> CMakeExtension:
    """Get setuptools extension object for this AMICI model package"""

    # Build shared object
    prefix_path = Path(_get_amici_path())
    AmiciBuildCMakeExtension.extend_cmake_prefix_path(str(prefix_path))

    # handle parallel building
    # Note: can be empty to use all hardware threads
    if (parallel_jobs := os.environ.get("AMICI_PARALLEL_COMPILE")) is not None:
        os.environ["CMAKE_BUILD_PARALLEL_LEVEL"] = parallel_jobs
    else:
        os.environ["CMAKE_BUILD_PARALLEL_LEVEL"] = "1"

    debug_build = os.getenv("ENABLE_AMICI_DEBUGGING", "").lower() in [
        "1",
        "true",
    ] or os.getenv("ENABLE_GCOV_COVERAGE", "").lower() in ["1", "true"]

    cmake_prefix_path = os.getenv("CMAKE_PREFIX_PATH", "").split(os.pathsep)
    cmake_prefix_path.append(prefix_path.as_posix())

    # If scipy_openblas64 is installed, we make its cmake configuration
    # available
    amici_distribution = importlib.metadata.distribution("amici")
    amici_dir = Path(amici_distribution.locate_file(""))
    # this path is created during the amici build if scipy_openblas64 is used
    openblas_cmake_dir = amici_dir / "lib" / "cmake" / "openblas"
    if openblas_cmake_dir.exists():
        cmake_prefix_path.append(str(openblas_cmake_dir))

    return CMakeExtension(
        name="model_ext",
        source_dir=os.getcwd(),
        install_prefix="model_events_py",
        cmake_configure_options=[
            "-DCMAKE_VERBOSE_MAKEFILE=ON",
            f"-DCMAKE_PREFIX_PATH='{';'.join(cmake_prefix_path)}'",
            "-DAMICI_PYTHON_BUILD_EXT_ONLY=ON",
            f"-DPython3_EXECUTABLE={Path(sys.executable).as_posix()}",
        ],
        cmake_build_type="Debug" if debug_build else "Release",
    )


# Change working directory to setup.py location
os.chdir(os.path.dirname(os.path.abspath(__file__)))

MODEL_EXT = get_extension()

CLASSIFIERS = [
    "Development Status :: 3 - Alpha",
    "Intended Audience :: Science/Research",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python",
    "Programming Language :: C++",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

CMDCLASS = {
    # for CMake-based builds
    "build_ext": AmiciBuildCMakeExtension,
}

# Install
setup(
    name="model_events_py",
    cmdclass=CMDCLASS,
    version="0.1.0",
    description="AMICI-generated module for model model_events_py",
    url="https://github.com/AMICI-dev/AMICI",
    author="model-author-todo",
    author_email="model-author-todo",
    ext_modules=[MODEL_EXT],
    packages=find_namespace_packages(),
    install_requires=["amici==0.33.0"],
    python_requires=">=3.11",
    package_data={},
    zip_safe=False,
    classifiers=CLASSIFIERS,
)
