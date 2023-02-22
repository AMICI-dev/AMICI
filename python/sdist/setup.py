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

import cmake_build_extension
from setuptools import setup

# Add containing directory to path, as we need some modules from the AMICI
# package already for installation
sys.path.insert(0, os.path.dirname(__file__))

from amici.custom_commands import (
    AmiciInstall, AmiciDevelop,
    AmiciInstallLib, AmiciSDist, AmiciBuildPy,
    AmiciBuildCMakeExtension)


def get_extensions():
    """Get required extensions for build_ext"""
    # TODO no_clibs option

    amici_base_dir = Path('amici')
    suitesparse_base_dir = amici_base_dir / 'ThirdParty' / 'SuiteSparse'
    sundials_base_dir = amici_base_dir / 'ThirdParty' / 'sundials'

    # C++ extensions
    klu_inc_dir = Path(__file__).parent / 'amici' / 'ThirdParty' / 'SuiteSparse' / 'include'
    klu_inc_dir = klu_inc_dir.as_posix()
    sundials = cmake_build_extension.CMakeExtension(
        name='sundials',
        install_prefix='amici',
        source_dir='amici/ThirdParty/sundials',
        cmake_configure_options=[
            "-DBUILD_ARKODE=OFF",
            "-DBUILD_CVODE=OFF",
            "-DBUILD_IDA=OFF",
            "-DBUILD_KINSOL=OFF",
            "-DBUILD_SHARED_LIBS=OFF",
            "-DBUILD_STATIC_LIBS=ON",
            "-DBUILD_NVECTOR_MANYVECTOR=OFF",
            "-DBUILD_SUNNONLINSOL_PETSCSNES=OFF",
            "-DEXAMPLES_ENABLE_C=OFF",
            "-DEXAMPLES_INSTALL=OFF",
            "-DENABLE_KLU=ON",
            "-DKLU_LIBRARY_DIR='${build_dir}/amici/lib'",
            f"-DKLU_INCLUDE_DIR='{klu_inc_dir}'",
        ]
    )

    # TODO
    # 'macros': [
    #     ("gsl_CONFIG_CONTRACT_VIOLATION_THROWS", None),
    #     ("gsl_CONFIG_NARROW_THROWS_ON_TRUNCATION", 1),
    # ],
    amici_module = cmake_build_extension.CMakeExtension(
        name='_amici',
        install_prefix='amici',
        source_dir='amici',
        cmake_configure_options=[
            '-DAMICI_PYTHON_EXT_ONLY=ON',
            f'-DPython3_EXECUTABLE={Path(sys.executable).as_posix()}',
        ]
    )

    ss_config = cmake_build_extension.CMakeExtension(
        name='SuiteSparse_config',
        install_prefix='amici',
        source_dir='amici/ThirdParty/SuiteSparse/SuiteSparse_config',
        cmake_configure_options=[
            "-DBLA_VENDOR=All",
            "-DENABLE_CUDA=FALSE",
            "-DNFORTRAN=TRUE",
            # "--trace-expand",
            # "--debug-output",
            # "--debug-find",
        ]
    )
    install_dir = (Path(__file__).parent / "amici").absolute()
    prefix_path = install_dir #/ "lib" / "cmake" / "SuiteSparse"
    AmiciBuildCMakeExtension.extend_cmake_prefix_path(str(prefix_path))
    amd = cmake_build_extension.CMakeExtension(
        name='amd',
        install_prefix='amici',
        source_dir='amici/ThirdParty/SuiteSparse/AMD',
        cmake_configure_options=[
            "-DNFORTRAN=TRUE",
        ]
    )
    btf = cmake_build_extension.CMakeExtension(
        name='btf',
        install_prefix='amici',
        source_dir='amici/ThirdParty/SuiteSparse/BTF',
        cmake_configure_options=[
            "-DNFORTRAN=TRUE",
        ]
    )
    colamd = cmake_build_extension.CMakeExtension(
        name='colamd',
        install_prefix='amici',
        source_dir='amici/ThirdParty/SuiteSparse/COLAMD',
        cmake_configure_options=[
            "-DNFORTRAN=TRUE",
        ]
    )

    klu = cmake_build_extension.CMakeExtension(
        name='klu',
        install_prefix='amici',
        source_dir='amici/ThirdParty/SuiteSparse/KLU',
        cmake_configure_options=[
            "-DNCHOLMOD=ON",
            "-DENABLE_CUDA=FALSE",
            "-DNFORTRAN=TRUE",
        ]
    )
    return [ss_config, amd, btf, colamd, klu, sundials, amici_module]


def main():
    # Readme as long package description to go on PyPi
    # (https://pypi.org/project/amici/)
    with open(os.path.join(os.path.dirname(__file__), "README.md"),
              "r", encoding="utf-8") as fh:
        long_description = fh.read()

    ext_modules = get_extensions()

    # Install
    setup(
        cmdclass={
            'install': AmiciInstall,
            'sdist': AmiciSDist,
            'build_ext': AmiciBuildCMakeExtension,
            'install_lib': AmiciInstallLib,
            'develop': AmiciDevelop,
            'build_py': AmiciBuildPy,
        },
        long_description=long_description,
        long_description_content_type="text/markdown",
        ext_modules=ext_modules,
    )


if __name__ == '__main__':
    main()
