"""Setuptools file for creating AMICI module

This file is based on setuptools alone and does not require CMake.
All sources are compiled anew.

This file expects to be run from within its directory.

Non-python-package requirements:
- swig3.0
- Optional: hdf5 libraries and headers
"""

import os
import sys
from pathlib import Path

import cmake_build_extension

# Add containing directory to path, as we need some modules from the AMICI
# package already for installation
sys.path.insert(0, os.path.dirname(__file__))

from setuptools import setup

from amici.custom_commands import (
    AmiciInstall, AmiciDevelop,
    AmiciInstallLib, AmiciSDist, AmiciBuildPy,
    AmiciBuildCMakeExtension)
from amici.setuptools import (
    get_blas_config,
    get_hdf5_config,
    add_coverage_flags_if_required,
    add_debug_flags_if_required,
    add_openmp_flags,
)


def main():
    # TODO no_clibs option
    # TODO to CMake
    # Extra compiler flags
    cxx_flags = []
    amici_module_linker_flags = []
    define_macros = []

    add_openmp_flags(cxx_flags=cxx_flags, ldflags=amici_module_linker_flags)

    blaspkgcfg = get_blas_config()
    amici_module_linker_flags.extend(blaspkgcfg['extra_link_args'])
    amici_module_linker_flags.extend(
        f'-l{lib}' for lib in blaspkgcfg['libraries'])
    define_macros.extend(blaspkgcfg['define_macros'])

    h5pkgcfg = get_hdf5_config()

    if h5pkgcfg['found']:
        # Manually add linker flags. The libraries passed to Extension will
        # end up in front of the clibs in the linker line and not after, where
        # they are required.
        print("HDF5 library found. Building AMICI with HDF5 support.")
        amici_module_linker_flags.extend(
            [f'-l{lib}' for lib in
             ['hdf5_hl_cpp', 'hdf5_hl', 'hdf5_cpp', 'hdf5']])
        define_macros.extend(h5pkgcfg['define_macros'])
    else:
        print("HDF5 library NOT found. Building AMICI WITHOUT HDF5 support.")
        define_macros.append(('AMICI_SWIG_WITHOUT_HDF5', None))

    add_coverage_flags_if_required(
        cxx_flags,
        amici_module_linker_flags,
    )

    add_debug_flags_if_required(
        cxx_flags,
        amici_module_linker_flags,
    )

    # compiler and linker flags for libamici
    if 'AMICI_CXXFLAGS' in os.environ:
        cxx_flags.extend(os.environ['AMICI_CXXFLAGS'].split(' '))
    if 'AMICI_LDFLAGS' in os.environ:
        amici_module_linker_flags.extend(
            os.environ['AMICI_LDFLAGS'].split(' '))

    amici_base_dir = Path('amici')
    suitesparse_base_dir = amici_base_dir / 'ThirdParty' / 'SuiteSparse'
    sundials_base_dir = amici_base_dir / 'ThirdParty' / 'sundials'

    # Readme as long package description to go on PyPi
    # (https://pypi.org/project/amici/)
    with open(os.path.join(os.path.dirname(__file__), "README.md"),
              "r", encoding="utf-8") as fh:
        long_description = fh.read()

    # C++ extensions
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
            f"-DKLU_INCLUDE_DIR='{Path(__file__).parent / 'amici' / 'ThirdParty' / 'SuiteSparse' }/include'",
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
        ]
    )

    # TODO _DNSTATIC
    ss_config = cmake_build_extension.CMakeExtension(
        name='SuiteSparse_config',
        install_prefix='amici',
        source_dir='amici/ThirdParty/SuiteSparse/SuiteSparse_config',
        cmake_configure_options=[
            "-DALLOW_64BIT_BLAS=ON",
            "-DBLA_VENDOR=All",
            "-DENABLE_CUDA=FALSE",
            "-DNFORTRAN=TRUE",
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

    # Monkey-patch extension (see
    # `custom_commands.set_compiler_specific_extension_options`)
    amici_module.extra_compile_args_mingw32 = ['-std=c++14']
    amici_module.extra_compile_args_unix = ['-std=c++14']
    amici_module.extra_compile_args_msvc = ['/std:c++14']

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
        ext_modules=[ss_config, amd, btf, colamd, klu, sundials, amici_module],
    )


if __name__ == '__main__':
    main()
