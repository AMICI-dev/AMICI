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

# Add containing directory to path, as we need some modules from the AMICI
# package already for installation
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
import setup_clibs  # Must run from within containing directory
from setuptools import find_packages, setup, Extension

from amici.custom_commands import (
    AmiciInstall, AmiciBuildCLib, AmiciDevelop,
    AmiciInstallLib, AmiciBuildExt, AmiciSDist)
from amici.setuptools import (
    get_blas_config,
    get_hdf5_config,
    add_coverage_flags_if_required,
    add_debug_flags_if_required,
    add_openmp_flags,
)

def main():
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

    extension_sources = [
        'amici/amici_wrap.cxx',  # swig interface
    ]

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

    libamici = setup_clibs.get_lib_amici(
        h5pkgcfg=h5pkgcfg, blaspkgcfg=blaspkgcfg,
        extra_compiler_flags=cxx_flags)
    libsundials = setup_clibs.get_lib_sundials(extra_compiler_flags=cxx_flags)
    libsuitesparse = setup_clibs.get_lib_suite_sparse(
        extra_compiler_flags=cxx_flags + ['-DDLONG']
    )

    # Readme as long package description to go on PyPi
    # (https://pypi.org/project/amici/)
    with open(os.path.join(os.path.dirname(__file__), "README.md"),
              "r", encoding="utf-8") as fh:
        long_description = fh.read()

    # Build shared object
    amici_module = Extension(
        name='amici._amici',
        sources=extension_sources,
        include_dirs=['amici/include',
                      'amici/ThirdParty/gsl/',
                      *libsundials[1]['include_dirs'],
                      *libsuitesparse[1]['include_dirs'],
                      *h5pkgcfg['include_dirs'],
                      *blaspkgcfg['include_dirs'],
                      np.get_include()
                      ],
        # Cannot use here, see above
        # libraries=[
        #    'hdf5_hl_cpp', 'hdf5_hl', 'hdf5_cpp', 'hdf5'
        # ],
        define_macros=define_macros,
        library_dirs=[
            *h5pkgcfg['library_dirs'],
            *blaspkgcfg['library_dirs'],
            'amici/libs',  # clib target directory
        ],
        extra_compile_args=cxx_flags,
        extra_link_args=amici_module_linker_flags
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
            'build_ext': AmiciBuildExt,
            'build_clib': AmiciBuildCLib,
            'install_lib': AmiciInstallLib,
            'develop': AmiciDevelop,
        },
        long_description=long_description,
        long_description_content_type="text/markdown",
        libraries=[libamici, libsundials, libsuitesparse],
        ext_modules=[amici_module],
    )


if __name__ == '__main__':
    main()
