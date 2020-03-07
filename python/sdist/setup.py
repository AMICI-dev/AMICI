"""Setuptools file for creating AMICI module

This file is based on setuptools alone and does not require CMake.
All sources are compiled anew.

This file expects to be run from within its directory.

Requires:
- swig3.0
- setuptools
- pkgconfig python+executables
- hdf5 libraries and headers
"""

from setuptools import find_packages, setup, Extension

import os
import sys
import subprocess
import setup_clibs  # Must run from within containing directory

# Add current directory to path, as we need some modules from the AMICI
# package already for installation
sys.path.insert(0, os.getcwd())

from amici import __version__
from amici.custom_commands import (
    my_install, my_build_clib, my_develop,
    my_install_lib, my_build_ext, my_sdist)
from amici.setuptools import (
    get_blas_config,
    get_hdf5_config,
    add_coverage_flags_if_required,
    add_debug_flags_if_required,
)


def try_install(package):
    """Try installing the given package using pip. Exit on error."""
    errno = subprocess.call([sys.executable, "-m", "pip", "install", package])
    if errno:
        print(f"Failed trying to install {package}. Please install manually.")
        raise SystemExit(errno)


try:
    # required for include directory
    import numpy as np
except ImportError:
    # We need numpy, but setup_requires fires too late
    try_install('numpy')
    # retry
    import numpy as np

# Python version check. We need >= 3.6 due to e.g. f-strings
if sys.version_info < (3, 6):
    sys.exit('amici requires at least Python version 3.6')


def main():
    # Extra compiler flags
    cxx_flags = []
    amici_module_linker_flags = []
    define_macros = []

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
        extension_sources = [
            'amici/amici_wrap.cxx',  # swig interface
        ]
        define_macros.extend(h5pkgcfg['define_macros'])
    else:
        print("HDF5 library NOT found. Building AMICI WITHOUT HDF5 support.")
        extension_sources = [
            'amici/amici_wrap_without_hdf5.cxx',  # swig interface
        ]

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
    with open("README.md", "r", encoding="utf-8") as fh:
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
        name='amici',
        cmdclass={
            'install': my_install,
            'sdist': my_sdist,
            'build_ext': my_build_ext,
            'build_clib': my_build_clib,
            'install_lib': my_install_lib,
            'develop': my_develop,
        },
        version=__version__,
        description='Advanced multi-language Interface to CVODES and IDAS',
        long_description=long_description,
        long_description_content_type="text/markdown",
        url='https://github.com/ICB-DCM/AMICI',
        author='Fabian Froehlich, Jan Hasenauer, Daniel Weindl and '
               'Paul Stapor',
        author_email='fabian_froehlich@hms.harvard.edu',
        license='BSD',
        libraries=[libamici, libsundials, libsuitesparse],
        ext_modules=[amici_module],
        py_modules=['amici/amici',  # the swig interface
                    'amici/amici_without_hdf5',  # the swig interface
                    ],
        packages=find_packages(),
        package_dir={'amici': 'amici'},
        entry_points={
            'console_scripts': [
                'amici_import_petab = amici.petab_import:main',
                # for backwards compatibility
                'amici_import_petab.py = amici.petab_import:main'
            ]
        },
        install_requires=['sympy',
                          'python-libsbml',
                          'h5py',
                          'pandas',
                          'pkgconfig',
                          'wurlitzer'],
        setup_requires=['setuptools>=40.6.3'],
        python_requires='>=3.6',
        extras_require={
            'petab': ['petab>=0.1.4']
        },
        package_data={
            'amici': ['amici/include/amici/*',
                      'src/*template*',
                      'swig/*',
                      'libs/*',
                      'amici.py',
                      'amici_without_hdf5.py',
                      'setup.py.template',
                      ],
        },
        zip_safe=False,
        include_package_data=True,
        exclude_package_data={
            '': ['README.txt'],
        },
        test_suite="tests",
        classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Operating System :: POSIX :: Linux',
            'Operating System :: MacOS :: MacOS X',
            'Programming Language :: Python',
            'Programming Language :: C++',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    )


if __name__ == '__main__':
    main()
