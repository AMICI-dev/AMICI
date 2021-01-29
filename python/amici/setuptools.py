"""
setuptools
----------
Helper functions for AMICI core and module package preparation
"""

import os
import sys
import shlex
import subprocess

from distutils import log
from .swig import find_swig, get_swig_version

try:
    import pkgconfig  # optional

    # pkgconfig python module might be installed without pkg-config binary
    # being available
    pkgconfig.exists('somePackageName')
except (ModuleNotFoundError, EnvironmentError):
    pkgconfig = None

from typing import Dict, List, Union, Tuple, Any

PackageInfo = Dict[str, List[Union[str, Tuple[str, Any]]]]


def get_blas_config() -> PackageInfo:
    """
    Find CBLAS-compatible BLAS

    :return:
        blas related package information
    """

    blaspkgcfg = {'include_dirs': [],
                  'library_dirs': [],
                  'libraries': [],
                  'define_macros': [],
                  'extra_compile_args': [],
                  'extra_link_args': []
                  }

    # Check environment variables
    if 'BLAS_CFLAGS' in os.environ:
        blaspkgcfg['extra_compile_args'].extend(
            shlex.split(os.environ['BLAS_CFLAGS'])
        )

    if 'BLAS_LIBS' in os.environ:
        blaspkgcfg['extra_link_args'].extend(
            shlex.split(os.environ['BLAS_LIBS'])
        )

    if 'BLAS_CFLAGS' in os.environ or 'BLAS_LIBS' in os.environ:
        # If options have been provided by the user, we don't try to detect
        # anything by ourselves
        return blaspkgcfg

    # Try environment modules
    # MKL
    if 'MKLROOT' in os.environ:
        if 'MKL_INC' in os.environ:
            blaspkgcfg['extra_compile_args'].extend(
                shlex.split(os.environ['MKL_INC'])
            )
        if 'MKL_LIB' in os.environ:
            blaspkgcfg['extra_link_args'].extend(
                shlex.split(os.environ['MKL_LIB'])
            )
        blaspkgcfg['define_macros'].append(('AMICI_BLAS_MKL', None), )
        return blaspkgcfg

    # Try pkgconfig
    if pkgconfig:
        for blas_name in ['cblas', 'openblas']:
            if pkgconfig.exists(blas_name):
                blaspkgcfg = pkgconfig.parse(blas_name)
                blaspkgcfg['extra_compile_args'] = [
                    pkgconfig.cflags(blas_name)
                ]
                blaspkgcfg['extra_link_args'] = [
                    pkgconfig.libs(blas_name)
                ]

                return blaspkgcfg

    # If none of the previous worked, fall back to libcblas in default paths
    blaspkgcfg['libraries'] = ['cblas']

    return blaspkgcfg


def get_hdf5_config() -> PackageInfo:
    """
    Find HDF5 include dir and libs

    :return:
        hdf5 related package information
    """

    h5pkgcfg = {'include_dirs': [],
                'library_dirs': [],
                'libraries': [],
                'define_macros': []
                }
    hdf5_include_dir_found = False
    hdf5_library_dir_found = False

    # try for hdf5 in standard locations
    hdf5_include_dir_hints = [
        '/usr/include/hdf5/serial',
        '/usr/local/include',
        '/usr/include',  # travis ubuntu xenial, centos
        '/usr/local/Cellar/hdf5/1.10.2_1/include'  # travis macOS
    ]
    hdf5_library_dir_hints = [
        '/usr/lib/x86_64-linux-gnu/',  # travis ubuntu xenial
        '/usr/lib/x86_64-linux-gnu/hdf5/serial',
        '/usr/local/lib',
        '/usr/lib64/',  # CentOS
        '/usr/local/Cellar/hdf5/1.10.2_1/lib'  # travis macOS
    ]

    # special treatment for conda environments
    # as the conda library dir is provided first, we should also check for
    # conda header files first
    if 'CONDA_DIR' in os.environ:
        hdf5_include_dir_hints.insert(
            0, os.path.join(os.environ['CONDA_DIR'], 'include'))
        hdf5_library_dir_hints.insert(
            0, os.path.join(os.environ['CONDA_DIR'], 'lib'))

    # Check for Environment Modules variables
    if 'HDF5_BASE' in os.environ:
        hdf5_include_dir_hints.insert(
            0, os.path.join(os.environ['HDF5_BASE'], 'include'))
        hdf5_library_dir_hints.insert(
            0, os.path.join(os.environ['HDF5_BASE'], 'lib'))

    for hdf5_include_dir_hint in hdf5_include_dir_hints:
        hdf5_include_dir_found = os.path.isfile(
            os.path.join(hdf5_include_dir_hint, 'hdf5.h'))
        if hdf5_include_dir_found:
            log.info('hdf5.h found in %s' % hdf5_include_dir_hint)
            h5pkgcfg['include_dirs'] = [hdf5_include_dir_hint]
            break

    for hdf5_library_dir_hint in hdf5_library_dir_hints:
        # check for static or shared library
        for lib_filename in ['libhdf5.a', 'libhdf5.so']:
            hdf5_library_dir_found = os.path.isfile(
                os.path.join(hdf5_library_dir_hint, lib_filename))
            if hdf5_library_dir_found:
                log.info(f'{lib_filename} found in {hdf5_library_dir_hint}')
                h5pkgcfg['library_dirs'] = [hdf5_library_dir_hint]
                break
        if hdf5_library_dir_found:
            # break to not override hdf5_library_dir_found
            break

    h5pkgcfg['found'] = hdf5_include_dir_found and hdf5_library_dir_found
    if h5pkgcfg['found']:
        return h5pkgcfg

    if pkgconfig:
        try:
            h5pkgcfg = pkgconfig.parse('hdf5')
        except pkgconfig.PackageNotFoundError:
            pass
        # NOTE: Cannot use pkgconfig.exists('hdf5f'), since this is true
        # although no libraries or include dirs are available
        h5pkgcfg['found'] = 'include_dirs' in h5pkgcfg \
                            and h5pkgcfg['include_dirs'] and \
                            'library_dirs' in h5pkgcfg \
                            and h5pkgcfg['library_dirs']

    return h5pkgcfg


def add_coverage_flags_if_required(cxx_flags: List[str],
                                   linker_flags: List[str]) -> None:
    """
    Add compiler and linker flags if gcov coverage requested

    :param cxx_flags:
        list of existing cxx flags

    :param linker_flags:
        list of existing linker flags
    """
    if 'ENABLE_GCOV_COVERAGE' in os.environ and \
            os.environ['ENABLE_GCOV_COVERAGE'].upper() == 'TRUE':
        log.info("ENABLE_GCOV_COVERAGE was set to TRUE."
                 " Building AMICI with coverage symbols.")
        cxx_flags.extend(['-g', '-O0', '--coverage'])
        linker_flags.extend(['--coverage', '-g'])


def add_debug_flags_if_required(cxx_flags: List[str],
                                linker_flags: List[str]) -> None:
    """
    Add compiler and linker debug flags if requested

    Arguments:
    :param cxx_flags:
        list of existing cxx flags

    :param linker_flags:
        list of existing linker flags
    """
    if 'ENABLE_AMICI_DEBUGGING' in os.environ \
            and os.environ['ENABLE_AMICI_DEBUGGING'] == 'TRUE':
        log.info("ENABLE_AMICI_DEBUGGING was set to TRUE."
                 " Building AMICI with debug symbols.")
        cxx_flags.extend(['-g', '-O0', '-UNDEBUG'])
        linker_flags.extend(['-g'])


def generate_swig_interface_files(swig_outdir: str = None,
                                  with_hdf5: bool = None) -> None:
    """
    Compile the swig python interface to amici
    """

    swig_exe = find_swig()
    swig_version = get_swig_version(swig_exe)

    swig_args = [
        '-c++',
        '-python',
        '-py3',
        '-threads',
        '-Wall',
        f'-Iamici{os.sep}swig',
        f'-Iamici{os.sep}include',
    ]

    log.info(f"Found SWIG version {swig_version}")

    # Are HDF5 includes available to generate the wrapper?
    if with_hdf5 is None:
        with_hdf5 = get_hdf5_config()['found']

    if not with_hdf5:
        swig_args.append('-DAMICI_SWIG_WITHOUT_HDF5')

    if swig_outdir is not None:
        swig_args.extend(['-outdir', swig_outdir])

    # Do we have -doxygen?
    if swig_version >= (4, 0, 0):
        swig_args.append('-doxygen')

    swig_cmd = [swig_exe,
                *swig_args,
                '-o', os.path.join("amici", "amici_wrap.cxx"),
                os.path.join("amici", "swig", "amici.i")]

    log.info(f"Running SWIG: {' '.join(swig_cmd)}")
    sp = subprocess.run(swig_cmd, stdout=subprocess.PIPE,
                        stderr=sys.stdout.buffer)
    if not sp.returncode == 0:
        raise AssertionError('Swigging AMICI failed:\n'
                             + sp.stdout.decode('utf-8'))


def add_openmp_flags(cxx_flags: List, ldflags: List) -> None:
    """Add OpenMP flags to lists for compiler/linker flags (in-place)"""

    # Enable OpenMP support for Linux / OSX:
    if sys.platform == 'linux':
        log.info("Adding OpenMP flags...")
        cxx_flags.insert(0, "-fopenmp")
        ldflags.insert(0, "-fopenmp")
    elif sys.platform == 'darwin':
        if os.path.exists('/usr/local/lib/libomp.a'):
            log.info("Adding OpenMP flags...")
            cxx_flags[0:0] = ["-Xpreprocessor", "-fopenmp"]
            ldflags[0:0] = ["-Xpreprocessor", "-fopenmp", "-lomp"]
        else:
            log.info("Not adding OpenMP flags, because /usr/local/lib/libomp.a"
                     " does not exist. To enable, run `brew install libomp` "
                     "or add flags manually.")
