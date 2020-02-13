"""Helper functions for AMICI core and module package preparation"""

import os
import sys
import shlex
import subprocess
import shutil

from .swig import find_swig, get_swig_version

try:
    import pkgconfig # optional
    # pkgconfig python module might be installed without pkg-config binary
    # being available
    pkgconfig.exists('somePackageName')
except (ModuleNotFoundError, EnvironmentError):
    pkgconfig = None


def getBlasConfig():
    """
    Find CBLAS-compatible BLAS

    Arguments:

    Returns:

    Raises:

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
        blaspkgcfg['define_macros'].append(('AMICI_BLAS_MKL', None),)
        return blaspkgcfg


    # Try pkgconfig
    if pkgconfig:
        if pkgconfig.exists('cblas'):
            blaspkgcfg = pkgconfig.parse('cblas')
            blaspkgcfg['extra_compile_args'] = [pkgconfig.cflags('cblas')]
            blaspkgcfg['extra_link_args'] = [pkgconfig.libs('cblas')]

            return blaspkgcfg

    # If none of the previous worked, fall back to libcblas in default paths
    blaspkgcfg['libraries'] = ['cblas']

    return blaspkgcfg


def getHdf5Config():
    """Find HDF5 include dir and libs

    Arguments:

    Returns:

    Raises:

    """
    if pkgconfig:
        h5pkgcfg = {}
        try:
            h5pkgcfg = pkgconfig.parse('hdf5')
        except pkgconfig.PackageNotFoundError:
            pass
        # NOTE: Cannot use pkgconfig.exists('hdf5f'), since this is true
        # although no libraries or include dirs are available
        h5pkgcfg['found'] = 'include_dirs' in h5pkgcfg \
                            and h5pkgcfg['include_dirs']
        if h5pkgcfg['found']:
            return h5pkgcfg

    h5pkgcfg = {'include_dirs': [],
                'library_dirs': [],
                'libraries': [],
                'define_macros': []
                }

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
            print('hdf5.h found in %s' % hdf5_include_dir_hint)
            h5pkgcfg['include_dirs'] = [hdf5_include_dir_hint]
            break

    for hdf5_library_dir_hint in hdf5_library_dir_hints:
        # check for static or shared library
        for lib_filename in ['libhdf5.a', 'libhdf5.so']:
            hdf5_library_dir_found = os.path.isfile(
                os.path.join(hdf5_library_dir_hint, lib_filename))
            if hdf5_library_dir_found:
                print(f'{lib_filename} found in {hdf5_library_dir_hint}')
                h5pkgcfg['library_dirs'] = [hdf5_library_dir_hint]
                break
        if hdf5_library_dir_found:
            # break to not override hdf5_library_dir_found
            break
    h5pkgcfg['found'] = hdf5_include_dir_found and hdf5_library_dir_found

    return h5pkgcfg


def addCoverageFlagsIfRequired(cxx_flags, linker_flags):
    """Add compiler and linker flags if gcov coverage requested

    Arguments:
    cxx_flags: list of existing cxx flags
    linker_flags: list of existing linker flags

    Returns:

    Raises:

    """
    if 'ENABLE_GCOV_COVERAGE' in os.environ and \
            os.environ['ENABLE_GCOV_COVERAGE'] == 'TRUE':
        print("ENABLE_GCOV_COVERAGE was set to TRUE."
              " Building AMICI with coverage symbols.")
        cxx_flags.extend(['-g', '-O0',  '--coverage'])
        linker_flags.extend(['--coverage','-g'])


def addDebugFlagsIfRequired(cxx_flags, linker_flags):
    """Add compiler and linker debug flags if requested

    Arguments:
    cxx_flags: list of existing cxx flags
    linker_flags: list of existing linker flags
    force: flag to force debug mode ignoring env settings

    Returns:

    Raises:

    """
    if 'ENABLE_AMICI_DEBUGGING' in os.environ \
            and os.environ['ENABLE_AMICI_DEBUGGING'] == 'TRUE':
        print("ENABLE_AMICI_DEBUGGING was set to TRUE."
              " Building AMICI with debug symbols.")
        cxx_flags.extend(['-g', '-O0'])
        linker_flags.extend(['-g'])


def generateSwigInterfaceFiles():
    """Compile the swig python interface to amici
    """
    swig_outdir = '%s/amici' % os.path.abspath(os.getcwd())
    swig_exe = find_swig()
    swig_version = get_swig_version(swig_exe)

    print(f"Found SWIG version {swig_version}")

    # Swig AMICI interface without HDF5 dependency
    swig_cmd = [swig_exe,
                '-c++',
                '-python',
                '-threads',
                '-Iamici/swig', '-Iamici/include',
                '-DAMICI_SWIG_WITHOUT_HDF5',
                '-outdir', swig_outdir,
                '-o', 'amici/amici_wrap_without_hdf5.cxx',
                'amici/swig/amici.i']

    # Do we have -doxygen?
    if swig_version >= (4, 0, 0):
        swig_cmd.insert(1, '-doxygen')

    print(f"Running SWIG: {' '.join(swig_cmd)}")
    sp = subprocess.run(swig_cmd, stdout=subprocess.PIPE,
                        stderr=sys.stdout.buffer)
    if not sp.returncode == 0:
        raise AssertionError('Swigging AMICI failed:\n'
                             + sp.stdout.decode('utf-8'))
    shutil.move(os.path.join(swig_outdir, 'amici.py'),
                os.path.join(swig_outdir, 'amici_without_hdf5.py'))

    # Swig AMICI interface with HDF5 dependency
    swig_cmd = [swig_exe,
                '-c++',
                '-python',
                '-threads',
                '-Iamici/swig', '-Iamici/include',
                '-outdir', swig_outdir,
                '-o', 'amici/amici_wrap.cxx',
                'amici/swig/amici.i']
    if swig_version >= (4, 0, 0):
        swig_cmd.insert(1, '-doxygen')
    print(f"Running SWIG: {' '.join(swig_cmd)}")
    sp = subprocess.run(swig_cmd, stdout=subprocess.PIPE,
                        stderr=sys.stdout.buffer)
    if not sp.returncode == 0:
        raise AssertionError('Swigging AMICI failed:\n'
                             + sp.stdout.decode('utf-8'))
