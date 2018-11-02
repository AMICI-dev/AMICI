"""Helper functions for AMICI core and module package preparation"""

import os
import platform
import shlex
import sys
import subprocess
import shutil

try:
    import pkgconfig # optional
    # pkgconfig python module might be installed without pkg-config binary being available
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

    if platform.system() in ['Linux', 'Darwin']:
        blaspkgcfg['libraries'] = ['cblas']

    if pkgconfig:
        if pkgconfig.exists('cblas'):
            blaspkgcfg = pkgconfig.parse('cblas')
            blaspkgcfg['extra_compile_args'] = [pkgconfig.cflags('cblas')]
            blaspkgcfg['extra_link_args'] = [pkgconfig.libs('cblas')]

    if 'BLAS_CFLAGS' in os.environ:
        blaspkgcfg['extra_compile_args'].extend(
            shlex.split(os.environ['BLAS_CFLAGS'])
        )

    if 'BLAS_LIBS' in os.environ:
        blaspkgcfg['extra_link_args'].extend(
            shlex.split(os.environ['BLAS_LIBS'])
        )

    return blaspkgcfg


def getHdf5Config():
    """Find HDF5 include dir and libs

    Arguments:

    Returns:

    Raises:

    """
    if pkgconfig:
        h5pkgcfg = pkgconfig.parse('hdf5')
        # NOTE: Cannot use pkgconfig.exists('hdf5f'), since this is true
        # althoughno libraries or include dirs are available
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
        '/usr/include',  # travis ubuntu xenial
        '/usr/local/Cellar/hdf5/1.10.2_1/include'  # travis macOS
    ]
    hdf5_library_dir_hints = [
        '/usr/lib/x86_64-linux-gnu/',  # travis ubuntu xenial
        '/usr/lib/x86_64-linux-gnu/hdf5/serial',
        '/usr/local/lib',
        '/usr/local/Cellar/hdf5/1.10.2_1/lib'  # travis macOS
    ]

    for hdf5_include_dir_hint in hdf5_include_dir_hints:
        hdf5_include_dir_found = os.path.isfile(
            os.path.join(hdf5_include_dir_hint, 'hdf5.h'))
        if hdf5_include_dir_found:
            print('hdf5.h found in %s' % hdf5_include_dir_hint)
            h5pkgcfg['include_dirs'] = [hdf5_include_dir_hint]
            break

    for hdf5_library_dir_hint in hdf5_library_dir_hints:
        hdf5_library_dir_found = os.path.isfile(
            os.path.join(hdf5_library_dir_hint, 'libhdf5.a'))
        if hdf5_library_dir_found:
            print('libhdf5.a found in %s' % hdf5_library_dir_hint)
            h5pkgcfg['library_dirs'] = [hdf5_library_dir_hint]
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
    swig_cmd = findSwig()
    sp = subprocess.run([swig_cmd,
                         '-c++',
                         '-python',
                         '-Iamici/swig', '-Iamici/include',
                         '-DAMICI_SWIG_WITHOUT_HDF5',
                         '-outdir', swig_outdir,
                         '-o', 'amici/amici_wrap_without_hdf5.cxx',
                         'amici/swig/amici.i'])
    assert (sp.returncode == 0)
    shutil.move(os.path.join(swig_outdir, 'amici.py'),
                os.path.join(swig_outdir, 'amici_without_hdf5.py'))
    sp = subprocess.run([swig_cmd,
                         '-c++',
                         '-python',
                         '-Iamici/swig', '-Iamici/include',
                         '-outdir', swig_outdir,
                         '-o', 'amici/amici_wrap.cxx',
                         'amici/swig/amici.i'])
    assert (sp.returncode == 0)


def findSwig():
    """Get name of SWIG executable

    We need version 3.0.
    Probably we should try some default paths and names, but this should do the trick for now.
    Debian/Ubuntu systems have swig3.0 ('swig' is older versions), OSX has swig 3.0 as 'swig'."""
    if sys.platform != 'linux':
        return 'swig'
    return 'swig3.0'