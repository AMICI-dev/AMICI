"""Helper functions for AMICI core and module package preparation"""

import os
import platform

try:
    import pkgconfig # optional
except ModuleNotFoundError:
    pkgconfig = None


def getBlasConfig():
    """Find cblas

    Arguments:

    Returns:

    Raises:

    """

    # Find cblas
    blaspkgcfg = {'include_dirs': [],
                  'library_dirs': [],
                  'libraries': [],
                  'define_macros': []
                  }
    
    if platform.system() == 'Linux':
        blaspkgcfg['libraries'] = ['cblas']
    
    if 'BLAS_INCDIR' in os.environ:
        blaspkgcfg['include_dirs'].extend(os.environ['BLAS_INCDIR'].split(' '))
    
    if 'BLAS_LIB' in os.environ:
        blaspkgcfg['libraries'].extend(os.environ['BLAS_LIB'].split(' '))

    return blaspkgcfg


def getHdf5Config():
    """Find HDF5 include dir and libs
    
    Arguments:

    Returns:

    Raises:

    """
    if pkgconfig:
        h5pkgcfg = pkgconfig.parse('hdf5')
        # NOTE: Cannot use pkgconfig.exists('hdf5f'), since this is true although
        # no libraries or include dirs are available
        h5pkgcfg['found'] = 'include_dirs' in h5pkgcfg and h5pkgcfg['include_dirs']
        if h5pkgcfg['found']:
            return h5pkgcfg
        
    h5pkgcfg = {'include_dirs': [],
                'library_dirs': [],
                'libraries': [],
                'define_macros': []
                }
    
    # try for hdf5 in standard locations
    hdf5_include_dir_hints = ['/usr/include/hdf5/serial',
                              '/usr/local/include',
                              '/usr/include', # travis ubuntu xenial
                              '/usr/local/Cellar/hdf5/1.10.2_1/include'] # travis macOS
    hdf5_library_dir_hints = ['/usr/lib/x86_64-linux-gnu/', # travis ubuntu xenial
                              '/usr/lib/x86_64-linux-gnu/hdf5/serial',
                              '/usr/local/lib',
                              '/usr/local/Cellar/hdf5/1.10.2_1/lib'] # travis macOS

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
    cxx_flags: list
    linker_flags: list

    Returns:

    Raises:

    """
    if 'ENABLE_GCOV_COVERAGE' in os.environ and os.environ['ENABLE_GCOV_COVERAGE'] == 'TRUE':
        print("ENABLE_GCOV_COVERAGE was set to TRUE. Building AMICI with coverage symbols.")
        cxx_flags.extend(['-g', '-O0',  '--coverage'])
        linker_flags.extend(['--coverage','-g'])


def addDebugFlagsIfRequired(cxx_flags, linker_flags):
    """Add compiler and linker debug flags if requested
    
    Arguments:
    cxx_flags: list
    linker_flags: list

    Returns:

    Raises:

    """
    if 'ENABLE_AMICI_DEBUGGING' in os.environ and os.environ['ENABLE_AMICI_DEBUGGING'] == 'TRUE':
        print("ENABLE_AMICI_DEBUGGING was set to TRUE. Building AMICI with debug symbols.")
        cxx_flags.extend(['-g', '-O0'])
        linker_flags.extend(['-g'])

    