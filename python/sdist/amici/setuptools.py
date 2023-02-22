"""
setuptools
----------
Helper functions for AMICI core and module package preparation
"""

import os
import sys

try:
    import pkgconfig  # optional

    # pkgconfig python module might be installed without pkg-config binary
    # being available
    pkgconfig.exists('somePackageName')
except (ModuleNotFoundError, EnvironmentError):
    pkgconfig = None

from typing import Dict, List, Union, Tuple, Any

PackageInfo = Dict[str, List[Union[str, Tuple[str, Any]]]]


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
        print("ENABLE_GCOV_COVERAGE was set to TRUE."
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
        print("ENABLE_AMICI_DEBUGGING was set to TRUE."
              " Building AMICI with debug symbols.")
        cxx_flags.extend(['-g', '-O0', '-UNDEBUG'])
        if sys.platform != "win32":
            # these options are incompatible with MSVC, but there is no easy
            # way to detect which compiler will be used. so we just skip them
            # altogether on windows.
            cxx_flags.extend(['-Werror', '-Wno-error=deprecated-declarations'])

        linker_flags.extend(['-g'])

