"""Package-level entrypoint"""

import sys

from . import (
    __version__,
    compiledWithOpenMP,
    has_clibs,
    hdf5_enabled,
    CpuTimer,
)


def print_info():
    """Displays information on the current AMICI installation.

    Useful for verifying package installation of submitting bug reports"""
    features = []

    if has_clibs:
        features.append("extensions")

    if compiledWithOpenMP():
        features.append("OpenMP")

    if hdf5_enabled:
        features.append("HDF5")

    if CpuTimer.uses_thread_clock:
        features.append("thread_clock")

    print(
        f"AMICI ({sys.platform}) version {__version__} ({','.join(features)})"
    )


if __name__ == "__main__":
    print_info()
