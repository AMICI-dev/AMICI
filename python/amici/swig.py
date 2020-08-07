"""Functions for downloading/building/finding SWIG"""

from typing import Tuple
import os
import subprocess
import re


def find_swig() -> str:
    """Get name and version of SWIG executable

    We need version >=3.0. Probably we should try some default paths and names,
    but this should do the trick for now.

    Debian/Ubuntu systems have swig3.0 ('swig' is older versions),
    OSX has swig 3.0 as 'swig'.
    """

    candidates = ['swig4.0', 'swig3.0', 'swig']
    # Environment variable has priority
    if 'SWIG' in os.environ:
        candidates.insert(0, os.environ['SWIG'])

    for candidate in candidates:
        if swig_works(candidate):
            return candidate

    raise RuntimeError(
        "Unable to find SWIG executable with default names. "
        "Ensure you have SWIG installed, e.g. by "
        "`sudo apt install swig` or `brew install swig`. "
        "As non-root user, you can install SWIG using "
        "https://github.com/AMICI-dev/AMICI/blob/master/scripts/"
        "downloadAndBuildSwig.sh, or by following the "
        "instructions at http://www.swig.org/Doc4.0/"
        "SWIGDocumentation.html#Preface_installation. "
        "If was not found despite being installed, set the SWIG"
        " environment variable to the full path of the correct "
        "executable."
    )


def swig_works(swig: str, verbose: bool = True) -> bool:
    """Test if `swig` looks like a working SWIG executable."""

    try:
        # For python3.6 compatibility no `capture_output=True`
        result = subprocess.run([swig, '-version'],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    except (FileNotFoundError, PermissionError):
        if verbose:
            print(f'Testing SWIG executable {swig}... FAILED.')
        return False

    if verbose:
        if result.returncode == 0:
            print(f'Testing SWIG executable {swig}... SUCCEEDED.')
        else:
            print(f'Testing SWIG executable {swig}... FAILED.')

    return result.returncode == 0


def get_swig_version(swig_exe: str) -> Tuple:
    """Determine version of the given SWIG executable

    Returns:
        Version tuple
    """
    result = subprocess.run([swig_exe, '-version'],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    assert result.returncode == 0

    version = re.sub(r'(?s).*Version\s+([\S]+).*', r'\1',
                     result.stdout.decode('utf-8'))

    return tuple(int(x) for x in version.split('.'))
