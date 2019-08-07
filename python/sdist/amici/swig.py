"""Functions for downloading/building/finding SWIG"""

import sys
import os
import subprocess


def find_swig():
    """Get name of SWIG executable

    We need version >=3.0. Probably we should try some default paths and names,
    but this should do the trick for now.

    Debian/Ubuntu systems have swig3.0 ('swig' is older versions),
    OSX has swig 3.0 as 'swig'."""

    candidates = ['swig4.0', 'swig3.0', 'swig']
    # Environment variable has priority
    if 'SWIG' in os.environ:
        candidates.insert(0, os.environ['SWIG'])

    for candidate in candidates:
        if swig_works(candidate):
            return candidate

    raise RuntimeError("Unable to find SWIG executable with default names. "
                       "Ensure you have SWIG installed, e.g. by "
                       "`sudo apt install swig3.0` or `brew install swig`. "
                       "As non-root user, you can install SWIG using "
                       "https://github.com/ICB-DCM/AMICI/blob/master/scripts/"
                       "downloadAndBuildSwig.sh, or by following the "
                       "instructions at http://www.swig.org/Doc4.0/"
                       "SWIGDocumentation.html#Preface_installation. "
                       "If was not found despite being installed, set the SWIG"
                       " environment variable to the full path of the correct "
                       "executable.")


def swig_works(swig, verbose = True):
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
            print(result.stdout.decode('utf-8'))
        else:
            print(f'Testing SWIG executable {swig}... FAILED.')

    return result.returncode == 0
