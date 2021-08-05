#!/bin/bash
# Test python model wrapping inside virtual environment

# set up some useful path variables
script_path=$(dirname $BASH_SOURCE)
amici_path=$(cd "$script_path"/.. && pwd)
set -e

# check for presence of BioNetGen package
if [[ -z "${BNGPATH}" ]]; then
    export BNGPATH=${amici_path}/ThirdParty/BioNetGen-2.5.2
fi

# required packages for testing
PACKAGE_LIST=(scipy h5py pytest pytest-cov)

# satisfied prerequisites for testing
cd "${amici_path}"/python/tests
if [ -f "${amici_path}/build/venv/bin/activate" ]; then 
   source "${amici_path}/build/venv/bin/activate"
   pip install ${PACKAGE_LIST[*]}
else
   pip install --user ${PACKAGE_LIST[*]}
fi

# actually run the tests but exclude PEtab tests (run separetely)
pytest --ignore-glob=*petab*
