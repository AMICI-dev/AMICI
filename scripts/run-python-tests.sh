#!/bin/bash
# Test python model wrapping inside virtual environment

script_path=$(dirname $BASH_SOURCE)
amici_path=$(cd "$script_path"/.. && pwd)

set -e

if [[ -z "${BNGPATH}" ]]; then
    export BNGPATH=${amici_path}/ThirdParty/BioNetGen-2.3.2
fi

cd "${amici_path}"/python/tests
source "${amici_path}"/build/venv/bin/activate
pip install scipy h5py pytest pytest-cov

# PEtab tests are run separately
pytest --ignore-glob=*petab*
