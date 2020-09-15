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
pip install scipy h5py pytest pytest-cov pytest-valgrind

# PEtab tests are run separately
PYTHONMALLOC=malloc valgrind --show-leak-kinds=definite --error-exitcode=1 --leak-check=full --log-file=/tmp/valgrind-output python -m pytest -vv --valgrind --valgrind-log=/tmp/valgrind-output --ignore-glob=*petab*
