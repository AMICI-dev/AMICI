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
pip install scipy h5py pytest

if [[ ! -e "valgrind-python.supp" ]]; then
    wget -O valgrind-python.supp http://svn.python.org/projects/python/trunk/Misc/valgrind-python.supp
fi

# PEtab tests are run separately
PYTHONMALLOC=malloc valgrind --suppressions=valgrind-python.supp --show-leak-kinds=definite --error-exitcode=1 --leak-check=full python -m pytest -vv --ignore-glob=*petab*
