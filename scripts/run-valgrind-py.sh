#!/bin/bash
# Test python model wrapping inside virtual environment

script_path=$(dirname $BASH_SOURCE)
amici_path=$(cd "$script_path"/.. && pwd)

set -e

if [[ -z "${BNGPATH}" ]]; then
    export BNGPATH=${amici_path}/ThirdParty/BioNetGen-2.7.0
fi

cd "${amici_path}"/python/tests
source "${amici_path}"/build/venv/bin/activate
pip install scipy h5py pytest

# PEtab tests are run separately
python -m pytest -vv --ignore-glob=*petab* \
  --ignore-glob=test_sbml_import_special_functions.py

#PYTHONMALLOC=malloc valgrind \
#  --suppressions=valgrind-python.supp \
#  --show-leak-kinds=definite \
#  --errors-for-leak-kinds=definite \
#  --error-exitcode=1 \
#  --leak-check=full \
#  --gen-suppressions=all \
#  -v \
#  python -m pytest -vv --ignore-glob=*petab* \
#    -k "not test_sbml2amici_observable_dependent_error"
