#!/bin/bash
# Check code coverage via codecov

script_path=$(dirname $BASH_SOURCE)
amici_path=$(cd "$script_path"/.. && pwd)

source "${amici_path}"/build/venv/bin/activate
pip install coverage pytest pytest-cov

if [[ -z "${BNGPATH}" ]]; then
    export BNGPATH="${amici_path}"/ThirdParty/BioNetGen-2.7.0
fi

pytest \
  --ignore-glob=*petab* \
  --cov=amici \
  --cov-report=xml:"${amici_path}"/coverage_py.xml \
  --cov-append \
  "${amici_path}"/python/tests

ret=$?
if [[ $ret != 0 ]]; then exit $ret; fi
