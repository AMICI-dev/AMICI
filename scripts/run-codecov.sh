#!/bin/bash
# Check code coverage via codecov

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

source ${AMICI_PATH}/build/venv/bin/activate
pip install coverage

if [[ -z "${BNGPATH}" ]]; then
    export BNGPATH=${AMICI_PATH}/ThirdParty/BioNetGen-2.3.2
fi

python ./tests/testCoverage.py
ret=$?
if [[ $ret != 0 ]]; then exit $ret; fi

rm -rf ./test_model_steadystate_scaled
