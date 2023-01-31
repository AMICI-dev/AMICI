#!/bin/bash
# generate code documentation via sphinx for upload to rtd

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

python3 -m venv ${AMICI_PATH}/doc-venv --clear
source ${AMICI_PATH}/doc-venv/bin/activate
python -m pip install --upgrade --no-cache-dir pip
(cd ${AMICI_PATH}/ && python -m pip install --exists-action=w --no-cache-dir -r documentation/rtd_requirements.txt)
(cd ${AMICI_PATH}/ && python -m pip install --exists-action=w --no-cache-dir -r documentation/rtd_requirements2.txt)

${AMICI_PATH}/scripts/run-sphinx-hasenv.sh

ret=$?
if [[ $ret != 0 ]]; then exit $ret; fi
