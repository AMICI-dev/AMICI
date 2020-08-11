#!/bin/bash
# generate code documentation via sphinx and upload to rtd

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

cd ${AMICI_PATH}/documentation
python3 -m venv ${AMICI_PATH}/doc-venv --clear
source ${AMICI_PATH}/doc-venv/bin/activate
pip3 install -r ${AMICI_PATH}/documentation/rtd_requirements.txt

sphinx-build -T -E -W -D language=en -b html . _build
ret=$?
if [[ $ret != 0 ]]; then exit $ret; fi
