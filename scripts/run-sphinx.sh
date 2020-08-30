#!/bin/bash
# generate code documentation via sphinx for upload to rtd

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

cd ${AMICI_PATH}/documentation
python3 -m venv ${AMICI_PATH}/doc-venv --clear
source ${AMICI_PATH}/doc-venv/bin/activate
python -m pip install --upgrade --no-cache-dir pip
python -m pip install git+https://github.com/readthedocs/readthedocs-sphinx-ext
python -m pip install --exists-action=w --no-cache-dir -r ${AMICI_PATH}/documentation/rtd_requirements.txt

rm -rf ${AMICI_PATH}/documentation/generated
sphinx-build -T -E --keep-going -b readthedocs -d _build/doctrees-readthedocs -D language=en . _build/html
ret=$?
if [[ $ret != 0 ]]; then exit $ret; fi
