#!/bin/bash
# generate code documentation via sphinx and upload to rtd

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

cd ${AMICI_PATH}
python3 -m venv ${AMICI_PATH}/doc-venv --clear
source ${AMICI_PATH}/doc-venv/bin/activate
python -m pip install --upgrade --no-cache-dir pip
python -m pip install --upgrade --no-cache-dir Pygments==2.3.1
setuptools==41.0.1 docutils==0.14 mock==1.0.1 pillow==5.4.1 alabaster==0.7.12  commonmark==0.8.1 recommonmark==0.5.0 sphinx<2 sphinx-rtd-theme<0.5 readthedocs-sphinx-ext<1.1
python -m pip install --exists-action=w --no-cache-dir -r ${AMICI_PATH}/documentation/rtd_requirements.txt

rm -rf ${AMICI_PATH}/documentation/generated
sphinx-build -T -E -W --keep-going -b html -d _build/doctrees-readthedocs -D language=en -c documentation . _build/html
ret=$?
if [[ $ret != 0 ]]; then exit $ret; fi
