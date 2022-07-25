#!/bin/bash
# Build the sphinx documentation in an environment prepared
# as in run-sphinx.sh already

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

source ${AMICI_PATH}/doc-venv/bin/activate

pip3 install git+https://github.com/PEtab-dev/libpetab-python.git@model_pysb

cd ${AMICI_PATH}/documentation

rm -rf ${AMICI_PATH}/documentation/generated

sphinx-build -T -E -W --keep-going -b readthedocs -d _build/doctrees-readthedocs -D language=en . _build/html

ret=$?
if [[ $ret != 0 ]]; then exit $ret; fi
