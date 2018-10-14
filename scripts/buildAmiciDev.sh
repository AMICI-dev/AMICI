#!/bin/bash
#
# Build libamici
#
set -e

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

cd ${AMICI_PATH}/build
python3 -m venv ${AMICI_PATH}/build/venvDev --clear
source ${AMICI_PATH}/build/venvDev/bin/activate
pip3 install --upgrade pip setuptools pkgconfig wheel numpy scipy matplotlib jupyter
export ENABLE_AMICI_DEBUGGING=TRUE
pip3 install --verbose -e ${AMICI_PATH}/python/sdist
