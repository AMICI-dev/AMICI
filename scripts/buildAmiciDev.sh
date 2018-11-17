#!/bin/bash
#
# Build libamici
#
set -e

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

cd ${AMICI_PATH}/build
python3 -m venv ${AMICI_PATH}/build/venvDev --clear
rm -f ${AMICI_PATH}/python/sdist/amici/*.so
rm -f ${AMICI_PATH}/python/sdist/amici/*.cxx
rm -f ${AMICI_PATH}/python/sdist/amici/*.so
rm -f ${AMICI_PATH}/python/sdist/amici/amici.py
rm -f ${AMICI_PATH}/python/sdist/amici/amici_without_hdf5.py
source ${AMICI_PATH}/build/venvDev/bin/activate
pip3 install --upgrade pip setuptools pkgconfig wheel numpy scipy matplotlib jupyter pysb
export ENABLE_AMICI_DEBUGGING=TRUE
pip3 install --verbose -e ${AMICI_PATH}/python/sdist
