#!/bin/bash
#
# Build libamici
#
set -e

AMICI_PATH="`dirname \"$BASH_SOURCE\"`"
AMICI_PATH="`( cd \"$AMICI_PATH/..\" && pwd )`"

mkdir -p ${AMICI_PATH}/build
cd ${AMICI_PATH}/build
CPPUTEST_BUILD_DIR=${AMICI_PATH}/ThirdParty/cpputest-master/build/
CppUTest_DIR=${CPPUTEST_BUILD_DIR} cmake -DCMAKE_BUILD_TYPE=Debug ..
make

# Disabled until cmake package is made compatible with updated setup.py
#make python-wheel
#pip3 install --user --prefix= `ls -t ${AMICI_PATH}/build/python/amici-*.whl | head -1`

make python-sdist
set -x
python3 -m venv ${AMICI_PATH}/build/venv --clear
source ${AMICI_PATH}/build/venv/bin/activate
pip3 install --upgrade pip
pip3 install pkgconfig wheel
pip3 install --verbose `ls -t ${AMICI_PATH}/build/python/amici-*.tar.gz | head -1`
