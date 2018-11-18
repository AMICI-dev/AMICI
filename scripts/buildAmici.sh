#!/bin/bash
#
# Build libamici
#
set -e

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

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

rm -f ${AMICI_PATH}/python/sdist/amici/*.cxx
rm -f ${AMICI_PATH}/python/sdist/amici/*.so
rm -f ${AMICI_PATH}/python/sdist/amici/amici.py
rm -f ${AMICI_PATH}/python/sdist/amici/amici_without_hdf5.py

# test install from archive
python3 -m venv ${AMICI_PATH}/build/venvArchive --clear
source ${AMICI_PATH}/build/venvArchive/bin/activate
pip3 install --upgrade pip setuptools pkgconfig wheel numpy
pip3 install $(ls -t ${AMICI_PATH}/build/python/amici-*.tar.gz | head -1)
deactivate

# test install from setup.py
python3 -m venv ${AMICI_PATH}/build/venv --clear
source ${AMICI_PATH}/build/venv/bin/activate
pip3 install --upgrade pip setuptools pkgconfig wheel numpy scipy matplotlib pysb
pip3 install --verbose -e ${AMICI_PATH}/python/sdist
