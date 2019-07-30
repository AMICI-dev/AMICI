#!/bin/bash
#
# Build libamici
#
set -e

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

# Disabled until cmake package is made compatible with updated setup.py
#make python-wheel
#pip3 install --user --prefix= `ls -t ${AMICI_PATH}/build/python/amici-*.whl | head -1`

rm -f ${AMICI_PATH}/python/sdist/amici/*.cxx
rm -f ${AMICI_PATH}/python/sdist/amici/*.so
rm -f ${AMICI_PATH}/python/sdist/amici/amici.py
rm -f ${AMICI_PATH}/python/sdist/amici/amici_without_hdf5.py

# test install from archive
python3 -m venv ${AMICI_PATH}/build/venvArchive --clear
source ${AMICI_PATH}/build/venvArchive/bin/activate
pip install --upgrade pip setuptools pkgconfig wheel
pip install $(ls -t ${AMICI_PATH}/build/python/amici-*.tar.gz | head -1)
deactivate
