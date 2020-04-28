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

# test install from setup.py
set +e
python3 -m venv ${AMICI_PATH}/build/venv --clear
# in case this fails (usually due to missing ensurepip, try getting pip
# manually
if [[ $? ]]; then
    set -e
    python3 -m venv ${AMICI_PATH}/build/venv --clear --without-pip
    source ${AMICI_PATH}/build/venv/bin/activate
    curl https://bootstrap.pypa.io/get-pip.py -o ${AMICI_PATH}/build/get-pip.py
    python3 ${AMICI_PATH}/build/get-pip.py
else
    set -e
    source ${AMICI_PATH}/build/venv/bin/activate
fi

# install wheel separately to prevent build_wheel fail in next step
pip install --upgrade wheel
pip install --upgrade pip setuptools pkgconfig scipy matplotlib pysb coverage pytest pytest-cov
pip install --verbose -e ${AMICI_PATH}/python/sdist[petab]
deactivate
