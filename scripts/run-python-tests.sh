#!/bin/bash
# Test python model wrapping inside virtual environment

AMICI_PATH="`dirname \"$0\"`"
AMICI_PATH="`( cd \"$AMICI_PATH/..\" && pwd )`"

cd ${AMICI_PATH}/tests
source ${AMICI_PATH}/build/venv/activate
pip3 install --upgrade pip
pip3 install --user scipy h5py
python3 testModels.py
