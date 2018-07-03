#!/bin/bash
# Test python model wrapping inside virtual environment

AMICI_PATH="`dirname \"$0\"`"
AMICI_PATH="`( cd \"$AMICI_PATH/..\" && pwd )`"

set -e
cd ${AMICI_PATH}/tests
source ${AMICI_PATH}/build/venv/bin/activate
pip3 install --upgrade pip
pip3 install --user scipy h5py
python3 testModels.py
