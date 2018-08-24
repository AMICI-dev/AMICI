#!/bin/bash
# Test python model wrapping inside virtual environment

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

set -e
cd ${AMICI_PATH}/tests
source ${AMICI_PATH}/build/venv/bin/activate
pip3 install scipy h5py
python3 testModels.py
