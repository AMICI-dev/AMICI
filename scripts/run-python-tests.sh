#!/bin/bash
# Test python model wrapping inside virtual environment

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

set -e
cd ${AMICI_PATH}/tests
source ${AMICI_PATH}/build/venv/bin/activate
pip3 install scipy h5py
python3 testModels.py
python3 testSBML.py
python3 testPYSB.py
python3 testCPP.py
python3 testPreequilibration.py
python3 testMisc.py
python3 testPandas.py
