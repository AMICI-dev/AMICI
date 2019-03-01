#!/bin/bash
# Test python model wrapping inside virtual environment

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

set -e
cd ${AMICI_PATH}/tests
source ${AMICI_PATH}/build/venv/bin/activate
pip install scipy h5py
pip install -U git+https://github.com/pysb/pysb
python testModels.py
python testSBML.py
python testPYSB.py
python testCPP.py
python testPreequilibration.py
python testMisc.py
python testPandas.py
