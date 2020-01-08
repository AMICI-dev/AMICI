#!/bin/bash
# Test python model wrapping inside virtual environment

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

set -e

if [[ -z "${BNGPATH}" ]]; then
    export BNGPATH=${AMICI_PATH}/ThirdParty/BioNetGen-2.3.2
fi

cd ${AMICI_PATH}/tests
source ${AMICI_PATH}/build/venv/bin/activate
pip install scipy h5py
python testModels.py
python testSBML.py
python testPYSB.py
python testCPP.py
python testPreequilibration.py
python testMisc.py
python testPandas.py

./test-petab-import.sh
