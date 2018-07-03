#!/bin/bash
# Test python model wrapping inside virtual environment

AMICI_PATH="`dirname \"$0\"`"
AMICI_PATH="`( cd \"$AMICI_PATH/..\" && pwd )`"

cd ${AMICI_PATH}/tests

${AMICI_PATH}/build/venv/bin/pip3 install --user scipy h5py
${AMICI_PATH}/build/venv/bin/python3 testModels.py
