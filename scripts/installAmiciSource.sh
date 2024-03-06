#!/bin/bash
# Create a virtual environment and perform an editable amici installation
set -e

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd "$SCRIPT_PATH/.." && pwd)

# Disabled until cmake package is made compatible with updated setup.py
#make python-wheel
#pip3 install --user --prefix= `ls -t ${AMICI_PATH}/build/python/amici-*.whl | head -1`

# test install from setup.py
source "${AMICI_PATH}/venv/bin/activate"
python -m pip install --upgrade pip wheel
python -m pip install --upgrade pip setuptools cmake_build_extension numpy
python -m pip install git+https://github.com/FFroehlich/pysb@fix_pattern_matching # pin to PR for SPM with compartments
AMICI_BUILD_TEMP="${AMICI_PATH}/python/sdist/build/temp" \
  python -m pip install --verbose -e "${AMICI_PATH}/python/sdist[petab,test,vis]" --no-build-isolation
deactivate
