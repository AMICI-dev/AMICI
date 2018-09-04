#!/bin/bash
# Check code coverage via codecov

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

source build/venv/bin/activate
pip install coverage
python ./tests/testCoverage.py
ret=$?
if [[ $ret != 0 ]]; then exit $ret; fi


lcov --compat-libtool --no-external --directory ${AMICI_PATH}/build/CMakeFiles/amici.dir/src --base-directory ${AMICI_PATH} --capture --output-file coverage.info

rm -rf ./test_model_steadystate_scaled
