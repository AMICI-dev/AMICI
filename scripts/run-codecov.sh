#!/bin/bash
# Check code coverage via codecov

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

lcov --base-directory ${AMICI_PATH} --directory ${AMICI_PATH} --zerocounters -q

cd ${AMICI_PATH}/build
ctest -V
rm ${AMICI_PATH}/tests/cpputest/writeResults.h5
cd ${AMICI_PATH}

lcov --compat-libtool --no-external --directory ${AMICI_PATH}/build/CMakeFiles/amici.dir/src --base-directory ${AMICI_PATH} --capture --output-file coverage.info


source build/venv/bin/activate
pip install coverage
python ./tests/testCoverage.py

rm -rf ./test_model_steadystate_scaled
