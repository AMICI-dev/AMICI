#!/bin/bash
SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

lcov --base-directory ${AMICI_PATH} --directory ${AMICI_PATH}/build/CMakeFiles/amici.dir/src --zerocounters -q

# run tests
cd ${AMICI_PATH}/build

ctest -V
mv ${AMICI_PATH}/tests/cpputest/writeResults.h5 ${AMICI_PATH}/tests/cpputest/writeResults.h5.bak
