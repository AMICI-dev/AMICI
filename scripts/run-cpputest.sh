#!/bin/bash
AMICI_PATH="`dirname \"$0\"`"
AMICI_PATH="`( cd \"$AMICI_PATH/..\" && pwd )`"

# run tests
cd ${AMICI_PATH}/build

ctest -V
mv ${AMICI_PATH}/tests/cpputest/writeResults.h5 ${AMICI_PATH}/tests/cpputest/writeResults.h5.bak
