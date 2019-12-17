#!/bin/bash
#
# Build libamici
#
set -e
CMAKE=${CMAKE:-cmake}
MAKE=${MAKE:-make}

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

mkdir -p ${AMICI_PATH}/build
cd ${AMICI_PATH}/build
CPPUTEST_BUILD_DIR=${AMICI_PATH}/ThirdParty/cpputest-master/build/
CppUTest_DIR=${CPPUTEST_BUILD_DIR} \
  ${CMAKE} -DCMAKE_BUILD_TYPE=Debug -DPython_EXECUTABLE=$(which python3) ..
${MAKE}

${MAKE} python-sdist
set -x
