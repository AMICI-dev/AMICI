#!/bin/bash
#
# Build amici tests
#
set -e

AMICI_PATH="`dirname \"$BASH_SOURCE\"`"
AMICI_PATH="`( cd \"$AMICI_PATH/..\" && pwd )`"

${AMICI_PATH}/scripts/buildSuiteSparse.sh
${AMICI_PATH}/scripts/buildSundials.sh
${AMICI_PATH}/scripts/buildAmici.sh

# Build test suite
CPPUTEST_BUILD_DIR=${AMICI_PATH}/ThirdParty/cpputest-master/
cd ${AMICI_PATH}/tests/cpputest/
mkdir -p build_xcode
cd build_xcode
cmake -G"Xcode" -DCMAKE_CXX_STANDARD=11 -DCMAKE_CXX_STANDARD_REQUIRED=ON -DCMAKE_BUILD_TYPE=Debug -DCPPUTEST_DIR=${CPPUTEST_BUILD_DIR} ..

cp ./../expectedResults.h5 ./expectedResults.h5
cp ./../writeResults.h5 ./writeResults.h5

mkdir -p ${AMICI_PATH}/build_xcode
cd ${AMICI_PATH}/build_xcode
cmake -G"Xcode" -DCMAKE_BUILD_TYPE=Debug ..



# Save list of testmodels
ENV_FILE=${AMICI_PATH}/scripts/env.sh
touch ${ENV_FILE} && echo export TESTMODELS="${TESTMODELS}" >> ${ENV_FILE}
