#!/bin/bash
#
# Build amici tests
#
set -e

AMICI_PATH="`dirname \"$BASH_SOURCE\"`"
AMICI_PATH="`( cd \"$AMICI_PATH/..\" && pwd )`"

# Build test models
TESTMODELS="model_dirac model_steadystate model_jakstat_adjoint model_jakstat_adjoint_o2 model_neuron model_neuron_o2 model_events model_nested_events model_robertson"
for MODEL in $TESTMODELS; do
    mkdir -p ${AMICI_PATH}/models/${MODEL}/build
    cd ${AMICI_PATH}/models/${MODEL}/build
    cmake -DCMAKE_CXX_STANDARD=11 -DCMAKE_CXX_STANDARD_REQUIRED=ON -DCMAKE_BUILD_TYPE=Debug ..
    make
done;


# Build test suite
CPPUTEST_BUILD_DIR=${AMICI_PATH}/ThirdParty/cpputest-master/
cd ${AMICI_PATH}/tests/cpputest/
mkdir -p build
cd build
cmake -DCMAKE_CXX_STANDARD=11 -DCMAKE_CXX_STANDARD_REQUIRED=ON -DCMAKE_BUILD_TYPE=Debug -DCPPUTEST_DIR=${CPPUTEST_BUILD_DIR} ..
make



# Save list of testmodels
ENV_FILE=${AMICI_PATH}/scripts/env.sh
touch ${ENV_FILE} && echo export TESTMODELS="${TESTMODELS}" >> ${ENV_FILE}
