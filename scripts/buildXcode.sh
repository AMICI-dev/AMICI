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
${AMICI_PATH}/scripts/buildCpputest.sh


cp ${AMICI_PATH}/tests/cpputest/expectedResults.h5 ./expectedResults.h5

mkdir -p ${AMICI_PATH}/build_xcode
cd ${AMICI_PATH}/build_xcode
CPPUTEST_BUILD_DIR=${AMICI_PATH}/ThirdParty/cpputest-master/build/
CppUTest_DIR=${CPPUTEST_BUILD_DIR} cmake -G"Xcode" -DCMAKE_BUILD_TYPE=Debug ..

for model in steadystate robertson neuron neuron_o2 jakstat_adjoint jakstat_adjoint_o2 dirac events nested_events
do
    cp ${AMICI_PATH}/build/tests/cpputest/external_model_${model}-prefix/src/external_model_${model}-build/libmodel_${model}.a ${AMICI_PATH}/build_xcode/tests/cpputest/external_model_${model}-prefix/src/external_model_${model}-build/libmodel_${model}.a
done
