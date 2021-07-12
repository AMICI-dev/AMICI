#!/bin/bash
#
# Build amici tests
#
set -e

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

CMAKE=${CMAKE:-cmake}

${AMICI_PATH}/scripts/buildSuiteSparse.sh
${AMICI_PATH}/scripts/buildSundials.sh
${AMICI_PATH}/scripts/buildAmici.sh


cp ${AMICI_PATH}/tests/cpp/expectedResults.h5 ./expectedResults.h5

mkdir -p ${AMICI_PATH}/build_xcode
cd ${AMICI_PATH}/build_xcode
${CMAKE} -G"Xcode" -DCMAKE_BUILD_TYPE=Debug ..

for model in steadystate robertson neuron neuron_o2 jakstat_adjoint jakstat_adjoint_o2 dirac events nested_events
do
    cp ${AMICI_PATH}/build/tests/cpp/external_model_${model}-prefix/src/external_model_${model}-build/libmodel_${model}.a \
       ${AMICI_PATH}/build_xcode/tests/cpp/external_model_${model}-prefix/src/external_model_${model}-build/libmodel_${model}.a
done
