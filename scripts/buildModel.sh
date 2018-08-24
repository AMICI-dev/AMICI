#!/bin/bash
#
# Build amici tests
#
set -e

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

mkdir -p ${AMICI_PATH}/models/$1/build
cd ${AMICI_PATH}/models/$1/build
cmake -DCMAKE_CXX_STANDARD=11 -DCMAKE_CXX_STANDARD_REQUIRED=ON -DCMAKE_BUILD_TYPE=Release ..
make

