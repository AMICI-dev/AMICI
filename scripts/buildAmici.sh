#!/bin/bash
#
# Build libamici
#
set -e

AMICI_PATH="`dirname \"$BASH_SOURCE\"`"
AMICI_PATH="`( cd \"$AMICI_PATH/..\" && pwd )`"

mkdir -p ${AMICI_PATH}/build
cd ${AMICI_PATH}/build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make

python3 ${AMICI_PATH}/build/swig/python/setup.py install --user
