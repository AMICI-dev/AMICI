#!/usr/bin/env bash
#
# Check that CMake build succeeds for one of the models generated during
# Python tests
#

set -e

SCRIPT_DIR=$(dirname $BASH_SOURCE)
MODEL_DIR=${SCRIPT_DIR}/test_model_presimulation/

cd ${MODEL_DIR}
mkdir -p build
cd build

cmake ..
make
