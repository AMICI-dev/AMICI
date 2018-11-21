#!/bin/bash
#
# Build AMICI along with dependencies and test suite
#
set -e

SCRIPT_PATH=$(dirname $BASH_SOURCE)
SCRIPT_PATH=$(cd $SCRIPT_PATH && pwd)

${SCRIPT_PATH}/buildSuiteSparse.sh
${SCRIPT_PATH}/buildSundials.sh
${SCRIPT_PATH}/buildCpputest.sh
${SCRIPT_PATH}/buildBNGL.sh
${SCRIPT_PATH}/buildAmici.sh
