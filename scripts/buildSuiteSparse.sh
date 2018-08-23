#!/bin/bash
#
# Build SuiteSpare
#
set -e

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

SUITESPARSE_ROOT="${AMICI_PATH}/ThirdParty/SuiteSparse"

for subdir in SuiteSparse_config BTF AMD CAMD COLAMD KLU
  do cd ${SUITESPARSE_ROOT}/${subdir} && make library
done
