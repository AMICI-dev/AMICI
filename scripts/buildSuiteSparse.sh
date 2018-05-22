#!/bin/bash
#
# Build SuiteSpare
#
set -e

AMICI_PATH="`dirname \"$BASH_SOURCE\"`"
AMICI_PATH="`( cd \"$AMICI_PATH/..\" && pwd )`"

SUITESPARSE_ROOT="${AMICI_PATH}/ThirdParty/SuiteSparse"

for subdir in SuiteSparse_config BTF AMD CAMD COLAMD KLU
  do cd ${SUITESPARSE_ROOT}/${subdir} && make library
done
