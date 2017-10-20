#!/bin/bash
#
# Build SuiteSpare
#
set -e

AMICI_PATH="`dirname \"$BASH_SOURCE\"`"
AMICI_PATH="`( cd \"$AMICI_PATH/..\" && pwd )`"

SUITESPARSE_ROOT="${AMICI_PATH}/SuiteSparse"

for subdir in SuiteSparse_config BTF CAMD COLAMD KLU
  do cd ${SUITESPARSE_ROOT}/${subdir} && make library
done

# Set DYLD_LIBRARY_PATH for osx
ENV_FILE=${AMICI_PATH}/scripts/env.sh
touch ${ENV_FILE} && echo export DYLD_LIBRARY_PATH=\"\${DYLD_LIBRARY_PATH}:${SUITESPARSE_ROOT}/lib\" >> ${ENV_FILE}
