#!/bin/bash
#
# Build SuiteSpare
#
set -e

script_path=$(dirname "$BASH_SOURCE")
amici_path=$(cd "$script_path/.." && pwd)

suitesparse_root="${amici_path}/ThirdParty/SuiteSparse"

export CFLAGS="-D'SuiteSparse_long=int64_t' -D'SuiteSparse_long_max=9223372036854775807' -D'SuiteSparse_long_idd=\"lld\"'"
for subdir in SuiteSparse_config BTF AMD CAMD COLAMD KLU
  do cd "${suitesparse_root}/${subdir}" && make library
done
