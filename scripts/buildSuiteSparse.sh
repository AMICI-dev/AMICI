#!/bin/bash
#
# Build SuiteSpare
#
set -e

script_path=$(dirname "$BASH_SOURCE")
amici_path=$(cd "$script_path/.." && pwd)

suitesparse_root="${amici_path}/ThirdParty/SuiteSparse"
export CMAKE_OPTIONS="-DBLA_VENDOR=All -DALLOW_64BIT_BLAS=ON -DENABLE_CUDA=FALSE -DNFORTRAN=TRUE"
for subdir in SuiteSparse_config BTF AMD CAMD COLAMD KLU
  do cd "${suitesparse_root}/${subdir}" && make local install
done
