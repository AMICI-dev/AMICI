#!/bin/bash
#
# Build SuiteSpare
#
set -e

script_path=$(dirname "$BASH_SOURCE")
amici_path=$(cd "$script_path/.." && pwd)

suitesparse_root="${amici_path}/ThirdParty/SuiteSparse"

cd "${suitesparse_root}/build"
cmake -DSUITESPARSE_ENABLE_PROJECTS="amd;btf;colamd;klu" \
  -DCMAKE_INSTALL_PREFIX="${suitesparse_root}/install" \
  -DCHOLMOD_CAMD=OFF \
  -DKLU_USE_CHOLMOD=OFF \
  -DSUITESPARSE_CONFIG_USE_OPENMP=OFF \
  -DSUITESPARSE_USE_CUDA=OFF \
  -DSUITESPARSE_USE_FORTRAN=OFF \
  -DSUITESPARSE_USE_PYTHON=OFF \
  -DSUITESPARSE_USE_OPENMP=OFF \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
  -DSUITESPARSE_USE_64BIT_BLAS=ON \
  -DBLAS_LIBRARIES=dummy \
  ..

cmake --build .
cmake --install .
