#!/bin/bash
#
# Build SuiteSpare
#
set -e

script_path=$(dirname "$BASH_SOURCE")
amici_path=$(cd "$script_path/.." && pwd)

suitesparse_root="${amici_path}/ThirdParty/SuiteSparse"

cd "${suitesparse_root}/build"

if [[ "${GITHUB_REPOSITORY:-}" = *"/AMICI" ]] ||
  [ "${ENABLE_AMICI_DEBUGGING:-}" = TRUE ]; then
  # Running on CI server
  build_type="Debug"
else
  build_type="RelWithDebInfo"
fi

cmake -DSUITESPARSE_ENABLE_PROJECTS="amd;btf;colamd;klu" \
  -DCMAKE_BUILD_TYPE="$build_type" \
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
  -DBUILD_SHARED_LIBS=OFF \
  -DBUILD_TESTING=OFF \
  ..

cmake --build .
cmake --install .
