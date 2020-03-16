#!/bin/bash
#
# Build Sundials
#
set -e

# shellcheck disable=SC2128
SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd "$SCRIPT_PATH"/.. && pwd)

SUITESPARSE_ROOT="${AMICI_PATH}/ThirdParty/SuiteSparse"
SUNDIALS_BUILD_PATH="${AMICI_PATH}/ThirdParty/sundials/build/"

CMAKE=${CMAKE:-cmake}
MAKE=${MAKE:-make}

if [[ $TRAVIS = true ]]; then
  # Running on CI server
  build_type="Debug"
else
  build_type="RelWithDebInfo"
fi

# enable SuperLUMT support if library exists
SuperLUMT=""
if [[ -f ${AMICI_PATH}/ThirdParty/SuperLU_MT_3.1/lib/libsuperlu_mt_PTHREAD.a ]]
then
    SuperLUMT="-DSUPERLUMT_ENABLE=ON -DBLAS_ENABLE=ON \
         -DSUPERLUMT_INCLUDE_DIR=${AMICI_PATH}/ThirdParty/SuperLU_MT_3.1/SRC/ \
         -DSUPERLUMT_LIBRARY_DIR=${AMICI_PATH}/ThirdParty/SuperLU_MT_3.1/lib/"
fi

mkdir -p "${SUNDIALS_BUILD_PATH}"
cd "${SUNDIALS_BUILD_PATH}"

${CMAKE} -DCMAKE_INSTALL_PREFIX="${SUNDIALS_BUILD_PATH}" \
  -DCMAKE_BUILD_TYPE=$build_type \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
  -DBUILD_ARKODE=OFF \
  -DBUILD_CVODE=OFF \
  -DBUILD_IDA=OFF \
  -DBUILD_KINSOL=OFF \
  -DBUILD_SHARED_LIBS=ON \
  -DBUILD_STATIC_LIBS=ON \
  -DEXAMPLES_ENABLE_C=OFF \
  -DEXAMPLES_INSTALL=OFF \
  -DKLU_ENABLE=ON \
  -DKLU_LIBRARY_DIR="${SUITESPARSE_ROOT}/lib" \
  -DKLU_INCLUDE_DIR="${SUITESPARSE_ROOT}/include" \
  "${SuperLUMT}" \
  ..

${MAKE}
${MAKE} install
