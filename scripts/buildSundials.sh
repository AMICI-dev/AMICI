#!/bin/bash
#
# Build Sundials
#
set -e

# shellcheck disable=SC2128
script_path=$(dirname "$BASH_SOURCE")
amici_path=$(cd "$script_path/.." && pwd)

suitesparse_root="${amici_path}/ThirdParty/SuiteSparse"
sundials_build_path="${amici_path}/ThirdParty/sundials/build/"

cmake=${CMAKE:-cmake}
make=${MAKE:-make}

if [[ "${GITHUB_REPOSITORY:-}" = *"/AMICI" ]] ||
  [ "${ENABLE_AMICI_DEBUGGING:-}" = TRUE ]; then
  # Running on CI server
  build_type="Debug"
else
  build_type="RelWithDebInfo"
fi

# enable SuperLUMT support if library exists
SuperLUMT=""
if [[ -f "${amici_path}/ThirdParty/SuperLU_MT_3.1/lib/libsuperlu_mt_PTHREAD.a" ]]
then
    SuperLUMT="-DSUNDIALS_ENABLE_SUPERLUMT=ON \
         -DSUPERLUMT_INCLUDE_DIR=\"${amici_path}/ThirdParty/SuperLU_MT_3.1/SRC/\" \
         -DSUPERLUMT_LIBRARY_DIR=\"${amici_path}/ThirdParty/SuperLU_MT_3.1/lib/\""
fi

mkdir -p "${sundials_build_path}"
cd "${sundials_build_path}"

${cmake} -DCMAKE_INSTALL_PREFIX="${sundials_build_path}" \
  -DCMAKE_BUILD_TYPE="$build_type" \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
  -DSUNDIALS_ENABLE_ARKODE=OFF \
  -DSUNDIALS_ENABLE_CVODE=OFF \
  -DSUNDIALS_ENABLE_IDA=OFF \
  -DSUNDIALS_ENABLE_KINSOL=OFF \
  -DBUILD_SHARED_LIBS=OFF \
  -DBUILD_STATIC_LIBS=ON \
  -DSUNDIALS_ENABLE_NVECTOR_MANYVECTOR=OFF \
  -DSUNDIALS_ENABLE_C_EXAMPLES=OFF \
  -DSUNDIALS_ENABLE_EXAMPLES_INSTALL=OFF \
  -DSUNDIALS_ENABLE_KLU=ON \
  -DCMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH};${suitesparse_root}/install" \
  ${SuperLUMT} \
  ..

${make}
${make} install
