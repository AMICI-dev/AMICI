#!/bin/bash
#
# Build libamici
#
set -e
CMAKE=${CMAKE:-cmake}
MAKE=${MAKE:-make}

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd "$SCRIPT_PATH"/.. && pwd)

mkdir -p "${AMICI_PATH}"/build
cd "${AMICI_PATH}"/build

CPPUTEST_BUILD_DIR=${AMICI_PATH}/ThirdParty/cpputest-master/build/

if [[ $TRAVIS = true ]]; then
  # Running on CI server
  build_type="Debug"
else
  build_type="RelWithDebInfo"
fi

CppUTest_DIR=${CPPUTEST_BUILD_DIR} \
  ${CMAKE} -DCMAKE_BUILD_TYPE=$build_type -DPython3_EXECUTABLE="$(command -v python3)" ..

${MAKE}

${MAKE} python-sdist
set -x
