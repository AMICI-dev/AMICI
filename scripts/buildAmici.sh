#!/bin/bash
#
# Build libamici
#
set -e
cmake=${CMAKE:-cmake}
make=${MAKE:-make}

script_path=$(dirname "$BASH_SOURCE")
amici_path=$(cd "$script_path/.." && pwd)
amici_build_dir="${amici_path}/build"
mkdir -p "${amici_build_dir}"
cd "${amici_build_dir}"

cpputest_build_dir="${amici_path}/ThirdParty/cpputest-master/build/"

if [[ $TRAVIS = true ]] || [[ $GITHUB_ACTIONS = true ]] || [[ $ENABLE_AMICI_DEBUGGING = TRUE ]];
then
  # Running on CI server
  build_type="Debug"
else
  build_type="RelWithDebInfo"
fi

CppUTest_DIR=${cpputest_build_dir} \
  ${cmake} \
    -DCMAKE_CXX_FLAGS="-Wall -Wextra -Werror" \
    -DCMAKE_BUILD_TYPE=$build_type \
    -DPython3_EXECUTABLE="$(command -v python3)" ..

# build, with or without sonarcloud wrapper
if [[ "$CI_SONARCLOUD" == "TRUE" ]]; then
  build-wrapper-linux-x86-64 \
    --out-dir "${amici_path}/bw-output" \
    cmake --build . --parallel
elif [[ "$TRAVIS" == "true" ]]; then
  cmake --build .
  ${make} python-sdist
else
  cmake --build . --parallel
  ${make} python-sdist
fi
