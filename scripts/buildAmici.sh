#!/bin/bash
#
# Build libamici
#
set -e
cmake=${CMAKE:-cmake}
make=${MAKE:-make}

script_path=$(dirname "$BASH_SOURCE")
amici_path=$(cd "$script_path/.." && pwd)

mkdir -p "${amici_path}/build"
cd "${amici_path}/build"

cpputest_build_dir="${amici_path}"/ThirdParty/cpputest-master/build/

if [[ $TRAVIS = true ]]; then
  # Running on CI server
  build_type="Debug"
else
  build_type="RelWithDebInfo"
fi

CppUTest_DIR=${cpputest_build_dir} \
  ${cmake} -DCMAKE_BUILD_TYPE=$build_type -DPython3_EXECUTABLE="$(command -v python3)" ..

${make}

${make} python-sdist
set -x
