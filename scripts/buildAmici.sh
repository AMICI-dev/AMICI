#!/bin/bash
#
# Build libamici
#
set -eou pipefail

cmake=${CMAKE:-cmake}
make=${MAKE:-make}
CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH:-}

script_path=$(dirname "$BASH_SOURCE")
amici_path=$(cd "$script_path/.." && pwd)
amici_build_dir="${amici_path}/build"
mkdir -p "${amici_build_dir}"
cd "${amici_build_dir}"

if [[ "${GITHUB_REPOSITORY:-}" = *"/AMICI" ]] ||
  [ "${ENABLE_AMICI_DEBUGGING:-}" = TRUE ] ||
  [ "${ENABLE_GCOV_COVERAGE:-}" = TRUE ]; then
  # Running on CI server
  build_type="Debug"
  # exceptions instead of terminate()
  extra_cxx_flags=";-Dgsl_CONFIG_CONTRACT_VIOLATION_THROWS;-Dgsl_CONFIG_NARROW_THROWS_ON_TRUNCATION=1;-Werror;-Wno-error=deprecated-declarations"
else
  build_type="RelWithDebInfo"
  extra_cxx_flags=""
fi

# required for build swig interface
venv_dir="${amici_path}/venv"
set +e
mkdir -p "${venv_dir}"
python3 -m venv "${venv_dir}" --clear
# in case this fails (usually due to missing ensurepip, try getting pip
# manually
if [[ $? ]]; then
    set -e
    python3 -m venv "${venv_dir}" --clear --without-pip
    source "${venv_dir}/bin/activate"
    get_pip=${amici_path}/get-pip.py
    curl "https://bootstrap.pypa.io/get-pip.py" -o "${get_pip}"
    python3 "${get_pip}"
    rm "${get_pip}"
else
    set -e
    source "${venv_dir}/bin/activate"
fi

# set python executable for cmake
export PYTHON_EXECUTABLE="${amici_path}/venv/bin/python"
python3 -m pip install numpy

suitesparse_root="${amici_path}/ThirdParty/SuiteSparse"

${cmake} \
  -Wdev -DAMICI_CXX_OPTIONS="-Wall;-Wextra${extra_cxx_flags}" \
  -DCMAKE_BUILD_TYPE=$build_type \
  -DCMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH};${suitesparse_root}/install" \
  -DPython3_EXECUTABLE="$(command -v python3)" \
   ..

# build, with or without sonarcloud wrapper
if [ "${CI_SONARCLOUD:-}" = "TRUE" ]; then
  build-wrapper-linux-x86-64 \
    --out-dir "${amici_path}/bw-output" \
    cmake --build . --parallel
elif [ "${GITHUB_ACTIONS:-}" = "true" ]; then
  cmake --build .
  ${make} python-sdist
else
  cmake --build . --parallel
  ${make} python-sdist
fi
