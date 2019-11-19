#!/usr/bin/env bash
# Download and build Doxygen (in case apt or homebrew version is buggy again)
set -e

SCRIPT_PATH=$(dirname "$BASH_SOURCE")
AMICI_PATH=$(cd "$SCRIPT_PATH"/.. && pwd)

cd "${AMICI_PATH}"/ThirdParty
if [[ ! -d ${SWIG_DIR} ]]; then
  git clone --depth 1 https://github.com/doxygen/doxygen.git
fi

cd doxygen
mkdir -p build
cd build
cmake -G "Unix Makefiles" ..
make -j2
make install
