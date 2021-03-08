#!/usr/bin/env bash
# Download and build Doxygen (in case apt or homebrew version is buggy again)
set -e

SCRIPT_PATH=$(dirname "$BASH_SOURCE")
AMICI_PATH=$(cd "$SCRIPT_PATH"/.. && pwd)

DOXYGEN_DIR="${AMICI_PATH}"/ThirdParty/doxygen
cd "${AMICI_PATH}"/ThirdParty
if [[ ! -d ${DOXYGEN_DIR} ]]; then
  # git clone --depth 1 https://github.com/doxygen/doxygen.git "${DOXYGEN_DIR}"
  git clone --single-branch --branch Release_1_9_1 --depth 1 https://github.com/doxygen/doxygen.git "${DOXYGEN_DIR}"
fi

cd "${DOXYGEN_DIR}"
mkdir -p build
cd build
cmake -G "Unix Makefiles" ..
make -j2
make install
