#!/bin/bash
# Check test suite with valgrind
# Note: Consider using ctest -T memcheck instead

set -euo pipefail

SCRIPT_PATH=$(dirname "$BASH_SOURCE")
AMICI_PATH=$(cd "$SCRIPT_PATH"/.. && pwd)

cd "${AMICI_PATH}"

cppcheck \
  "-i${AMICI_PATH}/src/doc" \
  "${AMICI_PATH}/src" \
  "-I${AMICI_PATH}/include/" \
  --enable=style \
  "--exitcode-suppressions=${AMICI_PATH}/.cppcheck-exitcode-suppressions"

