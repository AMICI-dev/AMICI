#!/bin/bash
set -eou pipefail

SCRIPT_PATH=$(dirname "$BASH_SOURCE")
AMICI_PATH=$(cd "$SCRIPT_PATH/.." && pwd)

if [[ "${ENABLE_GCOV_COVERAGE:-}" == TRUE ]]; then
  lcov --base-directory "${AMICI_PATH}" \
    --directory "${AMICI_PATH}/build/CMakeFiles/amici.dir/src" \
    --zerocounters -q
fi

# run tests
cd "${AMICI_PATH}/build"

ctest -V
ret=$?
if [[ $ret != 0 ]]; then exit $ret; fi
mv "${AMICI_PATH}/tests/cpp/writeResults.h5" \
  "${AMICI_PATH}/tests/cpp/writeResults.h5.bak"
