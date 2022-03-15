#!/bin/bash
# Check test suite with valgrind
# Note: Consider using ctest -T memcheck instead

SCRIPT_PATH=$(dirname "$BASH_SOURCE")
AMICI_PATH=$(cd "$SCRIPT_PATH/.." && pwd)

set -eou pipefail

# run tests
cd "${AMICI_PATH}/build/"
VALGRIND_OPTS="--leak-check=full --error-exitcode=1 --trace-children=yes --show-leak-kinds=definite"
valgrind ${VALGRIND_OPTS} ctest
