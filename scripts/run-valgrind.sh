#!/bin/bash
# Check test suite with valgrind
# Note: CppuTest memcheck should be disabled
# Note: Consider using ctest -T memcheck instead

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

set -e

# run tests
cd ${AMICI_PATH}/build/tests/cpputest/

VALGRIND_OPTS="--leak-check=full --error-exitcode=1 --trace-children=yes --show-leak-kinds=definite"
set -x
for MODEL in $(ctest -N | grep "Test[ ]*#" | grep -v unittests | sed --regexp-extended 's/ *Test[ ]*#[0-9]+: model_(.*)_test/\1/')
    do cd ${AMICI_PATH}/build/tests/cpputest/${MODEL}/ && valgrind ${VALGRIND_OPTS} ./model_${MODEL}_test
done
cd ${AMICI_PATH}/build/tests/cpputest/unittests/ && valgrind ${VALGRIND_OPTS} ./unittests
