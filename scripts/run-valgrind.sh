#!/bin/bash
# Check test suite with valgrind
# Note: CppuTest memcheck should be disabled
# Note: Consider using ctest -T memcheck instead

AMICI_PATH="`dirname \"$0\"`"
AMICI_PATH="`( cd \"$AMICI_PATH/..\" && pwd )`"

# read environment variables put there by build script
if [ -f ${AMICI_PATH}/scripts/env.sh ]; then
    . ${AMICI_PATH}/scripts/env.sh
fi

set -e 

# run tests
cd ${AMICI_PATH}/build/tests/cpputest/

VALGRIND_OPTS="--leak-check=full --error-exitcode=1 --trace-children=yes --show-leak-kinds=definite"
set -x
for MODEL in `ctest -N | grep "Test[ ]*#" | grep -v unittests | sed -E 's/ *Test[ ]*#[0-9]+: model_(.*)_test/\1/'`
    do cd ${AMICI_PATH}/build/tests/cpputest/${MODEL}/ && valgrind ${VALGRIND_OPTS} ./model_${MODEL}_test
done
cd ${AMICI_PATH}/tests/cpputest/build/unittests/ && valgrind ${VALGRIND_OPTS} ./unittests
rm ${AMICI_PATH}/tests/cpputest/writeResults.h5
