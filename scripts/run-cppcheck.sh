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

cd ${AMICI_PATH}

cppcheck ${AMICI_PATH}/src 2> cppcheck.txt
# check if warnings log was created
if [ -f cppcheck.txt  ]; then
    # check if warnings log is empty
    if [ -s cppcheck.txt ]; then
        exit 1
    else
        exit 0
    fi
else
    exit 1
fi
