#!/bin/bash
# Check test suite with valgrind
# Note: CppuTest memcheck should be disabled
# Note: Consider using ctest -T memcheck instead

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

cd ${AMICI_PATH}

cppcheck -i${AMICI_PATH}/src/doc ${AMICI_PATH}/src 2> cppcheck.txt
# check if error log was created
if [ -f cppcheck.txt  ]; then
    # check if error log is empty
    if [ -s cppcheck.txt ]; then
        echo "CPPCHECK failed:"
        cat cppcheck.txt
        rm cppcheck.txt
        exit 1
    else
        exit 0
    fi
else
    exit 1
fi
