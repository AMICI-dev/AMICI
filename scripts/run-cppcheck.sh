#!/bin/bash
# Check test suite with valgrind
# Note: CppuTest memcheck should be disabled
# Note: Consider using ctest -T memcheck instead

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

cd ${AMICI_PATH}

cppcheck -i${AMICI_PATH}/src/doc ${AMICI_PATH}/src  -I${AMICI_PATH}/include/ --enable=performance,portability,missingInclude,style 2> cppcheck_pre.txt

# suppress alloca warnings
grep -v "(warning) Obsolete function 'alloca' called." cppcheck_pre.txt > cppcheck_pre1.txt
rm cppcheck_pre.txt

# suppress header warnings for standard libraries
grep -v "Cppcheck cannot find all the include files" cppcheck_pre1.txt > cppcheck_pre2.txt
rm cppcheck_pre1.txt

grep -v "'AmiVectorArray' does not have a operator=" cppcheck_pre2.txt > cppcheck.txt
rm cppcheck_pre2.txt

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
