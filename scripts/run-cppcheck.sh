#!/bin/bash
# Check test suite with valgrind
# Note: CppuTest memcheck should be disabled
# Note: Consider using ctest -T memcheck instead

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

cd ${AMICI_PATH}

cppcheck -i${AMICI_PATH}/src/doc ${AMICI_PATH}/src  -I${AMICI_PATH}/include/ --enable=style 2> cppcheck.txt

# suppress alloca warnings
grep -v "(warning) Obsolete function 'alloca' called." cppcheck.txt > cppcheck_tmp.txt
mv cppcheck_tmp.txt cppcheck.txt


# suppress header warnings for standard libraries
grep -v "Cppcheck cannot find all the include files" cppcheck.txt > cppcheck_tmp.txt
mv cppcheck_tmp.txt cppcheck.txt

grep -v "'AmiVectorArray' does not have a operator=" cppcheck.txt > cppcheck_tmp.txt
mv cppcheck_tmp.txt cppcheck.txt

grep -v "Member variable 'ExpData::nytrue_' is not initialized in the constructor" cppcheck.txt > cppcheck_tmp.txt
mv cppcheck_tmp.txt cppcheck.txt

grep -v "Member variable 'ExpData::nztrue_' is not initialized in the constructor" cppcheck.txt > cppcheck_tmp.txt
mv cppcheck_tmp.txt cppcheck.txt

grep -v "Member variable 'ExpData::nmaxevent_' is not initialized in the constructor" cppcheck.txt > cppcheck_tmp.txt
mv cppcheck_tmp.txt cppcheck.txt

# check if error log was created
if [ -f cppcheck.txt  ]; then
    # check if error log is empty
    if [ -s cppcheck.txt ]; then
        echo "CPPCHECK failed:"
        cat cppcheck.txt
        rm cppcheck.txt
        exit 1
    else
        rm cppcheck.txt
        exit 0
    fi
else
    exit 1
fi
