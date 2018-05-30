#!/bin/bash
# Check test suite with valgrind
# Note: CppuTest memcheck should be disabled
# Note: Consider using ctest -T memcheck instead

AMICI_PATH="`dirname \"$0\"`"
AMICI_PATH="`( cd \"$AMICI_PATH/..\" && pwd )`"

set -e

# run tests
if [ ! -d "tests/sbml-test-suite" ]; then
	git clone http://github.com/sbmlteam/sbml-test-suite
	mv -f ./sbml-test-suite ./tests/sbml-test-suite
fi

python3 ./tests/testSBML.py

cat ./test.txt
rm ./test.txt
