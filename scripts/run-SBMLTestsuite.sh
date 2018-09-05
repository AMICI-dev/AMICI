#!/bin/bash

set -e

# run tests
if [ ! -d "tests/sbml-test-suite" ]; then
	git clone --depth=1 http://github.com/sbmlteam/sbml-test-suite
	mv -f ./sbml-test-suite ./tests/sbml-test-suite
fi

source build/venv/bin/activate
python ./tests/testSBMLSuite.py

cat ./testSuite.txt
