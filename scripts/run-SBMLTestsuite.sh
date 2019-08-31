#!/bin/bash
# Download and run SBML test suite
REPO_URL="https://github.com/sbmlteam/sbml-test-suite/"
set -e

# run tests
if [ ! -d "tests/sbml-test-suite" ]; then
    git clone --depth=1 ${REPO_URL}
    mv -f ./sbml-test-suite ./tests/sbml-test-suite
fi

source build/venv/bin/activate
if [[ -z "$*" ]]; then
  args="1-1780" # run all tests
else
  args="$@" # use user selection
fi
python ./tests/testSBMLSuite.py "${args}"
ret=$?
if [[ $ret != 0 ]]; then exit $ret; fi
