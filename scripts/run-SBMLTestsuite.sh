#!/bin/bash
# Download and run SBML test suite
REPO_URL="https://github.com/sbmlteam/sbml-test-suite/"
set -e

# run tests
if [[ ! -d "tests/sbml-test-suite" ]]; then
    git clone --depth=1 ${REPO_URL}
    mv -f ./sbml-test-suite ./tests/sbml-test-suite
fi

source build/venv/bin/activate
pip show pytest-xdist > /dev/null 2>&1 || pip install pytest-xdist
pip install coverage pytest-cov

if [[ -z "$*" ]]; then
  args="1-1780" # run all tests
else
  args="$@" # use user selection
fi

# delete old result directory and recreate
RESULT_DIR=tests/amici-semantic-results
if [[ -d "${RESULT_DIR}" ]]; then
  rm -rf "${RESULT_DIR}"
fi
mkdir "${RESULT_DIR}"

pytest ./tests/testSBMLSuite.py --cases="${args}" -rfsE -n auto \
  --cov=amici --cov-report=xml:"coverage_SBMLSuite.xml" --cov-append
