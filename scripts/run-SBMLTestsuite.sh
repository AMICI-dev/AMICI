#!/bin/bash
# Download and run SBML test suite. Use --jax to run JAX backend.
REPO_URL="https://github.com/sbmlteam/sbml-test-suite/"
set -e

# run tests
if [[ ! -d "tests/sbml/sbml-test-suite" ]]; then
    git clone --depth=1 ${REPO_URL}
    mv -f ./sbml-test-suite ./tests/sbml/sbml-test-suite
fi

source venv/bin/activate
export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}."
pip show pytest-xdist > /dev/null 2>&1 || pip install pytest-xdist
pip install coverage pytest-cov

BACKEND="python"
if [[ "$1" == "--jax" ]]; then
  BACKEND="jax"
  shift
fi

if [[ -z "$*" ]]; then
  cases=""
else
  cases="--cases=$@"
fi

if [[ "$BACKEND" == "jax" ]]; then
  RESULT_DIR=tests/sbml/amici-semantic-results-jax
  TEST_SCRIPT=./tests/sbml/testSBMLSuiteJax.py
  COVERAGE_FILE=coverage_SBMLSuite_jax.xml
else
  RESULT_DIR=tests/sbml/amici-semantic-results
  TEST_SCRIPT=./tests/sbml/testSBMLSuite.py
  COVERAGE_FILE=coverage_SBMLSuite.xml
fi

if [[ -d "${RESULT_DIR}" ]]; then
  rm -rf "${RESULT_DIR}"
fi
mkdir "${RESULT_DIR}"

pytest "$TEST_SCRIPT" $cases -rfsE -n auto \
  --cov=amici --cov-report=xml:"$COVERAGE_FILE" --cov-append
