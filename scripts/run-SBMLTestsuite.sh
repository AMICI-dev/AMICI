#!/bin/bash
# Download and run SBML test suite. Use --jax to run JAX backend.
# The test suite is cloned to $SBML_TEST_SUITE_ROOT if set,
# otherwise to tests/sbml/sbml-test-suite.
REPO_URL="https://github.com/sbmlteam/sbml-test-suite/"
set -e

SBML_TEST_SUITE_ROOT="${SBML_TEST_SUITE_ROOT:-tests/sbml/sbml-test-suite}"
export SBML_TEST_SUITE_ROOT

if [[ ! -d "${SBML_TEST_SUITE_ROOT}" ]]; then
    # only cases/semantic is needed
    git clone --depth=1 --filter=blob:none --sparse "${REPO_URL}" \
      "${SBML_TEST_SUITE_ROOT}"
    git -C "${SBML_TEST_SUITE_ROOT}" sparse-checkout set cases/semantic
fi

source venv/bin/activate

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

pytest "$TEST_SCRIPT" $cases -rfsE -n auto --dist=loadgroup \
  --cov=amici --cov-report=xml:"$COVERAGE_FILE" --cov-append
