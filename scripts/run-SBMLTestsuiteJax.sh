#!/bin/bash
# Download and run SBML test suite with JAX
REPO_URL="https://github.com/sbmlteam/sbml-test-suite/"
set -e

if [[ ! -d "tests/sbml-test-suite" ]]; then
    git clone --depth=1 ${REPO_URL}
    mv -f ./sbml-test-suite ./tests/sbml-test-suite
fi

source venv/bin/activate
export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}."
pip show pytest-xdist > /dev/null 2>&1 || pip install pytest-xdist
pip install coverage pytest-cov

if [[ -z "$*" ]]; then
  cases=""
else
  cases="--cases=$@"
fi

RESULT_DIR=tests/sbml/amici-semantic-results-jax
if [[ -d "${RESULT_DIR}" ]]; then
  rm -rf "${RESULT_DIR}"
fi
mkdir "${RESULT_DIR}"

pytest ./tests/sbml/testSBMLSuiteJax.py $cases -rfsE -n auto \
  --cov=amici --cov-report=xml:"coverage_SBMLSuite_jax.xml" --cov-append
