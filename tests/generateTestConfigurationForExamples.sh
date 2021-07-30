#!/bin/bash
# Generate AMICI configuration for test models

set -eou pipefail

# AMICI root directory
AMICI_PATH=$(dirname "$BASH_SOURCE")
AMICI_PATH=$( cd "$AMICI_PATH/.." && pwd )

# File with test configuration
TEST_FILE="${AMICI_PATH}/tests/cpp/testOptions.h5"

# Delete old config
rm "${TEST_FILE}"

cd "${AMICI_PATH}/tests/generateTestConfig"
./example_dirac.py "${TEST_FILE}"
./example_events.py "${TEST_FILE}"
./example_jakstat.py "${TEST_FILE}"
./example_nested_events.py "${TEST_FILE}"
./example_neuron.py "${TEST_FILE}"
./example_robertson.py "${TEST_FILE}"
./example_steadystate.py "${TEST_FILE}"
./example_calvetti.py "${TEST_FILE}"
