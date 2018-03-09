#!/bin/bash
AMICI_PATH="`dirname \"$BASH_SOURCE\"`"
AMICI_PATH="`( cd \"$AMICI_PATH/..\" && pwd )`"

cd ${AMICI_PATH}/tests

TEST_FILE="cpputest/testOptions.h5"
rm ${TEST_FILE}
./example_dirac.py ${TEST_FILE}
./example_events.py ${TEST_FILE}
./example_jakstat.py ${TEST_FILE}
./example_nested_events.py ${TEST_FILE}
./example_neuron.py ${TEST_FILE}
./example_robertson.py ${TEST_FILE}
./example_steadystate.py ${TEST_FILE}

