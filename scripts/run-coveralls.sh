#!/bin/bash
# Check code coverage via coveralls

AMICI_PATH="`dirname \"$0\"`"
AMICI_PATH="`( cd \"$AMICI_PATH/..\" && pwd )`"

# read environment variables put there by build script
if [ -f ${AMICI_PATH}/scripts/env.sh ]; then
    . ${AMICI_PATH}/scripts/env.sh
fi

lcov --base-directory ${AMICI_PATH} --directory ${AMICI_PATH} --zerocounters -q

for MODELSTR in $TESTMODELS; do
    MODEL=${MODELSTR#model_}
    cd ${AMICI_PATH}/tests/cpputest/build/${MODEL}
    ./model_${MODEL}_test
done

cd ${AMICI_PATH}/build/CMakeFiles/amici.dir/src

lcov --compat-libtool --no-external --directory . --base-directory ${AMICI_PATH} --capture --output-file coverage.info 

coveralls-lcov coverage.info

