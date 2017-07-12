#!/bin/bash
AMICI_PATH="`dirname \"$0\"`"
AMICI_PATH="`( cd \"$AMICI_PATH/..\" && pwd )`"

cd ${AMICI_PATH}/tests/cpputest/build

ctest -V
if [ $? -ne 0 ] ; then
exit 1
fi

cd ${AMICI_PATH}
