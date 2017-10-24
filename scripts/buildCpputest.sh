#!/bin/bash
#
# Build CppUTest
#
set -e

AMICI_PATH="`dirname \"$BASH_SOURCE\"`"
AMICI_PATH="`( cd \"$AMICI_PATH/..\" && pwd )`"

# Cpputest
mkdir -p ${AMICI_PATH}/ThirdParty
cd ${AMICI_PATH}/ThirdParty
CPPUTEST_BUILD_DIR=${AMICI_PATH}/ThirdParty/cpputest-master/

if [ ! -d "cpputest-master" ]; then
    if [ ! -e "cpputest-master.zip" ]; then
        wget -O cpputest-master.zip https://codeload.github.com/cpputest/cpputest/zip/master
    fi
    unzip cpputest-master.zip
    cd cpputest-master/ && ./autogen.sh && ./configure && make
fi
