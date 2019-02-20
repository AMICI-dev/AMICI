#!/bin/bash
#
# Build CppUTest
#
set -e

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

# Cpputest
mkdir -p ${AMICI_PATH}/ThirdParty
cd ${AMICI_PATH}/ThirdParty
export CPPUTEST_BUILD_DIR=${AMICI_PATH}/ThirdParty/cpputest-master/

if [ ! -d "cpputest-master" ]; then
    if [ ! -e "cpputest-master.zip" ]; then
        wget -q -O cpputest-master.zip https://codeload.github.com/cpputest/cpputest/zip/23a8ec83bd129b0afe4d9e38509e9bdb489fce55 # https://codeload.github.com/cpputest/cpputest/zip/master
    fi
    unzip -q cpputest-master.zip
    #cd cpputest-master/ && ./autogen.sh && ./configure && make
fi

cd cpputest-master
mkdir -p build
cd build
cmake -DTESTS=OFF -DBUILD_TESTING=OFF -DC++11=ON ..
make -j4

