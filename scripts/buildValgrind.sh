#!/bin/bash
#
# Build Valgrind
#
set -e

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

# Valgrind
mkdir -p ${AMICI_PATH}/ThirdParty
cd ${AMICI_PATH}/ThirdParty

if [ ! -d "valgrind-3.14.0" ]; then
    if [ ! -e "valgrind.tar.bz2" ]; then
        wget -q -O valgrind-3.14.0.tar.bz2 http://www.valgrind.org/downloads/valgrind-3.14.0.tar.bz2
    fi
    tar xjf valgrind-3.14.0.tar.bz2
fi

cd valgrind-3.14.0
./configure
make
sudo make install
cd ..
