#!/bin/bash
#
# Build CppUTest
#
set -e

script_path=$(dirname "$BASH_SOURCE")
amici_path=$(cd "$script_path/.." && pwd)

cmake=${CMAKE:-cmake}
make=${MAKE:-make}

# Cpputest
mkdir -p "${amici_path}/ThirdParty"
cd "${amici_path}/ThirdParty"
export CPPUTEST_BUILD_DIR="${amici_path}/ThirdParty/cpputest-master/"

if [ ! -d "cpputest-master" ]; then
    if [ ! -e "cpputest-master.zip" ]; then
        wget -q -O cpputest-master.zip https://codeload.github.com/cpputest/cpputest/zip/master
    fi
    unzip -q cpputest-master.zip
    #cd cpputest-master/ && ./autogen.sh && ./configure && make
fi

cd cpputest-master
mkdir -p build
cd build
${cmake} -DTESTS=OFF -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release \
  -DC++11=ON -DMEMORY_LEAK_DETECTION=OFF ..
${make} -j4

