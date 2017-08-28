#!/bin/bash
AMICI_PATH="`dirname \"$BASH_SOURCE\"`"
AMICI_PATH="`( cd \"$AMICI_PATH/..\" && pwd )`"

# Build dependencies

# SuiteSpare
SUITESPARSE_ROOT="${AMICI_PATH}/SuiteSparse"
cd ${SUITESPARSE_ROOT}/SuiteSparse_config
make library
if [ $? -ne 0 ] ; then
exit 1
fi


cd ${SUITESPARSE_ROOT}/AMD
make library
if [ $? -ne 0 ] ; then
exit 1
fi


cd ${SUITESPARSE_ROOT}/BTF
make library
if [ $? -ne 0 ] ; then
exit 1
fi


cd ${SUITESPARSE_ROOT}/CAMD
make library
if [ $? -ne 0 ] ; then
exit 1
fi


cd ${SUITESPARSE_ROOT}/COLAMD
make library
if [ $? -ne 0 ] ; then
exit 1
fi


cd ${SUITESPARSE_ROOT}/KLU
make library
if [ $? -ne 0 ] ; then
exit 1
fi


# Sundials
SUNDIALS_BUILD_PATH="${AMICI_PATH}/sundials/build/"
mkdir -p ${SUNDIALS_BUILD_PATH}
cd ${SUNDIALS_BUILD_PATH}

cmake -DCMAKE_INSTALL_PREFIX="${SUNDIALS_BUILD_PATH}" \
-DBUILD_ARKODE=OFF \
-DBUILD_CVODE=OFF \
-DBUILD_IDA=OFF \
-DBUILD_KINSOL=OFF \
-DBUILD_SHARED_LIBS=ON \
-DBUILD_STATIC_LIBS=ON \
-DEXAMPLES_ENABLE=OFF \
-DEXAMPLES_INSTALL=OFF \
-DKLU_ENABLE=ON \
-DKLU_LIBRARY_DIR="${SUITESPARSE_ROOT}/lib" \
-DKLU_INCLUDE_DIR="${SUITESPARSE_ROOT}/include" \
..

make
if [ $? -ne 0 ] ; then
exit 1
fi

make install
if [ $? -ne 0 ] ; then
exit 1
fi

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
    if [ $? -ne 0 ] ; then
        exit 1
    fi
fi

# done building dependencies

# Prepare tests
# libamici
mkdir -p ${AMICI_PATH}/build
cd ${AMICI_PATH}/build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
if [ $? -ne 0 ] ; then
    exit 1
fi


TESTMODELS="model_dirac model_steadystate model_jakstat_adjoint model_jakstat_adjoint_o2 model_neuron model_neuron_o2"
for MODEL in $TESTMODELS; do
mkdir -p ${AMICI_PATH}/models/${MODEL}/build
cd ${AMICI_PATH}/models/${MODEL}/build
cmake -DCMAKE_CXX_STANDARD=11 -DCMAKE_CXX_STANDARD_REQUIRED=ON -DCMAKE_BUILD_TYPE=Debug ..
make
if [ $? -ne 0 ] ; then
exit 1
fi
done;


# Build test suite

cd ${AMICI_PATH}/tests/cpputest/
mkdir -p build
cd build
cmake -DCMAKE_CXX_STANDARD=11 -DCMAKE_CXX_STANDARD_REQUIRED=ON -DCMAKE_BUILD_TYPE=Debug -DCPPUTEST_DIR=${CPPUTEST_BUILD_DIR} ..
make
if [ $? -ne 0 ] ; then
exit 1
fi

cd ${AMICI_PATH}

cat >scripts/env.sh <<EOL
export DYLD_LIBRARY_PATH=\"\${DYLD_LIBRARY_PATH}:${SUNDIALS_BUILD_PATH}/lib:${SUITESPARSE_ROOT}/lib\"
export TESTMODELS="${TESTMODELS}"
EOL
