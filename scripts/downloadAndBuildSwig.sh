#!/usr/bin/env bash
# Download and build SWIG
set -e

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

SWIG_URL='http://prdownloads.sourceforge.net/swig/swig-4.0.1.tar.gz'
SWIG_ARCHIVE='swig-4.0.1.tar.gz'
SWIG_DIR='swig-4.0.1'
PREFIX=${AMICI_PATH}/ThirdParty/${SWIG_DIR}/install
SWIG_BIN_DIR=${PREFIX}/bin

cd ${AMICI_PATH}/ThirdParty/

if [[ ! -d ${SWIG_DIR} ]]; then
    if [[ ! -f ${SWIG_ARCHIVE} ]]
        then wget ${SWIG_URL}
    fi

    tar -xzf ${SWIG_ARCHIVE}
fi

cd ${SWIG_DIR}
./configure \
  --prefix=${PREFIX} \
  --without-alllang \
  --with-python

make
make install

echo
echo "================"
echo "SWIG installation successful"
echo
echo "To use this version of SWIG, add directory ${SWIG_BIN_DIR} to your PATH,"
echo "e.g. adding the following line to your .bashrc:"
echo "    export PATH=${SWIG_BIN_DIR}:\$PATH"
echo "================"
