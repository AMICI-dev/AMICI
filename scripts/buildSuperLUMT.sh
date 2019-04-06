#!/usr/bin/env bash
#
# Build SuperLUMT
#
set -e

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd ${SCRIPT_PATH}/.. && pwd)

cd ${AMICI_PATH}/ThirdParty

if [[ ! -d SuperLU_MT_3.1 ]]; then
    if [[ ! -f superlu_mt_3.1.tar.gz ]]; then
        wget https://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_mt_3.1.tar.gz
    fi
    tar -xzf superlu_mt_3.1.tar.gz SuperLU_MT_3.1/
fi

cd SuperLU_MT_3.1/
cp MAKE_INC/make.pthread make.inc
# Add -fPIC
sed -ri 's/(CFLAGS\W*=\W*.*)(\#.*)/\1-fPIC \2/' make.inc
make superlulib
