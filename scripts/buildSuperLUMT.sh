#!/usr/bin/env bash
#
# Build SuperLUMT
#
set -e

script_path=$(dirname "$BASH_SOURCE")
amici_path=$(cd "${script_path}/.." && pwd)

cd "${amici_path}/ThirdParty"

if [[ ! -d "SuperLU_MT_3.1" ]]; then
    if [[ ! -f "superlu_mt_3.1.tar.gz" ]]; then
        wget "https://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_mt_3.1.tar.gz"
    fi
    tar -xzf "superlu_mt_3.1.tar.gz SuperLU_MT_3.1/"
fi

cd SuperLU_MT_3.1/
# interferes with std::queue
test -f SRC/queue && mv SRC/queue SRC/queue.bak
cp MAKE_INC/make.pthread make.inc
# Add -fPIC
sed -ri 's/(CFLAGS\W*=\W*.*)(\#.*)/\1-fPIC \2/' make.inc
# Use 64bit integers
sed -ri 's/# (CFLAGS\W*+=.*-D_LONGINT.*)/\1/' make.inc
sed -ri 's/# (FFLAGS\W*+=.*-fdefault-integer-8.*)/\1/' make.inc

make superlulib
