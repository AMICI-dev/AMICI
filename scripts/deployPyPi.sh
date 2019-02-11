#!/bin/bash
#
# Deploy amici to pypi
#
set -e

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

pip3 install twine
# in case we are running with pyenv, we need to update pyenv shims after installing packages with binaries
if [[ -n "${PYENV_VERSION}" ]]; then
    pyenv rehash
fi

# authentication via env variables set in travis
twine upload $(ls -t ${AMICI_PATH}/build/python/amici-*.tar.gz | head -1)
