#!/bin/bash
#
# Deploy amici to pypi
#
set -e

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

pip3 install twine

# authentication via env variables set in travis
twine upload $(ls -t ${AMICI_PATH}/build/python/amici-*.tar.gz | head -1)
