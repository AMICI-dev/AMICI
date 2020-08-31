#!/bin/bash
# Download and build mtocpp (Doxygen filter for Matlab)

set -euo pipefail
set -x

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd "$SCRIPT_PATH/.." && pwd)

MTOC_CONFIG_PATH=${AMICI_PATH}/matlab/mtoc/config

cd "${AMICI_PATH}/ThirdParty"

# download mtocpp?
if [ ! -d "mtocpp-master" ]; then
    # mtocpp requires ragel
    if ! command -v ragel &> /dev/null; then
      echo "ragel not found"

      if [ ! -d "ragel-6.10" ]; then
        if [ ! -e "ragel-6.10.tar.gz" ]; then
            echo "Downloading ragel ..."
            wget -O ragel-6.10.tar.gz http://www.colm.net/files/ragel/ragel-6.10.tar.gz
        fi
        # build ragel?
        tar -xzf ragel-6.10.tar.gz
        (cd ragel-6.10 && ./configure && make && make install)
      fi
    fi

    if [ ! -e "mtocpp-master.zip" ]; then
        echo "Downloading mtocpp ..."
        wget -O mtocpp-master.zip https://github.com/mdrohmann/mtocpp/archive/master.zip
    fi
    # build mtocpp?
    unzip mtocpp-master.zip
    mkdir ./mtocpp-master/build

    echo "Building mtocpp ..."
    cd ./mtocpp-master/build && cmake .. && make mtocpp mtocpp_post
    if [ $? -ne 0 ] ; then
        exit 1
    fi
fi

# generate filter
echo "$AMICI_PATH/ThirdParty/mtocpp-master/build/mtocpp \$1 ${MTOC_CONFIG_PATH}/mtocpp.conf" > "${MTOC_CONFIG_PATH}/mtocpp_filter.sh"
chmod +x "${MTOC_CONFIG_PATH}/mtocpp_filter.sh"
