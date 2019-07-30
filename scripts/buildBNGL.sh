#!/bin/bash
#
# Build BNGL (required for pysb)
#
set -e

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

mkdir -p ${AMICI_PATH}/ThirdParty
cd ${AMICI_PATH}/ThirdParty

if [ ! -d "BioNetGen-2.3.2" ]; then
    if [ ! -e "bionetgen.tar.gz" ]; then
        if [[ "$OSTYPE" == "linux-gnu" || "$OSTYPE" == "linux" ]]; then
            wget -q -O bionetgen.tar.gz https://bintray.com/jczech/bionetgen/download_file?file_path=BioNetGen-2.3.2-linux.tar.gz
        elif [[ "$OSTYPE" == "darwin"* ]]; then
            wget -q -O bionetgen.tar.gz https://bintray.com/jczech/bionetgen/download_file?file_path=BioNetGen-2.3.2-osx.tar.gz
        fi
    fi
    tar -xzf bionetgen.tar.gz
fi
