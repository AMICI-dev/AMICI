#!/bin/bash
#
# Build BNGL (required for pysb)
#
set -e

script_path=$(dirname "$BASH_SOURCE")
amici_path=$(cd "$script_path/.." && pwd)

mkdir -p "${amici_path}/ThirdParty"
cd "${amici_path}/ThirdParty"

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
