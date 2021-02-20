#!/bin/bash
#
# Build BNGL (required for pysb)
#
set -e

script_path=$(dirname "$BASH_SOURCE")
amici_path=$(cd "$script_path/.." && pwd)

mkdir -p "${amici_path}/ThirdParty"
cd "${amici_path}/ThirdParty"

if [ ! -d "BioNetGen-2.5.2" ]; then
    if [ ! -e "bionetgen.tar.gz" ]; then
        if [[ "$OSTYPE" == "linux-gnu" || "$OSTYPE" == "linux" ]]; then
            wget -q -O bionetgen.tar.gz https://github.com/RuleWorld/bionetgen/releases/download/BioNetGen-2.5.2/BioNetGen-2.5.2-linux.tgz
        elif [[ "$OSTYPE" == "darwin"* ]]; then
            wget -q -O bionetgen.tar.gz https://github.com/RuleWorld/bionetgen/releases/download/BioNetGen-2.5.2/BioNetGen-2.5.2-mac.tgz
        fi
    fi
    tar -xf bionetgen.tar.gz
fi
