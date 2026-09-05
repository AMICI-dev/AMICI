#!/bin/bash
#
# Build BNGL (required for pysb)
#
set -euo pipefail

BNG_VERSION="2.9.2"

script_path=$(dirname "$BASH_SOURCE")
amici_path=$(cd "$script_path/.." && pwd)

mkdir -p "${amici_path}/ThirdParty"
cd "${amici_path}/ThirdParty"

if [ ! -d "BioNetGen-${BNG_VERSION}" ]; then
    if [ ! -e "bionetgen.tar.gz" ]; then
        if [[ "$OSTYPE" == "linux-gnu" || "$OSTYPE" == "linux" ]]; then
            os="linux"
        elif [[ "$OSTYPE" == "darwin"* ]]; then
            os="mac"
        else
            echo "Unsupported OSTYPE for BioNetGen download: ${OSTYPE}" >&2
            exit 1
        fi
        wget -q -O bionetgen.tar.gz \
            "https://github.com/RuleWorld/bionetgen/releases/download/BioNetGen-${BNG_VERSION}/BioNetGen-${BNG_VERSION}-${os}.tar.gz"
    fi
    tar -xf bionetgen.tar.gz

    if [ ! -d "BioNetGen-${BNG_VERSION}" ]; then
        echo "Error: expected directory BioNetGen-${BNG_VERSION} not found after extracting bionetgen.tar.gz" >&2
        exit 1
    fi
fi
