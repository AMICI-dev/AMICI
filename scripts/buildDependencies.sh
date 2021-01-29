#!/bin/bash
#
# Build AMICI along with dependencies and test suite
#
set -e

script_path=$(dirname "$BASH_SOURCE")
script_path=$(cd "$script_path" && pwd)

"${script_path}/buildSuiteSparse.sh"
"${script_path}/buildSundials.sh"
"${script_path}/buildCpputest.sh"
"${script_path}/buildBNGL.sh"
