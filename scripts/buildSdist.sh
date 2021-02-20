#!/bin/bash
#
# Create AMICI sdist for PyPI deployment
#
set -e

script_path=$(dirname "$BASH_SOURCE")
amici_path=$(cd "$script_path/.." && pwd)

pip3 install build

cd "${amici_path}/python/sdist"
python -m build --sdist
