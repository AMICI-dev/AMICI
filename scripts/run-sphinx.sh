#!/bin/bash
# generate code documentation via sphinx and upload to rtd

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

cd ${AMICI_PATH}/documentation/
source ../build/venv/bin/activate
pip3 install sphinx nbsphinx recommonmark sphinx_rtd_theme petab sphinx-autodoc-typehints
make clean
make html

