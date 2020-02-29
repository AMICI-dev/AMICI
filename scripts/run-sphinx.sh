#!/bin/bash
# generate code documentation via sphinx and upload to rtd

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

cd ${AMICI_PATH}/documentation
source ../build/venv/bin/activate
pip3 install sphinx nbsphinx recommonmark sphinx_rtd_theme petab sphinx-autodoc-typehints

# add linebreaks that swig is missing, we can repeat this as many times as
# we want
AMICI_FILE=${AMICI_PATH}/python/sdist/amici/amici.py
sed -i -e -E $'s/([ ]+):(rtype|type|param|return)/\\\n\\1:\\2/g' $AMICI_FILE

sphinx-build -b html . _build

# cleanup
rm ${AMICI_FILE}
mv ${AMICI_FILE}-e ${AMICI_FILE}


