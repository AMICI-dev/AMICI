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
sed -i -e -E "s#std::vector< amici::realtype,std::allocator< amici::realtype > >#DoubleVector#g" $AMICI_FILE
sed -i -e -E "s#std::vector< double,std::allocator< double > >#DoubleVector#g" $AMICI_FILE
sed -i -e -E "s#std::vector< int,std::allocator< int > >#IntVector#g" $AMICI_FILE
sed -i -e -E "s#std::vector< amici::ParameterScaling,std::allocator< amici::ParameterScaling >#ScalingVector#g" $AMICI_FILE
sed -i -e -E "s#std::vector< std::string,std::allocator< std::string > >#StringVector#g" $AMICI_FILE
sed -i -e -E "s#std::vector< bool,std::allocator< bool > >#BoolVector#g" $AMICI_FILE
sed -i -e -E "s#std::map< std::string,amici::realtype,std::less< std::string >,std::allocator< std::pair< std::string const,amici::realtype > > >#StringDoubleMap#g" $AMICI_FILE
sed -i -e -E "s#std::vector< amici::ExpData *,std::allocator< amici::ExpData * > >#ExpDataVector#g" $AMICI_FILE
sed -i -e -E "s#std::vector< std::unique_ptr< amici::ReturnData >,std::allocator< std::unique_ptr< amici::ReturnData > > >#ReturnDataVector#g" $AMICI_FILE
sed -i -e -E "s#std::unique_ptr< amici::ExpData >#ExpData#g" $AMICI_FILE
sed -i -e -E "s#std::unique_ptr< amici::ReturnData >#ReturnData#g" $AMICI_FILE
sed -i -e -E "s#std::unique_ptr< amici::Solver >#Solver#g" $AMICI_FILE
sphinx-build -b html . _build


