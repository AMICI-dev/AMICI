AMICI_PATH="`dirname \"$0\"`"
AMICI_PATH="`( cd \"$AMICI_PATH/..\" && pwd )`"

# Prepare tests

TESTMODELS="model_dirac model_steadystate model_jakstat_adjoint model_jakstat_adjoint_o2 model_neuron model_neuron_o2"
for MODEL in $TESTMODELS; do 
	mkdir -p ${AMICI_PATH}/models/${MODEL}/build
	cd ${AMICI_PATH}/models/${MODEL}/build
	cmake -DCMAKE_CXX_STANDARD=11 -DCMAKE_CXX_STANDARD_REQUIRED=ON -DCMAKE_BUILD_TYPE=Debug ..
# make clean
	make
done;


# Build test suite

cd ${AMICI_PATH}/tests/cpputest/
mkdir -p build
cd build
cmake -DCMAKE_CXX_STANDARD=11 -DCMAKE_CXX_STANDARD_REQUIRED=ON -DCMAKE_BUILD_TYPE=Debug ..
# make clean
make

export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH}:${SUNDIALS_BUILD_PATH}/lib:${SUITESPARSE_ROOT}/lib"
