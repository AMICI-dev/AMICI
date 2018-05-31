#!/bin/bash
# Check code coverage via codecov

AMICI_PATH="`dirname \"$0\"`"
AMICI_PATH="`( cd \"$AMICI_PATH/..\" && pwd )`"

lcov --base-directory ${AMICI_PATH} --directory ${AMICI_PATH} --zerocounters -q

cd ${AMICI_PATH}/build
ctest -V
rm ${AMICI_PATH}/tests/cpputest/writeResults.h5
cd ${AMICI_PATH}

lcov --compat-libtool --no-external --directory ${AMICI_PATH}/build/CMakeFiles/amici.dir/src  --capture --output-file coverage.info

wget http://raw.githubusercontent.com/eriwen/lcov-to-cobertura-xml/master/lcov_cobertura/lcov_cobertura.py

python3 ./lcov_cobertura.py coverage.info --output coverage_cpp.xml --demangle

python3 ./tests/testCoverage.py

wget https://gist.githubusercontent.com/tgsoverly/ef975d5b430fbce1eb33/raw/a4836655814bf09ac34bd42a6dd99f37aea7265d/merge-xml-coverage.py

python ./merge-xml-coverage.py coverage_py.xml coverage_cpp.xml

# cleanup
rm lcov_cobertura.py
rm merge-xml-coverage.py
rm -rf ./test
rm coverage_py.xml
rm coverage_cpp.xml
rm coverage.info

bash <(curl -s https://codecov.io/bash) -X fix -f coverage-merged.xml  || echo 'Codecov failed to upload'

