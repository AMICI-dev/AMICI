# Continuous integration (CI) and tests

AMICI uses a continuous integration pipeline running on https://travis-ci.org/.
This includes the following steps:

- Checking existence and format of documentation
- Static code analysis (http://cppcheck.sourceforge.net/)
- Unit and integration tests
- Memory leak detection

More details are provided in the sections below.

The CI scripts and tests can be found in `tests/` and `scripts/`. Some of the
tests are integrated with CMake, see `make help` in the build directory. 


## C++ unit and integration tests

To run C++ tests, build AMICI with `make` or `scripts/buildAll.sh`,
then run `scripts/run-cpputest.sh`.


## Python unit and integration tests

To run Python tests, run `../scripts/run-python-tests.sh` from anywhere
(assumes build directory is `build/`) or run `make python-tests` in your build
directory.


## Matlab tests (not included in CI pipeline)

To execute the Matlab test suite, run `tests/testModels.m`.


## Model simulation integration tests

Many of our integration tests are model simulations. The simulation results
obtained from the Python and C++ are compared to results saved in an HDF5 file
(`tests/cpputest/expectedResults.h5`).
Settings and data for the test simulations are also specified in this file.

**Note:** The C++ code for the models is included in the repository under 
`models/`.
This code is to be updated whenever `amici::Model` changes.


### Regenerating C++ code of the test models

Regeneration of the model code must done whenever `amici::Model` or
the Matlab model import routines change.

This is done with
  
    tests/cpputest/wrapTestModels.m

**Note:** This is currently only possible from Matlab < R2018a. This should
change as soon as 1) all second-order sensitivity code is ported to C++/Python,
2) a non-SBML import exists for Python and 3) support for events has been added
for Python.
    
    
### Regenerating expected results

To update test results, run `make test` in the build directory,
replace `tests/cpputest/expectedResults.h5` by 
`tests/cpputest/writeResults.h5.bak` 
[ONLY DO THIS AFTER TRIPLE CHECKING CORRECTNESS OF RESULTS]
Before replacing the test results, confirm that only expected datasets have
changed, e.g. using 

    h5diff -v -r 1e-8 tests/cpputest/expectedResults.h5 tests/cpputest/writeResults.h5.bak | less


## Adding/Updating tests

To add new tests add a new corresponding python script (see, e.g.,
`./tests/generateTestConfig/example_dirac.py`) and add it to and run 
`tests/generateTestConfigurationForExamples.sh`.
Then regenerate the expected test results (see above).
