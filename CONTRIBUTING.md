# How to contribute

We are happy about contributions to AMICI in any form (new functionality, documentation, bug reports, ...).

## Making code changes

When making code changes:

* Check if you agree to release your contribution under the conditions provided in `LICENSE`
* Start a new branch from `develop`
* Implement your changes
* Submit a pull request to the `develop` branch
* Make sure your code is documented appropriately
  * Run `mtoc/makeDocumentation.m` to check completeness of your documentation
* Make sure your code is compatible with C++11, `gcc` and `clang`
* when adding new functionality, please also provide test cases (see `tests/cpputest/`)
* Write meaningful commit messages
* Run all tests to ensure nothing was broken
  * Run `tests/cpputest/wrapTestModels.m` followed by CI tests `scripts/buildAll.sh && scripts/run-cpputest.sh`
  * Run `tests/testModels.m`
  * Run `make python-tests` in `build`
* When all tests are passing and you think your code is ready to merge, request a code review

## Adding/Updating tests

To add new tests add a new corresponding python script (see, e.g.,  `tests/example_dirac.py`) and add it to and run `tests/generateTestConfigurationForExamples.sh`
To update test results replace  `tests/cpputest/expectedResults.h5` by `tests/cpputest/writeResults.h5.bak` [ONLY DO THIS AFTER TRIPLE CHECKING CORRECTNESS OF RESULTS]
Before replacing the test results, confirm that only expected datasets have changed, e.g. using 
`h5diff -v -r 1e-8 tests/cpputest/expectedResults.h5 tests/cpputest/writeResults.h5.bak | less`
