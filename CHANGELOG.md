# Changelog 

## v0.X Series

### v0.11.21 (2021-11-21)

Fixes:
 * Fixed a bug in recursion depth computation for model expressions. This may
   have resulted in incorrect sensitivities for models with expressions nested
   more than 2 levels. (#1595)
 * Fixed improper handling of Piecewise functions in PySB import which may have
   produced incorrect simulation results. (#1594)
 * Fixed changed googletest reference which broke the CMake-based build if
   tests were enabled (#1592)

New:
 * It's now possible to build AMICI using Ninja (#1593)


### v0.11.20 (2021-11-12)

New: 
 * Changed parameter mappings such that unassigned values have non-nan default values. This fixes erroneous evaluation of `llh` as `NaN` in some situations (#1574)
 * Added support for Python 3.10 (#1555)

Fixes:
 * Fixed a bug when simulation start time was not transferred when copying a solver instance (#1573)
 * Fixed a bug which led to incorrect sensitivies for models with multiple assignment rules or rate rules (#1584)

Other:
 * Update CI and documentation settings (#1569, #1527, #1572, #1575, #1579, #1580, #1589, #1581)
 * Extend set of validated benchmark models that is checked during CI (#1571, #1577)
 * Fixed string formatting in derivative checks (#1585)
 * Added helper methods to save and restore model instance-only settings (#1576)


### v0.11.19 (2021-10-13)

New:
* Added support for observable transformations (lin/log/log10) (#1567). Thereby supporting additional noise distributions in combination with least squares solvers.

Fixes:
* Fixed a bug when Newton sensitivity computation was activated despite specifying newton_steps == 0. The error occurs when simulation converges to a steadystate but simulation sensitivities are not converged according to convergence criteria. In that case simulation returned failure, but the newton rootfinding "finds" a steadystate even before the iteration check, leading to the erroneous computation of sensitivities via Newton/IFT. For singular jacobians this means the overall simulation still fails, but a different, more informative error message is displayed. (#1541)
* Fixed a bug where argument "outdir" in ODEExporter.__init__ would not be used (#1543)

Other:
* Improve checking support for SBML extensions (#1546)
* SBML import: Use more descriptive IDs for flux expressions (#1551)
* Optimized SUNMatrixWrapper functions (#1538)
* C++: Changed test suite from CppUTest to gtest (#1532)
* Add CITATION.cff (#1559)
* Updated documentation (#1563, #1554, #1536)
* Removed distutils dependency (#1557)
* Require sympy<1.9


### v0.11.18 (2021-07-12)

New:
* Allow specifying maximum integration time via
  `amici::Solver::setMaxTime()` (#1530)
* Py: Add `failfast` and `num_threads` argument to `simulate_petab`
  (#1528, #1524)
* Enable typehints / static type checking for AMICI-generated model modules
  (#1514) (`amici.ModelModule`, available with Python>=3.8)

Fixes:
* Python: Fix unused argument `generate_sensitivity_code` in `pysb2amici`
  (#1526)
* Python: Fix C(++) stdout redirection which could have let to deadlocks in
  exotic situations (#1529)
* Py: Fixed deprecation warning (#1525)
* Py: Proper typehints for SWIG interface (#1523), enabling better static type
  checking and IDE integration (available with Python>=3.9)
* C++: Fixed clang compiler warning (#1521)
* C++: Fix inherited variadic ctor in exception class (#1520)
* PEtab: More informative output for unhandled species overrides (#1519)
* Return SbmlImporter from PEtab model import (#1517)


### v0.11.17 (2021-05-30)

Fixes:
* Fix "maybe-uninitialized" compiler warning (#1495)
* Fix substitution of expressions in `drootdt_total (#1512)
* C++: Fix serialization and == operator (#1493)
* C++: Avoid `w` in `root` and `stau` headers (refactor) (#1503)

Documentation:
* Updated OpenBLAS Windows installation instructions (#1496)
* Updated how-to-cite to Bioinformatics paper (#1499)
* Updated list of papers using AMICI (#1509)

Other:
* Remove sllh computation from `petab_objective.simulate_petab` (#1498)
* Add __main__.py to Python package to provide info on AMICI installation
  via `python -m amici` (#1500)

### v0.11.16 (2021-04-13)

Fixes:
* Fixed serialization bug  (#1490)

New:
* Construction of condition specific plist for parameter mappings (#1487, #1488)
* **Add support for error residuals** (#1489)

### v0.11.15 (2021-03-31)

Fixes:
* Fixed initial state sensitivities in adjoint preequilibration (#1468)
* Fixed various model import / parameter mapping issues (#1469, #1473, #1475)

New:
* New AMICI releases will automatically trigger releases at
  https://biosimulators.org/simulators/amici/latest
* Transparent logo

### v0.11.14 (2021-03-16)

New features:
* **Python import now supports SBML Events** (#1443)
* Implement support for compilation without sensitivities (#1454)

Fixes:
  * Issue #1446: Check whether constant parameters are valid targets (#1450)
  * Issue #1422: Fix Steadystate solver failing if preequilibration starts in
    steadystate (#1461)
  * Issue #1401: Ensure diagnostics variables in ReturnData are always of
    expected length (#1438, #1447)
  * Fix FIM approximation for parameter dependent sigma (#1441)
  * Fix invalid SBML in PEtab/PySB import (#1433)
  * Fix: No context for inspect.getouterframes (#1432)

Documentation:
* Added this CHANGELOG
* Added feature request issue template (#1437)
* Updated reference list (#1430)
* Overhauled experimental conditions notebook (#1460)

CI:
* Test Matlab interface on GHA  (#1451)
* Include componentTags in SBML test suite output (#1462)
* Split SBML semantic test suite into multiple jobs (#1452)
* Fix Crauste ref val, fixes #1458 (#1459)

Misc:
* Various cleanup (#1465, #1448, #1455)
* Micro-optimize SUNMatrixWrapper::transpose (#1439)
* Remove constant triggers from roots in Heaviside (#1440)

### v0.11.13 (2021-02-20)

Breaking changes:
* AMICI requires Python>=3.7
* Updated package installation (PEP517/518): 
  Creating source distributions requires https://github.com/pypa/build (#1384)
  (but now handles all package building dependencies properly)

Features:
* More flexible state reinitialization (#1417)

Updated dependencies:
* Upgrade to sundials 5.7.0 (#1392)

Fixes:
* Python: account for heaviside functions in expressions (#1382)
* Python: allow loading of existing models in import_petab_problem (#1383)
* Python: Don't override user-provided compiler/linker flags (#1389)
* Python: PEtab import reinitialization fixes (#1417)
* Python: Fix PEtab observables for pysb models (#1390)
* Python: Substitute expressions in event condition expressions (#1404)
* Python: Unspecified initial states in PEtab conditions table default to SBML initial value (#1397)
* C++: Fix timepoint out of bounds access (#1402)
* C++: Fix exported CMake config (#1388)
* Fixed Dockerfile: add python3-venv (#1398, #1408)

Other:
* Slim exported swig interface (#1425)
* Updated documentation
    * Getting started tutorial (#1423)
    * List supported SBML test tags (#1428)
    * Add AMICI C++/Python/Matlab feature comparison (#1409)
    * ...
* Various minor CI improvements
* ...

### v0.11.12 (2021-01-26)

Features: 
* Add expression IDs and names to generated models (#1374)

Fixes:
* Raise minimum sympy version to 1.7.1 (Closes #1367)
* Fix species assignment rules in reactions (#1372)
* Fix id vector for DAEs (#1376)

Docs:
* Update how-to-cite (#1378)


### v0.11.11 (2020-12-15)

#### Python
* Restore support for species references (#1358)
* Add support for noise models in pysb (#1360)
* Proper handling of discontinuities in the ODE rhs (#1352)
* Fix directly calling AMICI from snakemake (#1348 )
* Extend mathml function support, particularly for numerical arguments (#1357)
* Restore support for sympy 1.6 (#1356)

#### C++
* Fix some compiler related warnings (#1349, #1362 )
* Fix a rare segfault for adjoint sensitivities (#1351)

#### CI
* Move windows tests to GHA (#1354)
* Pin breathe to 4.24.1

#### Docker
* Update ubuntu to 20.04

### v0.11.10 (2020-11-30)

Bugfix release that restores compatibility with sympy 1.7

### v0.11.9 (2020-11-29)

#### Python
* General improvements to SBML import (#1312, #1316, #1315, #1322 , #1324 #1333, #1329)
* Small bugfixes and improvements (#1321 )
* Improve derivative computation for instances of `power` (#1320 )

#### C++
* Fix FIM and residual computation for models with parameter dependent sigma. (#1338)
* Disable chi2/residual/FIM computation for non-gaussian objective functions. (#1338)
* Bugfix for integration failure during adjoint solve (#1327)

#### Doc
* Update references (#1331, #1336)

#### CI
* Update OpenBLAS for windows (#1334)

### v0.11.8 (2020-10-21)

#### Python
* Fix pysb-petab support (#1288)
* Fix ExpData constructor overloading (#1299)
* Fix support for positivity enforcement (#1306)
* **Refactor SBML import, adds support for parameter rate rules and initial assignments** (#1284, #1296, #1304)
* Improve model generation for models with many parameters (#1300)
* Add support for PEtab based synthetic data generation (#1283)

#### C++
* Make HDF5 an optional dependency (#1285)

#### Doc
* General Improvements to Documentation (#1289, #1291, #1292, #1293, #1294, #1286, #1277, #1281)

#### CI
* Add python 3.9 support test (#1282)
* Allow manual triggering of GitHub actions (#1287)
* Remove appveyor config (#1295)
* Update GHA env and path management (#1302)

### v0.11.7 (2020-09-22)

#### Python
* Improve and extend available objective functions (#1235)
* Fix processing of compartment definitions (#1223)
* Fix replacement of reserved symbols (#1265)
* Use Hierarchical Derivatives for Expressions (#1224, #1246)
* Fix duplicate running of swig (#1216)
* Overload python interface functions for amici.{Model,Solver,ExpData} and amici.{Model,Solver,ExpData}Ptr (#1271)

#### C++
* Fix and extend use of sparse matrix operations (#1230, #1240, #1244, #1247, #1271) 
* **Fix application of maximal number of steps**, MaxNumStep parameter now limit total number of steps, not number of steps between output times. (#1267)

#### Doc
* Move all Documentation to RTD (#1229, #1241)
* General Improvements to Documentation (#1225, #1227, #1219, #1228, #1232, #1233, #1234, #1237,  #1238, #1239, #1243, #1253, #1255, #1262)

#### CI
* Better check for doc building (#1226)
* Add more gradient checks (#1213)
* Update GHA to Ubuntu 20.04 (#1268)

### v0.11.6 (2020-08-20)

#### Python
* Bugfix for piecewise functions (#1199)
* Refactor swigging - generate one single wrapper (#1213)

#### C++
* Fix warnings: account for zero indexing in nan/inf error (#1112)

#### Doc
* Update Windows build instructions (#1200, #1202)
* Update README: Projects using AMICI (#1209)
* Add CODE_OF_CONDUCT.md (#1210)
* Update documentation for Python interface (#1208)

#### CI
* Create sdist on GHA using swig4.0.1 (#1204)  (Fixing broken pypi package)
* Fix links after repository move
* Speed-up swig build: disable all languages except python (#1211)
* Fix doc generation on readthedocs (#1196) 


### v0.11.5 (2020-08-07)

#### General
* Move repo to new organization (#1193)
* Update Bibliography

#### Python
* Fix bug for energyPySB models (#1191)

#### CI
* Fix release deployment (#1189)

### v0.11.4 (2020-08-06)

#### Python
* Skip unnecessary expressions in pysb models (#1185)
* MSVC compiler support (this time for real... #1152)

#### CI
* Implement MSVC tests (#1152)
* Rename and group GitHub actions (#1186)
* Fix release deployment (#1186)

### v0.11.3 (2020-08-06)

#### Python
* Fix simplification for pysb models (#1168)
* Pass verbosity flags to pysb network generation (#1173)
* Enable experimental pysb-petab support (#1175)
* Add installation instructions for Fedora (#1177)
* Implement support for SBML rate-references (#1180)

#### C++
* Refactoring (#1162, #1163)

#### CI
* Move majority of tests to Github Actions (#1166, #1160)
* Improve reporting of skipped tests in SBML testsuite (#1183)

### v0.11.2 (2020-07-17)

#### Python
* Speed up model import, compilation (#1123, #1112)
* Improve/Add steady-state solver documentation (#1102)
* Improve extension import (#1141)
* Bugfixes SBML import (#1135, #1134, #1145, #1154)
* Fixed issue that prevented simplification (#1158)

#### C++
* Bugfixes (#1121, #1125, #1131, #1132, #1136)
* Enable openMP by default (#1118)
* Improve memoy footprint for simulations with replicates (#1153)
* Improve steady-state solver and add option to to adjoint-steadystate hybrid (#1143, #1099, #1129, #1146)

#### CI
* Store build artifacts from github actions (#1138)

### v0.11.1 (2020-06-05)

#### Python
* Upgrade to sympy 1.6.0, which is now required minimum version  (#1098, #1103)
* Speed up model import 
  * Speed-up computation of sx0, reduce file size (#1109)
  * Replace terribly slow sympy.MutableDenseMatrix.is_zero_matrix by custom implementation (#1104)
* speedup dataframe creation in `get*AsDataFrame` (#1088)
* Allow caching edatas for simulate_petab (#1106)
* Fix wrong deprecation warning (Fixes #1093)
* Fix segmentation faults in NewtonSolver under certain conditions (#1089, #1090, #1097)
* fix wrong power function call in `unscale_parameter` (#1094)
* Fix MathML conversion (#1086)
* Fix deepcopy of SymPy objects (#1091)

#### Matlab
* handle empty rdata->{pre|post}eq_numlinsteps (Closes #1113), which previously made the matlab interface unusable
* Fix generation of compileMexFile.m for matlab compilation of python code (#1115)

#### C++
* Reduce memory requirements and speedup compilation of large models (#1105)
* Place generated code into own namespace (#937) (#1112)
* Fix several msvc compiler warnings (#1116) (Note that MSVC support is still experimental) **breaking change for users of C++ interface**
* Fix swig warning: ensure base class ContextManager is known before use (Fixes #1092) (#1101)

#### CI
* Don't install/run valgrind on travis CI (done with github actions… (#1111)


### v0.11.0 (2020-05-10)

Python:

- **Implement support for variable compartments (#1036)**
- Better handling of constant species (#1047)
- **Better handling of C++ enums, this makes `amici.SensitivityMethod_forward` available as `amici.SensitivityMethod.forward` (#1042)**
- Improve installation routines (#1055, #1056, #1058, #1076)
- Add option to reduce memory usage (#1044)
- **Fix handling of symbolic expressions in nested rules (#1081, 1069)**

Library:

- Update Sundials to 5.2.0 (#1039)
- Update SuiteSparse to 5.4.0 (#1040)
- Refactor use of ReturnData, now completely created post-hoc (#1002)
- **Fix propagation of reinitialization in ExpData constructor (#1041)**
- **Fix issue where InternalSensitivityParameter was sometimes not set (#1075)**
- **Fix or disable certain combinations of equilibraition, presimulation and adjoint sensitivity analysis**

CI:

- Move from Codacy to Sonarcloud (#1065)
- Run SBML Testsuite when appropriate (#1058)

### v0.10.21 (2020-04-04)

Library:
* Fix: Handle paths with blanks in build scripts
* Feature: Add function to write amici::Solver settings to HDF5 (#1023)
* Fix: typehints (#1018, #1022)
* Refactor: Move creation of parameter mapping for objective<->simulation to classes (#1020)

CI:
* Refactor: Cleanup and reorganize tests (#1026)
* Fix: benchmark problem test should fail on missing files (Closes #1015)


### v0.10.20 (2020-03-18)

* Fixed (re)initialization of sensitivities if ExpData::fixedParametersPreequilibration is set (#994)
* Fixed sensitivities for parameters in sigma expressions for Python/SBML in case provided expression was not just a single parameter ID
* Enable parallel compilation of model files from Python (#997) based on AMICI_PARALLEL_COMPILE enviroment variable
* Fixed computation of log-likelihood for log10-normal distributed noise
* Added `reinitializeFixedParameterInitialStates` to ExpData (#1000) (**breaking change**: overrides settings in `amici::Model`)
* Python model import now verifies that chosen model name is a valid identifier (Closes #928)
* Made w available in ReturnData (Closes #990) (#992)
* Fixed setting of log level when passing boolean values to verbose (#991)
* Documentation now on ReadTheDocs https://amici.readthedocs.io/en/
* Use proper state/observable names in plotting functions (#979)
* PEtab support:
  * Adapt to most recent PEtab (0.1.5)
  * Extended support for import of PEtab models
  * Added support for computing cost function based on PEtab problem
  * Implemented handling of species in condition table
  * petab_import.import_model now provides reproducible parameter list  (Closes #976)
  * Fix python import error in import_petab_problem: Add absolute paths to python path, invalidate caches and reload (#970)
  * Added example notebook
* CI: PEtab test suite integrated in CI workflow
* Added AMICI dockerfile and image deployment to dockerhub  (#948)
* Removed mention of 'mex' in warning/error ids (#968)
* More informative errors on SWIG interface import failures (#959)

### v0.10.19 (2020-02-13)

Python:
* Fix logo display on pypi
* Fix deadlocks in multithreaded python environments when using openMP parallelization

Matlab:
* Fix compilation errors due to switch to C++14

### v0.10.18 (2020-02-11)

General:
* AMICI now comes with a logo
* implement getName function for models
* Updated documentation / examples

Python:
* Enable MSVC compilation of Python extensions (#847)
* Always recompile clibs and extensions (Closes #700)
* Extended PEtab support (Running
* enable multithreading in swig (#938)
* Fixes pysb (#902) (#907)

C++
* Build optimized AMICI and sundials by default (Closes #934)

Matlab:
* Fix(matlab) Compile CalcMD5 on demand (Fixes #914)
* Don't pass empty include strings to mex
* Fix Matlab compilation error if AMICI or model path contains blanks

CI:
* Running additional test models

... and various minor fixes/updates

### v0.10.17 (2020-01-15)

- **added python 3.8 support, dropped python 3.6 support** (#898) 
- Added logging functionality (#900)
- Fixes PySB import (#879, #902)
- Fixes symbolic processing (#899)
- Improved build scripts (#894, 
- Improved petab support (#886, #888, #891)
- CI related fixes (#865, #896)

### v0.10.16 (2019-12-11)

* **Sparsify dwdp to reduce computation time for adjoints (#858)**
* Fix(matlab) update example name example_dae_events->example_calvetti (Closes #866)
* Fix nullptr deferencing for simulations with events when no measurements are provided (Fixes #866)
* Fix accessing empty vector during adjoint state event update (Closes #866)
* Fix pysb_import (fixes #878)


### v0.10.15 (2019-12-03)

Bugfix release due to incorrect sensitivities w.r.t. sigmas introduced in 0.10.14.

No other changes.

### v0.10.14 (2019-12-02)

**NOTE: For Python-imported SBML-models this release may compute incorrect sensitivities w.r.t. sigma. Bug introduced in 0.10.14, fixed in 0.10.15.**

Python: 

* Don't require use of ModelPtr.get to call ExpData(Model)
* Fix import in generated model Python package
* Setup AMICI standalone scripts as setuptools entrypoints
* Simplify symbolic sensitivity expressions during Python SBML import
        Fixes Infs in the Jacobian when using Hill-functions with states of 0.0.
* Extended Newton solver #848  
    The changes that allow performing Newton tests from the paper:    
    G. T. Lines, Ł. Paszkowski, L. Schmiester, D. Weindl, P. Stapor, and J. Hasenauer. Efficient computation of steady states in large-scale ODE models of biochemical reaction networks. accepted for Proceedings of the 8th IFAC Conference on Foundations of Systems Biology in Engineering (FOSBE), Valencia, Spain, October 2019.
* Use SWIG>=4.0 on travis to include PyDoc in sdist / pypi package (#841)
* **Fix choice of likelihood formula; failed if observable names were not equal to observable IDs**
* Fix(sbml-import) Compartment IDs in right-hand side of Rules are not replaced and lead to undefined identifiers in c++ files
* Fix invalid logging level
* Speed up sympy simplification (#871)

C++:

* Performance: Avoid unnecessary repeated function calls for SUNMatrixWrapper dimensions
*   Add AmiciApplication class as context for handling so far global settings.
    This allows for example setting custom logging functions for concurrent AMICI
    runs, e.g. in multi-thread applications (Closes #576).

Misc:

* Setup performance test on github actions (#853)
* Update documentation and FAQ for CBLAS requirement and others
* Update reference list

### v0.10.13 (2019-10-09)

* BREAKING CHANGE: Renaming {get|set}tNewtonPreequilibration to {get|set}Preequilibration (Closes #720)
* Make wurlitzer non-optional requirement for AMICI python package (Fixes missing AMICI errors when running from jupyter notebooks)
* Compute initial state for Model::getInitialStates if not already set (Fixes #818)
* Make swig generate pydoc comments from doxygen comments #830 (Closes #745) to provide Python docstrings for C++ wrapper functions
* feature(cmake) Add option to disable compiler optimizations for wrapfunctions.cpp (Fixes #828) (#829)
* Change SBML test suite to pytest to allow for parallel test execution… (#824)
* Fix(cmake): -E option is not available in all sed versions. Neither is the equivalent -r. Use --regexp-extended instead (Closes #826)
* Refactor(python) Move PEtab import code from command line script… (#825)
* Fix(core) Fix regular expressions for intel compiler (Closes #754) (#822)
* Update workflow figure to include PySB (Closes #799)
* Fix compiler warnings

### v0.10.12 (2019-09-28)

* Fix handling of species specified in PEtab condition table (#813)
* Fix some Visual C++ issues, update cppcheck handling, cleanup (VisualC++ still not fully supported)
* Minor fixups (#801)
* Create SBML test suite result files for upload to http://sbml.org/Facilities/Database/ (#798)


### v0.10.11 (2019-08-31)

* Fixed setting initial conditions for preequilibration (#784) 
* Fixed species->parameter conversion during PEtab import (#782) 
* Set correct Matlab include directories in CMake (#793)
* Extended and updated documentation (#785, #787)
* Fix various SBML import issues
* Run SBML test suite using github actions instead of travisCI (#789)

### v0.10.10 (2019-08-07)

* Simplify/fix AMICI installation
   * If available use environment modules to detect dependencies
 
  * Add SWIG installation script

* Update list of publication
* Update documentation
    * Update doc for SWIG build and custom SWIG location.
    * Add AMICI interface overview / workflow figure and show in README
    * Document environment variables for model/core compilation (Closes #737)

* Added handling of abs function, since there seem to be problems with case sensitivity (#713) Closes #770

Detaills:
    * cmake: Use package_ROOT environment variables
    * fix(cmake) Fix finding version.txt
    * cmake: Auto-detect loaded MKL environment module
    * cmake: Use new FindPython3 modules where possible
    * fix(python) Restore python3.6 compatibility
    * Inside venv, use pip instead of pip3 which should point to the correct version
    * fix(python) Workaround for missing ensurepip during venv creation [ci skip]
    * feature(python) Use MKL from environment modules to provide cblas
    * fix(python) Fix define_macros not being passed to setuptools for Extension
    * fix(python) Fix define_macros not being passed to setuptools for clibs
    * Do not always add 'cblas' library since users may want to override that by a cblas-compatible library with a different name (closes #736)   
    * Update HDF5 path hints; use shared library if static is not available.
    * Check for HDF5_BASE from environment module
    * Fix system-dependent sundials library directory (Fixes #749) (#750)
    * Handle OSTYPE==linux in scripts/buildBNGL.sh (Fixes #751)
    * Add SWIG download and build script
    * Improve finding swig executable and allow user override via SWIG environment variable
    * Provide installation hints if no SWIG found (Closes #724)
    * Allow overriding cmake executable with environment variables in build scripts (Closes #738)
 

### v0.10.9 (2019-07-24)

Fixup for missing version bump in v0.10.8 release. No code changes compared to v0.10.8.

### v0.10.8 (2019-07-24)

Changes in this release:

All:
- Updated / extended documentation
- Fix reuse of  `Solver` instances (#541)

C++:
-  Check for correct AMICI version for model in CMake 
- Add reporting of computation times (#699)

Python:
- Fix manifest file (#698)
- Fix initial amounts/concentrations in SBML import

... and various other minor fixes/improvements

### v0.10.7 (2019-05-01)

Python
* fix unset noise distribution when automatically generating observables in case None are passed (#691)

### v0.10.6 (2019-04-19)

C++
- Add SuperLUMT support (#681)
- Sparsified dJydy (#686)
- Enabled support of impulse-free events for DAE models (#687) - thanks to Sebastien Sten for providing a testcase for this

Python
- Enabled support for piecewise functions in SBML import (#662)
- Fix numeric type when constructing ExpData from Dataframes (#690)
- Fix dynamic override in PETab


### v0.10.5 (2019-04-08)

Bugfix release

Doc
- Update documentation of Windows installation

C++
- Fix missing source files in CMakeLists.txt (#658)
- Set CMake policies to prevent warnings (Closes #676) (#677)
- Start using gsl::span instead of raw pointers (#393) (#678) 

Python
- PySB parsing fix (#669)
- Fix failure to propagate BLAS_LIBS contents (#665) 
- Require setuptools at setup (#673) 
- Updated PEtab import to allow for different noise models


### v0.10.4 (2019-03-21)

Features / improvements:

* Implement ReturnData and ExpData wrappers as more efficient views (#657)
* Add list of references using AMICI (#659)
* Custom llh (normal/laplace, lin/log/log10) (#656)

Bugfixes:

+ Speedup and fix travis build

### v0.10.3 (2019-03-13)

Features / improvements:

- adds the option for early termination on integration failures for runAmiciSimulations
- improve runtime of `SUNMatrixWrapper::mutliply`
- expose finite difference routines in public API
- enable parallel compilation of clib source files

Bugfixes:

- fixed symbolic processing for unreleased pysb features

### v0.10.2 (2019-03-07)

Features / improvements:

- extended `ExpData` interface to allow for condition specific parameters, parameter scales, parameter lists, initial conditions and initial condition sensitivities.

Bugfixes:

- fixed output values of `ReturnData::x_ss` and `ReturnData::sx_ss`

### v0.10.1 (2019-03-04)

* travis-ci.com migration
* fix problem with has{variable} functions
* allow to import sbml model from string, not only file

### v0.10.0 (2019-03-01)

Features / improvements:

- updated sundials to 4.1.0
- updated SuiteSparse to 5.4.0
- added generic implementations for symbolic expressions that were sparse matrix vector products

Bugfixes:

- fixed return value of `rz` when no data is provided.

### v0.9.5 (2019-02-26)

Features / improvements:

- allow python installations without compilation of c++ extension
- improvements to ExpData <-> pandas.DataFrame interface
- allow generation of matlab models from python
- implement CLI interface for PEtab
- improve computation time for conservation laws from pysb import

Bugfixes:

- Fix sign in undamped Newton step.

Maintenance:

- use newer CI images 

### v0.9.4 (2019-02-11)

Minor fixes only:
- fix(core) Get solver diagnostics for first(last) timepoint (#588) (Closes #586)
- fix(ci) Fix autodeploy (Closes #589)

### v0.9.3 (2019-02-07)

**CRITICAL FIXES**
- **fix(python) fix symbolic computations for adjoint (#583)**

**Features**
- feature(python) Check for matching AMICI versions when importing model (Closes #556). Set exact AMICI version as model package requirement.
- feature(core) Add option to rethrow AMICI exception (Closes #552)
- feature(python) Redirect C/C++ output in stdout is redirected (e.g. in ipython notebooks) (Closes #456)

**Minor fixes**
- fix(python) Fix doc and rename sys_pipes to something more meaningful
- fix(ci) Fix premature exit of scripts/runNotebook.sh
- fix(deploy) Update pyenv shims to find twine



### v0.9.2 (2019-01-30)

Bugfixes:

- fixes a critical bug in the newton solver
- fixes multiple bugs in sbml import for degenerate models, empty stoichiometry assignments and conversion factors
- improved error messages for sbml import
- #560 
- #557 
- #559 


### v0.9.1 (2019-01-21)

Features / improvements:

- make pure steadystate results available as `rdata['x_ss']` and `rdata['sx_ss'] `
- add option to specify integration tolerances for the adjoint problem via `atolB` and `rtolB`

Bugfixes:

- improved conservation law identification to also account for constant species.
- fixed a bug where simulation results were written into results of the second newton solver attempt
- fixed an openMP related warning

Maintenance:

- attempt to fix automatic deploy to pypi via travis

### v0.9.0 (2019-01-18)

Features / improvements:

- Allow computation and application of conservation laws for pysb importet models. This enables use of NewtonSolver for preequilibration for models where it was previously not possible.
- Use `klu_refactor` in the sparse Newton solver instead of always using `klu_factor` and only perform symbolic factorization once (#421)
- Allow more detailed finiteness checks (#514)

Bugfixes:
 - #491 

Maintenance:
- Several improvements to travis log sizes and folding
- Use default copy constructor for Model by implementing class wrappers for sundials matrix types (#516)
- Reenable run of SBML testsuite


### v0.8.2 (2019-01-07)

Features / improvements:
* Speedup symbolic processing for ODE generation in python

Bugfixes:
* Fix(python) Add missing deepcopy (introduced in 6847ba675f2088854db583199b8754aaa6e01576)
* Fix(core) Prevent parameter scaling length mismatch (#488)
* Fix(python) Set distutils dependency to current version to fix </usr/lib/python3.6/distutils/dist.py:261: UserWarning: Unknown distribution option: 'long_description_content_type'>
* fix(python) add symlink to version.txt to be included in package

Backwards-compatibility:
* Replace 'newline' by literal to restore <R2016b compatibility (Fixes #493)

Maintenance:
* Remove obsolete swig library build via cmake and related file copying
* Provide issue template for bug reports
* Providing valid SBML document to import is not optional anymore
* Update documentation and tests
* Add python version check and raise required version to 3.6 to prevent cryptic error messages when encountering f-strings

### v0.8.1 (2018-11-25)

- [all] **critical** Fix long standing bugs in solving steadystate problems (including preequilibration) (#471)
- [all] Fix AmiVector constructor from std::vector (#471)
- [python] Reenable Solver and Model copy constructors
- Update documentation


### v0.8.0 (2018-11-25)

- replaced symengine by sympy for symbolic processing in python *which fixes several critical bugs* that were due to bugs in symengine (#467)


### v0.7.13 (2018-11-18)

- fixes a critical bug in objective function computation for models compiled using `sbml2amici` and `pysb2amici` that was introduced in v0.7.12
- fixes a critical bug in sensitivity computation when`model.reinitializeFixedParameterInitialStates` was set to true
- readds the python interface to the ExpData copy constructor that was inadvertently removed in 0.7.12 and streamlines the respective convenience wrapper to provide access to the full range of constructors.

### v0.7.12 (2018-11-17)

- fixed a critical bug in `amici.constructEdataFromDataFrame`
- enabled multithreaded simulation of multiple experiments (requires compilation with openMP)
- modularized sbml import and added pysb import

### v0.7.11 (2018-10-15)

- [python] Added numpy and python wrappers that provide a more user friendly python API
- [python] Enable import of SBML models with non-float assignment rules 
- [python] Enable handling of exceptions in python
- [python] Enable nativ python access to std::vector data-structures
- [core] Provide an API for more fine-grained control over sensitivity tolerances and steady-state tolerances
- [core] Provide an API to specify non-negative state variables (this feature is still preliminary)

### v0.7.10 (2018-08-29)

- Fixed python SBML import `log()` issues (#412)

### v0.7.9 (2018-08-24)

- fixes MATLAB compilation of models
- adds option to perform steady state sensitivity analysis via FSA
- condition dependent intitial conditions are now newly set after preequilibration is done

### v0.7.8 (2018-08-19)

- bugfixes for the ExpData interface
- created build configuration that enables debugging of c++ extensions on os x
- fixed python sbml import when stoichiometry is empty

### v0.7.7 (2018-08-17)

Fixes a couple of bugs just introduced in v0.7.6

### v0.7.6 (2018-08-13)

Important: **Use AMICI v0.7.7 due to https://github.com/ICB-DCM/AMICI/pull/403/commits/3a495d3db2fdbba70c2b0d52a3d4655c33c817a2**

Bug fixes:
- Python import: Fix log10 issues in observables (#382)
- Matlab: Fix broken model compilation (#392)
- Fixed simulation for models without observables (#390)
- Fixed potential matlab memory leaks (#392)

Breaking C++ API changes:
- Revised ExpData interface (#388)

### v0.7.5 (2018-07-30)

Features/enhancements:
- Add computation of residuals, residuals sensitivity, Fisher information matrix (#223)
- More efficient conversion of std::vector to numpy ndarray (#375)
- Allow specifying timepoints in ExpData (#370)

Minor fixes:
- Condition parameters in ExpData now only temporarily override Model parameters (#371)
- Ensure non-negative states for Newton solver (#378)

### v0.7.4 (2018-07-27)

Features/enhancements:
- Check SBML model validity (#343)
- Allow per-parameter setting of amioptions::pscale from matlab interface (#350)
- Documentation

Major fixes:
- Don't compile main.cpp into python model module which broke modules if amici was compiled without libhdf5 (#363)

Minor fixes:
- Fix compiler warnings (#353)
- Plotting, SBML example mode, ...

### v0.7.3 (2018-07-13)

Features:
- Added symbol names to python-wrapped models and make available via Model.getParameterNames(), model.getStateNames(), ...
- Extended Python interface example

Python package available via pypi: https://pypi.org/project/amici/0.7.3/

### v0.7.2 (2018-07-03)

Features:
- Python package: more flexible HDF5 library localization
- Extended CI: python tests, preequilibration tests, run in venv

Major bugfixes:
- Fix python sbml model import / compilation error (undefined function)
- Fix model preequilibration 

Minor fixes:
- Various fixes for mingw compilation of python source distribution
- Cmake compatibility with < 3.8 restored

### v0.7.1 (2018-06-12)

Features:
- Allow specifying sigma-parameters from Python interface

Major bugfixes:
- Fix dsigma_y/dp and downstream sensitivity errors

### v0.7.0 (2018-06-09)

- Major revision of documentation
- Improved Python interface
- More comprehensive Python interface example
- Fixed sensitivity computation in Python-generated models
- Various other bug fixes

WARNING:
- For models with sigma-parameters and dsigma_y/dp != 0, dsigma_y/dp was computed incorrectly. This propagates to all dependent sensitivities. This applies also to some older releases and has been fixed in v0.7.1.

### v0.6.0 (2018-05-22)

Implement experimental support for python via swig.
Python interface is now usable, but API will still receive some updates in the future.

WARNING: 
- There is a bug in sensitivity computation for  Python-generated models 
- Matlab C++ compilation will fail due to undefined M_PI
-> Please use v0.7.0

### v0.5.0 (2018-03-15)

Main new features are:

- Reimplemented support for DAE equations
- Added newton solver for steady state calculation and preequilibration
- Better caching for recompilation of models
- Blas support to  allow compilation of large models with many observables
- Improved SBML support
- Added c++ interface
- Added support for second order adjoint computation

- Rewrote large parts of the code as proper c++11 code to allow easier code maintanance
- Substantially extended testing in continuous integration to assure code quality

### v0.4.0 (2017-05-15)

* First citable version of AMICI (via zenodo integration).
* Better support for standalone compilation
* Improved SBML import scripts
* General Code Cleanup

### v0.3.0 (2016-09-05)

This update comes with many improvements, bug fixes and several new features. Most notably:

1) AMICI should now run on older versions of MATLAB (back to R2014a)
2) AMICI now compiles using older versions of Visual Studio
3) AMICI now also supports second order adjoint sensitivities (full (via the o2flag = 1 and as a vector product via o2flag = 2)
4) AMICI now supports more SBML, SBML v2 and rate rules


### 0.2.1 (2016-05-09)

Bugfix release. This release also includes some changes that should improve the performance on the new R2016a release of MATLAB.


### v0.2 (2016-03-17)

This update comes with many improvements to the computation time for both compilation and simulation. Moreover several new features were included:

1) Hessian Vector products for second order forward sensitivities
2) Correct treatment of parameter/state dependent discontinuities for both forward and adjoint sensitivities


### v0.1 (2015-11-05)

This is the initial release of the toolbox
