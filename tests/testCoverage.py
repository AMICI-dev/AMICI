#!/usr/bin/env python3

"""
Generate coverage reports for the testModels and testSBML scripts
exported format is cobertura xml
"""

import coverage
# only consider amici module and ignore the swig generated amici.py
cov = coverage.Coverage(
    source=['amici'],
    omit=['*/amici.py','*/amici_without_hdf5.py']
)

# ignore code blocks containing import statements
cov.exclude('import')
cov.set_option('report:exclude_lines', ['pragma: no cover', ' raise ',
                                        'except:'])

# coverage needs to be started before all other imports such that everything
# that happens during module import is also added to the report
cov.start()

import unittest
import sys

import testModels
import testSBML
import testPandas
import testPYSB
import testCPP
import testMisc
import testPreequilibration

test_modules = (
    testSBML, testModels, testPandas, testPYSB, testCPP, testMisc,
    testPreequilibration,
)

suite = unittest.TestSuite()
loader = unittest.TestLoader()
for test_module in test_modules:
    tests = loader.loadTestsFromModule(test_module)
    suite.addTests(tests)

testRunner = unittest.TextTestRunner(verbosity=0)
result = testRunner.run(suite)

cov.stop()
cov.xml_report(outfile='coverage_py.xml')

# propagate failure
sys.exit(not result.wasSuccessful())
