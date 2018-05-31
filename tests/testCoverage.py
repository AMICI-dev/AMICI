#!/usr/bin/env python3

import coverage
import unittest
import os
import sys
import amici

import testModels
import testSBML

cov = coverage.Coverage(source=['amici'],omit=['*/amici.py'])

cov.exclude('import')
cov.start()

print(testModels.TestAmiciPregeneratedModel())

suite = unittest.TestSuite()
suite.addTest(testModels.TestAmiciPregeneratedModel())
suite.addTest(testSBML.TestAmiciSBMLModel())
testRunner = unittest.TextTestRunner(verbosity=0)
testRunner.run(suite)

cov.stop()
cov.xml_report(outfile='coverage_py.xml')