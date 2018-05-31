#!/usr/bin/env python3

import sys
import h5py
import amici
import unittest
import importlib
import os
import re
import numpy as np

class TestAmiciSBMLModel(unittest.TestCase):
    '''
    TestCase class for tests that were pregenerated using the the matlab code generation routines and cmake
    build routines
    
    NOTE: requires having run scripts/buildTests.sh before to build the python modules for the test models
    '''

    expectedResultsFile = os.path.join(os.path.dirname(__file__), 'cpputest','expectedResults.h5')

    def runTest(self):
        '''
        test runner routine that loads data expectedResults.h5 hdf file and runs individual models/settings
        as subTests
        '''

        sbmlFile = os.path.join(os.path.dirname(__file__), '..', 'python', 'examples', 'example_steadystate', 'model_steadystate_scaled.sbml')
        sbmlImporter = amici.SbmlImporter(sbmlFile)
        sbml = sbmlImporter.sbml

        observables = amici.assignmentRules2observables(sbml, filter=lambda variableId: 
                                                        variableId.startswith('observable_') and not variableId.endswith('_sigma'))
        
        sbmlImporter.sbml2amici('test', 'test', 
                                observables=observables,
                                constantParameters=['k4'],
                                sigmas={'observable_x1withsigma': 'observable_x1withsigma_sigma'})

        sys.path.insert(0, 'test')
        import test as modelModule

        model = modelModule.getModel()
        model.setTimepoints(amici.DoubleVector(np.linspace(0, 60, 60))) 
        solver = model.getSolver()
        rdata = amici.runAmiciSimulation(model, solver)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestAmiciPregeneratedModel())
    unittest.main()
    