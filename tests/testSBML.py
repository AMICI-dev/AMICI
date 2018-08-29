#!/usr/bin/env python3

import sys
import amici
import unittest
import os
import numpy as np
import pickle

class TestAmiciSBMLModel(unittest.TestCase):
    '''
    TestCase class for testing SBML import and simulation from AMICI python interface
    '''

    expectedResultsFile = os.path.join(os.path.dirname(__file__), 'cpputest','expectedResults.h5')

    def runTest(self):
        '''
        Test SBML import and simulation from AMICI python interface
        '''

        sbmlFile = os.path.join(os.path.dirname(__file__), '..', 'python', 'examples', 'example_steadystate', 'model_steadystate_scaled.xml')
        sbmlImporter = amici.SbmlImporter(sbmlFile)
        sbml = sbmlImporter.sbml

        observables = amici.assignmentRules2observables(
            sbml,
            filter_function=lambda variableId:
                variableId.startswith('observable_') and
                not variableId.endswith('_sigma')
        )
        
        sbmlImporter.sbml2amici('test_model_steadystate_scaled',
                                'test_model_steadystate_scaled',
                                observables=observables,
                                constantParameters=['k0'],
                                sigmas={'observable_x1withsigma': 'observable_x1withsigma_sigma'})

        sys.path.insert(0, 'test_model_steadystate_scaled')
        import test_model_steadystate_scaled as modelModule

        model = modelModule.getModel()
        model.setTimepoints(amici.DoubleVector(np.linspace(0, 60, 60))) 
        solver = model.getSolver()
        rdata = amici.runAmiciSimulation(model, solver)
        edata = [amici.ExpData(rdata, 0.01, 0)]
        rdata = amici.runAmiciSimulations(model, solver, edata)

        df_res = amici.getResidualsAsDataFrame(model, edata, rdata)

        model.getParameterById('p1')
        model.setParameterById('p1',2.0)
        model.getParameterByName('p1')
        model.setParameterByName('p1', 2.0)
        model.getFixedParameterById('k0')
        model.setFixedParameterById('k0', 2.0)
        model.getFixedParameterByName('k0')
        model.setFixedParameterByName('k0', 2.0)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestAmiciSBMLModel())
    unittest.main()
    