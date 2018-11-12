#!/usr/bin/env python3

import sys
import amici
import unittest
import os
import pysb
import numpy as np

class TestAmiciPSBModel(unittest.TestCase):
    '''
    TestCase class for testing SBML import and simulation from AMICI python interface
    '''

    expectedResultsFile = os.path.join(os.path.dirname(__file__),
                                       'cpputest','expectedResults.h5')

    def setUp(self):
        self.resetdir = os.getcwd()
        if os.path.dirname(__file__) != '':
            os.chdir(os.path.dirname(__file__))

    def tearDown(self):
        os.chdir(self.resetdir)

    def runTest(self):
        self.test_presimulation()

    def test_presimulation(self):

        constantParameters = ['DRUG_0', 'KIN_0']

        # -------------- PYSB -----------------

        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..',
                                        'python', 'examples',
                                        'example_presimulation'))
        from createModel import model
        model.name = 'test_model_presimulation_pysb'
        amici.pysb2amici(model,
                         model.name,
                         verbose=False,
                         observables=['pPROT'],
                         constantParameters=constantParameters)
        sys.path.insert(0, model.name)
        import test_model_presimulation_pysb as modelModuleSBML
        model_pysb = modelModuleSBML.getModel()

        rdata_pysb = get_results(model_pysb)

        # -------------- SBML -----------------

        sbmlFile = os.path.join(os.path.dirname(__file__), '..', 'python',
                                'examples', 'example_presimulation',
                                'model_presimulation.xml')

        sbmlImporter = amici.SbmlImporter(sbmlFile)

        observables = amici.assignmentRules2observables(
            sbmlImporter.sbml,  # the libsbml model object
            filter_function=lambda variable: variable.getName() == 'pPROT'
        )
        outdir = 'test_model_presimulation_sbml'
        sbmlImporter.sbml2amici('test_model_presimulation_sbml',
                                outdir,
                                verbose=False,
                                observables=observables,
                                constantParameters=constantParameters)
        sys.path.insert(0, outdir)
        import test_model_presimulation_sbml as modelModuleSBML
        model_sbml = modelModuleSBML.getModel()

        rdata_sbml = get_results(model_sbml)







def get_results(model):
    solver = model.getSolver()
    model.setTimepoints(np.linspace(0, 60, 61))
    model.setReinitializeFixedParameterInitialStates(True)

    rdata_pre = amici.runAmiciSimulation(model, solver)
    edata = amici.ExpData(rdata_pre, 0.1, 0.0)
    edata.fixedParametersPreequilibration = [3, 0]
    edata.fixedParameters = [10, 2]
    edata.fixedParametersPresimulation = [10, 2]
    edata.fixedParametersPreequilibration = [3, 0]
    return amici.runAmiciSimulation(model, solver, edata),


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestAmiciSBMLModel())
    unittest.main()
    