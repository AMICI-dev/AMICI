#!/usr/bin/env python3

import sys
import amici
import unittest
import os
import numpy as np
from testModels import checkDerivatives

class TestAmiciSBMLModel(unittest.TestCase):
    '''
    TestCase class for testing SBML import and simulation from AMICI python interface
    '''

    expectedResultsFile = os.path.join(os.path.dirname(__file__),
                                       'cpputest', 'expectedResults.h5')

    def setUp(self):
        self.resetdir = os.getcwd()
        if os.path.dirname(__file__) != '':
            os.chdir(os.path.dirname(__file__))

    def tearDown(self):
        os.chdir(self.resetdir)

    def runTest(self):
        self.test_presimulation()
        self.test_steadystate_scaled()

    def test_presimulation(self):
        sbmlFile = os.path.join(os.path.dirname(__file__), '..', 'python',
                                'examples', 'example_presimulation',
                                'model_presimulation.xml')

        sbmlImporter = amici.SbmlImporter(sbmlFile)

        constantParameters = ['DRUG_0', 'KIN_0']

        observables = amici.assignmentRules2observables(
            sbmlImporter.sbml,  # the libsbml model object
            filter_function=lambda variable: variable.getName() == 'pPROT_obs'
        )
        outdir = 'test_model_presimulation'
        sbmlImporter.sbml2amici('test_model_presimulation',
                                outdir,
                                verbose=False,
                                observables=observables,
                                constantParameters=constantParameters)
        sys.path.insert(0, outdir)
        import test_model_presimulation as modelModule
        model = modelModule.getModel()
        solver = model.getSolver()
        solver.setNewtonMaxSteps(0)
        model.setTimepoints(np.linspace(0, 60, 61))
        model.setSteadyStateSensitivityMode(
            amici.SteadyStateSensitivityMode_simulationFSA
        )
        solver.setSensitivityOrder(amici.SensitivityOrder_first)
        model.setReinitializeFixedParameterInitialStates(True)

        rdata = amici.runAmiciSimulation(model, solver)
        edata = amici.ExpData(rdata, 0.1, 0.0)
        edata.fixedParameters = [10, 2]
        edata.fixedParametersPresimulation = [10, 2]
        edata.fixedParametersPreequilibration = [3, 0]
        self.assertIsInstance(
            amici.runAmiciSimulation(model, solver, edata),
            dict)

        solver.setRelativeTolerance(1e-12)
        solver.setAbsoluteTolerance(1e-12)
        checkDerivatives(model, solver, edata, epsilon=1e-6)

    def test_steadystate_scaled(self):
        '''
        Test SBML import and simulation from AMICI python interface
        '''

        sbmlFile = os.path.join(os.path.dirname(__file__), '..', 'python',
                                'examples', 'example_steadystate',
                                'model_steadystate_scaled.xml')
        sbmlImporter = amici.SbmlImporter(sbmlFile)

        observables = amici.assignmentRules2observables(
            sbmlImporter.sbml,
            filter_function=lambda variable:
                variable.getId().startswith('observable_') and
                not variable.getId().endswith('_sigma')
        )

        outdir = 'test_model_steadystate_scaled'
        sbmlImporter.sbml2amici('test_model_steadystate_scaled',
                                outdir,
                                observables=observables,
                                constantParameters=['k0'],
                                sigmas={'observable_x1withsigma':
                                            'observable_x1withsigma_sigma'})

        sys.path.insert(0, outdir)
        import test_model_steadystate_scaled as modelModule

        model = modelModule.getModel()
        model.setTimepoints(np.linspace(0, 60, 60))
        solver = model.getSolver()
        solver.setSensitivityOrder(amici.SensitivityOrder_first)
        rdata = amici.runAmiciSimulation(model, solver)
        edata = [amici.ExpData(rdata, 1, 0)]
        rdata = amici.runAmiciSimulations(model, solver, edata)

        # check roundtripping of DataFrame conversion
        df_edata = amici.getDataObservablesAsDataFrame(model, edata)
        edata_reconstructed = amici.getEdataFromDataFrame(model, df_edata)

        self.assertTrue(
            np.isclose(
                amici.edataToNumPyArrays(edata[0])
                ['observedData'],
                amici.edataToNumPyArrays(edata_reconstructed[0])
                ['observedData'],
            ).all()
        )
        self.assertTrue(
            np.isclose(
                amici.edataToNumPyArrays(edata[0])
                ['observedDataStdDev'],
                amici.edataToNumPyArrays(edata_reconstructed[0])
                ['observedDataStdDev'],
            ).all()
        )
        if len(edata[0].fixedParameters):
            self.assertListEqual(
                list(edata[0].fixedParameters),
                list(edata_reconstructed[0].fixedParameters),
            )
        else:
            self.assertListEqual(
                list(model.getFixedParameters()),
                list(edata_reconstructed[0].fixedParameters),
            )

        self.assertListEqual(
            list(edata[0].fixedParametersPreequilibration),
            list(edata_reconstructed[0].fixedParametersPreequilibration),
        )

        df_state = amici.getSimulationStatesAsDataFrame(model, edata, rdata)
        self.assertTrue(
            np.isclose(
                rdata[0]['x'],
                df_state[list(model.getStateIds())].values
            ).all()
        )
        df_obs = amici.getSimulationObservablesAsDataFrame(model, edata, rdata)
        self.assertTrue(
            np.isclose(
                rdata[0]['y'],
                df_obs[list(model.getObservableIds())].values
            ).all()
        )
        amici.getResidualsAsDataFrame(model, edata, rdata)

        solver.setRelativeTolerance(1e-12)
        solver.setAbsoluteTolerance(1e-12)
        checkDerivatives(model, solver, edata[0], atol=1e-3,
                         rtol=1e-3, epsilon=1e-7)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestAmiciSBMLModel())
    unittest.main()
    