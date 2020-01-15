#!/usr/bin/env python3

"""Various test cases for the AMICI Python interface"""

import sys
import amici
import unittest
import os
import copy
import numpy as np
from testModels import check_derivatives


class TestAmiciSBMLModel(unittest.TestCase):
    """
    TestCase class for testing SBML import and simulation from AMICI python
    interface
    """

    def setUp(self):
        self.default_path = copy.copy(sys.path)
        self.resetdir = os.getcwd()

        if os.path.dirname(__file__) != '':
            os.chdir(os.path.dirname(__file__))

    def tearDown(self):
        os.chdir(self.resetdir)
        sys.path = self.default_path

    def runTest(self):
        self.test_presimulation()
        self.test_steadystate_scaled()
        self.test_likelihoods()

    def test_presimulation(self):
        def assert_fun(x):
            return self.assertTrue(x)

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
            amici.ReturnDataView)

        solver.setRelativeTolerance(1e-12)
        solver.setAbsoluteTolerance(1e-12)
        check_derivatives(model, solver, edata, assert_fun, epsilon=1e-4)

    def test_steadystate_scaled(self):
        """
        Test SBML import and simulation from AMICI python interface
        """

        model_module = self.test_steadystate_import()

        self.steadystate_simulation(model_module=model_module)

        # Run some additional tests which need a working Model,
        # but don't need precomputed expectations.
        test_set_parameters_by_dict(model_module)

    def test_steadystate_import(self):
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
        import test_model_steadystate_scaled as model_module

        return model_module

    def steadystate_simulation(self, model_module):
        def assert_fun(x):
            return self.assertTrue(x)

        model = model_module.getModel()
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
                amici.ExpDataView(edata[0])
                ['observedData'],
                amici.ExpDataView(edata_reconstructed[0])
                ['observedData'],
            ).all()
        )
        self.assertTrue(
            np.isclose(
                amici.ExpDataView(edata[0])
                ['observedDataStdDev'],
                amici.ExpDataView(edata_reconstructed[0])
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
        check_derivatives(model, solver, edata[0], assert_fun, atol=1e-3,
                          rtol=1e-3, epsilon=1e-4)


    def test_likelihoods(self):
        """
        Test the custom noise distributions used to define cost functions.
        """
        def assert_fun(x):
            return self.assertTrue(x)

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

        # assign different noise models

        obs_keys = list(observables.keys())

        # exponentiate observable formulas
        obs1 = observables[obs_keys[1]]
        obs3 = observables[obs_keys[3]]
        obs1['formula'] = '10^(' + obs1['formula'] + ')'
        obs3['formula'] = 'exp(' + obs3['formula'] + ')'

        # customize noise distributions
        noise_distributions = {
            obs_keys[0]: 'normal',
            obs_keys[1]: 'log-normal',
            obs_keys[2]: 'laplace',
            obs_keys[3]: 'log10-laplace',
        }

        outdir = 'test_likelihoods'
        sbmlImporter.sbml2amici('test_likelihoods',
                                outdir,
                                observables=observables,
                                constantParameters=['k0'],
                                sigmas={'observable_x1withsigma':
                                        'observable_x1withsigma_sigma'},
                                noise_distributions=noise_distributions
                                )

        sys.path.insert(0, outdir)
        import test_likelihoods as modelModule

        model = modelModule.getModel()
        model.setTimepoints(np.linspace(0, 60, 60))
        solver = model.getSolver()
        solver.setSensitivityOrder(amici.SensitivityOrder_first)

        # run model once to create an edata
        rdata = amici.runAmiciSimulation(model, solver)
        edata = [amici.ExpData(rdata, 1, 0)]

        # just make all observables positive since some are logarithmic
        for ed in edata:
            y = ed.getObservedData()
            y = tuple([max(val, 1e-4) for val in y])
            ed.setObservedData(y)

        # and now run for real and also compute likelihood values
        rdata = amici.runAmiciSimulations(model, solver, edata)[0]

        # output for easy debugging
        for key in ['llh', 'sllh']:
            print(key, rdata[key])

        # it would be good to compute the expected llh+sllh by hand,
        # here, we only check if they make overall sense
        self.assertTrue(np.isfinite(rdata['llh']))
        self.assertTrue(np.all(np.isfinite(rdata['sllh'])))
        self.assertTrue(np.any(rdata['sllh']))


def test_set_parameters_by_dict(model_module):
    """Test setting parameter via id/name => value dicts"""

    model = model_module.getModel()
    old_parameter_values = model.getParameters()
    parameter_ids = model.getParameterIds()
    change_par_id = parameter_ids[-1]
    new_par_val = 0.1234
    old_par_val = model.getParameterById(change_par_id)

    assert model.getParameterById(change_par_id) != new_par_val
    model.setParameterById({change_par_id: new_par_val})
    assert model.getParameterById(change_par_id) == new_par_val
    # reset and check we are back to original
    model.setParameterById(change_par_id, old_par_val)
    assert model.getParameters() == old_parameter_values

    # Same for by-name
    parameter_names = model.getParameterNames()
    change_par_name = parameter_names[-1]
    model.setParameterByName({change_par_name: new_par_val})
    assert model.getParameterByName(change_par_name) == new_par_val
    model.setParameterByName(change_par_name, old_par_val)
    assert model.getParameters() == old_parameter_values


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestAmiciSBMLModel())
    unittest.main()
