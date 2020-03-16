#!/usr/bin/env python3

"""Run simulations with Matlab-AMICI pre-generated models and verify using
saved expectations."""

import sys
import h5py
import amici
import unittest
import importlib
import os
import copy
from amici.gradient_check import check_derivatives, check_results


class TestAmiciPregeneratedModel(unittest.TestCase):
    """
    TestCase class for tests that were pregenerated using the the matlab code
    generation routines and cmake build routines

    NOTE: requires having run `make python-tests` in /build/ before to build
    the python modules for the test models
    """

    expectedResultsFile = os.path.join(os.path.dirname(__file__),
                                       'cpputest', 'expectedResults.h5')

    def setUp(self):
        self.default_path = copy.copy(sys.path)
        self.resetdir = os.getcwd()

    def tearDown(self):
        os.chdir(self.resetdir)
        sys.path = self.default_path

    def runTest(self):
        """
        test runner routine that loads data expectedResults.h5 hdf file and
        runs individual models/settings as subTests
        """
        expected_results = h5py.File(self.expectedResultsFile, 'r')

        for subTest in expected_results.keys():
            for case in list(expected_results[subTest].keys()):
                if case.startswith('sensi2'):
                    model_name = subTest + '_o2'
                else:
                    model_name = subTest

                with self.subTest(modelName=model_name, caseName=case):
                    print(f'running {model_name}::{case}')

                    def assert_fun(x):
                        return self.assertTrue(x)

                    model_swig_folder = \
                        os.path.join(os.path.dirname(__file__), '..',
                                     'build', 'tests', 'cpputest',
                                     f'external_{model_name}-prefix',
                                     'src', f'external_{model_name}-build',
                                     'swig')
                    sys.path.insert(0, model_swig_folder)
                    test_model_module = importlib.import_module(model_name)
                    self.model = test_model_module.getModel()
                    self.solver = self.model.getSolver()
                    amici.readModelDataFromHDF5(
                        self.expectedResultsFile, self.model.get(),
                        f'/{subTest}/{case}/options'
                    )
                    amici.readSolverSettingsFromHDF5(
                        self.expectedResultsFile, self.solver.get(),
                        f'/{subTest}/{case}/options'
                    )

                    edata = None
                    if 'data' in expected_results[subTest][case].keys():
                        edata = amici.readSimulationExpData(
                            self.expectedResultsFile,
                            f'/{subTest}/{case}/data', self.model.get()
                        )
                    rdata = amici.runAmiciSimulation(self.model, self.solver,
                                                     edata)

                    check_derivative_opts = dict()

                    if model_name == 'model_nested_events':
                        check_derivative_opts['rtol'] = 1e-2
                    elif model_name == 'model_events':
                        check_derivative_opts['atol'] = 1e-3

                    if edata \
                            and self.solver.getSensitivityMethod() \
                            and self.solver.getSensitivityOrder() \
                            and len(self.model.getParameterList()) \
                            and not model_name.startswith('model_neuron') \
                            and not case.endswith('byhandpreeq'):
                        check_derivatives(self.model, self.solver, edata,
                                          assert_fun, **check_derivative_opts)

                    verify_simulation_opts = dict()

                    if model_name.startswith('model_neuron'):
                        verify_simulation_opts['atol'] = 1e-5
                        verify_simulation_opts['rtol'] = 1e-2

                    if model_name.startswith('model_robertson') and \
                            case == 'sensiforwardSPBCG':
                        verify_simulation_opts['atol'] = 1e-3
                        verify_simulation_opts['rtol'] = 1e-3

                    verify_simulation_results(
                        rdata, expected_results[subTest][case]['results'],
                        assert_fun, **verify_simulation_opts
                    )

                    if model_name == 'model_steadystate' and \
                            case == 'sensiforwarderrorint':
                        edata = amici.amici.ExpData(self.model.get())

                    if edata and model_name != 'model_neuron_o2' and not (
                        model_name == 'model_robertson' and
                        case == 'sensiforwardSPBCG'
                    ):
                        # Test runAmiciSimulations: ensure running twice
                        # with same ExpData yields same results
                        if isinstance(edata, amici.amici.ExpData):
                            edatas = [edata, edata]
                        else:
                            edatas = [edata.get(), edata.get()]

                        rdatas = amici.runAmiciSimulations(
                            self.model, self.solver, edatas, num_threads=2,
                            failfast=False
                        )
                        verify_simulation_results(
                            rdatas[0],
                            expected_results[subTest][case]['results'],
                            assert_fun, **verify_simulation_opts
                        )
                        verify_simulation_results(
                            rdatas[1],
                            expected_results[subTest][case]['results'],
                            assert_fun, **verify_simulation_opts
                        )

                    self.assertRaises(
                        RuntimeError,
                        self.model.getParameterByName,
                        'thisParameterDoesNotExist'
                    )


def verify_simulation_results(rdata, expected_results, assert_fun,
                              atol=1e-8, rtol=1e-4):
    """
    compares all fields of the simulation results in rdata against the
    expectedResults using the provided tolerances

    Arguments:
        rdata: simulation results as returned by amici.runAmiciSimulation
        expected_results: stored test results
        atol: absolute tolerance
        rtol: relative tolerance
    """

    if expected_results.attrs['status'][0] != 0:
        assert rdata['status'] == expected_results.attrs['status'][0]
        return

    for field in expected_results.keys():
        if field == 'diagnosis':
            for subfield in ['J', 'xdot']:
                check_results(rdata, subfield,
                              expected_results[field][subfield][()],
                              assert_fun, 1e-8, 1e-8)
        else:
            if field == 's2llh':
                check_results(rdata, field, expected_results[field][()],
                              assert_fun, 1e-4, 1e-3)
            else:
                check_results(rdata, field, expected_results[field][()],
                              assert_fun, atol, rtol)

    for attr in expected_results.attrs.keys():
        check_results(rdata, attr, expected_results.attrs[attr], assert_fun,
                      atol, rtol)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestAmiciPregeneratedModel())
    unittest.main()
