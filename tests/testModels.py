#!/usr/bin/env python3

import sys
import h5py
import amici
import unittest
import importlib
import os
import numpy as np
import copy


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

                    # todo: set higher tolerances in testcase options and
                    # regenerate results
                    if model_name == 'model_jakstat_adjoint_o2':
                        self.solver.setRelativeTolerance(1e-10)
                        self.solver.setAbsoluteTolerance(1e-10)

                    if model_name in [
                        'model_jakstat_adjoint', 'model_nested_events',
                        'model_steadystate'
                    ]:
                        self.solver.setRelativeTolerance(1e-12)
                        self.solver.setAbsoluteTolerance(1e-12)

                    if model_name in 'model_events':
                        epsilon = 1e-3
                    elif model_name == 'model_nested_events':
                        epsilon = 1e-4
                    else:
                        epsilon = 1e-5

                    if edata \
                            and self.solver.getSensitivityMethod() \
                            and self.solver.getSensitivityOrder() \
                            and len(self.model.getParameterList()) \
                            and not model_name.startswith('model_neuron') \
                            and not case.endswith('byhandpreeq'):
                        check_derivatives(self.model, self.solver, edata,
                                          assert_fun, epsilon=epsilon)

                    if model_name == 'model_neuron_o2':
                        self.solver.setRelativeTolerance(1e-14)
                        verify_simulation_results(
                            rdata, expected_results[subTest][case]['results'],
                            assert_fun,
                            atol=1e-5, rtol=1e-4
                        )
                    elif model_name == 'model_robertson':
                        verify_simulation_results(
                            rdata, expected_results[subTest][case]['results'],
                            assert_fun,
                            atol=1e-6, rtol=1e-2
                        )
                    else:
                        verify_simulation_results(
                            rdata, expected_results[subTest][case]['results'],
                            assert_fun,
                        )

                    if edata and model_name != 'model_neuron_o2':
                        # Test runAmiciSimulations: ensure running twice
                        # with same ExpData yields same results
                        edatas = [edata.get(), edata.get()]
                        rdatas = amici.runAmiciSimulations(
                            self.model, self.solver, edatas, num_threads=2
                        )
                        verify_simulation_results(
                            rdatas[0],
                            expected_results[subTest][case]['results'],
                            assert_fun,
                        )
                        verify_simulation_results(
                            rdatas[1],
                            expected_results[subTest][case]['results'],
                            assert_fun,
                        )

                    self.assertRaises(
                        RuntimeError,
                        self.model.getParameterByName,
                        'thisParameterDoesNotExist'
                    )


def check_close(result, expected, assert_fun, atol, rtol, field, ip=None):
    close = np.isclose(result, expected, atol=atol, rtol=rtol, equal_nan=True)

    if not close.all():
        if ip is None:
            index_str = ''
            check_type = 'Regression check  '
        else:
            index_str = f'at index ip={ip} '
            check_type = 'FD check '
        print(f'{check_type} failed for {field} {index_str}for '
              f'{close.sum()} indices:')
        adev = np.abs(result - expected)
        rdev = np.abs((result - expected)/(expected + atol))
        print(f'max(adev): {adev.max()}, max(rdev): {rdev.max()}')

    assert_fun(close.all())


def check_finite_difference(x0, model, solver, edata, ip, fields,
                            assert_fun, atol=1e-4, rtol=1e-4, epsilon=1e-5):
    old_sensitivity_order = solver.getSensitivityOrder()
    old_parameters = model.getParameters()
    old_plist = model.getParameterList()

    # sensitivity
    p = copy.deepcopy(x0)
    plist = [ip]

    solver.setSensitivityOrder(amici.SensitivityOrder_first)
    model.setParameters(p)
    model.setParameterList(plist)
    rdata = amici.runAmiciSimulation(model, solver, edata)

    # finite difference
    solver.setSensitivityOrder(amici.SensitivityOrder_first)

    # forward:
    p = copy.deepcopy(x0)
    p[ip] += epsilon/2
    model.setParameters(p)
    rdataf = amici.runAmiciSimulation(model, solver, edata)

    # backward:
    p = copy.deepcopy(x0)
    p[ip] -= epsilon/2
    model.setParameters(p)
    rdatab = amici.runAmiciSimulation(model, solver, edata)

    for field in fields:
        sensi_raw = rdata[f's{field}']
        fd = (rdataf[field]-rdatab[field])/epsilon
        if len(sensi_raw.shape) == 1:
            sensi = sensi_raw[0]
        elif len(sensi_raw.shape) == 2:
            sensi = sensi_raw[:, 0]
        elif len(sensi_raw.shape) == 3:
            sensi = sensi_raw[:, 0, :]
        else:
            assert_fun(False)  # not implemented

        check_close(sensi, fd, assert_fun, atol, rtol, field, ip=ip)

    solver.setSensitivityOrder(old_sensitivity_order)
    model.setParameters(old_parameters)
    model.setParameterList(old_plist)


def check_derivatives(model, solver, edata, assert_fun,
                      atol=1e-4, rtol=1e-4, epsilon=1e-5):
    """Finite differences check for likelihood gradient

    Arguments:
        model: amici model
        solver: amici solver
        edata: exp data
        atol: absolute tolerance
        rtol: relative tolerance
        epsilon: finite difference step-size
    """
    from scipy.optimize import check_grad

    p = np.array(model.getParameters())

    rdata = amici.runAmiciSimulation(model, solver, edata)

    fields = ['llh']

    leastsquares_applicable = \
        solver.getSensitivityMethod() == amici.SensitivityMethod_forward

    if 'ssigmay' in rdata.keys():
        if rdata['ssigmay'] is not None:
            if rdata['ssigmay'].any():
                leastsquares_applicable = False

    if leastsquares_applicable:
        fields += ['res', 'x', 'y']

        check_results(rdata, 'FIM',
                      np.dot(rdata['sres'].transpose(), rdata['sres']),
                      assert_fun,
                      1e-8, 1e-4)
        check_results(rdata, 'sllh',
                      -np.dot(rdata['res'].transpose(), rdata['sres']),
                      assert_fun,
                      1e-8, 1e-4)
    for ip in range(len(p)):
        check_finite_difference(p, model, solver, edata, ip, fields,
                                assert_fun, atol=atol, rtol=rtol,
                                epsilon=epsilon)


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
                              assert_fun, 0, 2)
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


def check_results(rdata, field, expected, assert_fun, atol, rtol):
    """
    checks whether rdata[field] agrees with expected according to provided
    tolerances
    
    Arguments:
        rdata: simulation results as returned by amici.runAmiciSimulation
        field: name of the field to check
        expected: expected test results 
        atol: absolute tolerance
        rtol: relative tolerance
    """

    result = rdata[field]
    if type(result) is float:
        result = np.array(result)

    check_close(result, expected, assert_fun, atol, rtol, field)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestAmiciPregeneratedModel())
    unittest.main()

