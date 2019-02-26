#!/usr/bin/env python3

import sys
import h5py
import amici
import unittest
import importlib
import os
import re
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
                    
                    if edata \
                            and self.solver.getSensitivityMethod() \
                            and self.solver.getSensitivityOrder() \
                            and len(self.model.getParameterList()) \
                            and not model_name.startswith('model_neuron') \
                            and not case.endswith('byhandpreeq'):
                        check_derivatives(self.model, self.solver, edata)

                    if model_name == 'model_neuron_o2':
                        self.solver.setRelativeTolerance(1e-14)
                        verify_simulation_results(
                            rdata, expected_results[subTest][case]['results'],
                            atol=1e-5, rtol=1e-4
                        )
                    elif model_name == 'model_robertson':
                        verify_simulation_results(
                            rdata, expected_results[subTest][case]['results'],
                            atol=1e-6, rtol=1e-2
                        )
                    else:
                        verify_simulation_results(
                            rdata, expected_results[subTest][case]['results']
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
                            expected_results[subTest][case]['results']
                        )
                        verify_simulation_results(
                            rdatas[1],
                            expected_results[subTest][case]['results']
                        )

                    self.assertRaises(
                        RuntimeError,
                        self.model.getParameterByName,
                        'thisParameterDoesNotExist'
                    )


def check_derivatives(model, solver, edata,
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

    def func(x0, symbol='llh', x0full=None, plist=None, verbose=False):
        """Function of which the gradient is to be checked"""
        if plist is None:
            plist = []
        p = copy.deepcopy(x0)
        if len(plist):
            p = copy.deepcopy(x0full)
            p[plist] = copy.deepcopy(x0)
        verbose and print('f: p=%s' % p)
        
        old_sensitivity_order = solver.getSensitivityOrder()
        old_parameters = model.getParameters()
        
        solver.setSensitivityOrder(amici.SensitivityOrder_none)
        model.setParameters(p)
        rdata = amici.runAmiciSimulation(model, solver, edata)
        
        solver.setSensitivityOrder(old_sensitivity_order)
        model.setParameters(old_parameters)

        res = np.sum(rdata[symbol])
        return res
    
    def grad(x0, symbol='llh', x0full=None, plist=None, verbose=False):
        """Gradient which is to be checked"""
        if plist is None:
            plist = []
        old_parameters = model.getParameters()
        old_plist = model.getParameterList()
        
        p = copy.deepcopy(x0)
        if len(plist):
            model.setParameterList(plist)
            p = copy.deepcopy(x0full)
            p[plist] = copy.deepcopy(x0)
        else:
            model.requireSensitivitiesForAllParameters()
        verbose and print('g: p=%s' % p)
        
        model.setParameters(p)

        rdata = amici.runAmiciSimulation(model, solver, edata)

        model.setParameters(old_parameters)
        model.setParameterList(old_plist)

        res = rdata['s%s' % symbol]
        if not isinstance(res, float):
            if len(res.shape) == 2:
                res = np.sum(res, axis=(0,))
            if len(res.shape) == 3:
                res = np.sum(res, axis=(0, 2))
        return res
    
    p = np.array(model.getParameters())

    rdata = amici.runAmiciSimulation(model, solver, edata)
      
    for ip in range(model.np()):
        plist = [ip]
        err_norm = check_grad(func, grad, copy.deepcopy(p[plist]), 'llh',
                              copy.deepcopy(p), [ip], epsilon=epsilon)
        print(f'sllh: p[{ip}]: |error|_2: {err_norm}')
        assert err_norm <= rtol * abs(rdata["sllh"][ip]) \
            or err_norm <= atol

    leastsquares_applicable = \
        solver.getSensitivityMethod() == amici.SensitivityMethod_forward

    if 'ssigmay' in rdata.keys():
        if rdata['ssigmay'] is not None:
            if rdata['ssigmay'].any():
                leastsquares_applicable = False

    if leastsquares_applicable:
        for ip in range(model.np()):
            plist = [ip]
            err_norm = check_grad(func, grad, copy.deepcopy(p[plist]), 'res',
                                  copy.deepcopy(p), [ip], epsilon=epsilon)
            print('sres: p[%d]: |error|_2: %f' % (ip, err_norm))
            assert err_norm < atol \
                or err_norm < rtol * abs(rdata['sres'][:, ip].sum())

        check_results(rdata, 'FIM',
                      np.dot(rdata['sres'].transpose(), rdata['sres']),
                      1e-8, 1e-4)
        check_results(rdata, 'sllh',
                      -np.dot(rdata['res'].transpose(), rdata['sres']),
                      1e-8, 1e-4)

        print()
        for ip in range(model.np()):
            plist = [ip]
            err_norm = check_grad(func, grad, copy.deepcopy(p[plist]), 'y',
                                  copy.deepcopy(p), [ip], epsilon=epsilon)
            print('sy: p[%d]: |error|_2: %f' % (ip, err_norm))
            assert err_norm < atol \
                or err_norm < rtol * abs(rdata['sy'][:, ip].sum())

        print()
        for ip in range(model.np()):
            plist = [ip]
            err_norm = check_grad(func, grad, copy.deepcopy(p[plist]), 'x',
                                  copy.deepcopy(p), [ip], epsilon=epsilon)
            print('sx: p[%d]: |error|_2: %f' % (ip, err_norm))
            assert err_norm < atol \
                or err_norm < rtol * abs(rdata['sx'][:, ip].sum())


def verify_simulation_results(rdata, expected_results, atol=1e-8, rtol=1e-4):
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
                              expected_results[field][subfield][()], 0, 2)
        else:
            if field == 's2llh':
                check_results(rdata, field,
                              expected_results[field][()], 1e-4, 1e-3)
            else:
                check_results(rdata, field,
                              expected_results[field][()], atol, rtol)

    for attr in expected_results.attrs.keys():
        check_results(rdata, attr, expected_results.attrs[attr], atol, rtol)


def check_results(rdata, field, expected, atol, rtol):
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

    adev = abs(result - expected)
    rdev = abs((result - expected)) / (abs(expected) + rtol)

    if np.any(np.isnan(expected)):
        if len(expected) > 1 :
            assert all(np.isnan(result[np.isnan(expected)]))
        else:  # subindexing fails for scalars
            assert np.isnan(result)
        adev = adev[~np.isnan(expected)]
        rdev = rdev[~np.isnan(expected)]

    if np.any(np.isinf(expected)):
        if len(expected) > 1 :
            assert all(np.isinf(result[np.isinf(expected)]))
        else:  # subindexing fails for scalars
            assert np.isinf(result)
        adev = adev[~np.isinf(expected)]
        rdev = rdev[~np.isinf(expected)]

    if not np.all(np.logical_or(rdev <= rtol, adev <= atol)):
        print('Failed to meet tolerances in ' + field + ':')
        print('adev:')
        print(adev[np.logical_and(rdev > rtol, adev > atol)])
        print('rdev:')
        print(rdev[np.logical_and(rdev > rtol, adev > atol)])
        assert np.all(np.logical_or(rdev <= rtol, adev <= atol))


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestAmiciPregeneratedModel())
    unittest.main()

