#!/usr/bin/env python3

import sys
import h5py
import amici
import unittest
import importlib
import os
import re
import numpy as np

class TestAmiciPregeneratedModel(unittest.TestCase):
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
        expectedResults = h5py.File(self.expectedResultsFile, 'r')

        for subTest in expectedResults.keys():
            for case in list(expectedResults[subTest].keys()):
                if re.search('^sensi2',case) != None:
                    modelName = subTest + '_o2'
                else:
                    modelName = subTest

                with self.subTest(modelName=modelName, caseName=case):
                    print('running subTest modelName = ' + modelName + ', caseName = ' + case)
                    modelSwigFolder = os.path.join(os.path.dirname(__file__), '..', 'build', 'tests', 
                                                   'cpputest', 'external_' + modelName + '-prefix', 
                                                   'src', 'external_' + modelName + '-build', 'swig')
                    sys.path.insert(0, modelSwigFolder)
                    testModelModule = importlib.import_module(modelName)
                    self.model = testModelModule.getModel()
                    self.solver = self.model.getSolver()
                    amici.readModelDataFromHDF5(self.expectedResultsFile,
                                                self.model.get(),
                                                "/" + subTest + "/" + case + "/options")
                    amici.readSolverSettingsFromHDF5(self.expectedResultsFile,
                                                 self.solver.get(),
                                                 "/" + subTest + "/" + case + "/options")

                    edata = None
                    if 'data' in expectedResults[subTest][case].keys():
                        edata = amici.readSimulationExpData(self.expectedResultsFile,
                                                            "/" + subTest + "/" + case + "/data",
                                                            self.model.get())
                    rdata = amici.runAmiciSimulation(self.model, self.solver, edata)
                    
                    if edata and self.solver.getSensitivityMethod() and self.solver.getSensitivityOrder():
                        if not modelName.startswith('model_neuron'):
                            checkGradient(self.model, self.solver, edata)

                    if modelName == 'model_neuron_o2':
                        verifySimulationResults(rdata, expectedResults[subTest][case]['results'],rtol=1e-3)
                    else:
                        verifySimulationResults(rdata, expectedResults[subTest][case]['results'])
                        

def checkGradient(model, solver, edata):
    """Finite differences check for likelihood gradient
    
    Arguments:
        model:
        solver:
        edata:
    """
    from scipy.optimize import check_grad

    def func(x0, symbol='llh', x0full=None, plist=[], verbose=False):
        """Function of which the gradient is to be checked"""
        p = x0
        if len(plist):
            p = x0full[:]
            p[plist] = x0
        verbose and print('f: p=%s' % p)
        
        old_sensitivity_order = solver.getSensitivityOrder()
        old_parameters = model.getParameters()
        
        solver.setSensitivityOrder(amici.AMICI_SENSI_ORDER_NONE)
        model.setParameters(amici.DoubleVector(p))
        rdata = amici.runAmiciSimulation(model, solver, edata)
        
        solver.setSensitivityOrder(old_sensitivity_order)
        model.setParameters(old_parameters)
        
        res = np.sum(rdata[symbol])
        return res
    
    def grad(x0, symbol='llh', x0full=None, plist=[], verbose=False):
        """Gradient which is to be checked"""
        old_parameters = model.getParameters()
        old_plist = model.getParameterList()
        
        p = x0
        if len(plist):
            model.setParameterList(amici.IntVector(plist))
            p = x0full[:]
            p[plist] = x0
        else:
            model.requireSensitivitiesForAllParameters()
        verbose and print('g: p=%s' % p)
        
        model.setParameters(amici.DoubleVector(p))
        rdata = amici.runAmiciSimulation(model, solver, edata)
        
        model.setParameters(old_parameters)
        model.setParameterList(old_plist)

        res = rdata['s%s' % symbol]
        if not isinstance(res, float):
            if len(res.shape) == 3:
                res = np.sum(res, axis=(0, 2))
        return res
    
    p = np.array(model.getParameters())
      
    for ip in range(model.np()):
        plist = [ip]
        err_norm = check_grad(func, grad, p[plist], 'llh', p, [ip])
        print('sllh: p[%d]: |error|_2: %f' % (ip, err_norm))
    
    '''
    print()
    for ip in range(model.np()):
        plist = [ip]
        err_norm = check_grad(func, grad, p[plist], 'y', p, [ip])
        print('sy: p[%d]: |error|_2: %f' % (ip, err_norm))
    
    print()
    for ip in range(model.np()):
        plist = [ip]
        err_norm = check_grad(func, grad, p[plist], 'x', p, [ip])
        print('sx: p[%d]: |error|_2: %f' % (ip, err_norm))
    
    '''
    
def verifySimulationResults(rdata, expectedResults, atol=1e-8, rtol=1e-4):
    '''
    compares all fields of the simulation results in rdata against the expectedResults using the provided
    tolerances
    
    Arguments:
        rdata: simulation results as returned by amici.runAmiciSimulation
        expectedResults: stored test results 
        atol: absolute tolerance
        rtol: relative tolerance
    '''

    if expectedResults.attrs['status'][0] != 0:
        assert rdata['status'] == expectedResults.attrs['status'][0]
        return

    for field in expectedResults.keys():
        if field == 'diagnosis':
           for subfield in expectedResults[field].keys():
               checkResults(rdata, subfield, expectedResults[field][subfield][()], 0, 2)
        else:
            if field == 's2llh':
                checkResults(rdata, field, expectedResults[field][()], 1e-4, 1e-3)
            else:
                checkResults(rdata, field, expectedResults[field][()], atol, rtol)

    for attr in expectedResults.attrs.keys():
        checkResults(rdata, attr, expectedResults.attrs[attr], atol, rtol)


def checkResults(rdata, field, expected, atol, rtol):
    '''
    checks whether rdata[field] agrees with expected according to provided tolerances
    
    Arguments:
        rdata: simulation results as returned by amici.runAmiciSimulation
        field: name of the field to check
        expected: expected test results 
        atol: absolute tolerance
        rtol: relative tolerance
    '''

    result = rdata[field]
    adev = abs(result-expected)
    rdev = abs((result-expected))/(abs(expected)+rtol)

    if np.any(np.isnan(expected)):
        if len(expected) > 1 :
            assert all(np.isnan(result[np.isnan(expected)]))
        else: # subindexing fails for scalars
            assert np.isnan(result)
        adev = adev[~np.isnan(expected)]
        rdev = rdev[~np.isnan(expected)]

    if np.any(np.isinf(expected)):
        if len(expected) > 1 :
            assert all(np.isinf(result[np.isinf(expected)]))
        else: # subindexing fails for scalars
            assert np.isinf(result)
        adev = adev[~np.isinf(expected)]
        rdev = rdev[~np.isinf(expected)]

    assert np.all(np.logical_or(rdev <= rtol, adev <= atol))



if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestAmiciPregeneratedModel())
    unittest.main()
    
