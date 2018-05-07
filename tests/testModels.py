#!/usr/bin/env python3

import sys
import h5py
import amici
import unittest
import importlib
import os
import numpy as np

class TestAmiciPregeneratedModel(unittest.TestCase):
    '''
    TestCase class for tests that were pregenerated using the the matlab code generation routines and cmake
    build routines
    '''

    expectedResultsFile = os.path.join(amici.amici_path,'tests','cpputest','expectedResults.h5')

    def runTest(self):


        expectedResults = h5py.File(self.expectedResultsFile, 'r')

        for subTest in expectedResults.keys():
            sys.path.insert(0, os.path.join(amici.amici_path, 'models', subTest, 'build', 'swig'))
            testModelModule = importlib.import_module(subTest)
            self.model = testModelModule.getModel()
            self.solver = self.model.getSolver()

            for case in list(expectedResults[subTest].keys()):
                with self.subTest(modelName=subTest, caseName=case):
                    amici.readModelDataFromHDF5(self.expectedResultsFile,
                                                self.model.get(),
                                                "/" + subTest + "/" + case + "/options")
                    amici.readSolverSettingsFromHDF5(self.expectedResultsFile,
                                                 self.solver.get(),
                                                 "/" + subTest + "/" + case + "/options")

                    if 'data' in expectedResults[subTest][case].keys():
                        edata = amici.readSimulationExpData(self.expectedResultsFile,
                                                            "/" + subTest + "/" + case + "/data",
                                                            self.model.get())
                        rdata = amici.runAmiciSimulation(self.solver.get(), edata.get(), self.model.get())
                    else:
                        rdata = amici.runAmiciSimulation(self.solver.get(), None, self.model.get())

                    verifySimulationResults(rdata, expectedResults[subTest][case]['results'])




def verifySimulationResults(rdata,expectedResults,atol=1e-8,rtol=1e-7):

    if expectedResults.attrs['status'][0] != 0:
        assert rdata.status == expectedResults.attrs['status'][0]
        return

    for field in expectedResults.keys():
        if field == 'diagnosis':
           for subfield in expectedResults[field].keys():
               checkResults(rdata, subfield, expectedResults[field][subfield][()],0,1)
        else:
            checkResults(rdata, field, expectedResults[field][()])

    for attr in expectedResults.attrs.keys():
        checkResults(rdata, attr, expectedResults.attrs[attr])


def checkResults(rdata, field, expected, atol=1e-9, rtol=1e-4):
    result = getFieldAsNumPyArray(rdata, field)
    adev = abs(result-expected)
    rdev = abs((result-expected))/(abs(expected)+rtol)

    if np.any(np.isnan(expected)):
        assert all(np.isnan(result[np.isnan(expected)]))
        adev = adev[~np.isnan(expected)]
        rdev = rdev[~np.isnan(expected)]

    assert np.all(np.logical_or(rdev <= rtol,adev <= atol))



def getFieldAsNumPyArray(rdata,field):

    if field == 't':
        field = 'ts'

    fieldDimensions = {'ts': [rdata.nt],
                       'x': [rdata.nt, rdata.nx],
                       'x0': [rdata.nx],
                       'sx': [rdata.nt, rdata.nplist, rdata.nx],
                       'sx0': [rdata.nx, rdata.nplist],
                       # observables
                       'y': [rdata.nt, rdata.ny],
                       'sigmay': [rdata.nt, rdata.ny],
                       'sy': [rdata.nt, rdata.nplist, rdata.ny],
                       'ssigmay': [rdata.nt, rdata.nplist, rdata.ny],
                       # event observables
                       'z': [rdata.nmaxevent, rdata.nz],
                       'rz': [rdata.nmaxevent, rdata.nz],
                       'sigmaz': [rdata.nmaxevent, rdata.nz],
                       'sz': [rdata.nmaxevent, rdata.nplist, rdata.nz],
                       'srz': [rdata.nmaxevent, rdata.nplist, rdata.nz],
                       'ssigmaz': [rdata.nmaxevent, rdata.nplist, rdata.nz],
                       # objective function
                       'sllh': [rdata.nplist],
                       's2llh': [rdata.np, rdata.nplist],
                       # diagnosis
                       'J': [rdata.nx,rdata.nx],
                       'xdot': [rdata.nx],
                       'newton_numlinsteps': [rdata.newton_maxsteps, 2],
                       'newton_numsteps': [2, 1],
                       'numsteps': [rdata.nt],
                       'numrhsevals': [rdata.nt],
                       'numerrtestfails': [rdata.nt],
                       'numnonlinsolvconvfails': [rdata.nt],
                       'order': [rdata.nt],
                       'numstepsB': [rdata.nt],
                       'numrhsevalsB': [rdata.nt],
                       'numerrtestfailsB': [rdata.nt],
                       'numnonlinsolvconvfailsB': [rdata.nt],
                      }
    if field in fieldDimensions.keys():
        if len(fieldDimensions[field]) == 1:
            return np.array(rdata.__getattr__(field))
        else:
            return np.array(rdata.__getattr__(field)).reshape(fieldDimensions[field])

    else:
        return float(rdata.__getattr__(field))


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestAmiciPregeneratedModel())
    unittest.main()








