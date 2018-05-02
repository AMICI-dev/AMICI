#!/usr/bin/env python3

from model import Model
import os
import sys
import importlib
import numpy as np
import amici

dirname = os.path.split(os.path.abspath(__file__))[0]
amici_path = os.path.split(dirname)[0]
test_path = os.path.join(amici_path,'tests','sbml-semantic-test-cases','cases','semantic')

def runTest(testId):
    try:
        current_test_path = os.path.join(test_path, testId)


        # settings
        settingsFile = os.path.join(current_test_path, testId + '-settings.txt')
        settings = {}
        with open(settingsFile) as f:
            for line in f:
                if not line == '\n':
                    (key, val) = line.split(':')
                    settings[key] = val
        ts = np.linspace(float(settings['start']),float(settings['start'])+float(settings['duration']),int(settings['steps'])+1)
        atol = float(settings['absolute'])
        rtol = float(settings['relative'])

        # results
        resultsFile = os.path.join(current_test_path, testId + '-results.csv')
        results = np.genfromtxt(resultsFile, delimiter=',')
        test_x = results[1:, 1:]

        # model
        sbmlFile = os.path.join(current_test_path, testId + '-sbml-l3v2.xml')
        wrapper = Model(sbmlFile,'SBMLTest' + testId)
        wrapper.wrapModel()
        sys.path.insert(0,os.path.join(wrapper.model_path,'build','swig'))
        mod = importlib.import_module(wrapper.modelname)

        model = mod.getModel()
        model.setTimepoints(mod.amici.DoubleVector(ts))
        solver = model.getSolver()
        #solver.setRelativeTolerance(rtol)
        #solver.setAbsoluteTolerance(atol)
        rdata = amici.runAmiciSimulation(solver.get(),None,model.get())
        simulated_x = np.array(rdata.x).reshape([len(ts),model.nx])


        adev = abs(simulated_x-test_x)
        adev = adev[~np.isnan(adev)]
        rdev = abs((simulated_x-test_x)/test_x)
        rdev = rdev[~np.isnan(rdev)]
        if(not np.all(adev<atol)):
            raise Exception('Absolute tolerance violated')

        if (not np.all(rdev < rtol)):
            raise Exception('Relative tolerance violated')


    except AssertionError:
        pass

def getTestStr(testId):
    testStr = str(testId)
    testStr = '0'*(5-len(testStr)) + testStr
    return testStr


for testId in range(17,1781):
    runTest(getTestStr(testId))

#model = Model('/Users/F.Froehlich/Downloads/Speedy_v3_r403445_v1.sbml','speedy')
#model.wrapModel()