#!/usr/bin/env python3
<<<<<<< HEAD

from model import Model
import traceback
=======
>>>>>>> origin/feature_swig
import os
import sys
import importlib
import numpy as np
import amici

dirname = os.path.split(os.path.abspath(__file__))[0]
amici_path = os.path.split(dirname)[0]
test_path = os.path.join(amici_path,'tests','sbml-semantic-test-cases','cases','semantic')


def runTest(testId, logfile):
    try:
        current_test_path = os.path.join(test_path, testId)

        # results
        resultsFile = os.path.join(current_test_path, testId + '-results.csv')
        results = np.genfromtxt(resultsFile, delimiter=',')

        # model
        sbmlFile = os.path.join(current_test_path, testId + '-sbml-l3v2.xml')

        # fallback l3v1
        if not os.path.isfile(sbmlFile):
            sbmlFile = os.path.join(current_test_path, testId + '-sbml-l3v1.xml')

        # fallback l2v5
        if not os.path.isfile(sbmlFile):
            sbmlFile = os.path.join(current_test_path, testId + '-sbml-l2v5.xml')

        wrapper = amici.SbmlImporter(sbmlFile)
        wrapper.sbml2amici('SBMLTest' + testId)

        # settings
        settingsFile = os.path.join(current_test_path, testId + '-settings.txt')
        settings = {}
        with open(settingsFile) as f:
            for line in f:
                if not line == '\n':
                    (key, val) = line.split(':')
                    settings[key] = val
        ts = np.linspace(float(settings['start']), float(settings['start']) + float(settings['duration']),
                         int(settings['steps']) + 1)
        atol = float(settings['absolute'])
        rtol = float(settings['relative'])

        sys.path.insert(0, os.path.join(wrapper.modelPath, 'build', 'swig'))
        mod = importlib.import_module(wrapper.modelName)

        model = mod.getModel()
        model.setTimepoints(mod.amici.DoubleVector(ts))
        solver = model.getSolver()
        solver.setMaxSteps(int(1e6))
        solver.setRelativeTolerance(rtol / 1000.0)
        solver.setAbsoluteTolerance(atol / 1000.0)
        rdata = amici.runAmiciSimulation(solver.get(), None, model.get())
        amountSpecies = settings['amount'].replace(' ', '').replace('\n', '').split(',')
        simulated_x = np.array(rdata.x).reshape([len(ts), model.nx])
        test_x = results[1:, [1+ wrapper.speciesIndex[variable]  for variable in settings['variables'].replace(' ', '').split(',') if variable in wrapper.speciesIndex.keys() ] ]

        for species in amountSpecies:
            if not species == '':
                volume = wrapper.speciesCompartment[wrapper.speciesIndex[species]].subs(wrapper.compartmentSymbols,
                                                                                        wrapper.compartmentVolume)
                simulated_x[:, wrapper.speciesIndex[species]] = simulated_x[:, wrapper.speciesIndex[species]] * volume
            pass

        adev = abs(simulated_x - test_x)
        adev[np.isnan(adev)] = True
        rdev = abs((simulated_x - test_x) / test_x)
        rdev[np.isnan(rdev)] = True
        if (not np.all(np.logical_or(adev < atol, rdev < rtol))):
            if (not np.all(adev < atol)):
                raise Exception('Absolute tolerance violated')

            if (not np.all(rdev < rtol)):
                raise Exception('Relative tolerance violated')

    except Exception as err:
        str = "Failed test " + testId + ": {0}".format(err)
        print(str)
        traceback.print_exc(10)
        logfile.write(str + '\n')
        return

def getTestStr(testId):
    testStr = str(testId)
    testStr = '0'*(5-len(testStr)) + testStr
    return testStr

'''
    currently failing due to https://github.com/symengine/symengine/issues/1444
    65
    121
    250
    253
    256
    259
    600
    803
'''
for testId in range(1,1781):
    with open("test.txt", "a") as logfile:
        runTest(getTestStr(testId), logfile)
