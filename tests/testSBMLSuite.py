#!/usr/bin/env python3
"""Run SBML Test Suite and verify simulation results [https://github.com/sbmlteam/sbml-test-suite/releases]"""
import traceback
import os
import sys
import importlib
import numpy as np
import amici

# directory with sbml semantic test cases
test_path = os.path.join(os.path.dirname(__file__), 'sbml-test-suite', 'cases', 'semantic')

def runTest(testId, logfile):
    try:
        current_test_path = os.path.join(test_path, testId)

        # results
        resultsFile = os.path.join(current_test_path, testId + '-results.csv')
        results = np.genfromtxt(resultsFile, delimiter=',')

        # model
        sbmlFile = findModelFile(current_test_path, testId)

        wrapper = amici.SbmlImporter(sbmlFile)

        modelDir = os.path.join(os.path.dirname(__file__),'SBMLTestModels',testId)
        if not os.path.exists(modelDir):
            os.makedirs(modelDir)
        wrapper.sbml2amici('SBMLTest' + testId, output_dir=modelDir)

        # settings
        settings = readSettingsFile(current_test_path, testId)
        ts = np.linspace(float(settings['start']), float(settings['start']) + float(settings['duration']),
                         int(settings['steps']) + 1)
        atol = float(settings['absolute'])
        rtol = float(settings['relative'])

        sys.path.insert(0, wrapper.modelPath)
        mod = importlib.import_module(wrapper.modelName)

        model = mod.getModel()
        model.setTimepoints(mod.amici.DoubleVector(ts))
        solver = model.getSolver()
        solver.setMaxSteps(int(1e6))
        solver.setRelativeTolerance(rtol / 1000.0)
        solver.setAbsoluteTolerance(atol / 1000.0)
        rdata = amici.runAmiciSimulation(model, solver)
        amountSpecies = settings['amount'].replace(' ', '').replace('\n', '').split(',')
        simulated_x = rdata['x']
        test_x = results[1:, [1+ wrapper.speciesIndex[variable]  for variable in settings['variables'].replace(' ', '').replace('\n', '').split(',') if variable in wrapper.speciesIndex.keys() ] ]

        for species in amountSpecies:
            if not species == '':
                volume = wrapper.speciesCompartment[wrapper.speciesIndex[species]].subs(wrapper.compartmentSymbols,
                                                                                        wrapper.compartmentVolume)
                simulated_x[:, wrapper.speciesIndex[species]] = simulated_x[:, wrapper.speciesIndex[species]] * volume

        adev = abs(simulated_x - test_x)
        adev[np.isnan(adev)] = True
        rdev = abs((simulated_x - test_x) / test_x)
        rdev[np.isnan(rdev)] = True
        if not np.all(np.logical_or(adev < atol, rdev < rtol)):
            if (not np.all(adev < atol)):
                raise Exception('Absolute tolerance violated')

            if (not np.all(rdev < rtol)):
                raise Exception('Relative tolerance violated')

    except amici.SBMLException as err:
        print("Did not run test " + testId + ": {0}".format(err))

    except Exception as err:
        str = "Failed test " + testId + ": {0}".format(err)
        traceback.print_exc(10)
        logfile.write(str + '\n')
        return

def findModelFile(current_test_path, testId):
    """Find model file for the given test (guess filename extension)"""
    sbmlFile = os.path.join(current_test_path, testId + '-sbml-l3v2.xml')

    # fallback l3v1
    if not os.path.isfile(sbmlFile):
        sbmlFile = os.path.join(current_test_path, testId + '-sbml-l3v1.xml')

    # fallback l2v5
    if not os.path.isfile(sbmlFile):
        sbmlFile = os.path.join(current_test_path, testId + '-sbml-l2v5.xml')
        
    return sbmlFile

def readSettingsFile(current_test_path, testId):
    """Read settings for the given test"""
    settingsFile = os.path.join(current_test_path, testId + '-settings.txt')
    settings = {}
    with open(settingsFile) as f:
        for line in f:
            if not line == '\n':
                (key, val) = line.split(':')
                settings[key] = val
    return settings


def getTestStr(testId):
    testStr = str(testId)
    testStr = '0'*(5-len(testStr)) + testStr
    return testStr

def main():
    for testId in range(1,1781):
        with open("test.txt", "w") as logfile:
            runTest(getTestStr(testId), logfile)

if __name__ == '__main__':
    main()