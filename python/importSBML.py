#!/usr/bin/env python3

from model import Model
import os
import sys
import importlib
sys.path.insert(0, "../build/swig/python/")
import amici

dirname = os.path.split(os.path.abspath(__file__))[0]
amici_path = os.path.split(dirname)[0]
test_path = os.path.join(amici_path,'tests','sbml-semantic-test-cases','cases','semantic')

def runTest(testId):
    try:
        sbmlfile = os.path.join(test_path, testId, testId + '-sbml-l3v2.xml')
        wrapper = Model(sbmlfile,'SBMLTest' + testId)
        wrapper.wrapModel()
        sys.path.insert(0,os.path.join(wrapper.model_path,'build','swig'))

        settings

        mod = importlib.import_module(wrapper.modelname)
        model = mod.getModel()
        solver = model.getSolver()
        rdata = amici.runAmiciSimulation(solver.get(),None,model.get())
        print(rdata)




    except AssertionError:
        None

def getTestStr(testId):
    testStr = str(testId)
    testStr = '0'*(5-len(testStr)) + testStr
    return testStr


for testId in range(1,1781):
    runTest(getTestStr(testId))

#model = Model('/Users/F.Froehlich/Downloads/Speedy_v3_r403445_v1.sbml','speedy')
#model.wrapModel()