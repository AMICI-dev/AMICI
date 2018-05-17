#!/usr/bin/env python3
import amici

import os
import sys
import numpy as np

os.chdir(os.path.dirname(__file__))

def createModule():
    sbmlImporter = amici.SbmlImporter('model_steadystate_scaled.sbml')
    sbml = sbmlImporter.sbml
    
    observables = amici.assignmentRules2observables(sbml, filter=lambda variableId: variableId.startswith('observable_'))
    
    print(observables)
    
    sbmlImporter.sbml2amici('test', 'test', observables=observables)


sys.path.insert(0, 'test')
import test as modelModule

model = modelModule.getModel()
model.setTimepoints(modelModule.amici.DoubleVector([0.0, 1, 2, 100])) 
solver = model.getSolver()
rdata = amici.runAmiciSimulation(solver, model)

print(dir(rdata))