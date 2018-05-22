#!/usr/bin/env python3
"""
Example using [model_steadystate_scaled.sbml] model to demonstrate and test SBML import and Amici Python interface.

To use python-generated model in matlab, run:
```
modelName = '';
modelDir = '';
amimodel.compileAndLinkModel(modelName, modelDir, [], [], [], []);
amimodel.generateMatlabWrapper(3, 6, 8, 1, 0, 0, [], [ modelDir '/simulate_test.m'], modelName, 'lin', 1, 1);
```
"""

import amici

import os
import sys
import numpy as np
import matplotlib.pyplot as plt


os.chdir(os.path.dirname(__file__))

def createModule():
    """Create Python module from SBML model"""
    sbmlImporter = amici.SbmlImporter('model_steadystate_scaled.sbml')
    sbml = sbmlImporter.sbml
    
    observables = amici.assignmentRules2observables(sbml, filter=lambda variableId: 
                                                    variableId.startswith('observable_') and not variableId.endswith('_sigma'))
    
    print(observables)
    
    sbmlImporter.sbml2amici('test', 'test', 
                            observables=observables,
                            constantParameters=['k4'],
                            sigmas={'observable_x1withsigma': 'observable_x1withsigma_sigma'})

def plotStateTrajectories(rdata):
    for ix in range(rdata['x'].shape[1]):
        plt.plot(rdata['t'], rdata['x'][:, ix], label='$x_%d$' % ix)
        plt.xlabel('$t$ (s)')
        plt.ylabel('$x_i(t)$ (mmol/ml)')
        plt.legend()
        plt.title('State trajectories')
    plt.show()
    
def plotObservableTrajectories(rdata):
    for iy in range(rdata['y'].shape[1]):
        plt.plot(rdata['t'], rdata['y'][:, iy], label='$y_%d$' % iy)
        plt.xlabel('$t$ (s)')
        plt.ylabel('$y_i(t)$ (AU)')
        plt.legend()
        plt.title('Observables')
    
    plt.show()

createModule()

sys.path.insert(0, 'test')
import test as modelModule

model = modelModule.getModel()
model.setTimepoints(amici.DoubleVector(np.linspace(0, 60, 60))) 
solver = model.getSolver()
rdata = amici.runAmiciSimulation(model, solver)

print(rdata)

plotStateTrajectories(rdata)
plotObservableTrajectories(rdata)

