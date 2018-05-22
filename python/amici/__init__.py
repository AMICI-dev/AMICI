"""The AMICI Python module

The AMICI Python module provides functionality for importing SBML models and turning them into C++ Python extensions.

Getting started:
```
# creating a extension module for an SBML model:
import amici
amiSbml = amici.SbmlImporter('mymodel.sbml')
amiSbml.sbml2amici('modelName', 'outputDirectory')

# using the created module (set python path)
import modelName
help(modelName)
```

Attributes:
-----------
    amici_path: absolute root path of the amici repository
    amiciSwigPath: absolute path of the amici swig directory
    amiciSrcPath: absolute path of the amici source directory
    amiciModulePath: absolute root path of the amici module
"""

import os
import numpy as np


def setAmiciPaths():
    """Set module-wide paths to components of the AMICI distribution"""

    global amici_path, amiciSwigPath, amiciSrcPath, amiciModulePath
    # determine package installation path, or, if used directly from git
    # repository, get repository root
    if os.path.exists(os.path.join(os.path.dirname(__file__), '..', '..', '.git')):
        amici_path = os.path.abspath(os.path.join(
            os.path.dirname(__file__), '..', '..'))
    else:
        amici_path = os.path.dirname(__file__)

    amiciSwigPath = os.path.join(amici_path, 'swig')
    amiciSrcPath = os.path.join(amici_path, 'src')
    amiciModulePath = os.path.dirname(__file__)


# If this file is inside the amici package, import swig interface,
# otherwise we are inside the git repository, then don't
dirname = os.path.dirname(__file__)
if os.path.isfile(os.path.join(dirname, 'amici.py')):
    import amici
    from .amici import *

setAmiciPaths()

from .sbml_import import *


def runAmiciSimulation(model, solver, edata=None):
    if edata:
        edata = edata.get()
    rdata = amici.runAmiciSimulation(solver.get(), edata, model.get())
    return rdataToNumPyArrays(rdata)


def rdataToNumPyArrays(rdata):
    npReturnData = {}
    fieldNames = ['t', 'x', 'x0', 'sx', 'sx0', 'y', 'sigmay', 'sy', 'ssigmay', 
                  'z', 'rz', 'sigmaz', 'sz', 'srz', 'ssigmaz', 'sllh', 's2llh', 
                  'J', 'xdot', 'status', 'llh', 'chi2',
                  'newton_numlinsteps', 'newton_numsteps', 
                  'numsteps', 'numrhsevals', 'numerrtestfails', 'numnonlinsolvconvfails', 
                  'order', 'numstepsB', 'numrhsevalsB', 'numerrtestfailsB', 'numnonlinsolvconvfailsB']

    for field in fieldNames:
        npReturnData[field] = getFieldAsNumPyArray(rdata, field)

    return npReturnData


def getFieldAsNumPyArray(rdata, field):

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
                       'J': [rdata.nx, rdata.nx],
                       'xdot': [rdata.nx],
                       'newton_numlinsteps': [rdata.newton_maxsteps, 2],
                       'newton_numsteps': [1, 2],
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
    if field == 't':
        field = 'ts'

    attr = getattr(rdata, field)
    if field in fieldDimensions.keys():
        if len(fieldDimensions[field]) == 1:
            return np.array(attr)
        elif len(attr) == 0:
            return None
        else:
            return np.array(attr).reshape(fieldDimensions[field])
    else:
        return float(attr)
