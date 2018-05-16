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

def setAmiciPaths():
    """Set module-wide paths to components of the AMICI distribution""" 
    
    global amici_path, amiciSwigPath, amiciSrcPath, amiciModulePath
    # determine package installation path, or, if used directly from git repository, get repository root
    if os.path.exists(os.path.join(os.path.dirname(__file__), '..', '..', '.git')):
        amici_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    else:
        amici_path = os.path.dirname(__file__)
    
    amiciSwigPath = os.path.join(amici_path, 'swig')
    amiciSrcPath = os.path.join(amici_path, 'src')
    amiciModulePath = os.path.dirname(__file__)

# If this file is inside the amici package, import swig interface, 
# otherwise we are inside the git repository, then don't
dirname = os.path.dirname(__file__)
if os.path.isfile(os.path.join(dirname, 'amici.py')):
    from .amici import *

setAmiciPaths()

from .sbml_import import *

