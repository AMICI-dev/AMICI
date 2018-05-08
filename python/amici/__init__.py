"""The AMICI Python module"""
import os

# determine package installation path, or, if used directly from git repository, get repository root
if os.path.exists(os.path.join(os.path.dirname(__file__), '..', '..', '.git')):
    amici_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
else:
    amici_path = os.path.dirname(__file__)

# if this file is inside the amici package, import swig interface, otherwise we are inside the git repository, then don't
dirname = os.path.dirname(__file__)
if os.path.isfile(os.path.join(dirname, 'amici.py')):
    #TODO: This should be more fine-grained
    from .amici import *

#TODO: This should be more fine-grained
from .sbml_import import *

