"""The AMICI Python module"""
import os

#TODO: work-around until all required files are included in the amici package
dirname = os.path.split(os.path.abspath(__file__))[0]
if os.path.isfile(os.path.join(dirname, 'amici.py')):
    #TODO: This should be more fine-grained
    from .amici import *

if os.path.isfile(os.path.join(dirname, 'paths.py')):
    from .paths import amici_path
else:
    amici_path = dirname + '/../../'

#TODO: This should be more fine-grained
from .sbml_import import *
