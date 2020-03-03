#!/usr/bin/env python3

"""Test SWIG interface using python/examples/test_model_presimulation_pysb

Test getters, setters, etc.
"""

import sys
import amici
import unittest
import os
import numbers
import pysb
import importlib
import copy

from amici.pysb_import import pysb2amici

class TestAmiciCPP(unittest.TestCase):
    """
    TestCase class for testing cpp API through swig
    """

    def setUp(self):
        self.resetdir = os.getcwd()
        self.default_path = copy.copy(sys.path)

        pysb.SelfExporter.cleanup()  # reset pysb
        pysb.SelfExporter.do_export = True

        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..',
                                        'python', 'examples',
                                        'example_presimulation'))
        if 'createModelPresimulation' in sys.modules:
            importlib.reload(sys.modules['createModelPresimulation'])
            model_module = sys.modules['createModelPresimulation']
        else:
            model_module = importlib.import_module('createModelPresimulation')
        model = copy.deepcopy(model_module.model)
        model.name = 'test_model_presimulation_pysb'
        pysb2amici(model,
                   model.name,
                   verbose=False,
                   observables=['pPROT_obs'],
                   constant_parameters=['DRUG_0', 'KIN_0'])
        sys.path.insert(0, model.name)
        import test_model_presimulation_pysb as modelModulePYSB
        self.model = modelModulePYSB.getModel()
        self.solver = self.model.getSolver()

    def tearDown(self):
        os.chdir(self.resetdir)
        sys.path = self.default_path

    def runTest(self):
        self.test_copy_constructors()
        self.test_version_number()

    def test_version_number(self):
        self.assertEqual(self.model.getAmiciVersion(), amici.__version__)
        self.assertEqual(self.model.getAmiciCommit(), amici.__commit__)

    def test_copy_constructors(self):
        # TODO: expand this to serialization
        for obj in [self.model, self.solver]:
            for attr in dir(obj):
                if attr.startswith('__') \
                        or attr == 'this' \
                        or is_callable_but_not_getter(obj, attr):
                    continue

                with self.subTest(obj=obj, attr=attr):
                    # objects will be initialized with default values so we
                    # can check if the clone routine does the right thing by
                    # modifying the value from default, cloning and checking
                    # if the change was carried over to the clone
                    val = get_val(obj, attr)

                    try:
                        modval = get_mod_val(val, attr)
                    except ValueError:
                        # happens for everything that is not bool or scalar
                        continue

                    try:
                        set_val(obj, attr, modval)
                    except AttributeError:
                        # some attributes cannot be set
                        continue

                    obj_clone = obj.clone()

                    self.assertEqual(
                        get_val(obj, attr),
                        get_val(obj_clone, attr)
                    )


def is_callable_but_not_getter(obj, attr):
    if not callable(getattr(obj, attr)):
        return False

    if attr.startswith('get'):
        return \
            'set' + attr[3:] not in dir(obj) \
            or attr.endswith('ById') \
            or attr.endswith('ByName')
    else:
        return True


def get_val(obj, attr):
    if callable(getattr(obj, attr)):
        return getattr(obj, attr)()
    else:
        return getattr(obj, attr)


def get_mod_val(val, attr):
    if attr == 'getStabilityLimitFlag':
        return val - 1
    elif isinstance(val, bool):
        return not val
    elif isinstance(val, numbers.Number):
        return val + 1

    raise ValueError('Cannot modify value')


def set_val(obj, attr, val):
    if callable(getattr(obj, attr)):
        getattr(obj, 'set' + attr[3:])(
            val
        )
    else:
        setattr(obj, attr, val)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestAmiciCPP())
    unittest.main()
