#!/usr/bin/env python3

import sys
import amici
import unittest
import os
import numbers

class TestAmiciCPP(unittest.TestCase):
    '''
    TestCase class for testing SBML import and simulation from AMICI python interface
    '''

    expectedResultsFile = os.path.join(os.path.dirname(__file__),
                                       'cpputest', 'expectedResults.h5')

    def setUp(self):
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..',
                                        'python', 'examples',
                                        'example_presimulation'))
        from createModel import model
        model.name = 'test_model_presimulation_pysb'
        amici.pysb2amici(model,
                         model.name,
                         verbose=False,
                         observables=['pPROT_obs'],
                         constant_parameters=['DRUG_0', 'KIN_0'])
        sys.path.insert(0, model.name)
        import test_model_presimulation_pysb as modelModulePYSB
        self.model = modelModulePYSB.getModel()
        self.solver = self.model.getSolver()

    def tearDown(self):
        pass

    def runTest(self):
        self.testCopyConstructors()

    def testCopyConstructors(self):
        for obj in [self.model, self.solver]:
            for attr in dir(obj):
                with self.subTest(obj=obj, attr=attr):
                    if callable(getattr(obj, attr)):
                        # only use setter and getter functions
                        if attr.startswith('get') \
                                and 'set' + attr[3:] in dir(obj):

                            if attr.endswith('ById') \
                                    or attr.endswith('ByName'):
                                # ignore parameter/fixedParameter getters
                                # and setters
                                continue

                            # modify the value
                            val = getattr(obj, attr)()
                            if attr == 'getStabilityLimitFlag':
                                modval = val - 1.0
                            else:
                                try:
                                    modval = get_mod_val(val)
                                except Exception:
                                    continue
                                getattr(obj, 'set' + attr[3:])(
                                    modval
                                )

                            # clone
                            obj_clone = obj.clone()

                            # check if change is also present in clone
                            self.assertEqual(
                                getattr(obj, attr)(),
                                getattr(obj_clone, attr)()
                            )
                        else:
                            pass
                    elif not attr.startswith('__') and attr != 'this':

                        # modify the value
                        try:
                            modval = get_mod_val(getattr(obj, attr))
                        except Exception:
                            continue

                        try:
                            setattr(obj, attr, modval)
                        except AttributeError:
                            # some attributes cant be set
                            pass

                        # clone
                        obj_clone = obj.clone()

                        # check if change is also present in clone
                        self.assertEqual(
                            getattr(obj, attr),
                            getattr(obj_clone, attr)
                        )

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestAmiciCPP())
    unittest.main()


def get_mod_val(val):
    if isinstance(val, bool):
        modval = not val
    elif isinstance(val, numbers.Number):
        modval = val + 1
    else:
        raise Exception('Cannot modify value')
    return modval