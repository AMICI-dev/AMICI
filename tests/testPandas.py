#!/usr/bin/env python3

"""Tests for ``amici.pandas``"""

import sys
import amici
import unittest
import os
import copy
import numpy as np
import itertools

from petab.sbml import assignment_rules_to_dict


class TestAmiciPandasImportExport(unittest.TestCase):
    """
    TestCase class for testing csv import using pandas
    """

    def setUp(self):
        self.default_path = copy.copy(sys.path)
        self.resetdir = os.getcwd()

        if os.path.dirname(__file__) != '':
            os.chdir(os.path.dirname(__file__))

        sbmlFile = os.path.join(os.path.dirname(__file__), '..', 'python',
                                'examples', 'example_presimulation',
                                'model_presimulation.xml')

        sbmlImporter = amici.SbmlImporter(sbmlFile)

        constantParameters = ['DRUG_0', 'KIN_0']

        observables = assignment_rules_to_dict(
            sbmlImporter.sbml,  # the libsbml model object
            filter_function=lambda variable: variable.getName() == 'pPROT',
            remove=True
        )
        outdir = 'test_model_presimulation'
        sbmlImporter.sbml2amici('test_model_presimulation',
                                outdir,
                                verbose=False,
                                observables=observables,
                                constantParameters=constantParameters)
        sys.path.insert(0, outdir)
        import test_model_presimulation as modelModule
        self.model = modelModule.getModel()
        self.model.setTimepoints(np.linspace(0, 60, 61))
        self.solver = self.model.getSolver()
        rdata = amici.runAmiciSimulation(self.model, self.solver)
        self.edata = [amici.ExpData(rdata, 0.01, 0)]
        # test copy constructor
        self.edata_copy = amici.ExpData(self.edata[0])

    def tearDown(self):
        os.chdir(self.resetdir)
        sys.path = self.default_path

    def runTest(self):
        self.tests_presimulation()

    def tests_presimulation(self):
        self.model.getFixedParameterNames()
        combos = itertools.product(
            [(10, 5), (5, 10), ()],
            repeat=3
        )
        cases = dict()
        for icombo, combo in enumerate(combos):
            cases[f'{icombo}'] = {
                'fixedParameters': combo[0],
                'fixedParametersPreequilibration': combo[1],
                'fixedParametersPresimulation': combo[2],
            }

        for case in cases:
            with self.subTest(**cases[case]):
                for fp in cases[case]:
                    setattr(self.edata[0], fp, cases[case][fp])

                df_edata = amici.getDataObservablesAsDataFrame(
                    self.model,
                    self.edata
                )
                edata_reconstructed = amici.getEdataFromDataFrame(
                    self.model,
                    df_edata
                )

                for fp in ['fixedParameters', 'fixedParametersPreequilibration',
                           'fixedParametersPresimulation']:

                    if fp != 'fixedParameters' or cases[case][fp] is not ():
                        self.assertTupleEqual(
                            getattr(self.edata[0], fp),
                            getattr(edata_reconstructed[0], fp),
                        )

                        self.assertTupleEqual(
                            cases[case][fp],
                            getattr(edata_reconstructed[0], fp),
                        )

                        self.assertTupleEqual(
                            getattr(self.edata[0], fp),
                            cases[case][fp],
                        )
                    else:
                        self.assertTupleEqual(
                            self.model.getFixedParameters(),
                            getattr(edata_reconstructed[0], fp),
                        )

                        self.assertTupleEqual(
                            self.model.getFixedParameters(),
                            getattr(edata_reconstructed[0], fp),
                        )

                        self.assertTupleEqual(
                            getattr(self.edata[0], fp),
                            cases[case][fp],
                        )


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestAmiciPandasImportExport())
    unittest.main()
