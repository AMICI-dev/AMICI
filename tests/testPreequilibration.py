#!/usr/bin/env python3

import sys
import amici
import unittest
import os
import numpy as np
import pysb
import importlib
import copy
import itertools
from testPYSB import get_data


class TestAmiciPreequilibration(unittest.TestCase):
    '''
    TestCase class for testing preequilibration
    '''

    expectedResultsFile = os.path.join(os.path.dirname(__file__),
                                       'cpputest', 'expectedResults.h5')

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
        amici.pysb2amici(model,
                         model.name,
                         verbose=False,
                         observables=['pPROT_obs'],
                         constant_parameters=['DRUG_0', 'KIN_0'])
        sys.path.insert(0, model.name)
        import test_model_presimulation_pysb as modelModulePYSB

        self.model = modelModulePYSB.getModel()

        self.solver = self.model.getSolver()
        self.solver.setSensitivityOrder(amici.SensitivityOrder_first)
        self.solver.setSensitivityMethod(amici.SensitivityMethod_forward)

        self.edata = get_data(self.model)
        self.edata.fixedParametersPresimulation = ()

        self.edata_preeq = amici.ExpData(self.edata)
        self.edata_preeq.setTimepoints([0])

        self.edata_sim = amici.ExpData(self.edata)
        self.edata_sim.fixedParametersPreequilibration = ()

        self.pscales = [
            amici.ParameterScaling_log10, amici.ParameterScaling_ln,
            amici.ParameterScaling_none,
            amici.parameterScalingFromIntVector([
                amici.ParameterScaling_log10, amici.ParameterScaling_ln,
                amici.ParameterScaling_none, amici.ParameterScaling_log10,
                amici.ParameterScaling_ln, amici.ParameterScaling_none
            ])
        ]

        self.plists = [
            [3, 1, 2, 4], [0, 1, 2, 3, 4, 5], [5, 3, 2, 0, 4, 1],
            [1, 2, 3, 4, 5], [1, 1, 1],
        ]

    def tearDown(self):
        os.chdir(self.resetdir)
        sys.path = self.default_path

    def runTest(self):
        self.test_manual_preequilibration()
        self.test_parameter_reordering()

    def test_manual_preequilibration(self):

        settings = itertools.product(self.pscales, self.plists)

        for pscale, plist in settings:

            self.model.setInitialStates([])
            self.model.setInitialStateSensitivities([])
            self.model.setParameterList(plist)
            self.model.setParameterScale(pscale)

            # combined
            rdata_auto = amici.runAmiciSimulation(
                self.model, self.solver, self.edata
            )

            # manual preqquilibration
            rdata_preeq = amici.runAmiciSimulation(
                self.model, self.solver, self.edata_preeq
            )

            # manual reinitialization + simulation
            self.model.setInitialStates(rdata_preeq['x0'])
            self.model.setInitialStateSensitivities(
                rdata_preeq['sx0'].flatten()
            )
            rdata_sim = amici.runAmiciSimulation(
                self.model, self.solver, self.edata_sim
            )

            for variable in ['x', 'sx']:
                with self.subTest(pscale=pscale, plist=plist,
                                  variable=variable):
                    self.assertTrue(np.isclose(
                        rdata_auto[variable][0, :],
                        rdata_preeq[variable][0, :],
                        1e-6, 1e-6
                    ).all())
                    self.assertTrue(np.isclose(
                        rdata_auto[variable],
                        rdata_sim[variable],
                        1e-6, 1e-6
                    ).all())

    def test_parameter_reordering(self):
        rdata_ordered = amici.runAmiciSimulation(
            self.model, self.solver, self.edata
        )

        for plist in self.plists:
            with self.subTest(plist=plist):

                self.model.setParameterList(plist)
                rdata_reordered = amici.runAmiciSimulation(
                    self.model, self.solver, self.edata
                )

                for ip, p_index in enumerate(plist):
                    self.assertTrue(np.isclose(
                        rdata_ordered['sx'][:, p_index, :],
                        rdata_reordered['sx'][:, ip, :],
                        1e-6, 1e-6
                    ).all())

    def test_parameter_in_expdata(self):
        rdata = amici.runAmiciSimulation(
            self.model, self.solver, self.edata
        )

        # set ExpData plist
        self.edata.plist = self.model.getParameterList()
        # perturb model parameter list
        self.model.setParameterList([
            i for i in reversed(self.model.getParameterList())
        ])

        # set ExpData parameters
        self.edata.parameters = self.model.getParameters()
        # perturb model parameters
        self.model.setParameters(tuple(
            p * 2 for p in self.model.getParameters()
        ))

        # set ExpData pscale
        self.edata.pscale = self.model.getParameterScale()
        # perturb model pscale, needs to be done after getting parameters,
        # otherwise we will mess up parameter value
        self.model.setParameterScale(amici.parameterScalingFromIntVector([
            amici.ParameterScaling_log10
            if scaling == amici.ParameterScaling_none
            else amici.ParameterScaling_none
            for scaling in self.model.getParameterScale()
        ]))

        self.edata.x0 = rdata['x_ss']
        self.edata.sx0 = rdata['sx_ss'].flatten()

        # perturb model initial states
        self.model.setInitialStates(rdata['x_ss'] * 4)
        self.model.setInitialStateSensitivities(rdata['sx_ss'].flatten() / 2)

        rdata_edata = amici.runAmiciSimulation(
            self.model, self.solver, self.edata
        )
        for variable in ['x', 'sx']:
            with self.subTest(variable=variable):
                self.assertTrue(np.isclose(
                    rdata[variable][0, :],
                    rdata_edata[variable][0, :],
                    1e-6, 1e-6
                ).all())


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestAmiciPreequilibration())
    unittest.main()

