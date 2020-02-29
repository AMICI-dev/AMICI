#!/usr/bin/env python3

import sys
import amici
import unittest
import os
import numpy as np
import pysb
from amici.pysb_import import pysb2amici
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
        outdir = model.name
        pysb2amici(model,
                   outdir,
                   verbose=False,
                   observables=['pPROT_obs'],
                   constant_parameters=['DRUG_0', 'KIN_0'])
        sys.path.insert(0, outdir)
        modelModulePYSB = importlib.import_module(outdir)

        self.model = modelModulePYSB.getModel()
        self.model.setReinitializeFixedParameterInitialStates(True)

        self.solver = self.model.getSolver()
        self.solver.setSensitivityOrder(amici.SensitivityOrder_first)
        self.solver.setSensitivityMethod(amici.SensitivityMethod_forward)

        self.edata = get_data(self.model)
        self.edata.t_presim = 2
        self.edata.fixedParameters = [10, 2]
        self.edata.fixedParametersPresimulation = [3, 2]
        self.edata.fixedParametersPreequilibration = [3, 0]
        self.edata.setTimepoints([1, 5])

        self.edata_preeq = amici.ExpData(self.edata)
        self.edata_preeq.t_presim = 0
        self.edata_preeq.setTimepoints([np.infty])
        self.edata_preeq.fixedParameters = \
            self.edata.fixedParametersPreequilibration
        self.edata_preeq.fixedParametersPresimulation = ()
        self.edata_preeq.fixedParametersPreequilibration = ()

        self.edata_presim = amici.ExpData(self.edata)
        self.edata_presim.t_presim = 0
        self.edata_presim.setTimepoints([self.edata.t_presim])
        self.edata_presim.fixedParameters = \
            self.edata.fixedParametersPresimulation
        self.edata_presim.fixedParametersPresimulation = ()
        self.edata_presim.fixedParametersPreequilibration = ()

        self.edata_sim = amici.ExpData(self.edata)
        self.edata_sim.t_presim = 0
        self.edata_sim.setTimepoints(self.edata.getTimepoints())
        self.edata_sim.fixedParameters = \
            self.edata.fixedParameters
        self.edata_sim.fixedParametersPresimulation = ()
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

            # manual reinitialization + presimulation
            x0 = rdata_preeq['x'][0, :]
            x0[1] = self.edata_presim.fixedParameters[0]
            x0[2] = self.edata_presim.fixedParameters[1]
            sx0 = rdata_preeq['sx'][0, :, :]
            sx0[:, 1] = 0
            sx0[:, 2] = 0
            self.model.setInitialStates(x0)
            self.model.setInitialStateSensitivities(
                sx0.flatten()
            )
            rdata_presim = amici.runAmiciSimulation(
                self.model, self.solver, self.edata_presim
            )

            # manual reinitialization + simulation
            x0 = rdata_presim['x'][0, :]
            x0[1] = self.edata_sim.fixedParameters[0]
            x0[2] = self.edata_sim.fixedParameters[1]
            sx0 = rdata_presim['sx'][0, :, :]
            sx0[:, 1] = 0
            sx0[:, 2] = 0
            self.model.setInitialStates(x0)
            self.model.setInitialStateSensitivities(
                sx0.flatten()
            )
            rdata_sim = amici.runAmiciSimulation(
                self.model, self.solver, self.edata_sim
            )

            for variable in ['x', 'sx']:
                with self.subTest(pscale=pscale, plist=plist,
                                  variable=variable):
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

