#!/usr/bin/env python3

import sys
import amici
import unittest
import os
import platform
import importlib
import copy
import numpy as np
from pysb.simulator import ScipyOdeSimulator

import pysb.examples

class TestAmiciPYSBModel(unittest.TestCase):
    '''
    TestCase class for testing SBML import and simulation from AMICI python
    interface
    '''

    expectedResultsFile = os.path.join(os.path.dirname(__file__),
                                       'cpputest', 'expectedResults.h5')

    def setUp(self):
        self.default_path = copy.copy(sys.path)
        self.resetdir = os.getcwd()

        if os.path.dirname(__file__) != '':
            os.chdir(os.path.dirname(__file__))

    def tearDown(self):
        os.chdir(self.resetdir)
        sys.path = self.default_path

    def runTest(self):
        self.test_compare_to_sbml_import()
        self.test_compare_to_pysb_simulation()

    def test_compare_to_sbml_import(self):
        constant_parameters = ['DRUG_0', 'KIN_0']

        # -------------- PYSB -----------------
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
                         constant_parameters=constant_parameters)
        sys.path.insert(0, model.name)
        import test_model_presimulation_pysb as modelModulePYSB
        model_pysb = modelModulePYSB.getModel()

        edata = get_data(model_pysb)
        rdata_pysb = get_results(model_pysb, edata)

        # -------------- SBML -----------------

        sbmlFile = os.path.join(os.path.dirname(__file__), '..', 'python',
                                'examples', 'example_presimulation',
                                'model_presimulation.xml')

        sbmlImporter = amici.SbmlImporter(sbmlFile)

        observables = amici.assignmentRules2observables(
            sbmlImporter.sbml,  # the libsbml model object
            filter_function=lambda variable: variable.getName() == 'pPROT_obs'
        )
        outdir = 'test_model_presimulation_sbml'
        sbmlImporter.sbml2amici('test_model_presimulation_sbml',
                                outdir,
                                verbose=False,
                                observables=observables,
                                constantParameters=constant_parameters)
        sys.path.insert(0, outdir)
        import test_model_presimulation_sbml as modelModuleSBML
        model_sbml = modelModuleSBML.getModel()

        rdata_sbml = get_results(model_sbml, edata)

        for field in rdata_pysb:
            if field not in ['ptr', 't_steadystate', 'numsteps',
                             'newton_numsteps', 'numrhsevals',
                             'numerrtestfails', 'order', 'J', 'xdot',
                             'wrms_steadystate', 'newton_cpu_time',
                             'cpu_time', 'cpu_timeB', ]:
                with self.subTest(field=field):
                    if rdata_pysb[field] is None:
                        self.assertIsNone(
                            rdata_sbml[field],
                        )
                    elif rdata_sbml[field] is None:
                        self.assertIsNone(
                            rdata_pysb[field],
                        )
                    else:
                        self.assertTrue(np.isclose(
                            rdata_sbml[field],
                            rdata_pysb[field],
                            atol=1e-6, rtol=1e-6
                        ).all())

    def test_compare_to_pysb_simulation(self):

        sys.path.insert(0, os.path.dirname(pysb.examples.__file__))

        pysb_models = [
            'tyson_oscillator', 'robertson', 'expression_observables',
            'bax_pore_sequential', 'bax_pore', 'bngwiki_egfr_simple',
            'bngwiki_enzymatic_cycle_mm', 'bngwiki_simple', 'earm_1_0',
            'earm_1_3', 'move_connected', 'michment', 'kinase_cascade',
            'hello_pysb', 'fricker_2010_apoptosis', 'explicit',
            'fixed_initial',
        ]

        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..',
                                        'tests', 'pysb_test_models'))

        custom_models = [
            'bngwiki_egfr_simple_deletemolecules',
        ]

        atol = 1e-8
        rtol = 1e-8

        for example in pysb_models + custom_models:
            with self.subTest(example=example):

                if example == 'earm_1_3' \
                        and platform.sys.version_info[0] == 3 \
                        and platform.sys.version_info[1] < 7:
                    continue

                # load example

                pysb.SelfExporter.cleanup()  # reset pysb
                pysb.SelfExporter.do_export = True

                module = importlib.import_module(example)
                pysb_model = module.model
                pysb_model.name = pysb_model.name.replace('pysb.examples.', '')
                # avoid naming clash for custom pysb models
                pysb_model.name += '_amici'

                # pysb part

                tspan = np.linspace(0, 100, 101)
                sim = ScipyOdeSimulator(
                    pysb_model,
                    tspan=tspan,
                    integrator_options={'rtol': rtol, 'atol': atol},
                    compiler='python'
                )
                pysb_simres = sim.run()

                # amici part

                outdir = pysb_model.name

                if pysb_model.name in ['move_connected_amici']:
                    self.assertRaises(
                        Exception,
                        amici.pysb2amici,
                        *[pysb_model, outdir],
                        **{'verbose': False, 'compute_conservation_laws': True}
                    )
                    compute_conservation_laws = False
                else:
                    compute_conservation_laws = True

                amici.pysb2amici(
                    pysb_model,
                    outdir,
                    verbose=False,
                    compute_conservation_laws=compute_conservation_laws
                )
                sys.path.insert(0, outdir)

                amici_model_module = importlib.import_module(pysb_model.name)

                model_pysb = amici_model_module.getModel()

                model_pysb.setTimepoints(tspan)

                solver = model_pysb.getSolver()
                solver.setMaxSteps(int(1e5))
                solver.setAbsoluteTolerance(atol)
                solver.setRelativeTolerance(rtol)
                rdata = amici.runAmiciSimulation(model_pysb, solver)

                # check agreement of species simulation

                self.assertTrue(np.isclose(
                    rdata['x'],
                    pysb_simres.species,
                    1e-4, 1e-4
                ).all())


def get_data(model):
    solver = model.getSolver()
    model.setTimepoints(np.linspace(0, 60, 61))
    model.setReinitializeFixedParameterInitialStates(True)
    model.setSteadyStateSensitivityMode(
        amici.SteadyStateSensitivityMode_simulationFSA
    )

    rdata = amici.runAmiciSimulation(model, solver)
    edata = amici.ExpData(rdata, 0.1, 0.0)
    edata.fixedParametersPreequilibration = [3, 0]
    edata.fixedParameters = [10, 2]
    edata.fixedParametersPresimulation = [10, 2]
    edata.fixedParametersPreequilibration = [3, 0]
    return edata


def get_results(model, edata):
    solver = model.getSolver()
    solver.setSensitivityOrder(1)
    model.setTimepoints(np.linspace(0, 60, 61))
    model.setReinitializeFixedParameterInitialStates(True)
    model.setSteadyStateSensitivityMode(
        amici.SteadyStateSensitivityMode_simulationFSA
    )
    return amici.runAmiciSimulation(model, solver, edata)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestAmiciPYSBModel())
    unittest.main()
