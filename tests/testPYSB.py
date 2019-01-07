#!/usr/bin/env python3

import sys
import amici
import unittest
import os
import importlib
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
        self.resetdir = os.getcwd()
        if os.path.dirname(__file__) != '':
            os.chdir(os.path.dirname(__file__))

    def tearDown(self):
        os.chdir(self.resetdir)

    def runTest(self):
        self.test_compare_to_sbml_import()
        self.test_compare_to_pysb_simulation()

    def test_compare_to_sbml_import(self):
        constant_parameters = ['DRUG_0', 'KIN_0']

        # -------------- PYSB -----------------

        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..',
                                        'tests', 'pysb_test_models'))
        from createModel import model
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
                             'numerrtestfails', 'order', 'J', 'xdot']:
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
        ]

        pysb_models = [
            'move_connected', 'fricker_2010_apoptosis',
        ]

        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..',
                                        'tests', 'pysb_test_models'))

        custom_models = [
            'bngwiki_egfr_simple_deletemolecules',
        ]

        for example in pysb_models:
            pysb.core.SelfExporter.cleanup()
            module = importlib.import_module(example)
            pysb_model = module.model
            pysb_model.name = pysb_model.name.replace('pysb.examples.', '')
            # avoid naming clash for custom pysb models
            pysb_model.name += '_amici'
            with self.subTest(example=pysb_model.name):
                # pysb part

                tspan = np.linspace(0, 100, 101)
                sim = ScipyOdeSimulator(
                    pysb_model,
                    tspan=tspan,
                    integrator_options={'rtol': 1e-8, 'atol': 1e-8},
                    compiler='python'
                )
                pysb_simres = sim.run()

                # amici part

                outdir = pysb_model.name
                amici.pysb2amici(pysb_model,
                                 outdir,
                                 verbose=False,
                                 compute_conservation_laws=True)
                sys.path.insert(0, outdir)

                amici_model_module = importlib.import_module(pysb_model.name)

                model_pysb = amici_model_module.getModel()

                model_pysb.setTimepoints(tspan)

                solver = model_pysb.getSolver()
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
