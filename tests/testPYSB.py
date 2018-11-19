#!/usr/bin/env python3

import sys
import amici
import unittest
import os
import importlib
import numpy as np
from pysb.simulator import ScipyOdeSimulator

from pysb.examples import tyson_oscillator, robertson, \
    expression_observables, bax_pore_sequential, bax_pore, \
    bngwiki_egfr_simple


class TestAmiciPYSBModel(unittest.TestCase):
    '''
    TestCase class for testing SBML import and simulation from AMICI python interface
    '''

    expectedResultsFile = os.path.join(os.path.dirname(__file__),
                                       'cpputest','expectedResults.h5')

    def setUp(self):
        self.resetdir = os.getcwd()
        if os.path.dirname(__file__) != '':
            os.chdir(os.path.dirname(__file__))

    def tearDown(self):
        os.chdir(self.resetdir)

    def runTest(self):
        self.compare_to_sbml_import()
        self.compare_to_pysb_simulation()

    def compare_to_sbml_import(self):
        constant_parameters = ['DRUG_0', 'KIN_0']

        # -------------- PYSB -----------------

        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..',
                                        'python', 'examples',
                                        'example_presimulation'))
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
            if field is not 'ptr':
                if rdata_pysb[field] is None:
                    self.assertIsNone(
                        rdata_sbml[field],
                        msg=f'checked field was {field} in sbml results'
                    )
                elif rdata_sbml[field] is None:
                    self.assertIsNone(
                        rdata_pysb[field],
                        msg=f'checked field was {field} in pysb results'
                    )
                else:
                    self.assertTrue(np.isclose(
                        rdata_sbml[field],
                        rdata_pysb[field],
                        atol=1e-8, rtol=1e-8
                    ).all(), msg=f'disagreement in field {field}')

    def compare_to_pysb_simulation(self):
        examples = [tyson_oscillator.model, robertson.model,
                  expression_observables.model,
                  bax_pore_sequential.model, bax_pore.model,
                  bngwiki_egfr_simple.model]
        for example in examples:
            example.name = example.name.replace('pysb.examples.','')
            with self.subTest(example=example.name):
                ## pysb part

                tspan = np.linspace(0, 100, 101)
                sim = ScipyOdeSimulator(
                    example,
                    tspan=tspan,
                    integrator_options={'rtol': 1e-8, 'atol': 1e-8},
                    compiler='python'
                )
                pysb_simres = sim.run()

                ## amici part

                amici.pysb2amici(example,
                                 example.name,
                                 verbose=False)
                sys.path.insert(0, example.name)
                amici_model_module = importlib.import_module(example.name)
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

    rdata = amici.runAmiciSimulation(model, solver)
    edata = amici.ExpData(rdata, 0.1, 0.0)
    edata.fixedParametersPreequilibration = [3, 0]
    edata.fixedParameters = [10, 2]
    edata.fixedParametersPresimulation = [10, 2]
    edata.fixedParametersPreequilibration = [3, 0]
    return edata


def get_results(model, edata):
    solver = model.getSolver()
    model.setTimepoints(np.linspace(0, 60, 61))
    model.setReinitializeFixedParameterInitialStates(True)
    return amici.runAmiciSimulation(model, solver, edata)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestAmiciPSBModel())
    unittest.main()