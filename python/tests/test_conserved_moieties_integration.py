import os
import unittest

import libsbml
import numpy as np


# from petab.petab_util import folder_base

class MyTest(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(MyTest, self).__init__(*args, **kwargs)
        self.model_setup = {}
        # SBML model we want to import
        sbml_file = '/home/stephan/Code/pyPESTO/doc/example/tmp/benchmark-models/Benchmark-Models/Blasi_CellSystems2016/model_Blasi_CellSystems2016.xml'

        # Name of the models that will also be the name of the python module
        model_name = 'model_Blasi_CellSystems2016'
        model_reduced_name = model_name
        # Directories to which the generated model code is written
        model_output_dir = model_name
        model_reduced_output_dir = model_reduced_name

        # Read the model and give some output
        sbml_reader = libsbml.SBMLReader()
        sbml_doc = sbml_reader.readSBML(sbml_file)
        sbml_model = sbml_doc.getModel()

        for reaction in sbml_model.getListOfReactions():
            reactants = ' + '.join(['%s %s' % (
            int(r.getStoichiometry()) if r.getStoichiometry() > 1 else '',
            r.getSpecies()) for r in reaction.getListOfReactants()])
            products = ' + '.join(['%s %s' % (
            int(r.getStoichiometry()) if r.getStoichiometry() > 1 else '',
            r.getSpecies()) for r in reaction.getListOfProducts()])
            reversible = '<' if reaction.getReversible() else ''
            print('%3s: %10s %1s->%10s\t\t[%s]' % (reaction.getId(),
                                                   reactants,
                                                   reversible,
                                                   products,
                                                   libsbml.formulaToL3String(
                                                       reaction.getKineticLaw().getMath())))

        # specify observables and constant parameters
        constantParameters = ['d']
        observables = {
            'observable_0ac': {'name': '', 'formula': 'x_0ac'},
            'observable_4ac': {'name': '', 'formula': 'x_4ac'},
            'observable_k12': {'name': '', 'formula': 'x_k12'},
            'observable_k12k16': {'name': '', 'formula': 'x_k12k16'},
            'observable_k16': {'name': '', 'formula': 'x_k16'},
            'observable_k5': {'name': '', 'formula': 'x_k5'},
        }

        sigmas = {'observable_0ac': 0.1,
                  'observable_4ac': 0.1,
                  'observable_k12': 0.1,
                  'observable_k12k16': 0.1,
                  'observable_k16': 0.1,
                  'observable_k5': 0.1}

        self.model_setup = {"model_reduced_name": model_reduced_name,
                            "model_reduced_output_dir": model_reduced_output_dir,
                            "constantParameters": constantParameters,
                            "sigmas": sigmas,
                            "observables": observables}

    def helper(self):
        # import the model
        # Create an SbmlImporter instance for our SBML model
        import amici
        sbml_file = '/home/stephan/Code/pyPESTO/doc/example/tmp/benchmark-models/Benchmark-Models/Blasi_CellSystems2016/model_Blasi_CellSystems2016.xml'
        sbml_importer = amici.SbmlImporter(sbml_file)
        sbml_importer.sbml2amici(self.model_setup["model_reduced_name"],
                                 self.model_setup["model_reduced_output_dir"],
                                 observables=self.model_setup["observables"],
                                 constantParameters=self.model_setup[
                                     "constantParameters"],
                                 sigmas=self.model_setup["sigmas"],
                                 compute_conservation_laws=self.model_setup[
                                     "compute_conservation_laws"])

        sbml_importer.sbml2amici(self.model_setup["model_reduced_name"],
                                 self.model_setup["model_reduced_output_dir"],
                                 observables=self.model_setup["observables"],
                                 constantParameters=self.model_setup[
                                     "constantParameters"],
                                 sigmas=self.model_setup["sigmas"],
                                 compute_conservation_laws=self.model_setup[
                                     "compute_conservation_laws"])

        # import the models and run some test simulations
        model_reduced_module = amici.import_model_module(
            self.model_setup['model_reduced_name'],
            os.path.abspath(self.model_setup['model_reduced_output_dir']))
        model_reduced = model_reduced_module.getModel()

        # simulate model with conservation laws
        model_reduced.setTimepoints(np.linspace(0, 2, 100))
        solver_reduced = model_reduced.getSolver()
        rdata_reduced = amici.runAmiciSimulation(model_reduced, solver_reduced)

        # Call postequilibration by setting an infinity timepoint
        model_reduced.setTimepoints(np.full(1, np.inf))

        # set the solver
        solver = model_reduced.getSolver()
        solver.setNewtonMaxSteps(10)
        solver.setMaxSteps(1000)
        rdata = amici.runAmiciSimulation(model_reduced, solver)

        # np.set_printoptions(threshold=8, edgeitems=2)
        for key, value in rdata.items():
            print('%12s: ' % key, value)

        return rdata

    def test_fail(self):
        """ Without computing conservation laws and removing them the test should always fail """
        self.model_setup['compute_conservation_laws'] = False
        rdata = self.helper()
        status = int(rdata["status"])
        self.assertEqual(status, 0,
                         "Newton solver required to fail for singular Jacobian")

    def test_success(self):
        """ With computing conservation laws and removing them the test should always succeed """
        self.model_setup['compute_conservation_laws'] = True
        rdata = self.helper()
        status = int(rdata["status"])
        self.assertEqual(status, 1,
                         "Newton solver still fails due to singular Jacobian")
