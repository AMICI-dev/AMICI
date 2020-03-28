#!/usr/bin/env python3

"""Miscellaneous AMICI Python interface tests"""

import amici
import sys
import copy
import os
import unittest
import importlib


class TestAmiciHDF5(unittest.TestCase):
    """TestCase class various AMICI Python interface functions"""

    def setUp(self):
        self.default_path = copy.copy(sys.path)
        self.resetdir = os.getcwd()

        sbml_file = os.path.join(os.path.dirname(__file__), '..', 'python',
                                 'examples', 'example_presimulation',
                                 'model_presimulation.xml')

        sbml_importer = amici.SbmlImporter(sbml_file)

        constant_parameters = ['DRUG_0', 'KIN_0']

        outdir = 'test_model_presimulation'
        sbml_importer.sbml2amici('test_model_presimulation',
                                 outdir,
                                 verbose=False,
                                 constant_parameters=constant_parameters)
        sys.path.insert(0, outdir)
        model_module = importlib.import_module('test_model_presimulation')
        self.model = model_module.getModel()

    def tearDown(self):
        os.chdir(self.resetdir)
        sys.path = self.default_path

    def test_solver_hdf5_roundtrip(self):
        solver = self.model.getSolver()

        # change to non-default values
        for attr in dir(solver):
            if not attr.startswith('set'):
                continue

            val = getattr(solver, attr.replace('set', 'get'))()

            if isinstance(val, bool):
                cval = not val
            elif attr == 'setStabilityLimitFlag':
                cval = 0
            else:
                cval = val + 1

            getattr(solver, attr)(
                cval
            )

        hdf5file = 'solverSettings.hdf5'

        amici.writeSolverSettingsToHDF5(solver.get(), hdf5file,
                                        'ssettings')

        new_solver = self.model.getSolver()

        # check that we changed everything
        for attr in dir(solver):
            if not attr.startswith('set'):
                continue
            self.assertNotEqual(
                getattr(solver, attr.replace('set', 'get'))(),
                getattr(new_solver, attr.replace('set', 'get'))()
            )

        amici.readSolverSettingsFromHDF5(hdf5file, new_solver.get(),
                                         'ssettings')

        # check that reading in settings worked
        for attr in dir(solver):
            if not attr.startswith('set'):
                continue

            self.assertAlmostEqual(
                getattr(solver, attr.replace('set', 'get'))(),
                getattr(new_solver, attr.replace('set', 'get'))()
            )

        os.remove(hdf5file)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestAmiciHDF5())
    unittest.main()
