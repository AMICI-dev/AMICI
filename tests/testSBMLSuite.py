#!/usr/bin/env python3
"""
Run SBML Test Suite and verify simulation results
[https://github.com/sbmlteam/sbml-test-suite/releases]

Usage:
    testSBMLSuite.py SELECTION
        SELECTION can be e.g.: `1`, `1,3`, or `-3,4,6-7` to select specific
        test cases or 1-1780 to run all.
"""

import re
import os
import sys
import importlib
import numpy as np
import sympy as sp
import amici
import unittest
import copy
from typing import List

# directory with sbml semantic test cases
test_path = os.path.join(os.path.dirname(__file__), 'sbml-test-suite', 'cases',
                         'semantic')


ALL_TESTS = set(range(1, 1781))


class TestAmiciSBMLTestSuite(unittest.TestCase):
    SBML_TEST_IDS = ALL_TESTS

    def setUp(self):
        self.resetdir = os.getcwd()
        self.default_path = copy.copy(sys.path)

    def tearDown(self):
        os.chdir(self.resetdir)
        sys.path = self.default_path

    def runTest(self):
        self.test_sbml_testsuite()

    def test_sbml_testsuite(self):
        for testId in self.SBML_TEST_IDS:
            if testId != 1395:  # we skip this test due to NaNs in the
                # jacobian
                with self.subTest(testId=testId):
                    self.run_sbml_testsuite_case(get_test_str(testId))

    def run_sbml_testsuite_case(self, test_id):
        try:
            current_test_path = os.path.join(test_path, test_id)

            # results
            results_file = os.path.join(current_test_path,
                                        test_id + '-results.csv')
            results = np.genfromtxt(results_file, delimiter=',')

            model, solver, wrapper = compile_model(current_test_path, test_id)
            atol, rtol = apply_settings(current_test_path, test_id, solver,
                                        model)

            rdata = amici.runAmiciSimulation(model, solver)

            amount_species, variables_species = get_amount_and_variables(
                current_test_path, test_id
            )

            simulated_x = rdata['x']
            test_x = results[1:, [
                                     1 + wrapper.speciesIndex[variable]
                                     for variable in variables_species
                                     if variable in wrapper.speciesIndex.keys()
                                 ]]

            for species in amount_species:
                if not species == '':
                    symvolume = wrapper.speciesCompartment[
                        wrapper.speciesIndex[species]
                    ]
                    volume = symvolume.subs({
                        comp: vol
                        for comp, vol in zip(
                            wrapper.compartmentSymbols,
                            wrapper.compartmentVolume
                        )
                    })
                    volume = volume.subs({
                        sp.Symbol(name): value
                        for name, value in zip(
                            model.getParameterIds(),
                            model.getParameters()
                        )
                    })

                    # required for 525-527, 530 as k is renamed to amici_k
                    volume = volume.subs({
                        sp.Symbol(name): value
                        for name, value in zip(
                        model.getParameterNames(),
                        model.getParameters()
                    )
                    })

                    simulated_x[:, wrapper.speciesIndex[species]] = \
                        simulated_x[:, wrapper.speciesIndex[species]] * volume

            self.assertTrue(
                np.isclose(simulated_x, test_x, atol, rtol).all()
            )
            print(f'TestCase {test_id} passed.')

        except amici.sbml_import.SBMLException as err:
            print(f'TestCase {test_id} was skipped: {err}')


def get_amount_and_variables(current_test_path, test_id):
    settings = read_settings_file(current_test_path, test_id)

    amount_species = settings['amount'] \
        .replace(' ', '') \
        .replace('\n', '') \
        .split(',')
    variables_species = settings['variables'] \
        .replace(' ', '') \
        .replace('\n', '') \
        .split(',')

    return amount_species, variables_species


def apply_settings(current_test_path, test_id, solver, model):
    settings = read_settings_file(current_test_path, test_id)

    ts = np.linspace(float(settings['start']),
                     float(settings['start'])
                     + float(settings['duration']),
                     int(settings['steps']) + 1)
    atol = float(settings['absolute'])
    rtol = float(settings['relative'])

    model.setTimepoints(ts)
    solver.setMaxSteps(int(1e6))
    solver.setRelativeTolerance(rtol / 1000.0)
    solver.setAbsoluteTolerance(atol / 1000.0)

    return atol, rtol


def compile_model(path, test_id):
    sbml_file = find_model_file(path, test_id)

    wrapper = amici.SbmlImporter(sbml_file)

    model_dir = os.path.join(os.path.dirname(__file__), 'SBMLTestModels',
                             test_id)
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)

    model_name = 'SBMLTest' + test_id
    wrapper.sbml2amici(model_name, output_dir=model_dir)

    # settings
    sys.path.insert(0, model_dir)
    model_module = importlib.import_module(model_name)

    model = model_module.getModel()
    solver = model.getSolver()

    return model, solver, wrapper


def find_model_file(current_test_path, testId):
    """Find model file for the given test (guess filename extension)"""
    sbmlFile = os.path.join(current_test_path, testId + '-sbml-l3v2.xml')

    # fallback l3v1
    if not os.path.isfile(sbmlFile):
        sbmlFile = os.path.join(current_test_path, testId + '-sbml-l3v1.xml')

    # fallback l2v5
    if not os.path.isfile(sbmlFile):
        sbmlFile = os.path.join(current_test_path, testId + '-sbml-l2v5.xml')

    return sbmlFile


def read_settings_file(current_test_path, test_id):
    """Read settings for the given test"""
    settings_file = os.path.join(current_test_path, test_id + '-settings.txt')
    settings = {}
    with open(settings_file) as f:
        for line in f:
            if not line == '\n':
                (key, val) = line.split(':')
                settings[key] = val
    return settings


def get_test_str(test_id):
    test_str = str(test_id)
    test_str = '0'*(5-len(test_str)) + test_str
    return test_str


def parse_selection(selection_str: str) -> List[int]:
    """
    Parse comma-separated list of integer ranges, return selected indices as
    integer list

    Valid input e.g.: 1 1,3 -3,4,6-7
    """
    indices = []
    for group in selection_str.split(','):
        if not re.match(r'^(?:-?\d+)|(?:\d+(?:-\d+))$', group):
            print("Invalid selection", group)
            sys.exit()
        spl = group.split('-')
        if len(spl) == 1:
            indices.append(int(spl[0]))
        elif len(spl) == 2:
            begin = int(spl[0]) if spl[0] else 0
            end = int(spl[1])
            indices.extend(range(begin, end + 1))
    return indices


if __name__ == '__main__':
    TestAmiciSBMLTestSuite.SBML_TEST_IDS = set(parse_selection(sys.argv.pop()))

    suite = unittest.TestSuite()
    suite.addTest(TestAmiciSBMLTestSuite())
    unittest.main()
