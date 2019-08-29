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

EXPECTED_TO_WORK = {
    *list(range(1, 26)), 27, *list(range(34, 38)), *list(range(42, 51)),
    *list(range(54, 66)), *list(range(75, 78)), *list(range(95, 104)),
    *list(range(107, 122)), *list(range(125, 129)), *list(range(132, 136)),
    *list(range(145, 151)), *list(range(186, 196)), 199,
    *list(range(202, 269)), *list(range(270, 276)), 277,
    *list(range(280, 286)), *list(range(462, 469)),
    *list(range(471, 477)), *list(range(480, 516)), *list(range(522, 528)),
    *list(range(529, 531)), *list(range(577, 607)), *list(range(616, 619)),
    *list(range(631, 634)), *list(range(676, 679)), 688, 697, 698, 706,
    *list(range(781, 789)),
    *list(range(794, 800)), *list(range(801, 827)), *list(range(830, 839)),
    *list(range(851, 870)), 877, 878, 882, *list(range(894, 898)),
    *list(range(969, 972)),
    *list(range(975, 978)), *list(range(1001, 1008)), *list(range(1009, 1014)),
    *list(range(1018, 1027)), *list(range(1030, 1039)),
    *list(range(1055, 1071)), *list(range(1077, 1083)),
    *list(range(1096, 1103)), 1104, 1107, *list(range(1109, 1112)), 1225, 1232,
    1245, 1246, 1307, 1341, 1342, *list(range(1420, 1435)), 1436, 1438, 1440,
    1442, 1465, 1552, 1554, 1555, 1557, 1574, *list(range(1631, 1641)),
    *list(range(1645, 1651)), 1722, *list(range(1724, 1727)),
    *list(range(1730, 1736)), *list(range(1739, 1742)), 1744, 1746, 1748, 1750,
    1752, 1753, 1760, 1762, 1764, 1767, 1773, 1774
}
ALL_TESTS = set(range(1, 1781))


class TestAmiciSBMLTestSuite(unittest.TestCase):
    SBML_TEST_IDS = ALL_TESTS
    SBML_TESTS_PASSED = set()

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
            TestAmiciSBMLTestSuite.SBML_TESTS_PASSED |= {int(test_id)}


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


def print_stats():
    """Print test statistics"""
    percent_passed = (len(TestAmiciSBMLTestSuite.SBML_TESTS_PASSED)
                     / len(TestAmiciSBMLTestSuite.SBML_TEST_IDS))
    known_to_fail = TestAmiciSBMLTestSuite.SBML_TEST_IDS - EXPECTED_TO_WORK
    failed = (TestAmiciSBMLTestSuite.SBML_TEST_IDS
              - TestAmiciSBMLTestSuite.SBML_TESTS_PASSED)
    unexpectedly_failed = failed - known_to_fail
    unexpectedly_succeeded = (TestAmiciSBMLTestSuite.SBML_TESTS_PASSED
                              - EXPECTED_TO_WORK)

    print()
    print('Tests executed:', len(TestAmiciSBMLTestSuite.SBML_TEST_IDS))
    print('Tests passed:', len(TestAmiciSBMLTestSuite.SBML_TESTS_PASSED),
          f'({percent_passed}% of tests run)')

    if unexpectedly_failed:
        print('ERROR:', len(unexpectedly_failed),
              'tests failed which were expected to work:', unexpectedly_failed)

    if unexpectedly_succeeded:
        print('Congratulations:', len(unexpectedly_succeeded),
              'tests passed which did not pass before:',
              unexpectedly_succeeded)
        print(f'Please update {__file__}.EXPECTED_TO_WORK accordingly.')


def get_return_code() -> int:
    """
    We return success if we pass at least all tests which were known to work
    before.

    Returns:
        exit code
    """
    known_to_fail = TestAmiciSBMLTestSuite.SBML_TEST_IDS - EXPECTED_TO_WORK
    failed = (TestAmiciSBMLTestSuite.SBML_TEST_IDS
              - TestAmiciSBMLTestSuite.SBML_TESTS_PASSED)
    unexpectedly_failed = failed - known_to_fail

    return len(unexpectedly_failed)


if __name__ == '__main__':
    TestAmiciSBMLTestSuite.SBML_TEST_IDS = set(parse_selection(sys.argv.pop()))

    suite = unittest.TestSuite()
    suite.addTest(TestAmiciSBMLTestSuite())
    unittest.main(buffer=False, exit=False)

    print_stats()

    sys.exit(get_return_code())
