#!/usr/bin/env python3
"""
Run SBML Test Suite and verify simulation results
[https://github.com/sbmlteam/sbml-test-suite/releases]

Usage:
    testSBMLSuite.py SELECTION
        SELECTION can be e.g.: `1`, `1,3`, or `-3,4,6-7` to select specific
        test cases or 1-1780 to run all.
"""

import os
import sys
import importlib
import pytest
import copy

import amici
import numpy as np
import sympy as sp
import pandas as pd

import libsbml


# directory with sbml semantic test cases
TEST_PATH = os.path.join(os.path.dirname(__file__), 'sbml-test-suite', 'cases',
                         'semantic')


@pytest.fixture(scope="session")
def result_path():
    # ensure directory for test results is empty
    upload_result_path = os.path.join(os.path.dirname(__file__),
                                      'amici-semantic-results')
    return upload_result_path


@pytest.fixture(scope="function", autouse=True)
def sbml_test_dir():
    # setup
    old_cwd = os.getcwd()
    old_path = copy.copy(sys.path)

    yield

    # teardown
    os.chdir(old_cwd)
    sys.path = old_path


def test_sbml_testsuite_case(test_number, result_path):

    test_id = format_test_id(test_number)

    try:
        current_test_path = os.path.join(TEST_PATH, test_id)

        # parse expected results
        results_file = os.path.join(current_test_path,
                                    test_id + '-results.csv')
        results = np.genfromtxt(results_file, delimiter=',')

        # setup model
        model, solver, wrapper = compile_model(current_test_path, test_id)
        settings = read_settings_file(current_test_path, test_id)

        atol, rtol = apply_settings(settings, solver, model)

        # simulate model
        rdata = amici.runAmiciSimulation(model, solver)

        # verify
        simulated_x = verify_results(settings, rdata, results, wrapper,
                                     model, atol, rtol)

        print(f'TestCase {test_id} passed.')

        # record results
        write_result_file(simulated_x, model, test_id, result_path)

    except amici.sbml_import.SBMLException as err:
        #assert False, f"Skipped {test_id}: {err}"
        print(f'TestCase {test_id} was skipped: {err}')


def verify_results(settings, rdata, results, wrapper,
                   model, atol, rtol):
    """Verify test results"""
    amount_species, variables_species = get_amount_and_variables(settings)

    # verify states
    simulated_x = rdata['x']
    sim_shared_cols = [variables_species.index(x_id) for x_id in wrapper.species_index.keys() if x_id in variables_species]
    expected_x = results[1:, [1+c for c in sim_shared_cols]]
    concentrations_to_amounts(amount_species, wrapper, model, simulated_x)
    #np.set_printoptions(linewidth=np.inf)
    assert np.isclose(simulated_x, expected_x, atol, rtol).all()

    # TODO: verify compartment volumes and parameters

    return simulated_x


def concentrations_to_amounts(amount_species, wrapper, model, simulated_x):
    """Convert AMICI simulated concentrations to amounts"""
    for species in amount_species:
        # Skip "species" that are actually compartments
        if not species == '' and species not in [
                sp.sstr(c) for c in wrapper.compartment_symbols]:
            symvolume = wrapper.species_compartment[
                wrapper.species_index[species]
            ]

            # Volumes are already reported for compartments with rate rules (not assignment rules)
            if symvolume in wrapper.compartment_rules and wrapper.compartment_rules[symvolume]['type_code'] == libsbml.SBML_RATE_RULE:
                simulated_x[:, wrapper.species_index[species]] = \
                    simulated_x[:, wrapper.species_index[species]] * \
                    simulated_x[:, wrapper.species_index[sp.sstr(symvolume)]]
                continue

            volume = symvolume.subs({
                comp: vol
                for comp, vol in zip(
                    wrapper.compartment_symbols,
                    wrapper.compartment_volume
                )
            })
            volume = volume.subs({
                sp.Symbol(name, real=True): value
                for name, value in zip(
                    model.getParameterIds(),
                    model.getParameters()
                )
            })

            # required for 525-527, 530 as k is renamed to amici_k
            volume = volume.subs({
                sp.Symbol(name, real=True): value
                for name, value in zip(
                    model.getParameterNames(),
                    model.getParameters()
                )
            })

            simulated_x[:, wrapper.species_index[species]] = \
                simulated_x[:, wrapper.species_index[species]] * volume


def write_result_file(simulated_x: np.array,
                      model: amici.Model,
                      test_id: str, result_path: str):
    """
    Create test result file for upload to
    http://sbml.org/Facilities/Database/Submission/Create

    Requires csv file with test ID in name and content of [time, Species, ...]
    """
    # TODO: only states are reported here, not compartments or parameters

    filename = os.path.join(result_path, f'{test_id}.csv')

    df = pd.DataFrame(simulated_x)
    df.columns = model.getStateIds()
    df.insert(0, 'time', model.getTimepoints())
    df.to_csv(filename, index=False)


def get_amount_and_variables(settings):
    """Read amount and species from settings file"""

    # species for which results are expected as amounts
    amount_species = settings['amount'] \
        .replace(' ', '') \
        .replace('\n', '') \
        .split(',')

    # IDs of all variables for which results are expected/provided
    variables_species = settings['variables'] \
        .replace(' ', '') \
        .replace('\n', '') \
        .split(',')

    return amount_species, variables_species


def apply_settings(settings, solver, model):
    """Apply model and solver settings as specified in the test case"""

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
    """Import the given test model to AMICI"""
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


def find_model_file(current_test_path: str, test_id: str):
    """Find model file for the given test (guess filename extension)"""

    sbml_file = os.path.join(current_test_path, test_id + '-sbml-l3v2.xml')

    # fallback l3v1
    if not os.path.isfile(sbml_file):
        sbml_file = os.path.join(current_test_path, test_id + '-sbml-l3v1.xml')

    # fallback l2v5
    if not os.path.isfile(sbml_file):
        sbml_file = os.path.join(current_test_path, test_id + '-sbml-l2v5.xml')

    return sbml_file


def read_settings_file(current_test_path: str, test_id: str):
    """Read settings for the given test"""
    settings_file = os.path.join(current_test_path, test_id + '-settings.txt')
    settings = {}
    with open(settings_file) as f:
        for line in f:
            if not line == '\n':
                (key, val) = line.split(':')
                settings[key] = val
    return settings


def format_test_id(test_id) -> str:
    """Format numeric to 0-padded string"""
    test_str = str(test_id)
    test_str = '0'*(5-len(test_str)) + test_str
    return test_str

