#!/usr/bin/env python3

"""Run simulations with Matlab-AMICI pre-generated models and verify using
saved expectations."""

import sys
import h5py
import amici
import unittest
import importlib
import os
import copy
from amici.gradient_check import check_derivatives, check_results
import pytest


expected_results_file = os.path.join(os.path.dirname(__file__), '..', '..',
                                     'tests', 'cpputest', 'expectedResults.h5')
expected_results = h5py.File(expected_results_file, 'r')

model_cases = [(sub_test, case)
               for sub_test in expected_results.keys()
               for case in list(expected_results[sub_test].keys())]


def assert_fun(x):
    assert x


@pytest.mark.parametrize("sub_test,case", model_cases)
def test_pregenerated_model(sub_test, case):
    """Tests models that were pregenerated using the the matlab code
    generation routines and cmake build routines.

    NOTE: requires having run `make python-tests` in /build/ before to build
    the python modules for the test models.
    """

    if case.startswith('sensi2'):
        model_name = sub_test + '_o2'
    else:
        model_name = sub_test

    model_swig_folder = \
        os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'tests',
                     'cpputest', f'external_{model_name}-prefix',
                     'src', f'external_{model_name}-build', 'swig')

    test_model_module = amici.import_model_module(
        module_name=model_name, module_path=model_swig_folder)
    model = test_model_module.getModel()
    solver = model.getSolver()
    amici.readModelDataFromHDF5(
        expected_results_file, model.get(),
        f'/{sub_test}/{case}/options'
    )
    amici.readSolverSettingsFromHDF5(
        expected_results_file, solver.get(),
        f'/{sub_test}/{case}/options'
    )

    edata = None
    if 'data' in expected_results[sub_test][case].keys():
        edata = amici.readSimulationExpData(
            expected_results_file,
            f'/{sub_test}/{case}/data', model.get()
        )
    rdata = amici.runAmiciSimulation(model, solver,
                                     edata)

    check_derivative_opts = dict()

    if model_name == 'model_nested_events':
        check_derivative_opts['rtol'] = 1e-2
    elif model_name == 'model_events':
        check_derivative_opts['atol'] = 1e-3

    if edata \
            and solver.getSensitivityMethod() \
            and solver.getSensitivityOrder() \
            and len(model.getParameterList()) \
            and not model_name.startswith('model_neuron') \
            and not case.endswith('byhandpreeq'):
        check_derivatives(model, solver, edata,
                          assert_fun, **check_derivative_opts)

    verify_simulation_opts = dict()

    if model_name.startswith('model_neuron'):
        verify_simulation_opts['atol'] = 1e-5
        verify_simulation_opts['rtol'] = 1e-2

    if model_name.startswith('model_robertson') and \
            case == 'sensiforwardSPBCG':
        verify_simulation_opts['atol'] = 1e-3
        verify_simulation_opts['rtol'] = 1e-3

    verify_simulation_results(
        rdata, expected_results[sub_test][case]['results'],
        assert_fun, **verify_simulation_opts
    )

    if model_name == 'model_steadystate' and \
            case == 'sensiforwarderrorint':
        edata = amici.amici.ExpData(model.get())

    if edata and model_name != 'model_neuron_o2' and not (
        model_name == 'model_robertson' and
        case == 'sensiforwardSPBCG'
    ):
        # Test runAmiciSimulations: ensure running twice
        # with same ExpData yields same results
        if isinstance(edata, amici.amici.ExpData):
            edatas = [edata, edata]
        else:
            edatas = [edata.get(), edata.get()]

        rdatas = amici.runAmiciSimulations(
            model, solver, edatas, num_threads=2,
            failfast=False
        )
        verify_simulation_results(
            rdatas[0],
            expected_results[sub_test][case]['results'],
            assert_fun, **verify_simulation_opts
        )
        verify_simulation_results(
            rdatas[1],
            expected_results[sub_test][case]['results'],
            assert_fun, **verify_simulation_opts
        )

    with pytest.raises(RuntimeError):
        model.getParameterByName('thisParameterDoesNotExist')


def verify_simulation_results(rdata, expected_results, assert_fun,
                              atol=1e-8, rtol=1e-4):
    """
    compares all fields of the simulation results in rdata against the
    expectedResults using the provided tolerances

    Arguments:
        rdata: simulation results as returned by amici.runAmiciSimulation
        expected_results: stored test results
        atol: absolute tolerance
        rtol: relative tolerance
    """

    if expected_results.attrs['status'][0] != 0:
        assert rdata['status'] == expected_results.attrs['status'][0]
        return

    for field in expected_results.keys():
        if field == 'diagnosis':
            for subfield in ['J', 'xdot']:
                check_results(rdata, subfield,
                              expected_results[field][subfield][()],
                              assert_fun, 1e-8, 1e-8)
        else:
            if field == 's2llh':
                check_results(rdata, field, expected_results[field][()],
                              assert_fun, 1e-4, 1e-3)
            else:
                check_results(rdata, field, expected_results[field][()],
                              assert_fun, atol, rtol)

    for attr in expected_results.attrs.keys():
        check_results(rdata, attr, expected_results.attrs[attr], assert_fun,
                      atol, rtol)
