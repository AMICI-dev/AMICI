#!/usr/bin/env python3

"""Run simulations with Matlab-AMICI pre-generated models and verify using
saved expectations."""

import os
from pathlib import Path

import amici
import h5py
import numpy as np
import pytest
from amici.gradient_check import check_derivatives, check_results

cpp_test_dir = Path(__file__).parents[2] / 'tests' / 'cpp'
options_file = str(cpp_test_dir / 'testOptions.h5')
expected_results_file = str(cpp_test_dir / 'expectedResults.h5')
expected_results = h5py.File(expected_results_file, 'r')

model_cases = [(sub_test, case)
               for sub_test in expected_results.keys()
               for case in list(expected_results[sub_test].keys())]


def assert_fun(x):
    assert x


@pytest.mark.skipif(os.environ.get('AMICI_SKIP_CMAKE_TESTS', '') == 'TRUE',
                    reason='skipping cmake based test')
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

    model_swig_folder = str(Path(__file__).parents[2] / 'build' / 'tests'
                            / 'cpp' / f'external_{model_name}-prefix' / 'src'
                            / f'external_{model_name}-build' / 'swig')

    test_model_module = amici.import_model_module(
        module_name=model_name, module_path=model_swig_folder)
    model = test_model_module.getModel()
    solver = model.getSolver()
    amici.readModelDataFromHDF5(
        options_file, model.get(),
        f'/{sub_test}/{case}/options'
    )
    amici.readSolverSettingsFromHDF5(
        options_file, solver.get(),
        f'/{sub_test}/{case}/options'
    )

    edata = None
    if 'data' in expected_results[sub_test][case].keys():
        edata = amici.readSimulationExpData(
            str(expected_results_file),
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
        check_derivatives(model, solver, edata, assert_fun,
                          **check_derivative_opts)

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
        **verify_simulation_opts
    )

    if model_name == 'model_steadystate' and \
            case == 'sensiforwarderrorint':
        edata = amici.amici.ExpData(model.get())

    # Test runAmiciSimulations: ensure running twice
    # with same ExpData yields same results
    if edata and model_name != 'model_neuron_o2' and not (
        model_name == 'model_robertson' and
        case == 'sensiforwardSPBCG'
    ):
        if isinstance(edata, amici.amici.ExpData):
            edatas = [edata, edata]
        else:
            edatas = [edata.get(), edata.get()]

        rdatas = amici.runAmiciSimulations(
            model, solver, edatas, num_threads=2,
            failfast=False
        )
        verify_simulation_results(
            rdatas[0], expected_results[sub_test][case]['results'],
            **verify_simulation_opts
        )
        verify_simulation_results(
            rdatas[1], expected_results[sub_test][case]['results'],
            **verify_simulation_opts
        )

    # test residuals mode
    if solver.getSensitivityMethod() == amici.SensitivityMethod.adjoint:
        with pytest.raises(RuntimeError):
            solver.setReturnDataReportingMode(amici.RDataReporting.residuals)
    else:
        solver.setReturnDataReportingMode(amici.RDataReporting.residuals)
        rdata = amici.runAmiciSimulation(model, solver, edata)
        verify_simulation_results(
            rdata, expected_results[sub_test][case]['results'],
            fields=['t', 'res', 'sres', 'y', 'sy', 'sigmay', 'ssigmay'],
            **verify_simulation_opts
        )
        with pytest.raises(RuntimeError):
            solver.setSensitivityMethod(amici.SensitivityMethod.adjoint)

    chi2_ref = rdata.chi2

    # test likelihood mode
    solver.setReturnDataReportingMode(amici.RDataReporting.likelihood)
    rdata = amici.runAmiciSimulation(model, solver, edata)
    verify_simulation_results(
        rdata, expected_results[sub_test][case]['results'],
        fields=['t', 'llh', 'sllh', 's2llh', 'FIM'], **verify_simulation_opts
    )

    # test sigma residuals

    if model_name == 'model_jakstat_adjoint' and \
            solver.getSensitivityMethod() != amici.SensitivityMethod.adjoint:
        model.setAddSigmaResiduals(True)
        solver.setReturnDataReportingMode(amici.RDataReporting.full)
        rdata = amici.runAmiciSimulation(model, solver, edata)
        # check whether activation changes chi2
        assert chi2_ref != rdata.chi2

        if edata and solver.getSensitivityMethod() and \
                solver.getSensitivityOrder() and len(model.getParameterList()):
            check_derivatives(model, solver, edata, assert_fun,
                              **check_derivative_opts)

        chi2_ref = rdata.chi2
        res_ref = rdata.res

        model.setMinimumSigmaResiduals(100)
        rdata = amici.runAmiciSimulation(model, solver, edata)
        # check whether changing the minimum changes res but not chi2
        assert np.isclose(chi2_ref, rdata.chi2)
        assert not np.allclose(res_ref, rdata.res)

        model.setMinimumSigmaResiduals(-10)
        rdata = amici.runAmiciSimulation(model, solver, edata)
        # check whether having a bad minimum results in nan chi2
        assert np.isnan(rdata.chi2)

    with pytest.raises(RuntimeError):
        model.getParameterByName('thisParameterDoesNotExist')


def verify_simulation_results(rdata, expected_results, fields=None,
                              atol=1e-8, rtol=1e-4):
    """
    compares all fields of the simulation results in rdata against the
    expectedResults using the provided tolerances

    :param rdata: simulation results as returned by amici.runAmiciSimulation
    :param expected_results: stored test results
    :param fields: subsetting of expected results to check
    :param atol: absolute tolerance
    :param rtol: relative tolerance
    """

    subfields = []
    if fields is None:
        attrs = expected_results.attrs.keys()
        fields = expected_results.keys()
        if 'diagnosis' in expected_results.keys():
            subfields = expected_results['diagnosis'].keys()

    else:
        attrs = [field for field in fields
                 if field in expected_results.attrs.keys()]
        if 'diagnosis' in expected_results.keys():
            subfields = [field for field in fields
                         if field in expected_results['diagnosis'].keys()]
        fields = [field for field in fields
                  if field in expected_results.keys()]

    if expected_results.attrs['status'][0] != 0:
        assert rdata['status'] == expected_results.attrs['status'][0]
        return

    for field in expected_results.keys():
        if field == 'diagnosis':
            for subfield in ['J', 'xdot']:
                if subfield not in subfields:
                    assert rdata[subfield] is None, field
                    continue
                check_results(rdata, subfield,
                              expected_results[field][subfield][()],
                              assert_fun, 1e-8, 1e-8)
        else:
            if field not in fields:
                assert rdata[field] is None, field
                continue
            if field == 's2llh':
                check_results(rdata, field, expected_results[field][()],
                              assert_fun, 1e-4, 1e-3)
            else:
                check_results(rdata, field, expected_results[field][()],
                              assert_fun, atol, rtol)

    for attr in attrs:
        check_results(rdata, attr, expected_results.attrs[attr], assert_fun,
                      atol, rtol)
