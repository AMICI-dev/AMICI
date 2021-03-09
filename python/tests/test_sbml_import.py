"""Tests related to amici.sbml_import"""

import os
import re
from urllib.request import urlopen
import shutil
from tempfile import TemporaryDirectory

import amici
import libsbml
import numpy as np
import pytest
from amici.gradient_check import check_derivatives
from amici.sbml_import import SbmlImporter


@pytest.fixture
def simple_sbml_model():
    """Some testmodel"""
    document = libsbml.SBMLDocument(3, 1)
    model = document.createModel()
    model.setTimeUnits("second")
    model.setExtentUnits("mole")
    model.setSubstanceUnits('mole')
    c1 = model.createCompartment()
    c1.setId('C1')
    model.addCompartment(c1)
    s1 = model.createSpecies()
    s1.setId('S1')
    s1.setCompartment('C1')
    model.addSpecies(s1)

    return document, model


def test_sbml2amici_no_observables(simple_sbml_model):
    """Test model generation works for model without observables"""
    sbml_doc, sbml_model = simple_sbml_model
    sbml_importer = SbmlImporter(sbml_source=sbml_model,
                                 from_file=False)

    with TemporaryDirectory() as tmpdir:
        sbml_importer.sbml2amici(model_name="test",
                                 output_dir=tmpdir,
                                 observables=None,
                                 compute_conservation_laws=False)


def test_nosensi(simple_sbml_model):
    sbml_doc, sbml_model = simple_sbml_model
    sbml_importer = SbmlImporter(sbml_source=sbml_model,
                                 from_file=False)

    with TemporaryDirectory() as tmpdir:
        sbml_importer.sbml2amici(model_name="test",
                                 output_dir=tmpdir,
                                 observables=None,
                                 compute_conservation_laws=False,
                                 generate_sensitivity_code=False)

        model_module = amici.import_model_module(module_name='test',
                                                 module_path=tmpdir)

        model = model_module.getModel()
        solver = model.getSolver()
        solver.setSensitivityOrder(amici.SensitivityOrder.first)
        with pytest.raises(RuntimeError):
            amici.runAmiciSimulation(model, solver)


def assert_fun(x):
    assert x


@pytest.fixture
def model_steadystate_module():
    sbml_file = os.path.join(os.path.dirname(__file__), '..',
                             'examples', 'example_steadystate',
                             'model_steadystate_scaled.xml')
    sbml_importer = amici.SbmlImporter(sbml_file)

    observables = amici.assignmentRules2observables(
        sbml_importer.sbml,
        filter_function=lambda variable:
        variable.getId().startswith('observable_') and
        not variable.getId().endswith('_sigma')
    )

    outdir = 'test_model_steadystate_scaled'
    module_name = 'test_model_steadystate_scaled'
    sbml_importer.sbml2amici(
        model_name=module_name,
        output_dir=outdir,
        observables=observables,
        constant_parameters=['k0'],
        sigmas={'observable_x1withsigma': 'observable_x1withsigma_sigma'})

    model_module = amici.import_model_module(module_name=module_name,
                                             module_path=outdir)
    yield model_module

    shutil.rmtree(outdir, ignore_errors=True)


@pytest.fixture
def model_units_module():
    sbml_file = os.path.join(os.path.dirname(__file__), '..',
                             'examples', 'example_units',
                             'model_units.xml')
    sbml_importer = amici.SbmlImporter(sbml_file)

    outdir = 'test_model_units'
    module_name = 'test_model_units'
    sbml_importer.sbml2amici(model_name=module_name,
                             output_dir=outdir)

    model_module = amici.import_model_module(module_name=module_name,
                                             module_path=outdir)
    yield model_module

    shutil.rmtree(outdir, ignore_errors=True)


def test_presimulation(sbml_example_presimulation_module):
    """Test 'presimulation' test model"""
    model = sbml_example_presimulation_module.getModel()
    solver = model.getSolver()
    solver.setNewtonMaxSteps(0)
    model.setTimepoints(np.linspace(0, 60, 61))
    model.setSteadyStateSensitivityMode(
        amici.SteadyStateSensitivityMode.simulationFSA
    )
    solver.setSensitivityOrder(amici.SensitivityOrder.first)
    model.setReinitializeFixedParameterInitialStates(True)

    rdata = amici.runAmiciSimulation(model, solver)
    edata = amici.ExpData(rdata, 0.1, 0.0)
    edata.fixedParameters = [10, 2]
    edata.fixedParametersPresimulation = [10, 2]
    edata.fixedParametersPreequilibration = [3, 0]
    assert isinstance(
        amici.runAmiciSimulation(model, solver, edata),
        amici.ReturnDataView)

    solver.setRelativeTolerance(1e-12)
    solver.setAbsoluteTolerance(1e-12)
    check_derivatives(model, solver, edata, assert_fun, epsilon=1e-4)


def test_steadystate_simulation(model_steadystate_module):
    model = model_steadystate_module.getModel()
    model.setTimepoints(np.linspace(0, 60, 60))
    solver = model.getSolver()
    solver.setSensitivityOrder(amici.SensitivityOrder.first)
    rdata = amici.runAmiciSimulation(model, solver)
    edata = [amici.ExpData(rdata, 1, 0)]
    rdata = amici.runAmiciSimulations(model, solver, edata)

    # check roundtripping of DataFrame conversion
    df_edata = amici.getDataObservablesAsDataFrame(model, edata)
    edata_reconstructed = amici.getEdataFromDataFrame(model, df_edata)

    assert np.isclose(
        amici.ExpDataView(edata[0])['observedData'],
        amici.ExpDataView(edata_reconstructed[0])['observedData']).all()

    assert np.isclose(
        amici.ExpDataView(edata[0])['observedDataStdDev'],
        amici.ExpDataView(edata_reconstructed[0])['observedDataStdDev']).all()

    if len(edata[0].fixedParameters):
        assert list(edata[0].fixedParameters) \
               == list(edata_reconstructed[0].fixedParameters)

    else:
        assert list(model.getFixedParameters()) \
               == list(edata_reconstructed[0].fixedParameters)

    assert list(edata[0].fixedParametersPreequilibration) == \
        list(edata_reconstructed[0].fixedParametersPreequilibration)

    df_state = amici.getSimulationStatesAsDataFrame(model, edata, rdata)
    assert np.isclose(rdata[0]['x'],
                      df_state[list(model.getStateIds())].values).all()

    df_obs = amici.getSimulationObservablesAsDataFrame(model, edata, rdata)
    assert np.isclose(rdata[0]['y'],
                      df_obs[list(model.getObservableIds())].values).all()
    amici.getResidualsAsDataFrame(model, edata, rdata)

    solver.setRelativeTolerance(1e-12)
    solver.setAbsoluteTolerance(1e-12)
    check_derivatives(model, solver, edata[0], assert_fun, atol=1e-3,
                      rtol=1e-3, epsilon=1e-4)

    # Run some additional tests which need a working Model,
    # but don't need precomputed expectations.
    _test_set_parameters_by_dict(model_steadystate_module)


@pytest.fixture
def model_test_likelihoods():
    """Test model for various likelihood functions."""
    # load sbml model
    sbml_file = os.path.join(os.path.dirname(__file__), '..',
                             'examples', 'example_steadystate',
                             'model_steadystate_scaled.xml')
    sbml_importer = amici.SbmlImporter(sbml_file)

    # define observables
    observables = {
        'o1': {'formula': 'x1'},
        'o2': {'formula': '10^x1'},
        'o3': {'formula': '10^x1'},
        'o4': {'formula': 'x1'},
        'o5': {'formula': '10^x1'},
        'o6': {'formula': '10^x1'},
        'o7': {'formula': 'x1'}
    }

    # define different noise models
    noise_distributions = {
        'o1': 'normal', 'o2': 'log-normal', 'o3': 'log10-normal',
        'o4': 'laplace', 'o5': 'log-laplace', 'o6': 'log10-laplace',
        'o7': lambda str_symbol: f'Abs({str_symbol} - m{str_symbol}) '
                                 f'/ sigma{str_symbol}',
    }

    module_name = 'test_likelihoods'
    outdir = 'test_likelihoods'
    sbml_importer.sbml2amici(
        model_name=module_name,
        output_dir=outdir,
        observables=observables,
        constant_parameters=['k0'],
        noise_distributions=noise_distributions,
    )

    yield amici.import_model_module(module_name=module_name,
                                    module_path=outdir)

    shutil.rmtree(outdir, ignore_errors=True)


def test_likelihoods(model_test_likelihoods):
    """Test the custom noise distributions used to define cost functions."""
    model = model_test_likelihoods.getModel()
    model.setTimepoints(np.linspace(0, 60, 60))
    solver = model.getSolver()
    solver.setSensitivityOrder(amici.SensitivityOrder.first)

    # run model once to create an edata

    rdata = amici.runAmiciSimulation(model, solver)
    sigmas = rdata['y'].max(axis=0) * 0.05
    edata = amici.ExpData(rdata, sigmas, [])
    # just make all observables positive since some are logarithmic
    while min(edata.getObservedData()) < 0:
        edata = amici.ExpData(rdata, sigmas, [])

    # and now run for real and also compute likelihood values
    rdata = amici.runAmiciSimulations(model, solver, [edata])[0]

    # check if the values make overall sense
    assert np.isfinite(rdata['llh'])
    assert np.all(np.isfinite(rdata['sllh']))
    assert np.any(rdata['sllh'])

    rdata_df = amici.getSimulationObservablesAsDataFrame(
        model, edata, rdata, by_id=True)
    edata_df = amici.getDataObservablesAsDataFrame(
        model, edata, by_id=True)

    # check correct likelihood value
    llh_exp = - sum([
        normal_nllh(edata_df['o1'], rdata_df['o1'], sigmas[0]),
        log_normal_nllh(edata_df['o2'], rdata_df['o2'], sigmas[1]),
        log10_normal_nllh(edata_df['o3'], rdata_df['o3'], sigmas[2]),
        laplace_nllh(edata_df['o4'], rdata_df['o4'], sigmas[3]),
        log_laplace_nllh(edata_df['o5'], rdata_df['o5'], sigmas[4]),
        log10_laplace_nllh(edata_df['o6'], rdata_df['o6'], sigmas[5]),
        custom_nllh(edata_df['o7'], rdata_df['o7'], sigmas[6]),
    ])
    assert np.isclose(rdata['llh'], llh_exp)

    # check gradient
    for sensi_method in [amici.SensitivityMethod.forward,
                         amici.SensitivityMethod.adjoint]:
        solver = model.getSolver()
        solver.setSensitivityMethod(sensi_method)
        solver.setSensitivityOrder(amici.SensitivityOrder.first)
        solver.setRelativeTolerance(1e-12)
        solver.setAbsoluteTolerance(1e-12)
        check_derivatives(
            model, solver, edata, assert_fun, atol=1e-2, rtol=1e-2,
            epsilon=1e-5, check_least_squares=False
        )


def test_likelihoods_error():
    """Test whether wrong inputs lead to expected errors."""
    sbml_file = os.path.join(os.path.dirname(__file__), '..',
                             'examples', 'example_steadystate',
                             'model_steadystate_scaled.xml')
    sbml_importer = amici.SbmlImporter(sbml_file)

    # define observables
    observables = {'o1': {'formula': 'x1'}}

    # define different noise models
    noise_distributions = {'o1': 'nörmal'}

    module_name = 'test_likelihoods'
    outdir = 'test_likelihoods'
    with pytest.raises(ValueError):
        sbml_importer.sbml2amici(
            model_name=module_name,
            output_dir=outdir,
            observables=observables,
            constant_parameters=['k0'],
            noise_distributions=noise_distributions,
        )


def test_units(model_units_module):
    """
    Test whether SBML import works for models using sbml:units annotations.
    """
    model = model_units_module.getModel()
    model.setTimepoints(np.linspace(0, 1, 101))
    solver = model.getSolver()

    rdata = amici.runAmiciSimulation(model, solver)
    assert rdata['status'] == amici.AMICI_SUCCESS


@pytest.mark.skipif(os.name == 'nt',
                    reason='Avoid `CERTIFICATE_VERIFY_FAILED` error')
def test_sympy_exp_monkeypatch():
    """
    This model contains a removeable discontinuity at t=0 that requires
    monkeypatching sympy.Pow._eval_derivative in order to be able to compute
    non-nan sensitivities
    """
    url = 'https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000529.2?' \
          'filename=BIOMD0000000529_url.xml'
    importer = amici.SbmlImporter(urlopen(url).read().decode('utf-8'),
                                  from_file=False)
    module_name = 'BIOMD0000000529'

    with TemporaryDirectory() as outdir:
        importer.sbml2amici(module_name, outdir)
        model_module = amici.import_model_module(module_name=module_name,
                                                 module_path=outdir)

        model = model_module.getModel()
        model.setTimepoints(np.linspace(0, 8, 250))
        model.requireSensitivitiesForAllParameters()
        model.setAlwaysCheckFinite(True)
        model.setParameterScale(amici.parameterScalingFromIntVector([
            amici.ParameterScaling.none
            if re.match(r'n[0-9]+$', par_id)
            else amici.ParameterScaling.log10
            for par_id in model.getParameterIds()
        ]))

        solver = model.getSolver()
        solver.setSensitivityMethod(amici.SensitivityMethod.forward)
        solver.setSensitivityOrder(amici.SensitivityOrder.first)

        rdata = amici.runAmiciSimulation(model, solver)

        # print sensitivity-related results
        assert rdata['status'] == amici.AMICI_SUCCESS
        check_derivatives(model, solver, None, assert_fun, atol=1e-2, rtol=1e-2,
                          epsilon=1e-3)


def normal_nllh(m, y, sigma):
    return sum(.5*(np.log(2*np.pi*sigma**2) + ((y-m)/sigma)**2))


def log_normal_nllh(m, y, sigma):
    return sum(.5*(np.log(2*np.pi*sigma**2*m**2)
                   + ((np.log(y)-np.log(m))/sigma)**2))


def log10_normal_nllh(m, y, sigma):
    return sum(.5*(np.log(2*np.pi*sigma**2*m**2*np.log(10)**2)
                   + ((np.log10(y) - np.log10(m))/sigma)**2))


def laplace_nllh(m, y, sigma):
    return sum(np.log(2*sigma) + np.abs(y-m)/sigma)


def log_laplace_nllh(m, y, sigma):
    return sum(np.log(2*sigma*m) + np.abs(np.log(y)-np.log(m))/sigma)


def log10_laplace_nllh(m, y, sigma):
    return sum(np.log(2*sigma*m*np.log(10))
               + np.abs(np.log10(y)-np.log10(m))/sigma)


def custom_nllh(m, y, sigma):
    return sum(np.abs(m-y)/sigma)


def _test_set_parameters_by_dict(model_module):
    """Test setting parameter via id/name => value dicts"""
    model = model_module.getModel()
    old_parameter_values = model.getParameters()
    parameter_ids = model.getParameterIds()
    change_par_id = parameter_ids[-1]
    new_par_val = 0.1234
    old_par_val = model.getParameterById(change_par_id)

    assert model.getParameterById(change_par_id) != new_par_val
    model.setParameterById({change_par_id: new_par_val})
    assert model.getParameterById(change_par_id) == new_par_val
    # reset and check we are back to original
    model.setParameterById(change_par_id, old_par_val)
    assert model.getParameters() == old_parameter_values

    # Same for by-name
    parameter_names = model.getParameterNames()
    change_par_name = parameter_names[-1]
    model.setParameterByName({change_par_name: new_par_val})
    assert model.getParameterByName(change_par_name) == new_par_val
    model.setParameterByName(change_par_name, old_par_val)
    assert model.getParameters() == old_parameter_values
