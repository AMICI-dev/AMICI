"""Tests related to amici.sbml_import"""

import os
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
                                 observables=None)


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


def test_presimulation(sbml_example_presimulation_module):
    """Test 'presimulation' test model"""

    model = sbml_example_presimulation_module.getModel()
    solver = model.getSolver()
    solver.setNewtonMaxSteps(0)
    model.setTimepoints(np.linspace(0, 60, 61))
    model.setSteadyStateSensitivityMode(
        amici.SteadyStateSensitivityMode_simulationFSA
    )
    solver.setSensitivityOrder(amici.SensitivityOrder_first)
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
    solver.setSensitivityOrder(amici.SensitivityOrder_first)
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

    # assign different noise models

    obs_keys = list(observables.keys())

    # exponentiate observable formulas
    obs1 = observables[obs_keys[1]]
    obs3 = observables[obs_keys[3]]
    obs1['formula'] = '10^(' + obs1['formula'] + ')'
    obs3['formula'] = 'exp(' + obs3['formula'] + ')'

    # customize noise distributions
    noise_distributions = {
        obs_keys[0]: 'normal',
        obs_keys[1]: 'log-normal',
        obs_keys[2]: 'laplace',
        obs_keys[3]: 'log10-laplace',
    }

    module_name = 'test_likelihoods'
    outdir = 'test_likelihoods'
    sbml_importer.sbml2amici(
        model_name=module_name,
        output_dir=outdir,
        observables=observables,
        constant_parameters=['k0'],
        sigmas={'observable_x1withsigma':
                    'observable_x1withsigma_sigma'},
        noise_distributions=noise_distributions
    )

    yield amici.import_model_module(module_name=module_name,
                                    module_path=outdir)

    shutil.rmtree(outdir, ignore_errors=True)


def test_likelihoods(model_test_likelihoods):
    """Test the custom noise distributions used to define cost functions."""

    model = model_test_likelihoods.getModel()
    model.setTimepoints(np.linspace(0, 60, 60))
    solver = model.getSolver()
    solver.setSensitivityOrder(amici.SensitivityOrder_first)

    # run model once to create an edata
    rdata = amici.runAmiciSimulation(model, solver)
    edata = [amici.ExpData(rdata, 1, 0)]

    # just make all observables positive since some are logarithmic
    for ed in edata:
        y = ed.getObservedData()
        y = tuple([max(val, 1e-4) for val in y])
        ed.setObservedData(y)

    # and now run for real and also compute likelihood values
    rdata = amici.runAmiciSimulations(model, solver, edata)[0]

    # output for easy debugging
    for key in ['llh', 'sllh']:
        print(key, rdata[key])

    # it would be good to compute the expected llh+sllh by hand,
    # here, we only check if they make overall sense
    assert np.isfinite(rdata['llh'])
    assert np.all(np.isfinite(rdata['sllh']))
    assert np.any(rdata['sllh'])


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
