import importlib
import amici
import os
import sys
import pytest
import numpy as np
import warnings


@pytest.fixture
def edata_fixture():
    """edata is generated to test preequilibration"""
    edata = amici.ExpData(2, 0, 0,
                          np.array([0., 0.1, 0.2, 0.5, 1., 2., 5., 10.]))
    edata.setObservedData([1.5] * 16)
    edata.fixedParameters = np.array([3., 5.])
    edata.fixedParametersPreequilibration = np.array([0., 5.])
    edata.reinitializeFixedParameterInitialStates = True

    return edata


def generate_models():
    # SBML model we want to import
    sbml_file = os.path.join(os.path.dirname(__file__), '..',
                             'examples', 'example_constant_species',
                             'model_constant_species.xml')
    sbml_importer = amici.SbmlImporter(sbml_file)

    # Name of the model that will also be the name of the python module
    model_name = model_output_dir = 'model_constant_species'
    model_name_cl = model_output_dir_cl = 'model_constant_species_cl'

    # Define constants, observables, sigmas
    constant_parameters = ['init_substrate', 'init_enzyme']
    observables = {
        'observable_product': {'name': '', 'formula': 'product'},
        'observable_substrate': {'name': '', 'formula': 'substrate'},
    }
    sigmas = {'observable_product': 1.0, 'observable_substrate': 1.0}

    # wrap models with and without conservations laws
    sbml_importer.sbml2amici(model_name_cl,
                             model_output_dir_cl,
                             observables=observables,
                             constant_parameters=constant_parameters,
                             sigmas=sigmas)
    sbml_importer.sbml2amici(model_name,
                             model_output_dir,
                             observables=observables,
                             constant_parameters=constant_parameters,
                             sigmas=sigmas,
                             compute_conservation_laws=False)

    # load both models
    sys.path.insert(0, os.path.abspath(model_name))
    sys.path.insert(0, os.path.abspath(model_name_cl))
    model_without_cl_module = importlib.import_module(model_name)
    model_with_cl_module = importlib.import_module(model_name_cl)

    # get the models and return
    model_without_cl = model_without_cl_module.getModel()
    model_with_cl = model_with_cl_module.getModel()
    return model_with_cl, model_without_cl


def get_results(model, edata=None, sensi_order=0):
    # get the solver, set the properties
    solver = model.getSolver()
    solver.setNewtonMaxSteps(20)
    solver.setSensitivityOrder(sensi_order)
    if edata is None:
        model.setTimepoints(np.linspace(0, 5, 101))

    # return simulation results
    return amici.runAmiciSimulation(model, solver, edata)


def test_compare_conservation_laws_sbml(edata_fixture):
    # first, create the modek
    model_with_cl, model_without_cl = generate_models()

    # ----- compare simulations wo edata, sensi = 0, states ------------------
    # run simulations
    rdata_cl = get_results(model_with_cl)
    assert rdata_cl['status'] == amici.AMICI_SUCCESS
    rdata = get_results(model_without_cl)
    assert rdata['status'] == amici.AMICI_SUCCESS

    # compare state trajectories
    assert np.isclose(rdata['x'], rdata_cl['x']).all()

    # ----- compare simulations wo edata, sensi = 1, states and sensis -------
    # run simulations
    rdata_cl = get_results(model_with_cl, sensi_order=1)
    assert rdata_cl['status'] == amici.AMICI_SUCCESS
    rdata = get_results(model_without_cl, sensi_order=1)
    assert rdata['status'] == amici.AMICI_SUCCESS

    # compare state trajectories
    for field in ['x', 'sx']:
        assert np.isclose(rdata[field], rdata_cl[field]).all(), field

    # ----- compare simulations wo edata, sensi = 0, states and sensis -------
    model_without_cl.setSteadyStateSensitivityMode(
        amici.SteadyStateSensitivityMode.simulationFSA
    )

    # run simulations
    rdata_cl = get_results(model_with_cl, edata=edata_fixture)
    assert rdata_cl['status'] == amici.AMICI_SUCCESS
    rdata = get_results(model_without_cl, edata=edata_fixture)
    assert rdata['status'] == amici.AMICI_SUCCESS

    # compare preequilibrated states
    for field in ['x', 'x_ss', 'llh']:
        assert np.isclose(rdata[field], rdata_cl[field]).all(), field

    # ----- compare simulations wo edata, sensi = 1, states and sensis -------

    # run simulations
    rdata_cl = get_results(model_with_cl, edata=edata_fixture, sensi_order=1)
    assert rdata_cl['status'] == amici.AMICI_SUCCESS
    rdata = get_results(model_without_cl, edata=edata_fixture, sensi_order=1)
    assert rdata['status'] == amici.AMICI_SUCCESS

    # compare state sensitivities with edata and preequilibration
    for field in ['x', 'x_ss', 'sx', 'llh', 'sllh']:
        assert np.isclose(rdata[field], rdata_cl[field]).all(), field

    # ----- check failure st.st. sensi computation if run wo CLs -------------
    # check failure of steady state senistivity computation if run wo CLs
    model_without_cl.setSteadyStateSensitivityMode(
        amici.SteadyStateSensitivityMode.newtonOnly
    )
    with pytest.raises(RuntimeError):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            get_results(model_without_cl, edata=edata_fixture, sensi_order=1)
