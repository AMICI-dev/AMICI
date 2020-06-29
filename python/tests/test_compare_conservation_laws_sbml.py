import amici
import os
import sys
import pytest
import numpy as np
import warnings


@pytest.fixture
def edata_fixture():
    """edata is generated to test pre- and postequilibration"""
    edata_pre = amici.ExpData(2, 0, 0,
                          np.array([0., 0.1, 0.2, 0.5, 1., 2., 5., 10.]))
    edata_pre.setObservedData([1.5] * 16)
    edata_pre.fixedParameters = np.array([5., 20.])
    edata_pre.fixedParametersPreequilibration = np.array([0., 10.])
    edata_pre.reinitializeFixedParameterInitialStates = True

    # edata for postequilibration
    edata_post = amici.ExpData(2, 0, 0,
                          np.array([float('inf')] * 3))
    edata_post.setObservedData([0.75] * 6)
    edata_post.fixedParameters = np.array([7.5, 30.])

    # edata with both equilibrations
    edata_full = amici.ExpData(2, 0, 0,
        np.array([0., 0., 0., 1., 2., 2., 4., float('inf'), float('inf')]))
    edata_full.setObservedData([3.14] * 18)
    edata_full.fixedParameters = np.array([1., 2.])
    edata_full.fixedParametersPreequilibration = np.array([3., 4.])
    edata_full.reinitializeFixedParameterInitialStates = True

    return edata_pre, edata_post, edata_full


def generate_models():
    # SBML model we want to import
    sbml_file = os.path.join(os.path.dirname(__file__), '..',
                             'examples', 'example_constant_species',
                             'model_constant_species.xml')
    sbml_importer = amici.SbmlImporter(sbml_file)

    # Name of the model that will also be the name of the python module
    model_name =  model_output_dir ='model_constant_species'
    model_name_cl = model_output_dir_cl = 'model_constant_species_cl'

    # Define constants, observables, sigmas
    constant_parameters = ['synthesis_substrate', 'init_enzyme']
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
    model_without_cl_module = amici.import_model_module(model_name)
    model_with_cl_module = amici.import_model_module(model_name_cl)

    # get the models and return
    model_without_cl = model_without_cl_module.getModel()
    model_with_cl = model_with_cl_module.getModel()
    return model_with_cl, model_without_cl


def get_results(model, edata=None, sensi_order=0,
                sensi_meth=amici.SensitivityMethod.forward,
                sensi_meth_preeq=amici.SensitivityMethod.forward,
                reinitialize_states=False):

    # set model and data properties
    model.setReinitializeFixedParameterInitialStates(reinitialize_states)

    # get the solver, set the properties
    solver = model.getSolver()
    solver.setNewtonMaxSteps(20)
    solver.setSensitivityOrder(sensi_order)
    solver.setSensitivityMethodPreequilibration(sensi_meth_preeq)
    solver.setSensitivityMethod(sensi_meth)
    if edata is None:
        model.setTimepoints(np.linspace(0, 5, 101))
    else:
        edata.reinitializeFixedParameterInitialStates = reinitialize_states

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
    edata, _, _ = edata_fixture
    rdata_cl = get_results(model_with_cl, edata=edata)
    assert rdata_cl['status'] == amici.AMICI_SUCCESS
    rdata = get_results(model_without_cl, edata=edata)
    assert rdata['status'] == amici.AMICI_SUCCESS

    # compare preequilibrated states
    for field in ['x', 'x_ss', 'llh']:
        assert np.isclose(rdata[field], rdata_cl[field]).all(), field

    # ----- compare simulations wo edata, sensi = 1, states and sensis -------

    # run simulations
    rdata_cl = get_results(model_with_cl, edata=edata, sensi_order=1)
    assert rdata_cl['status'] == amici.AMICI_SUCCESS
    rdata = get_results(model_without_cl, edata=edata, sensi_order=1)
    assert rdata['status'] == amici.AMICI_SUCCESS
    # check that steady state computation succeeded only by sim in full model
    assert (rdata['preeq_status'] == np.array([-3, 1, 0])).all()
    # check that steady state computation succeeded by Newton in reduced model
    assert (rdata_cl['preeq_status'] == np.array([1, 0, 0])).all()

    # compare state sensitivities with edata and preequilibration
    for field in ['x', 'x_ss', 'sx', 'llh', 'sllh']:
        assert np.isclose(rdata[field], rdata_cl[field]).all(), field

    # ----- check failure st.st. sensi computation if run wo CLs -------------
    # check failure of steady state senistivity computation if run wo CLs
    model_without_cl.setSteadyStateSensitivityMode(
        amici.SteadyStateSensitivityMode.newtonOnly
    )
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        rdata = get_results(model_without_cl, edata=edata, sensi_order=1)
        assert rdata['status'] == amici.AMICI_ERROR

def test_adjoint_pre_and_post_equilibration(edata_fixture):
    # get the model
    model_module = amici.import_model_module('model_constant_species_cl')
    model = model_module.getModel()

    # check gradient with and without state reinitialization
    for edata in edata_fixture:
        for reinit in [False, True]:
            # run simulations
            rff = get_results(model, edata=edata, sensi_order=1,
                              sensi_meth=amici.SensitivityMethod.forward,
                              sensi_meth_preeq=amici.SensitivityMethod.forward,
                              reinitialize_states=reinit)
            rfa = get_results(model, edata=edata, sensi_order=1,
                              sensi_meth=amici.SensitivityMethod.adjoint,
                              sensi_meth_preeq=amici.SensitivityMethod.forward,
                              reinitialize_states=reinit)
            raa = get_results(model, edata=edata, sensi_order=1,
                              sensi_meth=amici.SensitivityMethod.adjoint,
                              sensi_meth_preeq=amici.SensitivityMethod.adjoint,
                              reinitialize_states=reinit)

            # assert all are close
            assert np.isclose(rff['sllh'], rfa['sllh']).all()
            assert np.isclose(rfa['sllh'], raa['sllh']).all()
            assert np.isclose(raa['sllh'], rff['sllh']).all()
