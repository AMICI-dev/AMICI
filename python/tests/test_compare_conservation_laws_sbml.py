import os
import warnings

import amici
import numpy as np
import pytest
from amici import MeasurementChannel as MC
from amici.sim.sundials import (
    AMICI_ERROR,
    AMICI_SUCCESS,
    ExpData,
    SensitivityMethod,
    SteadyStateComputationMode,
    SteadyStateSensitivityMode,
    SteadyStateStatus,
    get_model_settings,
    run_simulation,
    set_model_settings,
)
from amici.testing import skip_on_valgrind
from numpy.testing import assert_allclose


@pytest.fixture
def edata_fixture():
    """edata is generated to test pre- and postequilibration"""
    edata_pre = ExpData(
        2, 0, 0, np.array([0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0])
    )
    edata_pre.set_measurements([1.5] * 16)
    edata_pre.fixed_parameters = np.array([5.0, 20.0])
    edata_pre.fixed_parameters_pre_equilibration = np.array([0.0, 10.0])
    edata_pre.reinitialize_fixed_parameter_initial_states = True

    # edata for postequilibration
    edata_post = ExpData(2, 0, 0, np.array([float("inf")] * 3))
    edata_post.set_measurements([0.75] * 6)
    edata_post.fixed_parameters = np.array([7.5, 30.0])

    # edata with both equilibrations
    edata_full = ExpData(
        2,
        0,
        0,
        np.array(
            [0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 4.0, float("inf"), float("inf")]
        ),
    )
    edata_full.set_measurements([3.14] * 18)
    edata_full.fixed_parameters = np.array([1.0, 2.0])
    edata_full.fixed_parameters_pre_equilibration = np.array([3.0, 4.0])
    edata_full.reinitialize_fixed_parameter_initial_states = True

    return edata_pre, edata_post, edata_full


@pytest.fixture(scope="session")
def models():
    from amici.importers.antimony import antimony2sbml

    ant_model = """model *model_constant_species()
  const compartment compartment_ = 1;
  var species substrate in compartment_ = 0;
  const species $enzyme in compartment_;
  var species complex in compartment_ = 0;
  species product in compartment_ = 0;

  var substrate_product := complex + product + substrate;

  creation:  => substrate; compartment_*(synthesis_substrate + k_create);
  binding: substrate + $enzyme -> complex; compartment_*(k_bind*substrate*enzyme - k_unbind*complex);
  conversion: complex => $enzyme + product; compartment_*k_convert*complex;
  decay: product => ; compartment_*k_decay*product;

  enzyme = init_enzyme;
  substrate = init_substrate;

  const init_enzyme = 2;
  const init_substrate = 5;
  const synthesis_substrate = 0;
  const k_bind = 10;
  const k_unbind = 1;
  const k_convert = 10;
  const k_create = 2;
  const k_decay = 1;

  unit volume = 1e-3 litre;
  unit substance = 1e-3 mole;
end"""
    sbml_importer = amici.SbmlImporter(
        antimony2sbml(ant_model), from_file=False
    )

    # Name of the model that will also be the name of the python module
    model_name = output_dir = "model_constant_species"
    model_name_cl = output_dir_cl = "model_constant_species_cl"

    # Define constants, observables, sigmas
    fixed_parameters = ["synthesis_substrate", "init_enzyme"]
    observation_model = [
        MC("observable_product", formula="product", sigma=1.0, name=""),
        MC("observable_substrate", formula="substrate", sigma=1.0, name=""),
    ]

    # wrap models with and without conservations laws
    sbml_importer.sbml2amici(
        model_name_cl,
        output_dir_cl,
        fixed_parameters=fixed_parameters,
        observation_model=observation_model,
    )
    sbml_importer.sbml2amici(
        model_name,
        output_dir,
        fixed_parameters=fixed_parameters,
        observation_model=observation_model,
        compute_conservation_laws=False,
    )

    # load both models
    model_without_cl_module = amici.import_model_module(
        model_name, module_path=os.path.abspath(model_name)
    )
    model_with_cl_module = amici.import_model_module(
        model_name_cl, module_path=os.path.abspath(model_name_cl)
    )

    # get the models and return
    model_without_cl = model_without_cl_module.get_model()
    model_with_cl = model_with_cl_module.get_model()
    return model_with_cl, model_without_cl


def get_results(
    model,
    edata=None,
    sensi_order=0,
    sensi_meth=SensitivityMethod.forward,
    sensi_meth_preeq=SensitivityMethod.forward,
    stst_mode=SteadyStateComputationMode.integrateIfNewtonFails,
    stst_sensi_mode=SteadyStateSensitivityMode.newtonOnly,
    reinitialize_states=False,
):
    # set model and data properties
    model.set_reinitialize_fixed_parameter_initial_states(reinitialize_states)

    # get the solver, set the properties
    solver = model.create_solver()
    solver.set_newton_max_steps(20)
    solver.set_sensitivity_order(sensi_order)
    solver.set_sensitivity_method_pre_equilibration(sensi_meth_preeq)
    solver.set_sensitivity_method(sensi_meth)
    model.set_steady_state_sensitivity_mode(stst_sensi_mode)
    model.set_steady_state_computation_mode(stst_mode)
    if edata is None:
        model.set_timepoints(np.linspace(0, 5, 101))
    else:
        edata.reinitialize_fixed_parameter_initial_states = reinitialize_states

    # return simulation results
    return run_simulation(model, solver, edata)


@skip_on_valgrind
def test_compare_conservation_laws_sbml(models, edata_fixture):
    # first, create the model
    model_with_cl, model_without_cl = models

    assert model_with_cl.ncl() > 0
    assert model_without_cl.nx_rdata == model_with_cl.nx_rdata
    assert model_with_cl.nx_solver < model_without_cl.nx_solver
    assert len(model_with_cl.get_state_ids_solver()) == model_with_cl.nx_solver
    assert (
        len(model_without_cl.get_state_ids_solver())
        == model_without_cl.nx_solver
    )

    # ----- compare simulations wo edata, sensi = 0, states ------------------
    # run simulations
    rdata_cl = get_results(model_with_cl)
    assert rdata_cl["status"] == AMICI_SUCCESS
    rdata = get_results(model_without_cl)
    assert rdata["status"] == AMICI_SUCCESS

    # compare state trajectories
    assert_allclose(
        rdata["x"],
        rdata_cl["x"],
        rtol=1.0e-5,
        atol=1.0e-8,
        err_msg="rdata.x mismatch",
    )

    # ----- compare simulations wo edata, sensi = 1, states and sensis -------
    # run simulations
    rdata_cl = get_results(model_with_cl, sensi_order=1)
    assert rdata_cl["status"] == AMICI_SUCCESS
    rdata = get_results(model_without_cl, sensi_order=1)
    assert rdata["status"] == AMICI_SUCCESS

    # compare state trajectories
    for field in ["x", "sx"]:
        assert_allclose(
            rdata[field],
            rdata_cl[field],
            rtol=1.0e-5,
            atol=1.0e-8,
            err_msg=f"rdata.{field} mismatch",
        )

    # ----- compare simulations wo edata, sensi = 0, states and sensis -------

    # run simulations
    edata, _, _ = edata_fixture
    rdata_cl = get_results(model_with_cl, edata=edata)
    assert rdata_cl["status"] == AMICI_SUCCESS
    rdata = get_results(model_without_cl, edata=edata)
    assert rdata["status"] == AMICI_SUCCESS

    # compare preequilibrated states
    for field in ["x", "x_ss", "llh"]:
        assert_allclose(
            rdata[field],
            rdata_cl[field],
            rtol=1.0e-5,
            atol=1.0e-8,
            err_msg=f"rdata.{field} mismatch",
        )

    # ----- compare simulations wo edata, sensi = 1, states and sensis -------

    # run simulations
    rdata_cl = get_results(model_with_cl, edata=edata, sensi_order=1)
    assert rdata_cl["status"] == AMICI_SUCCESS
    rdata = get_results(
        model_without_cl,
        edata=edata,
        sensi_order=1,
        stst_sensi_mode=SteadyStateSensitivityMode.integrateIfNewtonFails,
    )
    assert rdata["status"] == AMICI_SUCCESS
    # check that steady state computation succeeded only by sim in full model
    assert rdata["preeq_status"] == [
        SteadyStateStatus.failed_factorization,
        SteadyStateStatus.success,
        SteadyStateStatus.not_run,
    ]
    # check that steady state computation succeeded by Newton in reduced model
    assert rdata_cl["preeq_status"] == [
        SteadyStateStatus.success,
        SteadyStateStatus.not_run,
        SteadyStateStatus.not_run,
    ]

    # compare state sensitivities with edata and preequilibration
    for field in ["x", "x_ss", "sx", "llh", "sllh"]:
        assert_allclose(
            rdata[field],
            rdata_cl[field],
            rtol=1.0e-5,
            atol=1.0e-8,
            err_msg=f"rdata.{field} mismatch",
        )

    # ----- check failure st.st. sensi computation if run wo CLs -------------
    # check failure of steady state sensitivity computation if run wo CLs
    model_without_cl.set_steady_state_sensitivity_mode(
        SteadyStateSensitivityMode.newtonOnly
    )
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        rdata = get_results(model_without_cl, edata=edata, sensi_order=1)
        assert rdata["status"] == AMICI_ERROR


def test_adjoint_pre_and_post_equilibration(models, edata_fixture):
    # get both models
    model_cl, model = models

    # check gradient with and without state reinitialization
    for edata in edata_fixture:
        for reinit in [False, True]:
            # compare different ways of preequilibration, full rank Jacobian
            # forward preequilibration, forward simulation
            rff_cl = get_results(
                model_cl,
                edata=edata,
                sensi_order=1,
                sensi_meth=SensitivityMethod.forward,
                sensi_meth_preeq=SensitivityMethod.forward,
                reinitialize_states=reinit,
            )
            # forward preequilibration, adjoint simulation
            rfa_cl = get_results(
                model_cl,
                edata=edata,
                sensi_order=1,
                sensi_meth=SensitivityMethod.adjoint,
                sensi_meth_preeq=SensitivityMethod.forward,
                reinitialize_states=reinit,
            )
            # adjoint preequilibration, adjoint simulation
            raa_cl = get_results(
                model_cl,
                edata=edata,
                sensi_order=1,
                sensi_meth=SensitivityMethod.adjoint,
                sensi_meth_preeq=SensitivityMethod.adjoint,
                reinitialize_states=reinit,
            )

            assert rff_cl.status == AMICI_SUCCESS
            assert rfa_cl.status == AMICI_SUCCESS
            assert raa_cl.status == AMICI_SUCCESS

            # assert all are close
            assert_allclose(
                rff_cl["sllh"], rfa_cl["sllh"], rtol=1.0e-5, atol=1.0e-8
            )
            assert_allclose(
                rfa_cl["sllh"], raa_cl["sllh"], rtol=1.0e-5, atol=1.0e-8
            )
            assert_allclose(
                raa_cl["sllh"], rff_cl["sllh"], rtol=1.0e-5, atol=1.0e-8
            )

            # compare fully adjoint approach to simulation with singular
            #  Jacobian
            raa = get_results(
                model,
                edata=edata,
                sensi_order=1,
                sensi_meth=SensitivityMethod.adjoint,
                sensi_meth_preeq=SensitivityMethod.adjoint,
                stst_sensi_mode=SteadyStateSensitivityMode.integrateIfNewtonFails,
                reinitialize_states=reinit,
            )
            assert raa.status == AMICI_SUCCESS

            # assert gradients are close (quadrature tolerances are laxer)
            assert_allclose(raa_cl["sllh"], raa["sllh"], 1e-5, 1e-5)


@skip_on_valgrind
def test_get_set_model_settings(models):
    """test amici.sim.sundials.(get|set)_model_settings cycles for models with and without
    conservation laws"""

    for model in models:
        set_model_settings(model, get_model_settings(model))
