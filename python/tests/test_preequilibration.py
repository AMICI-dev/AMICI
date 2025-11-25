"""Tests for pre- and post-equilibration"""

import inspect
import itertools

import amici
import numpy as np
import pytest
from amici.debugging import get_model_for_preeq
from amici.importers.antimony import antimony2amici
from amici.sim.sundials import (
    AMICI_ERROR,
    AMICI_SUCCESS,
    ExpData,
    LogSeverity,
    ParameterScaling,
    SensitivityMethod,
    SensitivityOrder,
    SteadyStateComputationMode,
    SteadyStateSensitivityMode,
    SteadyStateStatus,
    parameter_scaling_from_int_vector,
    run_simulation,
)
from amici.sim.sundials.gradient_check import check_derivatives
from amici.testing import (
    TemporaryDirectoryWinSafe as TemporaryDirectory,
)
from amici.testing import (
    skip_on_valgrind,
)
from numpy.testing import assert_allclose, assert_equal
from test_pysb import get_data

pytestmark = pytest.mark.filterwarnings(
    # https://github.com/AMICI-dev/AMICI/issues/18
    "ignore:Adjoint sensitivity analysis for models with discontinuous "
    "right hand sides .*:UserWarning",
)


@pytest.fixture
def preeq_fixture(pysb_example_presimulation_module):
    model = pysb_example_presimulation_module.get_model()
    model.set_reinitialize_fixed_parameter_initial_states(True)
    model.set_steady_state_computation_mode(
        SteadyStateComputationMode.integrateIfNewtonFails
    )
    model.set_steady_state_sensitivity_mode(
        SteadyStateSensitivityMode.integrateIfNewtonFails
    )
    solver = model.create_solver()
    solver.set_sensitivity_order(SensitivityOrder.first)
    solver.set_sensitivity_method(SensitivityMethod.forward)

    edata = get_data(model)
    edata.t_presim = 2
    edata.fixed_parameters = [10, 2]
    edata.fixed_parameters_presimulation = [3, 2]
    edata.fixed_parameters_pre_equilibration = [3, 0]
    edata.set_timepoints([1, 5])

    edata_preeq = ExpData(edata)
    edata_preeq.t_presim = 0
    edata_preeq.set_timepoints([np.inf])
    edata_preeq.fixed_parameters = edata.fixed_parameters_pre_equilibration
    edata_preeq.fixed_parameters_presimulation = ()
    edata_preeq.fixed_parameters_pre_equilibration = ()

    edata_presim = ExpData(edata)
    edata_presim.t_presim = 0
    edata_presim.set_timepoints([edata.t_presim])
    edata_presim.fixed_parameters = edata.fixed_parameters_presimulation
    edata_presim.fixed_parameters_presimulation = ()
    edata_presim.fixed_parameters_pre_equilibration = ()

    edata_sim = ExpData(edata)
    edata_sim.t_presim = 0
    edata_sim.set_timepoints(edata.get_timepoints())
    edata_sim.fixed_parameters = edata.fixed_parameters
    edata_sim.fixed_parameters_presimulation = ()
    edata_sim.fixed_parameters_pre_equilibration = ()

    pscales = [
        ParameterScaling.log10,
        ParameterScaling.ln,
        ParameterScaling.none,
        parameter_scaling_from_int_vector(
            [
                ParameterScaling.log10,
                ParameterScaling.ln,
                ParameterScaling.none,
                ParameterScaling.log10,
                ParameterScaling.ln,
                ParameterScaling.none,
            ]
        ),
    ]

    plists = [
        [3, 1, 2, 4],
        [0, 1, 2, 3, 4, 5],
        [5, 3, 2, 0, 4, 1],
        [1, 2, 3, 4, 5],
        [1, 1, 1],
    ]

    return (
        model,
        solver,
        edata,
        edata_preeq,
        edata_presim,
        edata_sim,
        pscales,
        plists,
    )


def test_manual_preequilibration(preeq_fixture):
    """Manual preequilibration"""

    (
        model,
        solver,
        edata,
        edata_preeq,
        edata_presim,
        edata_sim,
        pscales,
        plists,
    ) = preeq_fixture

    settings = itertools.product(pscales, plists)

    for pscale, plist in settings:
        model.set_initial_state([])
        model.set_initial_state_sensitivities([])
        model.set_parameter_list(plist)
        model.set_parameter_scale(pscale)

        # combined
        rdata_auto = run_simulation(model, solver, edata)
        assert rdata_auto.status == AMICI_SUCCESS

        # manual preequilibration
        rdata_preeq = run_simulation(model, solver, edata_preeq)
        assert rdata_preeq.status == AMICI_SUCCESS

        # manual reinitialization + presimulation
        x0 = rdata_preeq["x"][0, :]
        x0[1] = edata_presim.fixed_parameters[0]
        x0[2] = edata_presim.fixed_parameters[1]
        sx0 = rdata_preeq["sx"][0, :, :]
        sx0[:, 1] = 0
        sx0[:, 2] = 0
        model.set_initial_state(x0)
        model.set_initial_state_sensitivities(sx0.flatten())
        rdata_presim = run_simulation(model, solver, edata_presim)
        assert rdata_presim.status == AMICI_SUCCESS

        # manual reinitialization + simulation
        x0 = rdata_presim["x"][0, :]
        x0[1] = edata_sim.fixed_parameters[0]
        x0[2] = edata_sim.fixed_parameters[1]
        sx0 = rdata_presim["sx"][0, :, :]
        sx0[:, 1] = 0
        sx0[:, 2] = 0
        model.set_initial_state(x0)
        model.set_initial_state_sensitivities(sx0.flatten())
        rdata_sim = run_simulation(model, solver, edata_sim)
        assert rdata_sim.status == AMICI_SUCCESS

        for variable in ["x", "sx"]:
            assert_allclose(
                rdata_auto[variable],
                rdata_sim[variable],
                atol=1e-6,
                rtol=1e-6,
                err_msg=str(
                    dict(pscale=pscale, plist=plist, variable=variable)
                ),
            )


def test_parameter_reordering(preeq_fixture):
    """Test parameter reordering"""

    (
        model,
        solver,
        edata,
        edata_preeq,
        edata_presim,
        edata_sim,
        pscales,
        plists,
    ) = preeq_fixture

    rdata_ordered = run_simulation(model, solver, edata)

    for plist in plists:
        model.set_parameter_list(plist)
        rdata_reordered = run_simulation(model, solver, edata)

        for ip, p_index in enumerate(plist):
            assert_allclose(
                rdata_ordered["sx"][:, p_index, :],
                rdata_reordered["sx"][:, ip, :],
                atol=1e-6,
                rtol=1e-6,
                err_msg=str(dict(variable="sx", plist=plist, p_index=p_index)),
            )


def test_data_replicates(preeq_fixture):
    """Test data replicates"""

    (
        model,
        solver,
        edata,
        edata_preeq,
        edata_presim,
        edata_sim,
        pscales,
        plists,
    ) = preeq_fixture

    sensi_meth = SensitivityMethod.forward
    solver.set_sensitivity_method(sensi_meth)

    # add infty timepoint
    y = edata.get_observed_data()
    stdy = edata.get_observed_data_std_dev()
    ts = np.hstack([*edata.get_timepoints(), np.inf])
    edata.set_timepoints(sorted(ts))
    edata.set_observed_data(np.hstack([y, y[0]]))
    edata.set_observed_data_std_dev(np.hstack([stdy, stdy[0]]))
    rdata_single = run_simulation(model, solver, edata)

    # duplicate data and timepoints
    y = edata.get_observed_data()
    stdy = edata.get_observed_data_std_dev()
    ts = np.hstack([*edata.get_timepoints(), *edata.get_timepoints()])
    idx = np.argsort(ts)
    edata.set_timepoints(sorted(ts))
    edata.set_observed_data(np.hstack([y, y])[idx])
    edata.set_observed_data_std_dev(np.hstack([stdy, stdy])[idx])

    rdata_double = run_simulation(model, solver, edata)

    for variable in ["llh", "sllh"]:
        assert_allclose(
            2 * rdata_single[variable],
            rdata_double[variable],
            atol=1e-6,
            rtol=1e-6,
            err_msg=str(dict(variable=variable, sensi_meth=sensi_meth)),
        )


def test_parameter_in_expdata(preeq_fixture):
    """Test parameter in ExpData"""

    (
        model,
        solver,
        edata,
        edata_preeq,
        edata_presim,
        edata_sim,
        pscales,
        plists,
    ) = preeq_fixture

    rdata = run_simulation(model, solver, edata)

    # get initial states will compute initial states if nothing is set,
    # this needs go first as we need unmodified model. Also set to
    # preequilibration fixpars first as this is where initial states would be
    # computed otherwise
    model.set_fixed_parameters(edata.fixed_parameters_pre_equilibration)
    edata.x0 = model.get_initial_state()
    edata.sx0 = model.get_initial_state_sensitivities()

    # perturb model initial states
    model.set_initial_state(rdata["x_ss"] * 4)
    model.set_initial_state_sensitivities(rdata["sx_ss"].flatten() / 2)

    # set ExpData plist
    edata.plist = model.get_parameter_list()
    # perturb model parameter list
    model.set_parameter_list([i for i in reversed(model.get_parameter_list())])

    # set ExpData parameters
    edata.free_parameters = model.get_free_parameters()
    # perturb model parameters
    model.set_free_parameters(
        tuple(p * 2 for p in model.get_free_parameters())
    )

    # set ExpData pscale
    edata.pscale = model.get_parameter_scale()
    # perturb model pscale, needs to be done after getting parameters,
    # otherwise we will mess up parameter value
    model.set_parameter_scale(
        parameter_scaling_from_int_vector(
            [
                ParameterScaling.log10
                if scaling == ParameterScaling.none
                else ParameterScaling.none
                for scaling in model.get_parameter_scale()
            ]
        )
    )

    rdata_edata = run_simulation(model, solver, edata)
    for variable in ["x", "sx"]:
        assert_allclose(
            rdata[variable][0, :],
            rdata_edata[variable][0, :],
            atol=1e-6,
            rtol=1e-6,
            err_msg=str(dict(variable=variable)),
        )


def test_equilibration_methods_with_adjoints(preeq_fixture):
    """Test different combinations of equilibration and simulation
    sensitivity methods"""

    (
        model,
        solver,
        edata,
        edata_preeq,
        edata_presim,
        edata_sim,
        pscales,
        plists,
    ) = preeq_fixture

    # we don't want presim
    edata.t_presim = 0.0
    edata.fixed_parameters_presimulation = ()

    # add infty timepoint
    y = edata.get_observed_data()
    stdy = edata.get_observed_data_std_dev()
    ts = np.hstack([*edata.get_timepoints(), np.inf])
    edata.set_timepoints(sorted(ts))
    edata.set_observed_data(np.hstack([y, y[0]]))
    edata.set_observed_data_std_dev(np.hstack([stdy, stdy[0]]))

    rdatas = {}
    equil_meths = [
        SteadyStateSensitivityMode.newtonOnly,
        SteadyStateSensitivityMode.integrationOnly,
        SteadyStateSensitivityMode.integrateIfNewtonFails,
    ]
    sensi_meths = [
        SensitivityMethod.forward,
        SensitivityMethod.adjoint,
    ]
    settings = itertools.product(equil_meths, sensi_meths)

    for setting in settings:
        # unpack, solver settings
        equil_meth, sensi_meth = setting
        model.set_steady_state_sensitivity_mode(equil_meth)
        solver.set_sensitivity_method(sensi_meth)
        solver.set_newton_max_steps(0)

        # add rdatas
        rdatas[setting] = run_simulation(model, solver, edata)
        # assert successful simulation

        assert rdatas[setting]["status"] == AMICI_SUCCESS

    for setting1, setting2 in itertools.product(settings, settings):
        # assert correctness of result
        for variable in ["llh", "sllh"]:
            assert_allclose(
                rdatas[setting1][variable],
                rdatas[setting2][variable],
                atol=1e-6,
                rtol=1e-6,
                err_msg=str(
                    dict(
                        variable=variable, setting1=setting1, setting2=setting2
                    )
                ),
            )


def test_newton_solver_equilibration(preeq_fixture):
    """Test newton solver for equilibration"""

    (
        model,
        solver,
        edata,
        edata_preeq,
        edata_presim,
        edata_sim,
        pscales,
        plists,
    ) = preeq_fixture

    # we don't want presim
    edata.t_presim = 0.0
    edata.fixed_parameters_presimulation = ()

    # add infty timepoint
    y = edata.get_observed_data()
    stdy = edata.get_observed_data_std_dev()
    ts = np.hstack([*edata.get_timepoints(), np.inf])
    edata.set_timepoints(sorted(ts))
    edata.set_observed_data(np.hstack([y, y[0]]))
    edata.set_observed_data_std_dev(np.hstack([stdy, stdy[0]]))

    rdatas = {}
    settings = [
        SteadyStateSensitivityMode.integrationOnly,
        SteadyStateSensitivityMode.newtonOnly,
    ]

    solver.set_newton_step_steady_state_check(True)
    solver.set_relative_tolerance_steady_state(1e-12)

    for equil_meth in settings:
        # set sensi method
        sensi_meth = SensitivityMethod.forward
        solver.set_sensitivity_method(sensi_meth)
        model.set_steady_state_sensitivity_mode(equil_meth)
        if equil_meth == SteadyStateSensitivityMode.newtonOnly:
            solver.set_newton_max_steps(10)

        # add rdatas
        rdatas[equil_meth] = run_simulation(model, solver, edata)

        # assert successful simulation
        assert rdatas[equil_meth]["status"] == AMICI_SUCCESS

    # assert correct results
    for variable in ["llh", "sllh", "sx0", "sx_ss", "x_ss"]:
        assert_allclose(
            rdatas[settings[0]][variable],
            rdatas[settings[1]][variable],
            atol=1e-5,
            rtol=1e-5,
            err_msg=str(dict(variable=variable)),
        )


def test_newton_steadystate_check(preeq_fixture):
    """Test NewtonStepSteadyStateCheck solver flag"""

    (
        model,
        solver,
        edata,
        edata_preeq,
        edata_presim,
        edata_sim,
        pscales,
        plists,
    ) = preeq_fixture

    # we don't want presim
    edata.t_presim = 0.0
    edata.fixed_parameters_presimulation = ()

    # add infty timepoint
    y = edata.get_observed_data()
    stdy = edata.get_observed_data_std_dev()
    ts = np.hstack([*edata.get_timepoints(), np.inf])
    edata.set_timepoints(sorted(ts))
    edata.set_observed_data(np.hstack([y, y[0]]))
    edata.set_observed_data_std_dev(np.hstack([stdy, stdy[0]]))

    # set sensi method
    sensi_meth = SensitivityMethod.forward
    solver.set_sensitivity_method(sensi_meth)

    solver.set_newton_max_steps(100)

    rdatas = {}
    for newton_check in [True, False]:
        solver.set_newton_step_steady_state_check(newton_check)

        # add rdatas
        rdatas[newton_check] = run_simulation(model, solver, edata)

        # assert successful simulation
        assert rdatas[newton_check]["status"] == AMICI_SUCCESS

    # assert correct results
    for variable in ["x_ss", "llh", "sx0", "sx_ss", "sllh"]:
        assert_allclose(
            rdatas[True][variable],
            rdatas[False][variable],
            atol=1e-6,
            rtol=1e-6,
            err_msg=str(dict(variable=variable, sensi_meth=sensi_meth)),
        )


def test_steadystate_computation_mode(preeq_fixture):
    """Test newtonOnly and integrationOnly steady-state computation modes"""
    (
        model,
        solver,
        edata,
        edata_preeq,
        edata_presim,
        edata_sim,
        pscales,
        plists,
    ) = preeq_fixture

    sensi_meth = SensitivityMethod.forward
    solver.set_sensitivity_order(SensitivityOrder.first)
    solver.set_sensitivity_method_pre_equilibration(sensi_meth)
    solver.set_newton_max_steps(10)

    rdatas = {}
    stst_computation_modes = [
        SteadyStateComputationMode.integrationOnly,
        SteadyStateComputationMode.newtonOnly,
    ]
    for mode in stst_computation_modes:
        model.set_steady_state_computation_mode(mode)
        rdatas[mode] = run_simulation(model, solver, edata)

        # assert successful simulation
        assert rdatas[mode]["status"] == AMICI_SUCCESS
    assert rdatas[SteadyStateComputationMode.integrationOnly][
        "preeq_status"
    ] == [
        SteadyStateStatus.not_run,
        SteadyStateStatus.success,
        SteadyStateStatus.not_run,
    ]

    assert (
        rdatas[SteadyStateComputationMode.integrationOnly]["preeq_numsteps"][0]
        == 0
    )

    assert rdatas[SteadyStateComputationMode.newtonOnly]["preeq_status"] == [
        SteadyStateStatus.success,
        SteadyStateStatus.not_run,
        SteadyStateStatus.not_run,
    ]
    assert (
        rdatas[SteadyStateComputationMode.newtonOnly]["preeq_numsteps"][0] > 0
    )

    # assert correct results
    for variable in ["llh", "sllh", "sx0", "sx_ss", "x_ss"]:
        assert_allclose(
            rdatas[stst_computation_modes[0]][variable],
            rdatas[stst_computation_modes[1]][variable],
            atol=1e-5,
            rtol=1e-5,
            err_msg=str(dict(variable=variable, sensi_meth=sensi_meth)),
        )


def test_simulation_errors(preeq_fixture):
    (
        model,
        solver,
        edata,
        edata_preeq,
        edata_presim,
        edata_sim,
        pscales,
        plists,
    ) = preeq_fixture

    solver.set_sensitivity_order(SensitivityOrder.first)
    solver.set_sensitivity_method_pre_equilibration(SensitivityMethod.forward)
    model.set_steady_state_sensitivity_mode(
        SteadyStateSensitivityMode.integrationOnly
    )
    solver.set_max_steps(1)

    # exceeded maxsteps
    # preeq & posteq
    for e in [edata, edata_preeq]:
        rdata = run_simulation(model, solver, e)
        assert rdata["status"] != AMICI_SUCCESS
        assert rdata._swigptr.messages[0].severity == LogSeverity.debug
        assert rdata._swigptr.messages[0].identifier == "EQUILIBRATION_FAILURE"
        assert (
            "exceeded maximum number of integration steps"
            in rdata._swigptr.messages[0].message
        )
        assert rdata._swigptr.messages[1].severity == LogSeverity.error
        assert rdata._swigptr.messages[1].identifier == "OTHER"


def test_t_overflow():
    """Test that equilibration fails with an informative error message
    upon time overflow."""
    module_name = inspect.stack()[0].function
    ant_str = f"""
    model {module_name}
        # Constant increase so the solver will take large steps;
        #  small enough to let `t` to overflow before `x`.
        dxx_dt = 1e-16
        xx' = dxx_dt
        xx = 0
    end
    """
    with TemporaryDirectory(prefix=module_name) as outdir:
        model = antimony2amici(
            ant_str,
            model_name=module_name,
            output_dir=outdir,
            fixed_parameters=["dxx_dt"],
        )
        model.set_steady_state_computation_mode(
            SteadyStateComputationMode.integrationOnly
        )

        # require simulation until forever
        solver = model.create_solver()
        solver.set_relative_tolerance_steady_state(0)
        solver.set_absolute_tolerance_steady_state(0)

        edata_preeq = ExpData(model)
        edata_preeq.set_timepoints([np.inf])
        edata_preeq.fixed_parameters_pre_equilibration = [1e-16]
        edata_posteq = ExpData(model)
        edata_posteq.set_timepoints([float("inf")])

        for edata in (edata_preeq, edata_posteq):
            rdata = run_simulation(model, solver, edata=edata)
            assert rdata.status == AMICI_ERROR
            for msg in rdata.messages:
                if "exceedingly long simulation time" in msg.message:
                    assert msg.identifier == "EQUILIBRATION_FAILURE"
                    break
            else:
                assert False, list(rdata.messages)


def test_get_model_for_preeq(preeq_fixture):
    (
        model,
        solver,
        edata,
        edata_preeq,
        edata_presim,
        edata_sim,
        pscales,
        plists,
    ) = preeq_fixture
    model.set_steady_state_sensitivity_mode(
        SteadyStateSensitivityMode.integrationOnly
    )
    model_preeq = get_model_for_preeq(model, edata)
    # the exactly same settings are used, so results should match exactly
    rdata1 = run_simulation(model_preeq, solver)
    rdata2 = run_simulation(model, solver, edata_preeq)
    assert_equal(
        rdata1.x,
        rdata2.x,
    )
    assert_equal(
        rdata1.sx,
        rdata2.sx,
    )


@skip_on_valgrind
def test_partial_eq():
    """Check that partial equilibration is possible."""
    from amici.importers.antimony import antimony2amici

    ant_str = """
    model test_partial_eq
        explodes = 1
        explodes' = explodes
        A = 1
        B = 0
        R: A -> B; k*A - k*B
        k = 1
    end
    """
    module_name = "test_partial_eq"
    with TemporaryDirectory(prefix=module_name) as outdir:
        antimony2amici(
            ant_str,
            model_name=module_name,
            output_dir=outdir,
        )
        model_module = amici.import_model_module(
            module_name=module_name, module_path=outdir
        )
        amici_model = model_module.get_model()
        amici_model.set_timepoints([np.inf])
        amici_solver = amici_model.create_solver()
        amici_solver.set_relative_tolerance_steady_state(1e-12)

        # equilibration of `explodes` will fail
        rdata = run_simulation(amici_model, amici_solver)
        assert rdata.status == AMICI_ERROR
        assert rdata.messages[0].identifier == "EQUILIBRATION_FAILURE"

        # excluding `explodes` should enable equilibration
        amici_model.set_steadystate_mask(
            [
                0 if state_id == "explodes" else 1
                for state_id in amici_model.get_state_ids_solver()
            ]
        )
        rdata = run_simulation(amici_model, amici_solver)
        assert rdata.status == AMICI_SUCCESS
        assert_allclose(
            rdata.by_id("A"),
            0.5,
            atol=amici_solver.get_absolute_tolerance_steady_state(),
            rtol=amici_solver.get_relative_tolerance_steady_state(),
        )
        assert_allclose(
            rdata.by_id("B"),
            0.5,
            atol=amici_solver.get_absolute_tolerance_steady_state(),
            rtol=amici_solver.get_relative_tolerance_steady_state(),
        )
        assert rdata.t_last < 100


def test_preequilibration_t0(tempdir):
    """Test that preequilibration uses the correct initial time."""
    from amici.importers.antimony import antimony2amici

    ant_str = """
    model test_preequilibration_t0
        preeq_indicator = 0
        t0_preeq = time * preeq_indicator
        # we need some state variable for simulation to work
        T = 0
        # at some large value of T, we will "reach steady state" due to
        # the tiny relative change
        T' = 1
    end
    """
    module_name = "test_preequilibration_t0"
    antimony2amici(
        ant_str,
        model_name=module_name,
        output_dir=tempdir,
        fixed_parameters=["preeq_indicator"],
    )
    model_module = amici.import_model_module(
        module_name=module_name, module_path=tempdir
    )
    amici_model = model_module.get_model()
    edata = ExpData(amici_model)
    edata.set_timepoints([0.0, 10_000.0])
    edata.fixed_parameters_pre_equilibration = [1.0]
    edata.fixed_parameters = [0.0]
    amici_model.set_t0_preeq(-10_000.0)
    amici_model.set_t0(-2.0)
    amici_solver = amici_model.create_solver()
    amici_solver.set_relative_tolerance_steady_state(1e-5)
    amici_model.set_steady_state_computation_mode(
        SteadyStateComputationMode.integrationOnly
    )

    rdata = run_simulation(amici_model, amici_solver, edata)
    assert rdata.status == AMICI_SUCCESS
    assert set(rdata.by_id("t0_preeq")) == {-10_000.0}
    idx_time_integral = amici_model.get_state_ids().index("T")
    assert np.isclose(
        rdata.x_ss[idx_time_integral], rdata.preeq_t - amici_model.t0_preeq()
    )


def test_preequilibration_events(tempdir):
    """Test that events are handled correctly during preequilibration."""
    from amici.importers.antimony import antimony2amici

    ant_str = """
    model test_preequilibration_events
        some_time = 0
        some_time' = 1
        target1 = 0
        target2 = 0
        target3 = 0
        target4 = 0
        is_preeq = 0
        bolus1 = 1
        bolus2 = 1
        bolus3 = 1
        bolus4 = 1
        # E1 & E2 will only trigger during pre-equilibration
        # (initialValue is only used once at the start of pre-equilibration)
        E1: at some_time >= 0, t0 = false: target1 = target1 + bolus1
        E2: at time >= 0, t0 = false: target2 = target2 + bolus2
        # requires early time point
        # https://github.com/AMICI-dev/AMICI/issues/2804
        trigger_time2 = 1e-3
        # E3 will trigger only during preequilibration
        # (some_time is not reset and trigger initial value is `true`)
        E3: at some_time >= trigger_time2: target3 = target3 + bolus3
        # will trigger during preequilibration and main simulation
        E4: at time >= trigger_time2: target4 = target4 + bolus4
    end
    """
    module_name = "test_preequilibration_events"
    antimony2amici(
        ant_str,
        model_name=module_name,
        output_dir=tempdir,
        fixed_parameters=["is_preeq"],
    )
    model_module = amici.import_model_module(
        module_name=module_name, module_path=tempdir
    )
    amici_model = model_module.get_model()
    target1_idx = amici_model.get_state_ids().index("target1")
    target2_idx = amici_model.get_state_ids().index("target2")
    target3_idx = amici_model.get_state_ids().index("target3")
    target4_idx = amici_model.get_state_ids().index("target4")

    amici_model.set_timepoints([0, 11])
    amici_solver = amici_model.create_solver()
    amici_solver.set_newton_max_steps(10**5)

    edata = ExpData(amici_model)
    edata.set_timepoints([0, 11])
    edata.fixed_parameters_pre_equilibration = [1.0]
    edata.fixed_parameters = [0.0]

    # Integration-only preequilibration should handle all events
    amici_model = model_module.get_model()
    amici_model.set_steady_state_sensitivity_mode(
        SteadyStateSensitivityMode.integrationOnly
    )
    amici_model.set_steady_state_computation_mode(
        SteadyStateSensitivityMode.integrationOnly
    )
    rdata = run_simulation(amici_model, amici_solver, edata)
    assert rdata.status == AMICI_SUCCESS
    assert rdata.preeq_t > 1e-3  # verifies that integration was done
    assert rdata.x_ss[target1_idx] == 1
    assert rdata.x_ss[target2_idx] == 1
    assert rdata.x_ss[target3_idx] == 1
    assert rdata.x_ss[target4_idx] == 1
    assert np.all(rdata.x[:, target1_idx] == [1, 1])
    assert np.all(rdata.x[:, target2_idx] == [1, 1])
    assert np.all(rdata.x[:, target3_idx] == [1, 1])
    assert np.all(rdata.x[:, target4_idx] == [1, 2])

    edata = ExpData(rdata, 1.0, 1.0, 1)
    edata.fixed_parameters_pre_equilibration = [1.0]
    edata.fixed_parameters = [0.0]

    amici_solver.set_sensitivity_order(SensitivityOrder.first)

    for sensi_meth, sensi_meth_preeq in (
        (SensitivityMethod.forward, SensitivityMethod.forward),
        (SensitivityMethod.adjoint, SensitivityMethod.forward),
        (SensitivityMethod.adjoint, SensitivityMethod.adjoint),
    ):
        amici_solver.set_sensitivity_method(sensi_meth)
        amici_solver.set_sensitivity_method_pre_equilibration(sensi_meth_preeq)

        # amici_model.requireSensitivitiesForAllParameters()
        # FIXME: finite differences w.r.t. trigger time are off
        #   need different epsilon for trigger time
        amici_model.set_parameter_list(
            [
                i
                for i, p in enumerate(amici_model.get_free_parameter_ids())
                if p != "trigger_time2"
            ]
        )
        rdata = run_simulation(amici_model, amici_solver, edata)
        assert rdata.status == AMICI_SUCCESS

        check_derivatives(
            amici_model,
            amici_solver,
            edata=edata,
            atol=1e-6,
            rtol=1e-6,
            epsilon=1e-8,
        )
