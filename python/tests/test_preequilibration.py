"""Tests for pre- and post-equilibration"""

import inspect
import itertools

import amici
import numpy as np
import pytest
from amici.adapters.fiddy import run_simulation_to_function_and_derivative
from amici.importers.antimony import antimony2amici
from amici.importers.utils import MeasurementChannel
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
    get_model_for_preeq,
    parameter_scaling_from_int_vector,
    run_simulation,
)
from amici.testing import (
    TemporaryDirectoryWinSafe as TemporaryDirectory,
)
from amici.testing import (
    skip_on_valgrind,
)
from fiddy import check_gradient, check_jacobian
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
    y = edata.get_measurements()
    stdy = edata.get_noise_scales()
    ts = np.hstack([*edata.get_timepoints(), np.inf])
    edata.set_timepoints(sorted(ts))
    edata.set_measurements(np.hstack([y, y[0]]))
    edata.set_noise_scales(np.hstack([stdy, stdy[0]]))
    rdata_single = run_simulation(model, solver, edata)

    # duplicate data and timepoints
    y = edata.get_measurements()
    stdy = edata.get_noise_scales()
    ts = np.hstack([*edata.get_timepoints(), *edata.get_timepoints()])
    idx = np.argsort(ts)
    edata.set_timepoints(sorted(ts))
    edata.set_measurements(np.hstack([y, y])[idx])
    edata.set_noise_scales(np.hstack([stdy, stdy])[idx])

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
    y = edata.get_measurements()
    stdy = edata.get_noise_scales()
    ts = np.hstack([*edata.get_timepoints(), np.inf])
    edata.set_timepoints(sorted(ts))
    edata.set_measurements(np.hstack([y, y[0]]))
    edata.set_noise_scales(np.hstack([stdy, stdy[0]]))

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
    y = edata.get_measurements()
    stdy = edata.get_noise_scales()
    ts = np.hstack([*edata.get_timepoints(), np.inf])
    edata.set_timepoints(sorted(ts))
    edata.set_measurements(np.hstack([y, y[0]]))
    edata.set_noise_scales(np.hstack([stdy, stdy[0]]))

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
    y = edata.get_measurements()
    stdy = edata.get_noise_scales()
    ts = np.hstack([*edata.get_timepoints(), np.inf])
    edata.set_timepoints(sorted(ts))
    edata.set_measurements(np.hstack([y, y[0]]))
    edata.set_noise_scales(np.hstack([stdy, stdy[0]]))

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


def test_newton_only_without_sensitivities(preeq_fixture):
    """Test that steady-state computation mode `newtonOnly` and steady-state
    sensitivity mode `integrationOnly` are only rejected with sensitivities."""
    model, solver, edata, *_ = preeq_fixture

    model.set_steady_state_computation_mode(
        SteadyStateComputationMode.newtonOnly
    )
    model.set_steady_state_sensitivity_mode(
        SteadyStateSensitivityMode.integrationOnly
    )
    solver.set_sensitivity_method(SensitivityMethod.forward)
    solver.set_newton_max_steps(10)

    solver.set_sensitivity_order(SensitivityOrder.none)
    rdata = run_simulation(model, solver, edata)
    assert rdata.status == AMICI_SUCCESS
    assert rdata.preeq_status[0] == SteadyStateStatus.success

    solver.set_sensitivity_order(SensitivityOrder.first)
    rdata = run_simulation(model, solver, edata)
    assert rdata.status != AMICI_SUCCESS


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
        dx_dt = 1e-16
        x' = dx_dt
        x = 0
    end
    """
    with TemporaryDirectory(prefix=module_name) as outdir:
        model = antimony2amici(
            ant_str,
            model_name=module_name,
            output_dir=outdir,
            fixed_parameters=["dx_dt"],
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
    edata.t_start = -2.0
    edata.t_start_preeq = -10_000.0
    amici_solver = amici_model.create_solver()
    amici_solver.set_relative_tolerance_steady_state(1e-5)
    amici_model.set_steady_state_computation_mode(
        SteadyStateComputationMode.integrationOnly
    )
    idx_time_integral = amici_model.get_state_ids().index("T")

    def check_preeq_start(t0_preeq):
        rdata = run_simulation(amici_model, amici_solver, edata)
        assert rdata.status == AMICI_SUCCESS
        assert set(rdata.by_id("t0_preeq")) == {t0_preeq}
        assert np.isclose(
            rdata.x_ss[idx_time_integral], rdata.preeq_t - t0_preeq
        )
        return rdata

    check_preeq_start(-10_000.0)

    # without ExpData.t_start_preeq, pre-equilibration starts at t_start
    edata.t_start_preeq = np.nan
    check_preeq_start(-2.0)

    # deprecated: Model.t0_preeq() is used if ExpData.t_start_preeq is NaN
    with pytest.warns(DeprecationWarning, match="set_t0_preeq"):
        amici_model.set_t0_preeq(-5_000.0)
    rdata = check_preeq_start(-5_000.0)
    assert any(msg.identifier == "DEPRECATION" for msg in rdata.messages)

    # ExpData.t_start_preeq takes precedence over Model.t0_preeq()
    edata.t_start_preeq = -10_000.0
    rdata = check_preeq_start(-10_000.0)
    assert not any(msg.identifier == "DEPRECATION" for msg in rdata.messages)
    assert amici_model.t0_preeq() == -5_000.0


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
    # Exclude the unbounded `some_time` clock -- not a meaningful sensitivity
    # target.
    observation_model = [
        MeasurementChannel(
            id_=f"obs_{target}",
            formula=target,
            noise_distribution="normal",
            sigma=1.0,
        )
        for target in ("target1", "target2", "target3", "target4")
    ]
    antimony2amici(
        ant_str,
        model_name=module_name,
        output_dir=tempdir,
        fixed_parameters=["is_preeq"],
        observation_model=observation_model,
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
        amici_model.require_sensitivities_for_all_parameters()

        rdata = run_simulation(amici_model, amici_solver, edata)
        assert rdata.status == AMICI_SUCCESS

        # `x`/`x0`/`x_ss` still expose `some_time`, unlike the observables.
        point = amici_model.get_free_parameters()
        # Adjoint sensitivities only provide `sllh`. `check_gradient`
        # differentiates the first output of `function`, so `function` must
        # not return anything but `llh` in that case.
        derivative_variables = (
            ["llh"]
            if sensi_meth == SensitivityMethod.adjoint
            else ["y", "res", "llh"]
        )
        function, derivative = run_simulation_to_function_and_derivative(
            amici_model=amici_model,
            amici_solver=amici_solver,
            amici_edata=edata,
            derivative_variables=derivative_variables,
        )
        expected = derivative(point)
        if list(expected) == ["llh"]:
            result = check_gradient(function, point, expected["llh"])
        else:
            result = check_jacobian(function, point, expected)
        result.assert_success(always_print=True)

        # sx0 is not computed (and not meaningful) when pre-equilibration
        # is run with adjoint sensitivities (see gh-1184)
        if sensi_meth_preeq == SensitivityMethod.adjoint:
            assert rdata.sx0 is None
        else:
            assert rdata.sx0 is not None


@pytest.fixture(scope="module")
def reinit_model_module(tmp_path_factory):
    """Model with fixed-parameter-dependent initial states and an event."""
    ant_str = """
    model test_preequilibration_reinit_adjoint_sensitivities
        kinit = 5
        tev = 1e6
        kdecay = 0.7
        ksyn = 2
        scale = 1.5
        c = 0.4
        x1 = kinit * scale
        x2 = kinit
        x1' = ksyn - kdecay * x1
        x2' = c * x1 - kdecay * x2
        E1: at time >= tev: x2 = x2 + 1
    end
    """
    module_name = "test_preequilibration_reinit_adjoint_sensitivities"
    output_dir = tmp_path_factory.mktemp(module_name)
    antimony2amici(
        ant_str,
        model_name=module_name,
        output_dir=output_dir,
        fixed_parameters=["kinit", "tev"],
    )
    return amici.import_model_module(
        module_name=module_name, module_path=output_dir
    )


@skip_on_valgrind
@pytest.mark.parametrize("preeq_event", [False, True])
@pytest.mark.parametrize("reinit_states", [None, ["x1"], ["x2"]])
def test_preequilibration_reinit_adjoint_sensitivities(
    reinit_model_module, reinit_states, preeq_event
):
    """Sensitivities with reinitialization after pre-equilibration.

    Adjoint sensitivities for the main simulation with forward (GH#3278) or
    adjoint (GH#1156) pre-equilibration, with reinitialization of all states
    or a subset of states. With an event during pre-equilibration, the
    pre-equilibration backward problem is integrated, otherwise the
    steady-state shortcut is used.
    """
    model = reinit_model_module.get_model()
    assert model.nx_reinit() == 2
    model.set_timepoints([0.0, 1.0, 2.0])
    if preeq_event:
        model.set_steady_state_computation_mode(
            SteadyStateComputationMode.integrationOnly
        )
    solver = model.create_solver()
    solver.set_sensitivity_order(SensitivityOrder.first)
    solver.set_relative_tolerance(1e-10)
    solver.set_absolute_tolerance(1e-14)
    solver.set_relative_tolerance_quadratures(1e-10)
    solver.set_absolute_tolerance_quadratures(1e-14)

    def fixed_parameters(**kwargs):
        return [kwargs[id_] for id_ in model.get_fixed_parameter_ids()]

    edata = ExpData(model)
    edata.fixed_parameters = fixed_parameters(kinit=10.0, tev=1e6)
    edata.fixed_parameters_pre_equilibration = fixed_parameters(
        kinit=3.0, tev=1.0 if preeq_event else 1e6
    )
    if reinit_states is None:
        edata.reinitialize_fixed_parameter_initial_states = True
    else:
        edata.reinitialization_state_idxs_sim = [
            model.get_state_ids().index(id_) for id_ in reinit_states
        ]
    edata.set_measurements(np.full(model.nt() * model.ny, 1.0))
    edata.set_noise_scales(np.full(model.nt() * model.ny, 0.1))

    sllh = {}
    for sensi_meth, sensi_meth_preeq in (
        (SensitivityMethod.forward, SensitivityMethod.forward),
        (SensitivityMethod.adjoint, SensitivityMethod.forward),
        (SensitivityMethod.adjoint, SensitivityMethod.adjoint),
    ):
        solver.set_sensitivity_method(sensi_meth)
        solver.set_sensitivity_method_pre_equilibration(sensi_meth_preeq)
        rdata = run_simulation(model, solver, edata)
        assert rdata.status == AMICI_SUCCESS, (sensi_meth, sensi_meth_preeq)
        if preeq_event:
            assert rdata.preeq_status[1] == SteadyStateStatus.success
        sllh[sensi_meth, sensi_meth_preeq] = np.array(rdata.sllh)

    expected = sllh[SensitivityMethod.forward, SensitivityMethod.forward]
    # `scale` only enters via the initial value of `x1`
    scale_idx = model.get_free_parameter_ids().index("scale")
    others = np.delete(expected, scale_idx)
    assert np.all(np.abs(others) > 1e-2)
    if reinit_states != ["x2"]:
        assert abs(expected[scale_idx]) > 1e-2
    for key, actual in sllh.items():
        assert_allclose(actual, expected, rtol=1e-6, atol=1e-5, err_msg=key)

    function, derivative = run_simulation_to_function_and_derivative(
        amici_model=model,
        amici_solver=solver,
        amici_edata=edata,
        derivative_variables=["llh"],
    )
    point = np.array(model.get_free_parameters())
    check_gradient(function, point, derivative(point)["llh"]).assert_success(
        always_print=True
    )


@skip_on_valgrind
@pytest.mark.parametrize("variant", ["log10_plist", "event_at_t0"])
def test_preequilibration_reinit_adjoint_sensitivities_variants(
    reinit_model_module, variant
):
    """Adjoint pre-equilibration with reinitialization: parameter scaling
    with a parameter subset, and an event at the start of the main
    simulation after the reinitialization."""
    model = reinit_model_module.get_model()
    model.set_timepoints([0.0, 1.0, 2.0])
    if variant == "log10_plist":
        model.set_parameter_scale(ParameterScaling.log10)
        plist = [
            model.get_free_parameter_ids().index(id_)
            for id_ in ("scale", "kdecay")
        ]
        model.set_parameter_list(plist)
    solver = model.create_solver()
    solver.set_sensitivity_order(SensitivityOrder.first)
    solver.set_relative_tolerance(1e-10)
    solver.set_absolute_tolerance(1e-14)
    solver.set_relative_tolerance_quadratures(1e-10)
    solver.set_absolute_tolerance_quadratures(1e-14)

    def fixed_parameters(**kwargs):
        return [kwargs[id_] for id_ in model.get_fixed_parameter_ids()]

    edata = ExpData(model)
    edata.fixed_parameters = fixed_parameters(
        kinit=10.0, tev=0.0 if variant == "event_at_t0" else 1e6
    )
    edata.fixed_parameters_pre_equilibration = fixed_parameters(
        kinit=3.0, tev=1e6
    )
    edata.reinitialize_fixed_parameter_initial_states = True
    edata.set_measurements(np.full(model.nt() * model.ny, 1.0))
    edata.set_noise_scales(np.full(model.nt() * model.ny, 0.1))

    sllh = {}
    for sensi_meth, sensi_meth_preeq in (
        (SensitivityMethod.forward, SensitivityMethod.forward),
        (SensitivityMethod.adjoint, SensitivityMethod.adjoint),
    ):
        solver.set_sensitivity_method(sensi_meth)
        solver.set_sensitivity_method_pre_equilibration(sensi_meth_preeq)
        rdata = run_simulation(model, solver, edata)
        assert rdata.status == AMICI_SUCCESS, (sensi_meth, sensi_meth_preeq)
        if variant == "event_at_t0":
            # the event fired at t0: x2 = kinit + 1
            x2_idx = model.get_state_ids().index("x2")
            assert_allclose(rdata.x[0, x2_idx], 11.0, rtol=1e-8)
        sllh[sensi_meth] = np.array(rdata.sllh)

    expected = sllh[SensitivityMethod.forward]
    assert np.all(np.abs(expected) > 1e-2)
    assert_allclose(
        sllh[SensitivityMethod.adjoint], expected, rtol=1e-6, atol=1e-5
    )


@skip_on_valgrind
@pytest.mark.parametrize(
    "period", ["preeq_switch", "preeq_no_switch", "presimulation"]
)
def test_adjoint_switch_differs_from_main_simulation(tempdir, period):
    """Adjoint sensitivities for switches differing from the main simulation.

    `r = 2` at the end of pre-equilibration or presimulation, and `r = 1`
    during the main simulation. The solution is linear in `k`, so
    `sllh = -sum(x^2) / k` for measurements 0 with sigma 1.
    """
    model_name = "test_adjoint_switch_differs_from_main_simulation"
    antimony2amici(
        """
        fp = 0
        k = 0.5
        x = 0; x' = k - r * x
        r := piecewise(2, time >= fp, 1)
        """,
        model_name=model_name,
        output_dir=tempdir,
        fixed_parameters=["fp"],
        observation_model=[MeasurementChannel("obs_x", formula="x")],
    )
    model = amici.import_model_module(model_name, tempdir).get_model()
    assert model.get_free_parameter_ids() == ("k",)
    model.set_steady_state_computation_mode(
        SteadyStateComputationMode.integrationOnly
    )
    model.set_steady_state_sensitivity_mode(
        SteadyStateSensitivityMode.integrationOnly
    )
    t = np.array([0.0, 0.5])
    model.set_timepoints(t)
    edata = ExpData(model)
    edata.set_measurements([0.0] * len(t))
    edata.set_noise_scales([1.0] * len(t))
    edata.fixed_parameters = [1.0]

    k = 0.5
    match period:
        case "preeq_switch":
            # the switch flips at t=1 during pre-equilibration
            edata.fixed_parameters_pre_equilibration = [1.0]
            x0 = k / 2
        case "preeq_no_switch":
            # the switch is on throughout pre-equilibration, i.e., without
            # discontinuities, which uses the steady-state backward problem
            edata.fixed_parameters_pre_equilibration = [-1.0]
            x0 = k / 2
        case "presimulation":
            # presimulation over [-2, 0], the switch flips at t=-1
            edata.t_presim = 2.0
            edata.fixed_parameters_presimulation = [-1.0]
            x0 = k / 2 + (k * (1 - np.exp(-1)) - k / 2) * np.exp(-2)
    x = k + (x0 - k) * np.exp(-t)
    expected_sllh = -np.sum(x**2) / k

    for sens_method in (SensitivityMethod.forward, SensitivityMethod.adjoint):
        solver = model.create_solver()
        solver.set_sensitivity_order(SensitivityOrder.first)
        solver.set_sensitivity_method(sens_method)
        solver.set_sensitivity_method_pre_equilibration(sens_method)
        solver.set_relative_tolerance(1e-12)
        solver.set_absolute_tolerance(1e-14)
        rdata = run_simulation(model, solver, edata)
        assert rdata.status == AMICI_SUCCESS
        assert_allclose(rdata.x[:, 0], x, rtol=1e-8)
        assert_allclose(
            rdata.sllh[0],
            expected_sllh,
            rtol=1e-6,
            err_msg=sens_method.name,
        )


@skip_on_valgrind
def test_event_assignment_switch_preequilibration_adjoint(tempdir):
    """Adjoint sensitivities of an event assignment during pre-equilibration
    that depends on a switch re-initialized for the main simulation.

    During pre-equilibration, `E` fires at t=2 with `r = 2`. The main
    simulation re-initializes the switch at t0 and flips it at t=1 again, so
    the backward pass of the main simulation ends with the switch value
    before t=1.
    """
    model_name = "test_event_assignment_switch_preeq_adjoint"
    antimony2amici(
        """
        fp = 0
        k = 0.5
        x = 1; x' = -k * x
        v = 0; v' = 0
        r := piecewise(2, time >= 1, 1)
        E: at time >= 2, fromTrigger=false: v = v + r * x
        """,
        model_name=model_name,
        output_dir=tempdir,
        # only to request pre-equilibration
        fixed_parameters=["fp"],
        observation_model=[MeasurementChannel("obs_v", formula="v")],
    )
    model = amici.import_model_module(model_name, tempdir).get_model()
    assert model.get_free_parameter_ids() == ("k",)
    model.set_steady_state_computation_mode(
        SteadyStateComputationMode.integrationOnly
    )
    model.set_steady_state_sensitivity_mode(
        SteadyStateSensitivityMode.integrationOnly
    )
    model.set_timepoints([0.0, 1.5])
    edata = ExpData(model)
    edata.set_measurements([0.0, 0.0])
    edata.set_noise_scales([1.0, 1.0])
    edata.fixed_parameters = [0.0]
    edata.fixed_parameters_pre_equilibration = [0.0]

    # v = r * x(2) after pre-equilibration, sllh = -sum(v * dv/dk)
    k = 0.5
    v = 2 * np.exp(-2 * k)
    dv_dk = -4 * np.exp(-2 * k)
    expected_sllh = -2 * v * dv_dk

    for sens_method in (SensitivityMethod.forward, SensitivityMethod.adjoint):
        solver = model.create_solver()
        solver.set_sensitivity_order(SensitivityOrder.first)
        solver.set_sensitivity_method(sens_method)
        solver.set_sensitivity_method_pre_equilibration(sens_method)
        solver.set_relative_tolerance(1e-12)
        solver.set_absolute_tolerance(1e-14)
        rdata = run_simulation(model, solver, edata)
        assert rdata.status == AMICI_SUCCESS
        assert_allclose(rdata.y[:, 0], [v, v], rtol=1e-8)
        assert_allclose(
            rdata.sllh[0], expected_sllh, rtol=1e-6, err_msg=sens_method.name
        )


@skip_on_valgrind
@pytest.mark.parametrize(
    "period, differing",
    [
        ("presimulation", "k"),
        ("preeq_simulation", "k"),
        ("preeq_steady_state", "k"),
        ("presimulation", "c0"),
        ("preeq_steady_state", "c0"),
    ],
)
def test_adjoint_condition_differs_from_main_simulation(
    tempdir, period, differing
):
    """Adjoint sensitivities and outputs for fixed parameters differing
    between pre-equilibration or presimulation and the main simulation.

    `k` enters the static expression `a`, `c0` the initial value of the
    constant species `C`, i.e., a conservation law total. Conservation laws
    are not supported with `piecewise`, so the switch `r`, which introduces
    a discontinuity during pre-equilibration, is only included for `k`.
    """
    switch = "piecewise(2, time >= tsw, 1)" if differing == "k" else "1"
    model_name = (
        f"test_adjoint_condition_differs_from_main_simulation_{differing}"
    )
    antimony2amici(
        f"""
        species x = 1
        species C
        C = c0
        J1: x -> ; a * C * x
        J2: -> x; q * r
        a := k * p
        r := {switch}
        k = 1
        c0 = 1
        tsw = 10
        p = 0.5
        q = 0.3
        """,
        model_name=model_name,
        output_dir=tempdir,
        fixed_parameters=["k", "c0", "tsw"],
        observation_model=[
            MeasurementChannel("obs_x", formula="a * x"),
            MeasurementChannel("obs_c", formula="C"),
        ],
    )
    model = amici.import_model_module(model_name, tempdir).get_model()
    assert set(model.get_free_parameter_ids()) == {"p", "q"}
    # `C` is eliminated by its conservation law
    assert model.nx_solver == (2 if differing == "k" else 1)
    model.set_steady_state_computation_mode(
        SteadyStateComputationMode.integrationOnly
    )
    model.set_steady_state_sensitivity_mode(
        SteadyStateSensitivityMode.integrationOnly
    )
    t = np.array([0.5, 1.0, 2.0])
    model.set_timepoints(t)
    edata = ExpData(model)
    edata.set_measurements([0.0] * len(t) * model.ny)
    edata.set_noise_scales([1.0] * len(t) * model.ny)

    def fixed_parameters(**kwargs):
        values = {"k": 1.0, "c0": 1.0, "tsw": 10.0} | kwargs
        return [values[id_] for id_ in model.get_fixed_parameter_ids()]

    edata.fixed_parameters = fixed_parameters()
    # reinitialize `C`, which is only constant with conservation laws
    edata.reinitialize_fixed_parameter_initial_states = differing == "c0"
    other = {"k": 5.0} if differing == "k" else {"c0": 3.0}
    match period:
        case "presimulation":
            edata.t_presim = 2.0
            edata.fixed_parameters_presimulation = fixed_parameters(**other)
        case "preeq_simulation":
            # the switch flips at t=1 during pre-equilibration
            edata.fixed_parameters_pre_equilibration = fixed_parameters(
                **other, tsw=1.0
            )
        case "preeq_steady_state":
            # without discontinuities, this uses the steady-state backward
            # problem
            edata.fixed_parameters_pre_equilibration = fixed_parameters(
                **other, tsw=-1.0
            )

    rdatas = {}
    for sens_method in (SensitivityMethod.forward, SensitivityMethod.adjoint):
        solver = model.create_solver()
        solver.set_sensitivity_order(SensitivityOrder.first)
        solver.set_sensitivity_method(sens_method)
        solver.set_sensitivity_method_pre_equilibration(sens_method)
        solver.set_relative_tolerance(1e-12)
        solver.set_absolute_tolerance(1e-14)
        rdata = run_simulation(model, solver, edata)
        assert rdata.status == AMICI_SUCCESS
        rdatas[sens_method] = rdata

    fsa = rdatas[SensitivityMethod.forward]
    asa = rdatas[SensitivityMethod.adjoint]
    # `C` has the value of the main simulation
    assert_allclose(fsa.y[:, 1], 1.0, rtol=1e-12)
    assert np.all(np.abs(fsa.sllh) > 1e-2)
    assert_allclose(asa.y, fsa.y, rtol=1e-10)
    assert_allclose(asa.llh, fsa.llh, rtol=1e-10)
    assert_allclose(asa.sllh, fsa.sllh, rtol=1e-6)


@skip_on_valgrind
@pytest.mark.parametrize("dae", [False, True])
def test_presimulation_reinitialization_after_preequilibration(tempdir, dae):
    """State reinitialization at the start of presimulation after
    pre-equilibration.

    `x` is reinitialized to `p * k` at the start of presimulation, so
    `x = a + (2 p - a) exp(-(t + 1))` and `dx/dp = 2 exp(-(t + 1))`.
    For the ODE, the reinitialization triggers an event at the start of
    presimulation. For the DAE, the algebraic state `z` has to be made
    consistent with the reinitialized `x`.
    """
    if dae:
        ant_str = """
        a = 1; p = 1.5; k = 1
        x = p * k; x' = a - z / 2
        z = 0; 0 = z - 2 * x
        """
    else:
        ant_str = """
        a = 1; p = 1.5; k = 1
        x = p * k; x' = a - x
        y = 0; y' = 0
        E: at x > 2: y = y + 1
        """
    module_name = f"test_presim_reinit_{'dae' if dae else 'ode'}"
    antimony2amici(
        ant_str,
        model_name=module_name,
        output_dir=tempdir,
        fixed_parameters=["k"],
    )
    model_module = amici.import_model_module(
        module_name=module_name, module_path=tempdir
    )
    model = model_module.get_model()
    model.set_parameter_scale(ParameterScaling.none)
    t = np.array([0.0, 1.0])
    model.set_timepoints(t)
    edata = ExpData(model)
    edata.fixed_parameters_pre_equilibration = [1.0]
    edata.fixed_parameters_presimulation = [2.0]
    edata.fixed_parameters = [2.0]
    edata.t_presim = 1.0
    ix = model.get_state_ids().index("x")
    edata.reinitialization_state_idxs_presim = [ix]

    solver = model.create_solver()
    solver.set_relative_tolerance(1e-12)
    solver.set_absolute_tolerance(1e-12)
    solver.set_sensitivity_order(SensitivityOrder.first)
    solver.set_sensitivity_method(SensitivityMethod.forward)
    rdata = run_simulation(model, solver, edata)
    assert rdata.status == AMICI_SUCCESS

    if not dae:
        # the event was triggered at the start of presimulation
        iy = model.get_state_ids().index("y")
        assert_allclose(rdata.x[:, iy], 1.0)
    ip = model.get_free_parameter_ids().index("p")
    assert_allclose(rdata.x[:, ix], 1 + 2 * np.exp(-(t + 1)), rtol=1e-8)
    assert_allclose(rdata.sx[:, ip, ix], 2 * np.exp(-(t + 1)), rtol=1e-6)
    if dae:
        iz = model.get_state_ids().index("z")
        assert_allclose(rdata.x[:, iz], 2 * rdata.x[:, ix], rtol=1e-8)
        assert_allclose(
            rdata.sx[:, ip, iz], 2 * rdata.sx[:, ip, ix], rtol=1e-6
        )


@skip_on_valgrind
def test_presimulation_without_preequilibration_custom_initial_state(tempdir):
    """Without pre-equilibration, a custom initial state is kept at the start
    of presimulation, and states are reinitialized only at the start of the
    main simulation."""
    model_name = "test_presim_custom_x0"
    antimony2amici(
        """
        k = 1
        x1 = k; x1' = 0
        x2 = 0; x2' = x1 - x2
        """,
        model_name=model_name,
        output_dir=tempdir,
        fixed_parameters=["k"],
    )
    model = amici.import_model_module(model_name, tempdir).get_model()
    model.set_timepoints([0.0])
    model.set_initial_state([5.0, 0.0])
    solver = model.create_solver()
    edata = ExpData(model)
    edata.fixed_parameters_presimulation = [2.0]
    edata.t_presim = 1.0
    edata.fixed_parameters = [3.0]
    edata.reinitialize_fixed_parameter_initial_states = True
    rdata = run_simulation(model, solver, edata)
    assert rdata.status == AMICI_SUCCESS
    # x1 = 5 during presimulation, reinitialized to k = 3 at t = 0
    assert_allclose(rdata.x[0], [3.0, 5.0 * (1 - np.exp(-1.0))], rtol=1e-6)


@pytest.fixture(scope="module")
def dae_posteq_event_module(tmp_path_factory):
    """A DAE model with a steady state and an event at `t_event`."""
    model_name = "dae_posteq_event"
    outdir = tmp_path_factory.mktemp(model_name)
    antimony2amici(
        """
        const k = 0.5, p = 2, q = 0.3, t_event = 5
        x = 1; y = 0.3
        x' = p - k * x + y
        0 = y - q * x
        E1: at time > t_event: x = x + 1
        """,
        observation_model=[MeasurementChannel("obs_y", formula="y")],
        fixed_parameters=["t_event"],
        model_name=model_name,
        output_dir=outdir,
    )
    return amici.import_model_module(model_name, outdir)


@skip_on_valgrind
@pytest.mark.parametrize(
    "with_event", [False, True], ids=["no_event", "event"]
)
def test_dae_posteq_forward_sensitivities_newton(
    dae_posteq_event_module, with_event
):
    """Forward sensitivities for post-equilibration of a DAE by simulation,
    with steady-state sensitivities from the Newton solver.

    Forward sensitivities are switched off during post-equilibration. With
    `with_event`, an event occurs during post-equilibration. Re-using the
    solver must give the same results as a fresh one.
    """
    model = dae_posteq_event_module.get_model()
    t_event = 5.0 if with_event else 1e10
    model.set_fixed_parameter_by_id("t_event", t_event)
    model.set_timepoints([1.0, np.inf])
    model.set_steady_state_computation_mode(
        SteadyStateComputationMode.integrationOnly
    )
    model.set_steady_state_sensitivity_mode(
        SteadyStateSensitivityMode.newtonOnly
    )
    solver = model.create_solver()
    solver.set_relative_tolerance(1e-12)
    solver.set_absolute_tolerance(1e-14)
    edata = ExpData(model)
    edata.set_measurements([0.5, 1.0])
    edata.set_noise_scales([1.0, 1.0])

    solver.set_sensitivity_order(SensitivityOrder.first)
    solver.set_sensitivity_method(SensitivityMethod.forward)
    rdata = run_simulation(model, solver, edata)
    assert rdata.status == AMICI_SUCCESS
    # the event only occurs during post-equilibration, if at all
    assert (rdata.posteq_t > t_event) == with_event

    # y_ss = p * q / (k - q), independent of the event
    k, p, q = 0.5, 2.0, 0.3
    sy_ss = [-p * q / (k - q) ** 2, q / (k - q), p * k / (k - q) ** 2]
    assert_allclose(rdata.sy[1, :, 0], sy_ss, rtol=1e-6)

    # the sensitivities and their error control are re-enabled
    rdata2 = run_simulation(model, solver, edata)
    rdata_fresh = run_simulation(model, solver.clone(), edata)
    assert rdata2.status == AMICI_SUCCESS
    assert rdata_fresh.status == AMICI_SUCCESS
    assert_equal(rdata2.numsteps, rdata_fresh.numsteps)
    assert_allclose(rdata2.sy, rdata_fresh.sy, rtol=1e-10)
    assert_allclose(rdata2.sy, rdata.sy, rtol=1e-10)


def test_expdata_from_model_reinitialization_presim(tempdir):
    """`ExpData(model)` keeps the model's reinitialization at presimulation
    start.

    Also without presimulation fixed parameters, in which case
    presimulation uses the main simulation fixed parameters.
    """
    ant_str = """
    k = 1
    x1 = k; x1' = 0
    x2 = 0; x2' = x1 - x2
    """
    module_name = "test_expdata_reinit_presim"
    antimony2amici(
        ant_str,
        model_name=module_name,
        output_dir=tempdir,
        fixed_parameters=["k"],
    )
    model = amici.import_model_module(module_name, tempdir).get_model()
    model.set_timepoints([0.0])
    solver = model.create_solver()
    k_sim = 3.0

    def simulate(edata, k_presim):
        edata.fixed_parameters_pre_equilibration = [1.0]
        if k_presim is not None:
            edata.fixed_parameters_presimulation = [k_presim]
        edata.t_presim = 1.0
        edata.fixed_parameters = [k_sim]
        rdata = run_simulation(model, solver, edata)
        assert rdata.status == AMICI_SUCCESS
        return rdata.x[0]

    model.set_reinitialize_fixed_parameter_initial_states(True)
    reinit_idxs = list(model.get_reinitialization_state_idxs())
    assert reinit_idxs
    edata = ExpData(model)
    assert list(edata.reinitialization_state_idxs_presim) == reinit_idxs
    assert list(edata.reinitialization_state_idxs_sim) == reinit_idxs

    for k_presim in (2.0, None):
        model.set_reinitialize_fixed_parameter_initial_states(True)
        x_model_setting = simulate(ExpData(model), k_presim)

        model.set_reinitialize_fixed_parameter_initial_states(False)
        edata = ExpData(model)
        edata.reinitialize_fixed_parameter_initial_states = True
        x_edata_setting = simulate(edata, k_presim)

        # x1 is reinitialized at the start of presimulation and of the main
        # simulation; x2 relaxes from its pre-equilibrium value 1 towards x1
        # during presimulation.
        k_pre = k_sim if k_presim is None else k_presim
        assert_allclose(
            x_edata_setting,
            [k_sim, k_pre - (k_pre - 1.0) * np.exp(-1.0)],
            rtol=1e-6,
        )
        assert_allclose(x_model_setting, x_edata_setting, rtol=1e-6)
