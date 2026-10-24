"""Tests for SBML events, including piecewise expressions."""

from copy import deepcopy

import numpy as np
import pytest
from amici import (
    MeasurementChannel as MC,
)
from amici import (
    import_model_module,
)
from amici.adapters.fiddy import run_simulation_to_function_and_derivative
from amici.importers.antimony import antimony2amici
from amici.sim.sundials import (
    AMICI_SUCCESS,
    ExpData,
    InterpolationType,
    SensitivityMethod,
    SensitivityOrder,
    SteadyStateComputationMode,
    SteadyStateSensitivityMode,
    run_simulation,
)
from amici.sim.sundials.gradient_check import check_derivatives
from amici.testing import skip_on_valgrind
from amici.testing.models import create_amici_model, create_sbml_model
from fiddy import check_gradient
from numpy.testing import assert_allclose
from util import (
    check_trajectories_with_adjoint_sensitivities,
    check_trajectories_with_forward_sensitivities,
    check_trajectories_without_sensitivities,
)

pytestmark = pytest.mark.filterwarnings(
    # https://github.com/AMICI-dev/AMICI/issues/18
    "ignore:Adjoint sensitivity analysis for models with discontinuous "
    "right hand sides .*:UserWarning",
)


@pytest.fixture(
    params=[
        pytest.param("events_plus_heavisides", marks=skip_on_valgrind),
        pytest.param(
            "piecewise_plus_event_simple_case", marks=skip_on_valgrind
        ),
        pytest.param(
            "piecewise_plus_event_semi_complicated", marks=skip_on_valgrind
        ),
        pytest.param(
            "piecewise_plus_event_trigger_depends_on_state",
            marks=skip_on_valgrind,
        ),
        pytest.param("nested_events", marks=skip_on_valgrind),
        pytest.param("event_state_dep_ddeltax_dtpx", marks=skip_on_valgrind),
    ]
)
def model(request):
    """Returns the requested AMICI model and analytical expressions."""
    (
        initial_assignments,
        parameters,
        rate_rules,
        species,
        events,
        timepoints,
        x_expected,
        sx_expected,
    ) = get_model_definition(request.param)

    # SBML model
    sbml_document, sbml_model = create_sbml_model(
        initial_assignments=initial_assignments,
        parameters=parameters,
        rate_rules=rate_rules,
        species=species,
        events=events,
        # uncomment `to_file` to save SBML model to file for inspection
        # to_file=sbml_test_models / (model_name + '.sbml'),
    )

    # AMICI model
    amici_model = create_amici_model(
        sbml_model=sbml_model,
        model_name=request.param,
    )
    amici_model.set_timepoints(timepoints)

    return amici_model, parameters, timepoints, x_expected, sx_expected


def get_model_definition(model_name):
    if model_name == "piecewise_plus_event_simple_case":
        return model_definition_piecewise_plus_event_simple_case()
    if model_name == "piecewise_plus_event_semi_complicated":
        return model_definition_piecewise_plus_event_semi_complicated()
    if model_name == "piecewise_plus_event_trigger_depends_on_state":
        return model_definition_piecewise_plus_event_trigger_depends_on_state()
    if model_name == "events_plus_heavisides":
        return model_definition_events_plus_heavisides()
    if model_name == "nested_events":
        return model_definition_nested_events()
    if model_name == "event_state_dep_ddeltax_dtpx":
        return model_definition_event_state_dep_ddeltax_dtpx()

    raise NotImplementedError(
        f"Model with name {model_name} is not implemented."
    )


def model_definition_events_plus_heavisides():
    """Test model for state- and parameter-dependent Heavisides.

    ODEs
    ----
    d/dt x_1:
        - {            0,    t <  delta
        - { -alpha * x_1,    t >= delta
    d/dt x_2:
        - beta * x_1 - gamma * x_2
    d/dt x_3:
        - {     -eta * x_3,    t <  zeta
        - { -eta * x_3 + 1,    t >= zeta

    Events:
    -------
    event_1:
        trigger: k1 - x_3
        bolus: [[ -x_3 / 2],
                [        0],
                [        0]]
    event_2:
        trigger: t - zeta
        bolus: [[        0],
                [        0],
                [ zeta / 3]]
    """
    # Model components
    species = ["x_1", "x_2", "x_3"]
    initial_assignments = {
        "x_1": "k1",
        "x_2": "k2",
        "x_3": "k3",
    }
    rate_rules = {
        "x_1": "piecewise( -alpha * x_1, time >= delta, 0)",
        "x_2": "beta * x_1 - gamma * x_2",
        "x_3": "-eta * x_3 + piecewise( 1, time >= zeta, 0)",
    }
    parameters = {
        "k1": 2,
        "k2": 0.01,
        "k3": 5,
        "alpha": 2,
        # FIXME: adjoint sensitivities w.r.t. beta are slightly off
        "beta": 3,
        "gamma": 2,
        "delta": 3,
        # FIXME: adjoint sensitivities w.r.t. eta are slightly off
        #  changing eta to e.g. 2.5 "fixes" python/tests/test_events.py::test_models[events_plus_heavisides]
        "eta": 1,
        "zeta": 5,
    }
    events = {
        "event_1": {
            "trigger": "x_3 < k1",
            "target": "x_1",
            "assignment": "x_1 - x_3 / 2",
        },
        "event_2": {
            "trigger": "time >= zeta",
            "target": "x_3",
            "assignment": "x_3 + zeta / 3",
        },
    }
    timepoints = np.linspace(0, 8, 400)

    # Analytical solution
    def x_expected(t, k1, k2, k3, alpha, beta, gamma, delta, eta, zeta):
        # The system reads dx/dt = Ax + b
        # x0 = (k1, k2, k3)
        x0 = np.array([[k1], [k2], [k3]])

        # gather event time points
        event_1_time = (np.log(k3) - np.log(k1)) / eta  # k1 > x3
        event_2_time = delta
        event_3_time = zeta

        def get_early_x(t):
            # compute dynamics
            if t < event_1_time:
                # Define A
                A = np.array([[0, 0, 0], [beta, -gamma, 0], [0, 0, -eta]])
                tmp_x = expm(t * A)
                return np.matmul(tmp_x, x0)

            elif t <= event_2_time:
                # "simulate" until first event
                A = np.array([[0, 0, 0], [beta, -gamma, 0], [0, 0, -eta]])
                tmp_x = expm(event_1_time * A)
                x1 = np.matmul(tmp_x, x0)
                # apply bolus
                delta_x = np.array([[float(-x1[2, 0] / 2)], [0], [0]])
                x1 += delta_x
                # "simulate" on
                tmp_x = expm((t - event_1_time) * A)
                return np.matmul(tmp_x, x1)

        if t < event_2_time:
            x = get_early_x(t).flatten()
        elif t < event_3_time:
            x2 = get_early_x(event_2_time)

            A = np.array([[-alpha, 0, 0], [beta, -gamma, 0], [0, 0, -eta]])
            tmp_x = expm((t - event_2_time) * A)
            x = np.matmul(tmp_x, x2).flatten()
        else:
            x2 = get_early_x(event_2_time)

            A = np.array([[-alpha, 0, 0], [beta, -gamma, 0], [0, 0, -eta]])
            tmp_x = expm((event_3_time - event_2_time) * A)
            x3 = np.matmul(tmp_x, x2)
            # apply bolus
            x3 += np.array([[0], [0], [zeta / 3]])

            hom_x = np.matmul(expm((t - event_3_time) * A), x3)
            inhom_x = [
                [0],
                [0],
                [-np.exp(-eta * (t - event_3_time)) / eta + 1 / eta],
            ]

            x = (hom_x + inhom_x).flatten()

        return np.array(x)

    def sx_expected(t, parameters):
        """get sx, w.r.t. parameters, via finite differences"""
        sx = []
        eps = 1e-6

        for ip in parameters:
            perturbed_params = deepcopy(parameters)
            perturbed_params[ip] += eps
            sx_p = x_expected(t, **perturbed_params)
            perturbed_params[ip] -= 2 * eps
            sx_m = x_expected(t, **perturbed_params)
            sx.append((sx_p - sx_m) / (2 * eps))

        return np.array(sx)

    return (
        initial_assignments,
        parameters,
        rate_rules,
        species,
        events,
        timepoints,
        x_expected,
        sx_expected,
    )


def model_definition_nested_events():
    """Test model for state- and parameter-dependent heavisides.

    ODEs
    ----
    d/dt x_1:
        inflow_1 - decay_1 * x1
    d/dt x_2:
        - decay_2 * x_2

    Events:
    -------
    event_1:
        trigger: x_1 > inflow_1 / decay_2
        bolus: [[          0],
                [ -1 / time]]
    event_2:
        trigger: x_2 > 0.5
        bolus: [[ bolus],
                [ bolus]]
    """
    # Model components
    species = ["x_1", "x_2"]
    initial_assignments = {
        "x_1": "k1",
        "x_2": "k2",
    }
    rate_rules = {
        "x_1": "inflow_1 - decay_1 * x_1",
        "x_2": "- decay_2 * x_2",
    }
    parameters = {
        "k1": 0,
        "k2": 0,
        "inflow_1": 4,
        "decay_1": 2,
        "decay_2": 5,
        "bolus": 0,  # for bolus != 0, nested event sensitivities are off!
    }
    events = {
        "event_1": {
            "trigger": "x_1 > inflow_1 / decay_2",
            "target": "x_2",
            "assignment": "x_2 - 1 / time",
        },
        "event_2": {
            "trigger": "x_2 < - 0.5",
            "target": ["x_1", "x_2"],
            "assignment": ["x_1 + bolus", "x_2 + bolus"],
        },
    }
    timepoints = np.linspace(0, 1, 101)

    # Analytical solution
    def x_expected(t, k1, k2, inflow_1, decay_1, decay_2, bolus):
        # gather temporary variables
        # event_time = x_1 > inflow_1 / decay_2
        equil = inflow_1 / decay_1
        tmp1 = inflow_1 / decay_2 - inflow_1 / decay_1
        tmp2 = k1 - inflow_1 / decay_1
        event_time = (-1 / decay_1) * np.log(tmp1 / tmp2)

        def get_early_x(t):
            # compute dynamics before event
            x_1 = equil * (1 - np.exp(-decay_1 * t)) + k1 * np.exp(
                -decay_1 * t
            )
            x_2 = k2 * np.exp(-decay_2 * t)
            return np.array([[x_1], [x_2]])

        if t < event_time:
            x = get_early_x(t).flatten()
        else:
            # compute state after event
            x_tau = get_early_x(event_time)
            tau_x1 = x_tau[0] + bolus
            tau_x2 = x_tau[1] - 1 / event_time + bolus

            # compute dynamics after event
            inhom = np.exp(decay_1 * event_time) * tau_x1
            x_1 = equil * (
                1 - np.exp(decay_1 * (event_time - t))
            ) + inhom * np.exp(-decay_1 * t)
            x_2 = tau_x2 * np.exp(decay_2 * event_time) * np.exp(-decay_2 * t)

            x = np.array([[x_1], [x_2]])

        return x.flatten()

    def sx_expected(t, parameters):
        """get sx, w.r.t. parameters, via finite differences"""
        sx = []
        eps = 1e-6

        for ip in parameters:
            perturbed_params = deepcopy(parameters)
            perturbed_params[ip] += eps
            sx_p = x_expected(t, **perturbed_params)
            perturbed_params[ip] -= 2 * eps
            sx_m = x_expected(t, **perturbed_params)
            sx.append((sx_p - sx_m) / (2 * eps))

        return np.array(sx)

    return (
        initial_assignments,
        parameters,
        rate_rules,
        species,
        events,
        timepoints,
        x_expected,
        sx_expected,
    )


def model_definition_piecewise_plus_event_simple_case():
    """Test model for boolean operations in a piecewise condition.

    ODEs
    ----
    d/dt x_1:
        - { 1,    (alpha <= t and t < beta)
        - { 0,    otherwise
    """
    # Model components
    species = ["x_1"]
    initial_assignments = {"x_1": "x_1_0"}
    rate_rules = {"x_1": "piecewise(1, (alpha < time && time < beta), 0)"}
    parameters = {
        "alpha": 2,
        "beta": 3,
        "gamma": 4.5,
        "x_1_0": 1,
    }
    timepoints = np.linspace(0.0, 5.0, 100)  # np.array((0.0, 4.0,))
    events = {
        "event_1": {
            "trigger": "time > alpha",
            "target": "x_1",
            "assignment": "gamma",
        },
        "event_2": {
            "trigger": "time > beta",
            "target": "x_1",
            "assignment": "x_1 + 2.5",
        },
    }

    # Analytical solution
    def x_expected(t, x_1_0, alpha, beta, gamma):
        t_event_1 = alpha
        t_event_2 = beta

        if t < t_event_1:
            x = x_1_0
        elif t < t_event_2:
            x = gamma + t - t_event_1
        else:
            x = gamma + t_event_2 - t_event_1 + 2.5

        return np.array((x,))

    def sx_expected(t, parameters):
        """get sx, w.r.t. parameters, via finite differences"""
        sx = []
        eps = 1e-6

        for ip in parameters:
            perturbed_params = deepcopy(parameters)
            perturbed_params[ip] += eps
            sx_p = np.array(x_expected(t, **perturbed_params))
            perturbed_params[ip] -= 2 * eps
            sx_m = np.array(x_expected(t, **perturbed_params))
            sx.append((sx_p - sx_m) / (2 * eps))

        return np.array(sx)

    return (
        initial_assignments,
        parameters,
        rate_rules,
        species,
        events,
        timepoints,
        x_expected,
        sx_expected,
    )


def model_definition_event_state_dep_ddeltax_dtpx():
    """Test model with state-dependent partial derivatives of update functions wrt parameters, time, and states."""
    # Model components
    species = ["x_1"]
    initial_assignments = {"x_1": "x_1_0"}
    rate_rules = {"x_1": "1"}
    parameters = {
        "alpha": 1.5,
        "beta": 2.5,
        "gamma": 3.5,
        "delta": 5.5,
        "x_1_0": 1,
    }
    timepoints = np.linspace(0.0, 5.0, 100)
    events = {
        # state-dependent ddeltaxdt
        "event_1": {
            "trigger": "time > alpha",
            "target": "x_1",
            "assignment": "x_1 * time",
        },
        # state-dependent ddeltaxdp
        "event_2": {
            "trigger": "time > beta",
            "target": "x_1",
            "assignment": "x_1 * delta",
        },
        # state-dependent ddeltaxdx
        "event_3": {
            "trigger": "time > gamma",
            "target": "x_1",
            "assignment": "2 * x_1 * x_1",
        },
    }

    # Analytical solution
    def x_expected(t, x_1_0, alpha, beta, gamma, delta):
        if t < alpha:
            # before first event triggered
            x = x_1_0 + t
        elif t < beta:
            # after first event triggered
            x = (x_1_0 + alpha) * alpha + (t - alpha)
        elif t < gamma:
            # after second event triggered
            x = ((x_1_0 + alpha) * alpha + (beta - alpha)) * delta + (t - beta)
        else:
            # after third event triggered
            x = (
                ((x_1_0 + alpha) * alpha + (beta - alpha)) * delta
                + (gamma - beta)
            ) ** 2 * 2 + (t - gamma)

        return np.array((x,))

    def sx_expected(t, parameters):
        """get sx, w.r.t. parameters, via finite differences"""
        sx = []
        eps = 1e-6

        for ip in parameters:
            perturbed_params = deepcopy(parameters)
            perturbed_params[ip] += eps
            sx_p = np.array(x_expected(t, **perturbed_params))
            perturbed_params[ip] -= 2 * eps
            sx_m = np.array(x_expected(t, **perturbed_params))
            sx.append((sx_p - sx_m) / (2 * eps))

        return np.array(sx)

    return (
        initial_assignments,
        parameters,
        rate_rules,
        species,
        events,
        timepoints,
        x_expected,
        sx_expected,
    )


def model_definition_piecewise_plus_event_semi_complicated():
    """Test model for boolean operations in a piecewise condition, discrete
    events and a non-vanishing quadrature for the adjoint state.
    """
    # Model components
    species = ["x_1", "x_2"]
    initial_assignments = {"x_1": "x_1_0", "x_2": "x_2_0"}
    rate_rules = {
        "x_1": "piecewise(delta * x_1, (alpha < time && time < beta), - x_1)",
        "x_2": "- eta * x_2",
    }
    parameters = {
        "alpha": 2,
        "beta": 3,
        "gamma": 4.5,
        "x_1_0": 1,
        "x_2_0": 5,
        "delta": 2.5,
        "eta": 1.4,
    }
    timepoints = np.linspace(0.0, 5.0, 100)
    events = {
        "event_1": {
            "trigger": "time > alpha / 2",
            "target": "x_1",
            "assignment": "gamma",
        },
        "event_2": {
            "trigger": "time > beta",
            "target": "x_1",
            "assignment": "x_1 + x_2",
        },
    }

    # Analytical solution
    def x_expected(t, x_1_0, x_2_0, alpha, beta, gamma, delta, eta):
        t_event_1 = alpha / 2
        t_event_2 = beta
        heaviside_1 = alpha

        x_2 = x_2_0 * np.exp(-eta * t)

        if t < t_event_1:
            x_1 = x_1_0 * np.exp(-t)
        elif t < heaviside_1:
            x_1 = gamma * np.exp(-(t - t_event_1))
        elif t < t_event_2:
            x_1_heaviside_1 = gamma * np.exp(-(heaviside_1 - t_event_1))
            x_1 = x_1_heaviside_1 * np.exp(delta * (t - heaviside_1))
        else:
            x_1_heaviside_1 = gamma * np.exp(-(heaviside_1 - t_event_1))
            x_1_at_event_2 = x_1_heaviside_1 * np.exp(
                delta * (t_event_2 - heaviside_1)
            )
            x_2_at_event_2 = x_2_0 * np.exp(-eta * t_event_2)
            x1_after_event_2 = x_1_at_event_2 + x_2_at_event_2
            x_1 = x1_after_event_2 * np.exp(-(t - t_event_2))

        return np.array((x_1, x_2))

    def sx_expected(t, parameters):
        """get sx, w.r.t. parameters, via finite differences"""
        sx = []
        eps = 1e-6

        for ip in parameters:
            perturbed_params = deepcopy(parameters)
            perturbed_params[ip] += eps
            sx_p = np.array(x_expected(t, **perturbed_params))
            perturbed_params[ip] -= 2 * eps
            sx_m = np.array(x_expected(t, **perturbed_params))
            sx.append((sx_p - sx_m) / (2 * eps))

        return np.array(sx)

    return (
        initial_assignments,
        parameters,
        rate_rules,
        species,
        events,
        timepoints,
        x_expected,
        sx_expected,
    )


def model_definition_piecewise_plus_event_trigger_depends_on_state():
    """Test model for boolean operations in a piecewise condition.

    ODEs
    ----
    d/dt x_1:
        - { 1,    (alpha <= t and t < beta)
        - { 0,    otherwise
    """
    # Model components
    species = ["x_1", "x_2"]
    initial_assignments = {"x_1": "x_1_0", "x_2": "x_2_0"}
    rate_rules = {
        "x_1": "piecewise(1, (alpha < time && time < beta), 0)",
        "x_2": "- x_2",
    }
    parameters = {
        "alpha": 2,
        "beta": 3,
        "gamma": 4.5,
        "x_1_0": 1,
        "x_2_0": 5,
    }
    timepoints = np.linspace(0.0, 5.0, 100)
    events = {
        "event_1": {
            "trigger": "x_1 > 1.4",
            "target": "x_1",
            "assignment": "x_1 + gamma",
        },
        "event_2": {
            "trigger": "time > beta",
            "target": "x_1",
            "assignment": "x_1 + x_2",
        },
    }

    # Analytical solution
    def x_expected(t, x_1_0, x_2_0, alpha, beta, gamma):
        heaviside_1 = alpha
        t_event_1 = alpha + 1.4 - x_1_0
        t_event_2 = beta
        # This should hold in order that the analytical solution is correct
        assert heaviside_1 < t_event_1

        # x_2 never gets perturbed
        x_2 = x_2_0 * np.exp(-t)

        if t < heaviside_1:
            x_1 = x_1_0
        elif t < t_event_1:
            x_1 = (t - heaviside_1) + x_1_0
        elif t < t_event_2:
            x_1 = gamma + (t - heaviside_1) + x_1_0
        else:
            x_2_at_event_2 = x_2_0 * np.exp(-t_event_2)
            x_1_at_event_2 = gamma + (t_event_2 - heaviside_1) + x_1_0
            x_1 = x_1_at_event_2 + x_2_at_event_2

        return np.array((x_1, x_2))

    def sx_expected(t, parameters):
        """get sx, w.r.t. parameters, via finite differences"""
        sx = []
        eps = 1e-6

        for ip in parameters:
            perturbed_params = deepcopy(parameters)
            perturbed_params[ip] += eps
            sx_p = np.array(x_expected(t, **perturbed_params))
            perturbed_params[ip] -= 2 * eps
            sx_m = np.array(x_expected(t, **perturbed_params))
            sx.append((sx_p - sx_m) / (2 * eps))

        return np.array(sx)

    return (
        initial_assignments,
        parameters,
        rate_rules,
        species,
        events,
        timepoints,
        x_expected,
        sx_expected,
    )


def test_models(model):
    amici_model, parameters, timepoints, x_expected, sx_expected = model

    result_expected_x = np.array(
        [x_expected(t, **parameters) for t in timepoints]
    )
    result_expected_sx = np.array(
        [sx_expected(t, parameters) for t in timepoints]
    )

    # assert correctness of trajectories
    check_trajectories_without_sensitivities(amici_model, result_expected_x)
    check_trajectories_with_forward_sensitivities(
        amici_model, result_expected_x, result_expected_sx
    )

    # For `nested_events`, `llh` is discontinuous in `bolus` at 0: for
    # `bolus < 0`, `event_1` triggers again right after its first occurrence.
    skip_fd = ("bolus",) if amici_model.get_name() == "nested_events" else ()
    check_trajectories_with_adjoint_sensitivities(amici_model, skip_fd)


def expm(x):
    """``expm`` wrapper

    Uses ``expm`` from ``mpmath``. *Something* changed in scipy's ``expm`` in
    version 1.9.0 breaking these tests"""
    from mpmath import expm

    return np.array(expm(x).tolist()).astype(float)


def test_handling_of_fixed_time_point_event_triggers(tempdir):
    """Test handling of events without solver-tracked root functions."""
    ant_model = """
    model test_events_time_based
        one = 1
        two = 2
        three = 3
        four = 4
        five = 5
        event_target = 0
        bolus = 1
        at (time > one): event_target = one
        at (time > two): event_target = event_target + bolus
        at (time > three): event_target = three
        at (time > four): event_target = four
        at (time > five): event_target = five
    end
    """
    module_name = "test_events_time_based"
    antimony2amici(
        ant_model,
        # test with constant parameters and non-constant parameters!
        fixed_parameters=["four"],
        model_name=module_name,
        output_dir=tempdir,
    )
    model_module = import_model_module(
        module_name=module_name, module_path=tempdir
    )
    amici_model = model_module.get_model()
    assert amici_model.ne == 5
    assert amici_model.ne_solver == 0

    amici_model.set_timepoints(np.linspace(0, 10, 20))
    amici_solver = amici_model.create_solver()
    rdata = run_simulation(amici_model, amici_solver)
    assert rdata.status == AMICI_SUCCESS
    assert (rdata.x[rdata.ts < 1] == 0).all()
    assert (rdata.x[(rdata.ts >= 1) & (rdata.ts < 2)] == 1).all()
    assert (rdata.x[(rdata.ts >= 2) & (rdata.ts < 3)] == 2).all()
    assert (rdata.x[(rdata.ts >= 3) & (rdata.ts < 4)] == 3).all()
    assert (rdata.x[(rdata.ts >= 4) & (rdata.ts < 5)] == 4).all()
    assert (rdata.x[(rdata.ts >= 5)] == 5).all()
    assert rdata.x[-1, :] == 5

    edata = ExpData(rdata, 1, 0, 0)

    for sens_meth in (
        SensitivityMethod.forward,
        SensitivityMethod.adjoint,
    ):
        amici_solver.set_sensitivity_method(sens_meth)
        amici_solver.set_sensitivity_order(SensitivityOrder.first)
        check_derivatives(amici_model, solver=amici_solver, edata=edata)


@skip_on_valgrind
def test_multiple_event_assignment_with_compartment(tempdir):
    """see https://github.com/AMICI-dev/AMICI/issues/2426"""
    ant_model = """
    model test_events_multiple_assignments
        compartment event_target = 1
        event_target' = 0
        species species_in_event_target in event_target = 1
        unrelated = 2

        # use different order of event assignments for the two events
        at (time > 5): unrelated = 4, event_target = 10
        at (time > 10): event_target = 1, unrelated = 2
    end
    """
    module_name = "test_multiple_event_assignment_with_compartment"
    antimony2amici(
        ant_model,
        model_name=module_name,
        output_dir=tempdir,
    )
    model_module = import_model_module(
        module_name=module_name, module_path=tempdir
    )
    amici_model = model_module.get_model()
    assert amici_model.ne == 2
    assert amici_model.ne_solver == 0
    assert amici_model.nx_rdata == 3
    amici_model.set_timepoints(np.linspace(0, 15, 16))
    amici_solver = amici_model.create_solver()
    rdata = run_simulation(amici_model, amici_solver)
    assert rdata.status == AMICI_SUCCESS
    idx_event_target = amici_model.get_state_ids().index("event_target")
    idx_unrelated = amici_model.get_state_ids().index("unrelated")
    idx_species_in_event_target = amici_model.get_state_ids().index(
        "species_in_event_target"
    )

    assert_allclose(
        rdata.x[(rdata.ts < 5) & (rdata.ts > 10), idx_event_target],
        1,
        rtol=0,
        atol=1e-15,
    )
    assert_allclose(
        rdata.x[(5 < rdata.ts) & (rdata.ts < 10), idx_event_target],
        10,
        rtol=0,
        atol=1e-15,
    )
    assert_allclose(
        rdata.x[(rdata.ts < 5) & (rdata.ts > 10), idx_unrelated],
        2,
        rtol=0,
        atol=1e-15,
    )
    assert_allclose(
        rdata.x[(5 < rdata.ts) & (rdata.ts < 10), idx_unrelated],
        4,
        rtol=0,
        atol=1e-15,
    )
    assert_allclose(
        rdata.x[(rdata.ts < 5) & (rdata.ts > 10), idx_species_in_event_target],
        1,
        rtol=0,
        atol=1e-15,
    )
    assert_allclose(
        rdata.x[(5 < rdata.ts) & (rdata.ts < 10), idx_species_in_event_target],
        0.1,
        rtol=0,
        atol=1e-15,
    )


@skip_on_valgrind
def test_event_assignment_with_expressions(tempdir):
    """Event assignments using (chains of) assignment rule variables."""
    model_name = "test_event_assignment_with_expressions"
    antimony2amici(
        """
        k = 0.5; p = 3
        x = 1; v = 0
        y := p * x
        z := y + 1
        x' = -k * x
        v' = 0
        E1: at time > 1: v = v + z
        """,
        model_name=model_name,
        output_dir=tempdir,
        observation_model=[
            MC("obs_x", formula="x"),
            MC("obs_v", formula="v"),
        ],
    )
    model = import_model_module(model_name, tempdir).get_model()
    model.set_timepoints(np.linspace(0.5, 4.5, 9))

    solver = model.create_solver()
    solver.set_relative_tolerance(1e-12)
    solver.set_absolute_tolerance(1e-12)
    rdata = run_simulation(model, solver)
    assert rdata.status == AMICI_SUCCESS
    k, p = model.get_free_parameters()
    v_expected = np.where(rdata.ts >= 1, p * np.exp(-k) + 1, 0)
    assert_allclose(rdata.by_id("v"), v_expected, rtol=1e-8)

    edata = ExpData(rdata, 1.0, 0.0)
    # fixed offset, so that the gradient is non-zero
    edata.set_measurements((np.asarray(rdata.y) + 0.5).flatten())
    function, _ = run_simulation_to_function_and_derivative(
        amici_model=model,
        amici_solver=solver,
        amici_edata=edata,
        derivative_variables=["llh"],
    )
    point = np.asarray(model.get_free_parameters())
    solver.set_sensitivity_order(SensitivityOrder.first)
    for sensi_meth in (
        SensitivityMethod.forward,
        SensitivityMethod.adjoint,
    ):
        solver.set_sensitivity_method(sensi_meth)
        rdata = run_simulation(model, solver, edata)
        assert rdata.status == AMICI_SUCCESS, sensi_meth
        assert np.all(np.abs(rdata.sllh) > 0.1), sensi_meth
        check_gradient(function, point, rdata.sllh).assert_success(
            always_print=True
        )


@skip_on_valgrind
def test_event_outputs_llh(tempdir):
    """Test the likelihood for multiple event outputs.

    Two outputs of an event that fires, and one output of an event that
    does not fire until the last output time, which contributes only via
    its root function value at that time.
    """
    model_name = "test_event_outputs_llh"
    antimony2amici(
        """
        k = 0.5; q = 0.1
        x = 1
        x' = -k * x + q
        E1: at (time > 1.3): x = x + 1
        E2: at (x > 2): x = 0
        """,
        model_name=model_name,
        output_dir=tempdir,
        observation_model=[
            MC("obs_x", formula="x"),
            MC("z1", formula="x", event_id="E1"),
            MC("z2", formula="2 * x", event_id="E1"),
            MC("z3", formula="x", event_id="E2"),
        ],
    )
    model = import_model_module(model_name, tempdir).get_model()
    model.set_timepoints(np.linspace(0.5, 4.5, 9))
    edata = ExpData(model)
    mz = np.full(model.n_max_event() * model.nztrue, np.nan)
    mz[:3] = [0.3, 0.6, 0.2]
    edata.set_event_measurements(mz)
    edata.set_event_noise_scales(np.ones_like(mz))
    rdata = run_simulation(model, model.create_solver(), edata)
    assert rdata.status == AMICI_SUCCESS

    z, rz = rdata.z[0], rdata.rz[0]
    assert np.all(np.abs(z - mz[:3]) > 0.1)
    assert np.abs(rz[2]) > 0.1
    res = np.append(z[:2] - mz[:2], rz[2])
    expected_llh = -np.sum(0.5 * np.log(2 * np.pi) + 0.5 * res**2)
    assert np.isclose(rdata.llh, expected_llh, rtol=1e-12)


@pytest.mark.filterwarnings(
    "ignore:Adjoint sensitivity analysis for models with discontinuous "
)
@skip_on_valgrind
def test_event_priorities(tempdir):
    """Test SBML event priorities."""

    model_name = "test_event_priorities"
    antimony2amici(
        r"""
        target1 = 1
        target2 = 2
        target3 = 3
        target3_rate = 0
        target3' = target3_rate
        trigger_time = 1

        # test time- and state-dependent triggers
        some_time = time
        some_time' = 1

        two = 2

        # three events with different priorities, where priorities
        #  don't match alphabetical order of IDs or anything the like
        E_two: \
            at some_time >= trigger_time, priority=22, fromTrigger=false:
            target2 = two * target1 + target2 - two;
        E_one: \
            at some_time >= trigger_time, priority=111, fromTrigger=false:
            target1 = 10 + time;
        E_three: \
            at some_time >= trigger_time, priority=3, fromTrigger=false:
            target3 = target1 + target2, target3_rate = 1;
        """,
        model_name=model_name,
        output_dir=tempdir,
    )

    model_module = import_model_module(model_name, tempdir)

    model = model_module.get_model()

    # check just after the trigger time,
    # the event does not fire at *exactly* 1
    model.set_timepoints([0, 1 + 1e-6, 2])

    solver = model.create_solver()
    solver.set_absolute_tolerance(1e-16)
    solver.set_relative_tolerance(1e-14)
    solver.set_sensitivity_order(SensitivityOrder.first)
    solver.set_sensitivity_method(SensitivityMethod.forward)

    rdata = run_simulation(model, solver)

    assert np.all(rdata.by_id("target1") == [1, 11, 11])
    assert np.all(rdata.by_id("target2") == [2, 22, 22])
    assert_allclose(rdata.by_id("target3"), [3, 33 + 1e-6, 33 + 1])

    # generate synthetic measurements
    edata = ExpData(rdata, 1, 0)

    # check forward sensitivities against finite differences
    check_derivatives(
        model,
        solver=solver,
        edata=edata,
        atol=1e-6,
        rtol=1e-6,
        # smaller than the offset from the trigger time
        epsilon=1e-8,
    )

    # TODO: test ASA after https://github.com/AMICI-dev/AMICI/pull/1539
    # FIXME: sensitivities w.r.t. the bolus and trigger parameter are totally off
    # solver.setSensitivityMethod(SensitivityMethod.adjoint)
    # edata.plist = []
    # model.requireSensitivitiesForAllParameters()
    # check_derivatives(
    #     model,
    #     solver=solver,
    #     edata=edata,
    #     atol=1e-6,
    #     rtol=1e-6,
    #     # smaller than the offset from the trigger time
    #     epsilon=1e-8,
    # )


@skip_on_valgrind
def test_switches_triggered_by_event_assignments(tempdir):
    """Sensitivities of switches that are only triggered by events.

    The `piecewise` conditions only depend on event-assigned variables,
    i.e., they don't change between events. They are only triggered by `E1`
    and share its event time sensitivity. The same holds for `E2`, whose
    trigger is either constant between events, depends on a state but is
    crossed only due to `E1`, or depends on time only via a switch.
    `v` doesn't depend on `k`, so `E2`'s own event time sensitivity
    w.r.t. `k` would be 0.
    """
    triggers = {
        "constant": "s1 >= 1",
        "state": "s1 + v >= 1.5",
        "switched_time": "w >= 0.5",
    }
    for trigger_id, trigger in triggers.items():
        model_name = (
            f"test_switches_triggered_by_event_assignments_{trigger_id}"
        )
        antimony2amici(
            f"""
            k = 1.5
            x = 0; x' = k
            v = 0; v' = 1
            s1 = 0; s2 = 0
            acc = 0
            w := piecewise(time, s1 >= 1, 0)
            z_on = 0; z_on' = piecewise(1, s2 >= 1, 0)
            z_off = 0; z_off' = piecewise(1, s2 < 1, 0)
            E1: at x >= 1: s1 = 1
            E2: at {trigger}: s2 = 1, acc = acc + time
            """,
            model_name=model_name,
            output_dir=tempdir,
            observation_model=[
                MC("obs_z_on", formula="z_on"),
                MC("obs_z_off", formula="z_off"),
                MC("obs_acc", formula="acc"),
            ],
        )
        model = import_model_module(model_name, tempdir).get_model()
        t_end = 3.0
        model.set_timepoints([t_end])
        ik = model.get_free_parameter_ids().index("k")
        k = model.get_free_parameters()[ik]

        # Both events fire at t1 = 1 / k
        t1 = 1 / k
        y = np.array([t_end - t1, t1, t1])
        dy_dk = np.array([1, -1, -1]) / k**2
        measurements = np.array([0.0, 1.0, 2.0])
        expected_sllh = -np.sum((y - measurements) * dy_dk)
        assert abs(expected_sllh) > 1

        edata = ExpData(model)
        edata.set_measurements(measurements)
        edata.set_noise_scales(np.ones_like(measurements))
        for sens_method in (
            SensitivityMethod.forward,
            SensitivityMethod.adjoint,
        ):
            solver = model.create_solver()
            solver.set_sensitivity_order(SensitivityOrder.first)
            solver.set_sensitivity_method(sens_method)
            solver.set_relative_tolerance(1e-12)
            solver.set_absolute_tolerance(1e-12)
            rdata = run_simulation(model, solver, edata=edata)
            err_msg = f"{trigger_id}, {sens_method}"
            assert rdata.status == AMICI_SUCCESS, err_msg
            assert_allclose(rdata.y[-1], y, rtol=1e-8, err_msg=err_msg)
            if sens_method == SensitivityMethod.forward:
                assert_allclose(
                    rdata.sy[-1, ik], dy_dk, rtol=1e-6, err_msg=err_msg
                )
            assert_allclose(
                rdata.sllh[ik], expected_sllh, rtol=1e-6, err_msg=err_msg
            )

    # Simultaneous events that each flip a switch
    model_name = "test_switches_simultaneous"
    antimony2amici(
        """
        k = 1.5
        x = 0; x' = k
        s = 0; u = 0
        z = 0; z' = piecewise(1, s >= 1, 0) + piecewise(3, u >= 1, 0)
        Ea: at x >= 1, priority=2: s = 1
        Eb: at x >= 1, priority=1: u = 1
        """,
        model_name=model_name,
        output_dir=tempdir,
        observation_model=[MC("obs_z", formula="z")],
    )
    model = import_model_module(model_name, tempdir).get_model()
    model.set_timepoints([t_end])
    y = np.array([4 * (t_end - t1)])
    dy_dk = np.array([4 / k**2])
    measurements = np.array([0.0])
    expected_sllh = -np.sum((y - measurements) * dy_dk)
    edata = ExpData(model)
    edata.set_measurements(measurements)
    edata.set_noise_scales(np.ones_like(measurements))
    for sens_method in (SensitivityMethod.forward, SensitivityMethod.adjoint):
        solver = model.create_solver()
        solver.set_sensitivity_order(SensitivityOrder.first)
        solver.set_sensitivity_method(sens_method)
        solver.set_relative_tolerance(1e-12)
        solver.set_absolute_tolerance(1e-12)
        rdata = run_simulation(model, solver, edata=edata)
        assert rdata.status == AMICI_SUCCESS
        assert_allclose(rdata.y[-1], y, rtol=1e-8)
        if sens_method == SensitivityMethod.forward:
            assert_allclose(rdata.sy[-1, 0], dy_dk, rtol=1e-6)
        assert_allclose(
            rdata.sllh[0], expected_sllh, rtol=1e-6, err_msg=str(sens_method)
        )

    # Only fixed-time events, i.e., the event time sensitivities are zero
    model_name = "test_switches_triggered_by_event_assignments_fixed"
    antimony2amici(
        """
        k = 1.5
        s = 0
        z = 0; z' = k * piecewise(1, s >= 1, 0)
        E1: at time >= 1: s = 1
        """,
        model_name=model_name,
        output_dir=tempdir,
        observation_model=[MC("obs_z", formula="z")],
    )
    model = import_model_module(model_name, tempdir).get_model()
    model.set_timepoints([t_end])
    y = np.array([k * (t_end - 1)])
    dy_dk = np.array([t_end - 1])
    measurements = np.array([0.0])
    expected_sllh = -np.sum((y - measurements) * dy_dk)
    edata = ExpData(model)
    edata.set_measurements(measurements)
    edata.set_noise_scales(np.ones_like(measurements))
    for sens_method in (SensitivityMethod.forward, SensitivityMethod.adjoint):
        solver = model.create_solver()
        solver.set_sensitivity_order(SensitivityOrder.first)
        solver.set_sensitivity_method(sens_method)
        rdata = run_simulation(model, solver, edata=edata)
        assert rdata.status == AMICI_SUCCESS
        assert_allclose(rdata.y[-1], y, rtol=1e-8)
        if sens_method == SensitivityMethod.forward:
            assert_allclose(rdata.sy[-1, 0], dy_dk, rtol=1e-6)
        assert_allclose(
            rdata.sllh[0], expected_sllh, rtol=1e-6, err_msg=str(sens_method)
        )


@skip_on_valgrind
@pytest.mark.parametrize("from_trigger", [True, False])
def test_event_assignment_switch_flipped_before_execution(
    tempdir, from_trigger
):
    """Sensitivities of event assignments depending on a switch that is
    flipped by a simultaneous event executed before.

    For `fromTrigger=true`, the assignment of `E2` uses the switch value from
    the trigger time, otherwise the one after `E1`.
    """
    ft = str(from_trigger).lower()
    model_name = f"test_switch_flipped_before_execution_{ft}"
    antimony2amici(
        f"""
        k = 1.5
        x = 0; x' = k
        s = 0; v = 0
        r := piecewise(2, s >= 1, 1)
        E1: at x >= 1, priority=2: s = 1
        E2: at x >= 1, priority=1, fromTrigger={ft}: v = v + r * time
        """,
        model_name=model_name,
        output_dir=tempdir,
        observation_model=[MC("obs_v", formula="v")],
    )
    model = import_model_module(model_name, tempdir).get_model()
    model.set_timepoints([2.0])
    k = model.get_free_parameters()[0]

    # both events fire at t1 = 1 / k
    r = 1 if from_trigger else 2
    v = r / k
    dv_dk = -r / k**2
    expected_sllh = -v * dv_dk

    edata = ExpData(model)
    edata.set_measurements([0.0])
    edata.set_noise_scales([1.0])
    for sens_method in (SensitivityMethod.forward, SensitivityMethod.adjoint):
        solver = model.create_solver()
        solver.set_sensitivity_order(SensitivityOrder.first)
        solver.set_sensitivity_method(sens_method)
        solver.set_relative_tolerance(1e-12)
        solver.set_absolute_tolerance(1e-12)
        rdata = run_simulation(model, solver, edata=edata)
        assert rdata.status == AMICI_SUCCESS
        assert_allclose(rdata.y[-1], [v], rtol=1e-8)
        if sens_method == SensitivityMethod.forward:
            assert_allclose(rdata.sy[-1, 0], [dv_dk], rtol=1e-6)
        assert_allclose(
            rdata.sllh[0], expected_sllh, rtol=1e-6, err_msg=str(sens_method)
        )


@skip_on_valgrind
def test_event_output_parameter_derivatives(tempdir):
    """Test forward sensitivities of event outputs that depend explicitly on
    the parameters, and of the root function value of an event that does
    not trigger."""
    from amici.sim.sundials import ParameterScaling

    model_name = "test_event_output_parameter_derivatives"
    model = antimony2amici(
        """
        k = 0.5; p = 0.1
        x = 1; y = 0
        x' = -k * x
        y' = 0
        E1: at time > 1.3: y = y + 1
        E2: at x < p: x = x + 1
        """,
        model_name=model_name,
        output_dir=tempdir,
        observation_model=[
            MC("obs_x", formula="x"),
            MC("z1", formula="p * x", event_id="E1"),
            MC("z2", formula="time", event_id="E2"),
        ],
    )
    model.set_parameter_scale(ParameterScaling.none)
    t = np.linspace(0.5, 2, 4)
    model.set_timepoints(t)
    assert model.get_free_parameter_ids() == ("k", "p")

    edata = ExpData(model)
    edata.set_measurements(np.full(len(t), np.nan))
    edata.set_noise_scales(np.ones(len(t)))
    mz = np.full(model.n_max_event() * model.nztrue, np.nan)
    mz[:2] = [1.0, 3.0]
    edata.set_event_measurements(mz)
    edata.set_event_noise_scales(np.ones_like(mz))

    solver = model.create_solver()
    solver.set_relative_tolerance(1e-12)
    solver.set_absolute_tolerance(1e-12)
    solver.set_sensitivity_order(SensitivityOrder.first)
    solver.set_sensitivity_method(SensitivityMethod.forward)
    rdata = run_simulation(model, solver, edata)
    assert rdata.status == AMICI_SUCCESS

    # E1 triggers at t1, E2 does not trigger until the last output time T
    k, p, t1, T = 0.5, 0.1, 1.3, t[-1]
    x1, xT = np.exp(-k * t1), np.exp(-k * T)
    z1 = p * x1
    dz1 = np.array([-p * t1 * x1, x1])
    rz2 = p - xT
    drz2 = np.array([T * xT, 1.0])
    assert_allclose(rdata.z[0], [z1, T], rtol=1e-8)
    assert_allclose(rdata.sz[0][:, 0], dz1, rtol=1e-6)
    assert_allclose(rdata.rz[0][1], rz2, rtol=1e-8)
    assert_allclose(rdata.srz[0][:, 1], drz2, rtol=1e-6)
    # Jz(z2) does not depend on the parameters (z2 = T)
    sllh_expected = -(z1 - mz[0]) * dz1 - rz2 * drz2
    assert np.all(np.abs(sllh_expected) > 0.1)
    assert_allclose(rdata.sllh, sllh_expected, rtol=1e-6)


@skip_on_valgrind
def test_random_event_ordering(tempdir):
    """For simultaneously executed events, the order of execution
    must be random."""

    model_name = "test_event_prio_rnd"
    antimony2amici(
        r"""
        target_rnd = 0
        target_first = 0
        target_last = 0
        # test time- and state-dependent triggers
        some_time = time
        some_time' = 1
        trigger_time = 1

        # {E1, E2, E3} must be executed in random order after E_first,
        # but before E_last
        E1: at some_time >= trigger_time, priority=1, fromTrigger=false:
            target_rnd = 1;
        E2: at some_time >= trigger_time, priority=1, fromTrigger=false:
            target_rnd = 2;
        E3: at some_time >= trigger_time, priority=1, fromTrigger=false:
            target_rnd = 3;
        E_first: \
            at some_time >= trigger_time, priority=10, fromTrigger=false:
            target_first = target_rnd + 2;
        E_last: \
            at some_time >= trigger_time, priority=-1, fromTrigger=false:
            target_last = target_rnd >= 1;
        """,
        model_name=model_name,
        output_dir=tempdir,
    )

    model_module = import_model_module(model_name, tempdir)

    model = model_module.get_model()
    model.set_timepoints([0, 2, 3])
    solver = model.create_solver()

    # the outcomes of the random assignment
    outcomes = []
    N = 1000
    for i in range(N):
        rdata = run_simulation(model, solver)
        assert np.all(rdata.by_id("target_first") == [0, 2, 2])
        assert np.all(rdata.by_id("target_last") == [0, 1, 1])
        traj = rdata.by_id("target_rnd")
        assert traj[0] == 0
        assert traj[1] == traj[2]
        # collect the random outputs
        outcomes.append(traj[2])

    assert set(outcomes) == {1, 2, 3}

    # check that the outcomes are about equally distributed
    # between 1, 2, and 3
    assert np.allclose(
        [outcomes.count(1), outcomes.count(2), outcomes.count(3)],
        [N / 3, N / 3, N / 3],
        rtol=0.25,
    )


@skip_on_valgrind
def test_event_uses_values_from_trigger_time(tempdir):
    """For simultaneously executed events, check that values from trigger
    times are used to compute the state update."""

    model_name = "test_event_vals_trig_time"
    antimony2amici(
        r"""
        some_time = time
        some_time' = 1
        trigger_time = 0.5
        target1 = 2
        target2 = 0
        one = 1
        three = 3

        E1: at some_time >= trigger_time, priority=10, fromTrigger=false:
            target1 = 10,
            target2 = target1 + one;

        E2: at some_time >= trigger_time, priority=-10, fromTrigger=true:
            # this must reset `target1` to its initial value!!
            target1 = target1,
            target2 = target2 + three
        """,
        model_name=model_name,
        output_dir=tempdir,
    )

    model_module = import_model_module(model_name, tempdir)

    model = model_module.get_model()
    model.set_timepoints([0, 1.1, 2])
    solver = model.create_solver()
    solver.set_sensitivity_order(SensitivityOrder.first)
    solver.set_sensitivity_method(SensitivityMethod.forward)

    rdata = run_simulation(model, solver)
    assert np.all(rdata.by_id("target1") == [2, 2, 2])
    assert np.all(rdata.by_id("target2") == [0, 3, 3])

    # generate synthetic measurements
    edata = ExpData(rdata, 1, 0)

    # check sensitivities against finite differences

    for sens_method in (
        SensitivityMethod.forward,
        SensitivityMethod.adjoint,
    ):
        if sens_method == SensitivityMethod.forward:
            # FIXME: forward sensitivities w.r.t. the bolus parameter
            #  of the first event are wrong
            model.set_parameter_list(
                [
                    ip
                    for ip, par in enumerate(model.get_free_parameter_ids())
                    if par not in ["one"]
                ]
            )
        elif sens_method == SensitivityMethod.adjoint:
            model.set_parameter_list(
                [
                    ip
                    for ip, par in enumerate(model.get_free_parameter_ids())
                    if par not in ["one"]
                ]
            )

        solver.set_sensitivity_method(sens_method)
        check_derivatives(
            model,
            solver=solver,
            edata=edata,
            atol=1e-6,
            rtol=1e-6,
            # smaller than the offset from the trigger time
            epsilon=1e-8,
        )


@pytest.fixture(scope="module", params=[False, True], ids=["ode", "dae"])
def state_dependent_event_model_module(request, tmp_path_factory):
    """A model with a state-dependent event.

    The DAE variant has an additional algebraic state `y = x`.
    """
    dae = request.param
    model_name = f"state_dependent_event_asa{'_dae' if dae else ''}"
    rhs = "y = 1; x' = -k * y + q; 0 = y - x" if dae else "x' = -k * x + q"
    outdir = tmp_path_factory.mktemp(model_name)
    antimony2amici(
        f"""
        k = 0.5; q = 0.1
        x = 1
        {rhs}
        E1: at (x < 0.5): x = x + 1
        """,
        model_name=model_name,
        output_dir=outdir,
        observation_model=[MC("obs_x", formula="x")],
    )
    return import_model_module(model_name, outdir)


@skip_on_valgrind
@pytest.mark.parametrize(
    "interpolation_type",
    [InterpolationType.hermite, InterpolationType.polynomial],
)
@pytest.mark.parametrize("t_end", [4.5, 5.0])
def test_adjoint_sensitivities_state_dependent_event(
    state_dependent_event_model_module, t_end, interpolation_type
):
    """Test adjoint sensitivities for a parameter whose right-hand side
    derivative jumps at a state-dependent event.

    For `t_end = 5`, the event fires twice.
    """
    model = state_dependent_event_model_module.get_model()
    t = np.linspace(0.5, t_end, 10)
    model.set_timepoints(t)
    edata = ExpData(model)
    edata.set_measurements(np.full(len(t), 0.5))
    edata.set_noise_scales(np.ones(len(t)))

    sllh = {}
    for sensi_meth in (SensitivityMethod.forward, SensitivityMethod.adjoint):
        solver = model.create_solver()
        solver.set_sensitivity_order(SensitivityOrder.first)
        solver.set_sensitivity_method(sensi_meth)
        solver.set_interpolation_type(interpolation_type)
        solver.set_relative_tolerance(1e-12)
        solver.set_absolute_tolerance(1e-12)
        rdata = run_simulation(model, solver, edata)
        assert rdata.status == AMICI_SUCCESS, sensi_meth
        sllh[sensi_meth] = rdata.sllh
    assert np.all(np.abs(sllh[SensitivityMethod.forward]) > 0.1)
    assert_allclose(
        sllh[SensitivityMethod.adjoint],
        sllh[SensitivityMethod.forward],
        rtol=1e-6,
    )


@skip_on_valgrind
def test_posteq_events_are_handled(tempdir):
    """Test that events are handled during post-equilibration."""

    model_name = "test_posteq_events_are_handled"
    antimony2amici(
        r"""
        some_time = 0
        some_time' = piecewise(1, (time < 10), 0)

        bolus = 1
        target_initial = 0
        target = target_initial
        E1: at time > 1: target = target + bolus
        E2: at some_time >= 2: target = target + bolus
        """,
        observation_model=[MC("obs_target", formula="target")],
        model_name=model_name,
        output_dir=tempdir,
        verbose=True,
    )

    model_module = import_model_module(model_name, tempdir)
    model = model_module.get_model()
    solver = model.create_solver()

    # test without post-equilibration
    model.set_timepoints([10])
    rdata = run_simulation(model, solver)
    assert rdata.status == AMICI_SUCCESS
    assert rdata.by_id("target").squeeze() == 2.0
    assert rdata.by_id("obs_target").squeeze() == 2.0

    # test with post-equilibration
    model.set_steady_state_computation_mode(
        SteadyStateComputationMode.integrationOnly
    )
    model.set_steady_state_sensitivity_mode(
        SteadyStateSensitivityMode.integrationOnly
    )
    model.set_timepoints([np.inf])
    rdata = run_simulation(model, solver)
    assert rdata.status == AMICI_SUCCESS
    assert rdata.by_id("target").squeeze() == 2.0
    assert rdata.by_id("obs_target").squeeze() == 2.0
    assert rdata.posteq_t == 10.0

    # check sensitivities against finite differences
    edata = ExpData(rdata, 1, 0, 0)
    for sens_method in (
        SensitivityMethod.forward,
        SensitivityMethod.adjoint,
    ):
        solver.set_sensitivity_order(SensitivityOrder.first)
        solver.set_sensitivity_method(sens_method)

        check_derivatives(
            model,
            solver=solver,
            edata=edata,
            atol=1e-12,
            rtol=1e-7,
            epsilon=1e-8,
        )


@skip_on_valgrind
def test_preeq_presim_preserve_heaviside_state(tempdir):
    """Test the state of event trigger functions is preserved after
    pre-equilibration and presimulation.

    I.e., the trigger.initialValue is only applied in the very beginning.
    """

    model_name = "test_preeq_presim_preserve_heaviside_state"
    antimony2amici(
        r"""
        some_time = 0
        some_time' = piecewise(1, (time <= 10), 0)

        target1 = 0
        target2 = 0

        # random parameter to be able to enable
        #  pre-equilibration and presimulation
        k1 = 42

        # this should only trigger at the beginning of the first period
        E1: at true, t0=false:
            target1 = target1 + 1;
        # this will only trigger at the beginning of the second period
        #  if the trigger is not re-initialized to `true`
        E2: at some_time >= 10 and time < 1, t0 = true:
            target2 = target2 + 1;
        """,
        fixed_parameters=["k1"],
        model_name=model_name,
        output_dir=tempdir,
    )

    model_module = import_model_module(model_name, tempdir)
    model = model_module.get_model()
    model.set_timepoints(np.linspace(0, 2, 3))
    model.set_steady_state_computation_mode(
        SteadyStateComputationMode.integrationOnly
    )
    solver = model.create_solver()

    # Only main simulation. E1 triggers, E2 does not.
    rdata = run_simulation(model, solver)
    assert rdata.status == AMICI_SUCCESS
    assert list(rdata.by_id("target1")) == [1.0, 1.0, 1.0]
    assert list(rdata.by_id("target2")) == [0.0, 0.0, 0.0]

    # Pre-equilibration + main simulation. Both E1 and E2 trigger.
    edata = ExpData(rdata, 1, 0, 0)
    edata.fixed_parameters_pre_equilibration = [1]
    rdata = run_simulation(model, solver, edata=edata)
    assert rdata.status == AMICI_SUCCESS
    assert list(rdata.by_id("target1")) == [1.0, 1.0, 1.0]
    assert list(rdata.by_id("target2")) == [1.0, 1.0, 1.0]

    # Pre-simulation + main simulation. Both E1 and E2 trigger.
    # E2 (which triggers shortly after t=0) has not fired yet at the first
    # main-simulation timepoint.
    edata = ExpData(rdata, 1, 0, 0)
    edata.fixed_parameters_presimulation = [1]
    edata.t_presim = 10
    rdata = run_simulation(model, solver, edata=edata)
    assert rdata.status == AMICI_SUCCESS
    assert list(rdata.by_id("target1")) == [1.0, 1.0, 1.0]
    assert list(rdata.by_id("target2")) == [0.0, 1.0, 1.0]

    # Pre-equilibration + pre-simulation + main simulation.
    # Both E1 and E2 trigger.
    edata = ExpData(rdata, 1, 0, 0)
    edata.fixed_parameters_pre_equilibration = [1]
    edata.fixed_parameters_presimulation = [1]
    edata.t_presim = 10
    rdata = run_simulation(model, solver, edata=edata)
    assert rdata.status == AMICI_SUCCESS
    assert list(rdata.by_id("target1")) == [1.0, 1.0, 1.0]
    assert list(rdata.by_id("target2")) == [1.0, 1.0, 1.0]


@skip_on_valgrind
def test_event_noise_parameter_derivative(tempdir):
    """Test the sensitivity of the likelihood with respect to a parameter of
    an event noise formula, with the noise scale from the model and from
    the measurements."""
    from amici.sim.sundials import ParameterScaling

    model_name = "test_event_noise_parameter_derivative"
    antimony2amici(
        """
        k = 0.5; sd = 0.3
        x = 1
        x' = -k * x
        E1: at time > 1: x = x + 1
        """,
        model_name=model_name,
        output_dir=tempdir,
        observation_model=[
            MC("obs_x", formula="x"),
            MC("z1", formula="x", event_id="E1", sigma="sd"),
        ],
    )
    model = import_model_module(model_name, tempdir).get_model()
    model.set_parameter_scale(ParameterScaling.none)
    t = np.linspace(0.5, 2, 4)
    model.set_timepoints(t)
    assert model.get_free_parameter_ids() == ("k", "sd")
    sd = 0.3

    for sigma_from_edata in (False, True):
        edata = ExpData(model)
        mz = np.full(model.n_max_event() * model.nztrue, np.nan)
        mz[0] = 1.0
        edata.set_event_measurements(mz)
        sigma_z = np.full_like(mz, np.nan)
        if sigma_from_edata:
            sigma_z[0] = sd
        edata.set_event_noise_scales(sigma_z)

        solver = model.create_solver()
        solver.set_sensitivity_order(SensitivityOrder.first)
        solver.set_sensitivity_method(SensitivityMethod.forward)
        rdata = run_simulation(model, solver, edata)
        assert rdata.status == AMICI_SUCCESS

        if sigma_from_edata:
            expected = 0.0
        else:
            res = rdata.z[0, 0] - mz[0]
            expected = -1 / sd + res**2 / sd**3
            assert abs(expected) > 1
        assert np.isclose(rdata.sllh[1], expected, atol=1e-6), sigma_from_edata
        # d sigma_z1 / d (k, sd)
        assert_allclose(
            rdata.ssigmaz[0, :, 0], [0, 0 if sigma_from_edata else 1]
        )


@skip_on_valgrind
def test_heaviside_reinit_after_preeq(tempdir):
    """Test that event triggers are re-evaluated from the re-initialized
    state after pre-equilibration / pre-simulation.

    A ``piecewise`` trigger reads a state that is re-initialized from a fixed
    parameter at the beginning of the main simulation. The trigger state must
    be computed from the re-initialized state, not from the
    pre-equilibration / pre-simulation state (gh-3316).
    """
    model_name = "test_heaviside_reinit_after_preeq"
    antimony2amici(
        r"""
        drug_level_0 = 0
        switch_time_0 = 0
        species drug_level = drug_level_0
        species switch_time = switch_time_0
        species readout_y = 0
        k_prod = 1
        k_deg = 1
        input := k_prod * drug_level * piecewise(0, time < switch_time, 1)
        readout_y' = input - k_deg * readout_y
        """,
        fixed_parameters=["drug_level_0", "switch_time_0"],
        model_name=model_name,
        output_dir=tempdir,
    )
    model_module = import_model_module(model_name, tempdir)
    model = model_module.get_model()
    model.set_reinitialize_fixed_parameter_initial_states(True)
    model.set_steady_state_computation_mode(
        SteadyStateComputationMode.integrationOnly
    )
    model.set_steady_state_sensitivity_mode(
        SteadyStateSensitivityMode.integrationOnly
    )
    timepoints = np.array([0.0, 1.0, 2.0, 3.0, 5.0])
    model.set_timepoints(timepoints)
    solver = model.create_solver()
    solver.set_absolute_tolerance(1e-14)
    solver.set_relative_tolerance(1e-10)

    # pre-equilibration: no drug, switch at t=0
    # main simulation: drug on, switched on at t=2
    edata = ExpData(model)
    edata.fixed_parameters = [1.0, 2.0]
    edata.fixed_parameters_pre_equilibration = [0.0, 0.0]
    edata.reinitialize_fixed_parameter_initial_states = True

    y_expected = np.where(timepoints < 2, 0.0, 1 - np.exp(-(timepoints - 2)))

    rdata = run_simulation(model, solver, edata=edata)
    assert rdata.status == AMICI_SUCCESS
    assert_allclose(rdata.by_id("readout_y"), y_expected, atol=1e-8)

    # pre-equilibration + pre-simulation + main simulation
    # pre-simulation on [-1, 0]: drug on, switched on at t=-2, i.e., the
    # input is on throughout the pre-simulation. With the pre-equilibration
    # state (switch at t=0), the trigger would evaluate to "off" at t=-1.
    # Similarly, at the start of the main simulation (switch at t=2), the
    # trigger is "off", while it would be "on" based on the pre-simulation
    # state (switch at t=-2).
    edata.t_presim = 1.0
    edata.fixed_parameters_presimulation = [1.0, -2.0]
    y_presim_end = 1 - np.exp(-1)
    y_expected_presim = y_presim_end * np.exp(-timepoints) + y_expected
    rdata = run_simulation(model, solver, edata=edata)
    assert rdata.status == AMICI_SUCCESS
    # state at the end of the pre-simulation
    # (= initial state of the main simulation)
    assert_allclose(
        rdata.x0[model.get_state_ids().index("readout_y")],
        y_presim_end,
        atol=1e-8,
    )
    assert_allclose(rdata.by_id("readout_y"), y_expected_presim, atol=1e-8)

    # forward sensitivities
    solver.set_sensitivity_order(SensitivityOrder.first)
    solver.set_sensitivity_method(SensitivityMethod.forward)
    check_derivatives(model, solver=solver, edata=edata)

    edata.t_presim = 0.0
    edata.fixed_parameters_presimulation = []
    check_derivatives(model, solver=solver, edata=edata)


@skip_on_valgrind
def test_event_measurement_without_event(tempdir):
    """Test the likelihood and its sensitivities for a measurement of an
    event that does not occur until the last output time.

    Only the root function value at the last output time enters the
    likelihood, not the event output.
    """
    model_name = "event_measurement_without_event"
    antimony2amici(
        """
        k = 0.5; q = 0.1
        x = 1
        x' = -k * x + q
        E1: at (x < 0.1): x = x + 1
        """,
        model_name=model_name,
        output_dir=tempdir,
        observation_model=[
            MC("obs_x", formula="x"),
            MC("z1", formula="2 * x", event_id="E1"),
        ],
    )
    model = import_model_module(model_name, tempdir).get_model()
    model.set_timepoints(np.linspace(0.5, 2, 4))
    edata = _event_output_problem(model, [0.3])

    solver = model.create_solver()
    solver.set_relative_tolerance(1e-12)
    solver.set_absolute_tolerance(1e-12)
    solver.set_sensitivity_order(SensitivityOrder.first)
    function, _ = run_simulation_to_function_and_derivative(
        amici_model=model,
        amici_solver=solver,
        amici_edata=edata,
        derivative_variables=["llh"],
    )
    point = np.asarray(model.get_free_parameters())

    for sensi_meth in (SensitivityMethod.forward, SensitivityMethod.adjoint):
        solver.set_sensitivity_method(sensi_meth)
        rdata = run_simulation(model, solver, edata)
        assert rdata.status == AMICI_SUCCESS, sensi_meth
        rz = rdata.rz[0, 0]
        assert abs(rz) > 0.1
        assert np.isclose(
            rdata.llh, -0.5 * np.log(2 * np.pi) - 0.5 * rz**2, rtol=1e-12
        )
        if sensi_meth == SensitivityMethod.forward:
            assert np.all(np.abs(rdata.sllh) > 0.01)
        check_gradient(function, point, rdata.sllh).assert_success(
            always_print=True
        )


@skip_on_valgrind
def test_gh2926(tempdir):
    """Two simultaneous events. Event `E1` changes the root function
    for the piecewise-switch from 0 to <0."""

    model_name = "test_gh2926"
    antimony2amici(
        r"""
        some_time = 0
        some_time' = 1

        threshold = 0

        x1 := piecewise(2, (some_time >= threshold), 1)
        E1: at some_time >= 0, t0 = false: threshold = 10
        """,
        model_name=model_name,
        output_dir=tempdir,
    )

    model_module = import_model_module(model_name, tempdir)
    model = model_module.get_model()
    # check output just after t=10, otherwise the solver stops at `10 - eps`
    #  and the event triggers only after the output was recorded
    model.set_timepoints([0, 5, 10.1, 11])
    solver = model.create_solver()

    rdata = run_simulation(model, solver)
    assert rdata.status == AMICI_SUCCESS
    assert rdata.by_id("x1").tolist() == [1.0, 1.0, 2.0, 2.0]


@skip_on_valgrind
def test_root_after_reinit_direction_is_checked(tempdir, caplog):
    """A root reported immediately after reinitialization must only be
    ignored if it matches the *direction* of a root recorded during the
    preceding discontinuity handling, not just the index.

    ``EA`` fires periodically and its assignment sets ``y`` to just
    below ``EB``'s threshold, causing a down-crossing of ``EB`` to be
    recorded as part of ``EA``'s own discontinuity handling.
    Within the post-reinitialization solver step, ``y`` increases,
    and `EB` immediately crosses back over the threshold,
    triggering ``EB``. If the post-reinit check skips all roots of a given
    event, ``EB`` would never fire, and ``y`` would never be reset.
    """
    model_name = "test_root_after_reinit_direction_is_checked"
    eps = 1e-15
    antimony2amici(
        rf"""
        x = 0
        x' = 1
        y = 5
        y' = 1

        n_fired_A = 0
        n_fired_B = 0

        EA: at x >= 1:
            n_fired_A = n_fired_A + 1,
            x = 0,
            y = 1 - {eps};

        EB: at y >= 1:
            n_fired_B = n_fired_B + 1,
            y = 5;
        """,
        model_name=model_name,
        output_dir=tempdir,
    )

    model_module = import_model_module(model_name, tempdir)
    model = model_module.get_model()
    n_periods = 5
    model.set_timepoints([n_periods + 0.5])
    solver = model.create_solver()

    rdata = run_simulation(model, solver)
    assert rdata.status == AMICI_SUCCESS
    assert rdata.by_id("n_fired_A")[-1] == n_periods
    assert rdata.by_id("n_fired_B")[-1] == n_periods
    assert "ROOT_AFTER_REINIT" not in caplog.text


@skip_on_valgrind
def test_root_after_reinit_ignores_just_fired_event(tempdir, caplog):
    """A root reported immediately after reinitialization, for the exact
    same index and direction as the event just fired, is treated as a
    spurious re-detection and ignored -- even though, read literally, the
    dynamics below really do cross the threshold, drop back below it, and
    cross again, i.e. two distinct trigger events did occur.

    Whether an immediate same-index, same-direction root reflects a
    genuine retrigger or floating-point noise is undecidable from inside
    the solver once the recurrence interval is below what the current
    tolerances can resolve.

    ``recovering`` is set by ``E1`` itself and never reset, so once ``x``
    is assigned a value just below the threshold, its derivative switches
    to an enormous rate, so the "immediate re-trigger" reliably falls
    within the first post-reinitialization solver step.

    ``E2`` is a second, independent event whose own threshold is crossed
    (for the first, genuine time) at the very same instant as ``E1``'s
    "spurious" repeat. It must still fire normally: suppression is scoped
    to the exact (index, direction) pair that was just handled, not to
    "anything reported in this step".
    """
    model_name = "test_root_after_reinit_ignores_just_fired_event"
    margin_below_threshold = 1e-13
    antimony2amici(
        rf"""
        x = 0
        y = 0
        z = 0
        recovering = 0
        recovering' = 0
        k = 1
        x' = k + 1e10 * recovering
        y' = 1e10 * recovering
        z' = k - z

        n_fired_1 = 0
        n_fired_2 = 0
        threshold = 1

        E1: at x >= threshold:
            n_fired_1 = n_fired_1 + 1,
            recovering = 1,
            x = {1.0 - margin_below_threshold},
            y = {1.0 - margin_below_threshold},
            z = z + 1;

        E2: at y >= threshold:
            n_fired_2 = n_fired_2 + 1,
            z = z + 2;
        """,
        model_name=model_name,
        output_dir=tempdir,
        observation_model=[MC("obs_z", formula="z")],
    )

    model_module = import_model_module(model_name, tempdir)
    model = model_module.get_model()
    t_end = 1.5
    model.set_timepoints([t_end])
    solver = model.create_solver()

    rdata = run_simulation(model, solver)
    assert rdata.status == AMICI_SUCCESS
    # E1's own immediate repeat is (deliberately) suppressed.
    assert rdata.by_id("n_fired_1")[-1] == 1
    # E2's genuine, simultaneous first crossing is still correctly handled.
    assert rdata.by_id("n_fired_2")[-1] == 1
    assert "ROOT_AFTER_REINIT" in caplog.text

    # Analytical solution: `z` jumps by 3 at t1 = threshold / k (1 from E1,
    # 2 from E2), so z(T) = k + (3 - k exp(-t1)) exp(-(T - t1)).
    k, threshold = model.get_free_parameters()
    assert list(model.get_free_parameter_ids()) == ["k", "threshold"]
    t1 = threshold / k
    z_end = k + (3 - k * np.exp(-t1)) * np.exp(-(t_end - t1))
    dz_dt1 = 3 * np.exp(-(t_end - t1))
    dz_dp = np.array(
        [1 - np.exp(-t_end) - dz_dt1 * threshold / k**2, dz_dt1 / k]
    )
    assert_allclose(rdata.by_id("z")[-1], z_end, rtol=1e-6)

    # Forward and adjoint sensitivities must both include E2's event-time
    # sensitivity, i.e., match the analytical gradient.
    measurement = 3.0
    edata = ExpData(rdata, 1, 0, 0)
    edata.set_measurements([measurement])
    expected_sllh = (measurement - z_end) * dz_dp
    for sens_method in (
        SensitivityMethod.forward,
        SensitivityMethod.adjoint,
    ):
        solver = model.create_solver()
        solver.set_sensitivity_order(SensitivityOrder.first)
        solver.set_sensitivity_method(sens_method)
        solver.set_relative_tolerance(1e-10)
        solver.set_absolute_tolerance(1e-14)
        rdata = run_simulation(model, solver, edata=edata)
        assert rdata.status == AMICI_SUCCESS
        assert_allclose(
            rdata.sllh, expected_sllh, rtol=1e-6, err_msg=str(sens_method)
        )


def test_root_after_reinit_ignores_all_simultaneous_events(tempdir):
    """Spurious re-detections are suppressed for every event of a
    simultaneous group, not only for the first one processed.

    ``E1`` and ``E2`` share the same trigger and fire together. ``E1``
    resets ``x`` to just below the threshold and switches on an enormous
    rate, so the solver re-reports *both* roots right after
    reinitialization. Both must be recognized as already handled.
    """
    model_name = "test_root_after_reinit_ignores_all_simultaneous_events"
    antimony2amici(
        """
        x = 0
        recovering = 0
        recovering' = 0
        x' = 1 + 1e10 * recovering

        n_fired_1 = 0
        n_fired_2 = 0

        E1: at x >= 1:
            n_fired_1 = n_fired_1 + 1,
            recovering = 1,
            x = 1 - 1e-13;

        E2: at x >= 1:
            n_fired_2 = n_fired_2 + 1;
        """,
        model_name=model_name,
        output_dir=tempdir,
    )

    model_module = import_model_module(model_name, tempdir)
    model = model_module.get_model()
    model.set_timepoints([1.5])

    rdata = run_simulation(model, model.create_solver())
    assert rdata.status == AMICI_SUCCESS
    assert rdata.by_id("n_fired_1")[-1] == 1
    assert rdata.by_id("n_fired_2")[-1] == 1


@skip_on_valgrind
@pytest.mark.parametrize(
    "from_trigger",
    [(True, True), (False, False), (False, True), (True, False)],
    ids=["true", "false", "false_true", "true_false"],
)
@pytest.mark.parametrize("e1_first", [True, False])
def test_simultaneous_events_sensitivities(tempdir, from_trigger, e1_first):
    """Sensitivities for simultaneous events, where the earlier one changes a
    state used by the later one.

    `from_trigger` is the `fromTrigger` attribute of `E1` and `E2`.
    """
    ft1, ft2 = from_trigger
    ft1_str, ft2_str = str(ft1).lower(), str(ft2).lower()
    model_name = (
        f"test_simultaneous_events_{int(ft1)}{int(ft2)}{int(e1_first)}"
    )
    prio1, prio2 = (2, 1) if e1_first else (1, 2)
    antimony2amici(
        f"""
        k = 0.5; p = 3
        x = 1; v = 0
        x' = -k * p * x
        v' = 0
        E1: at p * x < 0.9, priority={prio1}, fromTrigger={ft1_str}: x = 0.8 * x
        E2: at p * x < 0.9, priority={prio2}, fromTrigger={ft2_str}: v = v + 2 * p * x
        """,
        model_name=model_name,
        output_dir=tempdir,
        observation_model=[
            MC("obs_x", formula="x"),
            MC("obs_v", formula="v"),
        ],
    )
    model = import_model_module(model_name, tempdir).get_model()
    model.set_timepoints(np.linspace(0, 4, 9))

    solver = model.create_solver()
    solver.set_relative_tolerance(1e-12)
    solver.set_absolute_tolerance(1e-12)
    rdata = run_simulation(model, solver)
    assert rdata.status == AMICI_SUCCESS
    # E2 uses the value of x before E1, unless E1 is executed first and
    #  E2's assignment is evaluated at execution time
    x_e2 = 0.9 / 3 * (0.8 if e1_first and not ft2 else 1)
    assert_allclose(rdata.by_id("v")[-1], 2 * 3 * x_e2, rtol=1e-8)

    edata = ExpData(rdata, 1.0, 0.0)
    # fixed offset, so that the gradient is non-zero
    edata.set_measurements((np.asarray(rdata.y) + 0.5).flatten())
    function, _ = run_simulation_to_function_and_derivative(
        amici_model=model,
        amici_solver=solver,
        amici_edata=edata,
        derivative_variables=["llh"],
    )
    point = np.asarray(model.get_free_parameters())
    solver.set_sensitivity_order(SensitivityOrder.first)
    sllh = {}
    for sensi_meth in (
        SensitivityMethod.forward,
        SensitivityMethod.adjoint,
    ):
        solver.set_sensitivity_method(sensi_meth)
        rdata = run_simulation(model, solver, edata)
        assert rdata.status == AMICI_SUCCESS, sensi_meth
        assert np.all(np.abs(rdata.sllh) > 0.1), sensi_meth
        sllh[sensi_meth] = rdata.sllh
        # Explicit tolerance: the auto-derived one is too loose here, as the
        #  largest finite-difference steps cross the event.
        check_gradient(function, point, rdata.sllh, atol=1e-4).assert_success(
            always_print=True
        )
    assert_allclose(
        sllh[SensitivityMethod.adjoint],
        sllh[SensitivityMethod.forward],
        rtol=1e-6,
    )


@skip_on_valgrind
def test_simultaneous_events_heaviside_crossed_by_assignment(tempdir):
    """Adjoint sensitivities for simultaneous events, where the assignment of
    the first one crosses a `piecewise` switch before the second one is
    executed.

    The switch gets its own discontinuity, but shares the event time of the
    triggering events.
    """
    model_name = "test_simultaneous_events_heaviside_crossed"
    antimony2amici(
        """
        k = 0.5; p = 2.1
        x = 1; y = 0; z = 0
        x' = -k * x
        y' = piecewise(1, x > 2, 0)
        z' = 0
        E1: at time > p, priority=2: x = x + 2
        E2: at time > p, priority=1, fromTrigger=false: z = z + x
        """,
        model_name=model_name,
        output_dir=tempdir,
        observation_model=[
            MC("obs_x", formula="x"),
            MC("obs_y", formula="y"),
            MC("obs_z", formula="z"),
        ],
    )
    model = import_model_module(model_name, tempdir).get_model()
    model.set_timepoints(np.linspace(0, 6, 13))

    solver = model.create_solver()
    solver.set_relative_tolerance(1e-12)
    solver.set_absolute_tolerance(1e-12)
    rdata = run_simulation(model, solver)
    assert rdata.status == AMICI_SUCCESS
    # the switch is on right after the events
    assert rdata.by_id("y")[-1] > 0

    edata = ExpData(rdata, 1.0, 0.0)
    # fixed offset, so that the gradient is non-zero
    edata.set_measurements((np.asarray(rdata.y) + 0.5).flatten())
    solver.set_sensitivity_order(SensitivityOrder.first)
    sllh = {}
    for sensi_meth in (
        SensitivityMethod.forward,
        SensitivityMethod.adjoint,
    ):
        solver.set_sensitivity_method(sensi_meth)
        rdata = run_simulation(model, solver, edata)
        assert rdata.status == AMICI_SUCCESS, sensi_meth
        assert np.all(np.abs(rdata.sllh) > 0.1), sensi_meth
        sllh[sensi_meth] = rdata.sllh
    assert_allclose(
        sllh[SensitivityMethod.adjoint],
        sllh[SensitivityMethod.forward],
        rtol=1e-6,
    )


@skip_on_valgrind
@pytest.mark.parametrize("levels", [1, 2])
def test_trigger_time_values_after_secondary_event(tempdir, levels):
    """Sensitivities for an event using values from the trigger time, which
    is executed after secondary events were triggered.

    `E1` triggers `E3`, which, for two levels, triggers `E4` (otherwise,
    `E4` triggers later). `E2` is executed last and uses the value of `x`
    before `E1`, i.e., `v = 0.9 / p`.
    """
    model_name = f"test_trigger_time_values_after_secondary_event_{levels}"
    e3_x = ", x = x / 2" if levels == 2 else ""
    antimony2amici(
        f"""
        k = 0.5; p = 3
        x = 1; v = 0; w = 0; y = 0
        x' = -k * p * x
        v' = 0
        w' = 0
        y' = 0
        E1: at p * x < 0.9, priority=5, fromTrigger=false: x = x / 4
        E2: at p * x < 0.9, priority=1, fromTrigger=true: v = v + x
        E3: at x < 0.1, priority=4, fromTrigger=false: w = w + 1{e3_x}
        E4: at x < 0.05, priority=3, fromTrigger=false: y = y + 1
        """,
        model_name=model_name,
        output_dir=tempdir,
        observation_model=[
            MC("obs_x", formula="x"),
            MC("obs_v", formula="v"),
            MC("obs_w", formula="w"),
            MC("obs_y", formula="y"),
        ],
    )
    model = import_model_module(model_name, tempdir).get_model()
    model.set_timepoints(np.linspace(0, 4, 9))

    solver = model.create_solver()
    solver.set_relative_tolerance(1e-12)
    solver.set_absolute_tolerance(1e-12)
    rdata = run_simulation(model, solver)
    assert rdata.status == AMICI_SUCCESS
    assert_allclose(rdata.by_id("v")[-1], 0.9 / 3, rtol=1e-8)
    assert_allclose(rdata.by_id("w")[-1], 1.0)
    assert_allclose(rdata.by_id("y")[-1], 1.0)

    # dv/dk = 0, dv/dp = -0.9 / p^2
    solver.set_sensitivity_order(SensitivityOrder.first)
    solver.set_sensitivity_method(SensitivityMethod.forward)
    rdata = run_simulation(model, solver)
    assert rdata.status == AMICI_SUCCESS
    assert model.get_free_parameter_ids() == ("k", "p")
    ix_v = model.get_state_ids().index("v")
    assert_allclose(rdata.sx[-1, :, ix_v], [0.0, -0.1], rtol=1e-6, atol=1e-10)
    solver.set_sensitivity_order(SensitivityOrder.none)

    edata = ExpData(rdata, 1.0, 0.0)
    # fixed offset, so that the gradient is non-zero
    edata.set_measurements((np.asarray(rdata.y) + 0.5).flatten())
    function, _ = run_simulation_to_function_and_derivative(
        amici_model=model,
        amici_solver=solver,
        amici_edata=edata,
        derivative_variables=["llh"],
    )
    point = np.asarray(model.get_free_parameters())
    solver.set_sensitivity_order(SensitivityOrder.first)
    sllh = {}
    for sensi_meth in (
        SensitivityMethod.forward,
        SensitivityMethod.adjoint,
    ):
        solver.set_sensitivity_method(sensi_meth)
        rdata = run_simulation(model, solver, edata)
        assert rdata.status == AMICI_SUCCESS, sensi_meth
        assert np.all(np.abs(rdata.sllh) > 0.05), sensi_meth
        check_gradient(function, point, rdata.sllh, atol=1e-4).assert_success(
            always_print=True
        )
        sllh[sensi_meth] = rdata.sllh
    assert_allclose(
        sllh[SensitivityMethod.adjoint],
        sllh[SensitivityMethod.forward],
        rtol=1e-6,
    )


def _event_output_problem(model, event_measurements):
    """Experimental data with measurements only for the event outputs."""
    nt = model.nt()
    edata = ExpData(model)
    edata.set_measurements(np.full(nt * model.nytrue, np.nan))
    edata.set_noise_scales(np.ones(nt * model.nytrue))
    mz = np.full(model.n_max_event() * model.nztrue, np.nan)
    mz[: len(event_measurements)] = event_measurements
    edata.set_event_measurements(mz)
    edata.set_event_noise_scales(np.ones_like(mz))
    return edata


@skip_on_valgrind
def test_adjoint_sensitivities_event_outputs(tempdir):
    """Adjoint sensitivities for event outputs, with fewer events than
    `n_max_event` and with as many events as `n_max_event`."""
    model_name = "test_adjoint_sensitivities_event_outputs"
    antimony2amici(
        """
        k = 0.5; q = 0.1
        x = 1
        x' = -k * x + q
        E1: at (time > 1.3): x = x + 1
        """,
        model_name=model_name,
        output_dir=tempdir,
        observation_model=[
            MC("obs_x", formula="x"),
            MC("z1", formula="x", event_id="E1"),
        ],
    )
    model = import_model_module(model_name, tempdir).get_model()
    model.set_timepoints(np.linspace(0.5, 4.5, 9))
    assert model.n_max_event() > 1

    solver = model.create_solver()
    solver.set_relative_tolerance(1e-12)
    solver.set_absolute_tolerance(1e-12)
    solver.set_sensitivity_order(SensitivityOrder.first)

    # E1 occurs once
    for n_max_event in (model.n_max_event(), 1):
        model.set_n_max_event(n_max_event)
        edata = _event_output_problem(model, [0.3])
        function, _ = run_simulation_to_function_and_derivative(
            amici_model=model,
            amici_solver=solver,
            amici_edata=edata,
            derivative_variables=["llh"],
        )
        point = np.asarray(model.get_free_parameters())

        for sensi_meth in (
            SensitivityMethod.forward,
            SensitivityMethod.adjoint,
        ):
            solver.set_sensitivity_method(sensi_meth)
            rdata = run_simulation(model, solver, edata)
            msg = f"{sensi_meth=}, {n_max_event=}"
            assert rdata.status == AMICI_SUCCESS, msg
            if sensi_meth == SensitivityMethod.forward:
                assert np.all(np.abs(rdata.sllh) > 0.1), msg
            check_gradient(function, point, rdata.sllh).assert_success(
                always_print=True
            )


@skip_on_valgrind
def test_explicit_trigger_time_changed_by_events(tempdir):
    """Trigger times that depend on parameters changed by events.

    ``level``, ``switch_time``, and ``other_switch_time`` are targets of
    event assignments and are thus converted to state variables with zero
    time derivative. The trigger times of the switches of ``u`` and ``v``
    therefore only change at events, and the switches are handled at their
    precomputed trigger times instead of by root-finding (gh-3365).

    The structure mimics a PEtab v2 problem with time-dependent stimuli:
    The condition changes are applied by ``E0`` at ``t = 0``, after
    pre-equilibration without the stimuli. ``v`` is switched on before
    ``u``, which must not affect ``u``. ``E1`` postpones the switch of
    ``u``, switching ``u`` off and on again.

    The fast state ``x`` makes the solver reduce the step size below the
    floating point resolution of ``t`` right after the switch of ``u``. For
    root functions that are exactly zero at that time, root-finding can
    report a spurious reverse crossing.
    """
    model_name = "test_explicit_trigger_time_changed_by_events"
    antimony2amici(
        r"""
        preeq_indicator = 0
        level = 0
        level_on = 3
        switch_time = 0
        other_switch_time = 0
        t_on = 60
        t_other = 30
        t_off = 70
        t_on2 = 100
        k = 1
        x = 0
        u := level * piecewise(0, time < switch_time, 1)
        v := piecewise(0, time < other_switch_time, 1)
        x' = k * (u - x)
        E0: at (preeq_indicator < 0.5 && time >= 0), t0 = false:
            level = level_on, switch_time = t_on,
            other_switch_time = t_other;
        E1: at time >= t_off: switch_time = t_on2;
        """,
        fixed_parameters=["preeq_indicator"],
        model_name=model_name,
        output_dir=tempdir,
    )
    model = import_model_module(model_name, tempdir).get_model()
    # E0, E1, and two root functions for each of the piecewise switches,
    #  of which only the trigger of E0 requires root-finding
    assert model.ne == 6
    assert model.ne_solver == 1

    def create_edatas():
        """Create ExpData objects without and with pre-equilibration."""
        edata = ExpData(model)
        edata.fixed_parameters = [0.0]
        edata_preeq = ExpData(edata)
        edata_preeq.fixed_parameters_pre_equilibration = [1.0]
        return edata, edata_preeq

    # fast `x`, and tight tolerances, as in gh-3365
    model.set_free_parameter_by_id("k", 1e8)
    timepoints = np.array(
        [0, 29.99, 30, 30.01, 59.99, 60, 60.01, 61, 69.99, 70, 70.01]
        + [99.99, 100, 100.01]
    )
    model.set_timepoints(timepoints)
    u_expected = np.where(
        ((60 <= timepoints) & (timepoints < 70)) | (timepoints >= 100), 3, 0
    )
    v_expected = np.where(timepoints < 30, 0, 1)
    solver = model.create_solver()
    solver.set_absolute_tolerance(1e-16)
    solver.set_relative_tolerance(1e-8)
    for edata in create_edatas():
        rdata = run_simulation(model, solver, edata)
        assert rdata.status == AMICI_SUCCESS
        assert_allclose(rdata.by_id("u"), u_expected)
        assert_allclose(rdata.by_id("v"), v_expected)

    # analytical solution and sensitivities
    model.set_free_parameter_by_id("k", 1.0)
    timepoints = np.array([0, 45, 65, 69, 75, 99, 105, 120], dtype=float)
    model.set_timepoints(timepoints)
    x_70 = 3 * (1 - np.exp(-10))
    x_100 = x_70 * np.exp(-30)
    x_expected = np.select(
        [timepoints < 60, timepoints < 70, timepoints < 100],
        [
            0,
            3 * (1 - np.exp(-(timepoints - 60))),
            x_70 * np.exp(-(timepoints - 70)),
        ],
        3 + (x_100 - 3) * np.exp(-(timepoints - 100)),
    )
    solver = model.create_solver()
    solver.set_absolute_tolerance(1e-14)
    solver.set_relative_tolerance(1e-10)
    edata, edata_preeq = create_edatas()
    for cur_edata in (edata, edata_preeq):
        rdata = run_simulation(model, solver, cur_edata)
        assert rdata.status == AMICI_SUCCESS
        assert_allclose(rdata.by_id("x"), x_expected, rtol=1e-8, atol=1e-12)

    solver.set_sensitivity_order(SensitivityOrder.first)
    for sens_method in (SensitivityMethod.forward, SensitivityMethod.adjoint):
        solver.set_sensitivity_method(sens_method)
        check_derivatives(model, solver=solver, edata=edata)
    solver.set_sensitivity_method(SensitivityMethod.forward)
    check_derivatives(model, solver=solver, edata=edata_preeq)


@skip_on_valgrind
def test_explicit_trigger_time_preequilibration(tempdir):
    """Explicit trigger times that depend on parameters changed by events,
    during pre-equilibration by simulation and the main simulation.

    During pre-equilibration, ``u`` is switched on at ``switch_time``.
    The direction of this crossing is determined from the root function
    values at the start of the pre-equilibration.

    The trigger time of ``E_late`` is ``1 / rate``, i.e., infinite, until
    ``E_rate`` sets ``rate`` at the start of the main simulation.
    """
    model_name = "test_explicit_trigger_time_preequilibration"
    antimony2amici(
        r"""
        preeq_indicator = 0
        switch_time = 5
        rate = 0
        rate_main = 0.05
        n_late = 0
        k = 1
        x = 2
        u := piecewise(0, time < switch_time, 1)
        x' = u - k * x
        E_rate: at (preeq_indicator < 0.5 && time >= 0), t0 = false:
            rate = rate_main;
        E_late: at time * rate >= 1: n_late = n_late + 1;
        """,
        fixed_parameters=["preeq_indicator"],
        model_name=model_name,
        output_dir=tempdir,
    )
    model = import_model_module(model_name, tempdir).get_model()
    # E_rate, E_late, and two root functions for the piecewise switch of
    #  `u`, of which only the trigger of E_rate requires root-finding
    assert model.ne == 4
    assert model.ne_solver == 1
    model.set_steady_state_computation_mode(
        SteadyStateComputationMode.integrationOnly
    )
    model.set_steady_state_sensitivity_mode(
        SteadyStateSensitivityMode.integrationOnly
    )
    timepoints = np.array([0, 4, 6, 10, 19, 21, 30], dtype=float)
    model.set_timepoints(timepoints)
    edata = ExpData(model)
    edata.fixed_parameters = [0.0]
    edata.fixed_parameters_pre_equilibration = [1.0]

    solver = model.create_solver()
    solver.set_absolute_tolerance(1e-14)
    solver.set_relative_tolerance(1e-10)
    rdata = run_simulation(model, solver, edata)
    assert rdata.status == AMICI_SUCCESS

    # pre-equilibration: `u` switched on at t=5 -> x = 1 / k
    # main simulation: x(0) = 1, `u` switched on at t=5 again
    x_5 = np.exp(-5)
    x_expected = np.where(
        timepoints < 5,
        np.exp(-timepoints),
        1 + (x_5 - 1) * np.exp(-(timepoints - 5)),
    )
    assert_allclose(rdata.by_id("x"), x_expected, rtol=1e-6, atol=1e-8)
    # E_late triggers at 1 / rate_main = 20
    assert_allclose(
        rdata.by_id("n_late"), np.where(timepoints < 20, 0, 1), atol=0
    )

    solver.set_sensitivity_order(SensitivityOrder.first)
    solver.set_sensitivity_method(SensitivityMethod.forward)
    check_derivatives(model, solver=solver, edata=edata)


@skip_on_valgrind
def test_adjoint_sensitivities_outputs_of_multiple_events(tempdir):
    """Test adjoint sensitivities for measurements of the outputs of
    different events with the same occurrence index."""
    from amici.sim.sundials import ParameterScaling

    model_name = "test_adjoint_sensitivities_outputs_of_multiple_events"
    model = antimony2amici(
        """
        k = 0.5; x0 = 1
        x = x0; y = 0
        x' = -k * x
        y' = 0
        E1: at time > 1.3: y = y + 1
        E2: at x < 0.1: x = x + 1
        """,
        model_name=model_name,
        output_dir=tempdir,
        observation_model=[
            MC("obs_x", formula="x"),
            MC("z1", formula="x", event_id="E1"),
            MC("z2", formula="time", event_id="E2"),
        ],
    )
    model.set_parameter_scale(ParameterScaling.none)
    t = np.linspace(0.5, 2, 4)
    model.set_timepoints(t)
    assert model.get_free_parameter_ids() == ("k", "x0")

    # E1 triggers at t1, E2 does not trigger until the last output time T
    k, x0, t1, T = 0.5, 1.0, 1.3, t[-1]
    x1, xT = x0 * np.exp(-k * t1), x0 * np.exp(-k * T)
    dx1 = np.array([-t1 * x1, x1 / x0])
    dxT = np.array([-T * xT, xT / x0])

    # measurements of z1 only, of rz2 only, and of both
    for mz_set in ([1.0, np.nan], [np.nan, 3.0], [1.0, 3.0]):
        edata = ExpData(model)
        edata.set_measurements(np.full(len(t), np.nan))
        edata.set_noise_scales(np.ones(len(t)))
        mz = np.full(model.n_max_event() * model.nztrue, np.nan)
        mz[:2] = mz_set
        edata.set_event_measurements(mz)
        edata.set_event_noise_scales(np.ones_like(mz))

        # Jz(z2) does not depend on the parameters (z2 = T),
        # the root function value of E2 is +-(x - 0.1)
        sllh_expected = np.zeros(2)
        if not np.isnan(mz_set[0]):
            sllh_expected -= (x1 - mz_set[0]) * dx1
        if not np.isnan(mz_set[1]):
            sllh_expected -= (xT - 0.1) * dxT
        assert np.all(np.abs(sllh_expected) > 0.05)

        for method in (SensitivityMethod.forward, SensitivityMethod.adjoint):
            solver = model.create_solver()
            solver.set_relative_tolerance(1e-12)
            solver.set_absolute_tolerance(1e-12)
            solver.set_sensitivity_order(SensitivityOrder.first)
            solver.set_sensitivity_method(method)
            rdata = run_simulation(model, solver, edata)
            assert rdata.status == AMICI_SUCCESS
            assert_allclose(
                rdata.sllh,
                sllh_expected,
                rtol=1e-6,
                err_msg=f"{method.name}, {mz_set=}",
            )


@skip_on_valgrind
def test_adjoint_sensitivities_outputs_of_repeated_events(tempdir):
    """Test adjoint sensitivities for measurements of several outputs of
    several events that occur more often than `n_max_event`."""
    model_name = "test_adjoint_sensitivities_outputs_of_repeated_events"
    antimony2amici(
        """
        k = 0.1; q = 0.2
        x = 1; y = 0; c1 = 0; c2 = 0.5
        x' = -k * x
        y' = q * x
        c1' = 1
        c2' = 1
        E1: at c1 > 1: c1 = 0, x = x + 1
        E2: at c2 > 1: c2 = 0, y = y / 2
        """,
        model_name=model_name,
        output_dir=tempdir,
        observation_model=[
            MC("obs_x", formula="x"),
            MC("z1a", formula="y", event_id="E1"),
            MC("z1b", formula="x * y", event_id="E1"),
            MC("z2", formula="y", event_id="E2"),
        ],
    )
    model = import_model_module(model_name, tempdir).get_model()
    # E1 occurs at t = 1, 2, 3, 4, E2 at t = 0.5, 1.5, 2.5, 3.5
    model.set_timepoints(np.linspace(0.2, 4.2, 10))
    model.set_n_max_event(3)
    edata = _event_output_problem(
        model, np.arange(model.n_max_event() * model.nztrue) / 10
    )

    solver = model.create_solver()
    solver.set_relative_tolerance(1e-12)
    solver.set_absolute_tolerance(1e-12)
    solver.set_sensitivity_order(SensitivityOrder.first)

    rdatas = {}
    for sensi_meth in (SensitivityMethod.forward, SensitivityMethod.adjoint):
        solver.set_sensitivity_method(sensi_meth)
        rdatas[sensi_meth] = run_simulation(model, solver, edata)
        assert rdatas[sensi_meth].status == AMICI_SUCCESS

    rdata = rdatas[SensitivityMethod.forward]
    assert np.all(np.abs(rdata.sllh) > 1)
    assert_allclose(
        rdatas[SensitivityMethod.adjoint].sllh, rdata.sllh, rtol=1e-6
    )


@skip_on_valgrind
def test_event_with_w_dependent_trigger(tempdir):
    """Test sensitivities for events with trigger depending on
    cascading expressions in `w`."""

    model_name = "test_event_with_w_dependent_trigger"
    model = antimony2amici(
        r"""
        one = 1
        two = 2
        a := two^2 - 2 # 2
        b := a^2 - 1  # 3
        c := b * 2 + x / 10  # 6 + time / 10
        d := c + a + (a - 1) * time / 10 # -> d = 8 + time / 5

        x = 0
        x' = a / 2  # x' = 1
        target = 0

        E1: at x >= d:  # triggers at time = 8 + time / 5 <=> time = 10
            target = x + one;
        """,
        model_name=model_name,
        output_dir=tempdir,
    )

    model.set_timepoints([0, 5, 9, 11])
    solver = model.create_solver()
    solver.set_sensitivity_order(SensitivityOrder.first)
    solver.set_sensitivity_method(SensitivityMethod.forward)
    solver.set_relative_tolerance(1e-14)

    # generate synthetic measurements
    rdata = run_simulation(model, solver)
    assert rdata.status == AMICI_SUCCESS
    # check that event triggered correctly
    assert np.isclose(rdata.by_id("target")[-1], 11.0)
    edata = ExpData(rdata, 1.0, 0.0, 42)

    # check sensitivities against finite differences

    for sens_method in (
        SensitivityMethod.forward,
        SensitivityMethod.adjoint,
    ):
        solver.set_sensitivity_method(sens_method)
        check_derivatives(
            model,
            solver=solver,
            edata=edata,
            atol=1e-6,
            rtol=1e-6,
            epsilon=1e-8,
        )


@skip_on_valgrind
def test_dae_event_reinitialization(tempdir):
    """Simulation of a DAE with an event assignment to a differential state.

    After the event, the algebraic state and the derivatives have to be
    consistent with the updated differential state.
    """
    ant_model = """
    model dae_event_reinit
        const k = 0.5, p = 2, q = 0.3
        x = 1; y = 2
        x' = -k * y + q
        0 = y - p * x
        at (time >= 2.2): x = x + q
    end
    """
    module_name = "dae_event_reinit"
    antimony2amici(ant_model, model_name=module_name, output_dir=tempdir)
    model = import_model_module(module_name, tempdir).get_model()

    # analytical solution: x' = -k * p * x + q
    k, p, q = 0.5, 2.0, 0.3
    a = k * p
    x_inf = q / a
    t_event = 2.2
    x_pre_event = x_inf + (1.0 - x_inf) * np.exp(-a * t_event)

    # including an output at the event time, or within roundoff after it
    for t_output in (t_event, np.nextafter(t_event, np.inf)):
        t = np.sort(np.append(np.linspace(0, 5, 11), t_output))
        model.set_timepoints(t)
        x_expected = np.where(
            t < t_event,
            x_inf + (1.0 - x_inf) * np.exp(-a * t),
            x_inf + (x_pre_event + q - x_inf) * np.exp(-a * (t - t_event)),
        )
        y_expected = p * x_expected

        for rtol in (1e-6, 1e-12):
            msg = f"{t_output=}, {rtol=}"
            solver = model.create_solver()
            solver.set_relative_tolerance(rtol)
            solver.set_absolute_tolerance(rtol / 100)
            rdata = run_simulation(model, solver)
            assert rdata.status == AMICI_SUCCESS, msg
            assert_allclose(
                rdata.x[:, 0], x_expected, rtol=100 * rtol, err_msg=msg
            )
            assert_allclose(
                rdata.x[:, 1], y_expected, rtol=100 * rtol, err_msg=msg
            )


@skip_on_valgrind
def test_dae_event_forward_sensitivities(tempdir):
    """Forward sensitivities of a DAE with a state-dependent event.

    The event changes a differential state that the algebraic state depends
    on, so the algebraic state and the time derivative of the differential
    state jump at the event.
    """
    ant_model = """
    model dae_event_fsa
        const k = 0.5, p = 2, q = 0.3
        x = 1; y = 2
        x' = -k * y + q
        0 = y - p * x
        at (x < 0.5): x = x + 1
    end
    """
    module_name = "dae_event_fsa"
    antimony2amici(ant_model, model_name=module_name, output_dir=tempdir)
    model = import_model_module(module_name, tempdir).get_model()
    # output times before the second event at t ~ 3.04
    t = np.linspace(0.1, 2.5, 25)
    model.set_timepoints(t)

    def analytical_solution(k, p, q):
        # x' = -k * p * x + q, x(t_event) = 0.5
        a = k * p
        x_inf = q / a
        t_event = np.log((1.0 - x_inf) / (0.5 - x_inf)) / a
        x = np.where(
            t < t_event,
            x_inf + (1.0 - x_inf) * np.exp(-a * t),
            x_inf + (1.5 - x_inf) * np.exp(-a * (t - t_event)),
        )
        return np.stack([x, p * x], axis=1)

    pars = np.array([0.5, 2.0, 0.3])
    assert model.get_free_parameter_ids() == ("k", "p", "q")
    h = 1e-6
    sx_expected = np.stack(
        [
            (
                analytical_solution(*(pars + h * e))
                - analytical_solution(*(pars - h * e))
            )
            / (2 * h)
            for e in np.eye(len(pars))
        ],
        axis=1,
    )

    solver = model.create_solver()
    solver.set_sensitivity_order(SensitivityOrder.first)
    solver.set_sensitivity_method(SensitivityMethod.forward)
    solver.set_relative_tolerance(1e-12)
    solver.set_absolute_tolerance(1e-12)
    solver.set_relative_tolerance_fsa(1e-10)
    solver.set_absolute_tolerance_fsa(1e-10)
    rdata = run_simulation(model, solver)
    assert rdata.status == AMICI_SUCCESS
    assert_allclose(rdata.x, analytical_solution(*pars), rtol=1e-8)
    assert_allclose(rdata.sx, sx_expected, rtol=1e-6, atol=1e-8)


@skip_on_valgrind
def test_dae_algebraic_rule_piecewise(tempdir):
    """Forward sensitivities of a DAE with a state-dependent `piecewise` in
    an algebraic rule."""
    ant_model = """
    model dae_algebraic_rule_piecewise
        const k = 0.5, p = 3
        x = 1; y = 3
        x' = -k * y
        0 = y - piecewise(p * x, x > 0.3, x)
    end
    """
    module_name = "dae_algebraic_rule_piecewise"
    antimony2amici(ant_model, model_name=module_name, output_dir=tempdir)
    model = import_model_module(module_name, tempdir).get_model()
    # the switch is handled as an event
    assert model.ne > 0
    t = np.linspace(0, 4, 9)
    model.set_timepoints(t)
    assert model.get_free_parameter_ids() == ("k", "p")
    assert model.get_state_ids() == ("x", "y")

    # analytical solution: x' = -k * p * x until x = 0.3, then x' = -k * x
    k, p = 0.5, 3.0
    t_switch = np.log(1 / 0.3) / (k * p)
    assert t[0] < t_switch < t[-1]
    before = t < t_switch
    x = np.where(before, np.exp(-k * p * t), 0.3 * np.exp(-k * (t - t_switch)))
    dxdk = np.where(before, -p * t * x, -t * x)
    dxdp = np.where(before, -k * t * x, -k * t_switch * x / p)
    y = np.where(before, p * x, x)
    dydk = np.where(before, p * dxdk, dxdk)
    dydp = np.where(before, x + p * dxdp, dxdp)

    solver = model.create_solver()
    solver.set_sensitivity_order(SensitivityOrder.first)
    solver.set_sensitivity_method(SensitivityMethod.forward)
    solver.set_relative_tolerance(1e-12)
    solver.set_absolute_tolerance(1e-12)
    solver.set_relative_tolerance_fsa(1e-10)
    solver.set_absolute_tolerance_fsa(1e-10)
    rdata = run_simulation(model, solver)
    assert rdata.status == AMICI_SUCCESS
    assert_allclose(rdata.x, np.stack([x, y], axis=1), rtol=1e-8)
    sx_expected = np.stack(
        [np.stack([dxdk, dydk], axis=1), np.stack([dxdp, dydp], axis=1)],
        axis=1,
    )
    assert_allclose(rdata.sx, sx_expected, rtol=1e-6, atol=1e-8)


@skip_on_valgrind
@pytest.mark.parametrize(
    "trigger_state",
    ["differential", "algebraic"],
)
def test_dae_algebraic_rule_piecewise_adjoint(tempdir, trigger_state):
    """Adjoint sensitivities of a DAE with a state-dependent `piecewise` in
    an algebraic rule, triggered by a differential or an algebraic state."""
    if trigger_state == "differential":
        equations = """
            x = 1; z = 3
            x' = -k * z
            0 = z - piecewise(p * x, x > 0.3, x)
        """
    else:
        equations = """
            x = 1; y = 3; z = 3
            x' = -k * z
            0 = y - p * x
            0 = z - piecewise(y, y > 0.9, x)
        """
    ant_model = f"""
    model dae_algebraic_rule_piecewise_asa
        const k = 0.5, p = 3
        {equations}
    end
    """
    module_name = f"dae_algebraic_rule_piecewise_asa_{trigger_state}"
    antimony2amici(
        ant_model,
        model_name=module_name,
        output_dir=tempdir,
        observation_model=[MC("obs_x", formula="x")],
    )
    model = import_model_module(module_name, tempdir).get_model()
    assert model.ne > 0
    t = np.linspace(0, 4, 9)
    model.set_timepoints(t)
    assert model.get_free_parameter_ids() == ("k", "p")

    # analytical solution: x' = -k * p * x until x = x_switch,
    # then x' = -k * x
    k, p = 0.5, 3.0
    # x_switch and d log(x_switch) / dp
    x_switch, dlogx_switch_dp = (
        (0.3, 0.0) if trigger_state == "differential" else (0.9 / p, -1 / p)
    )
    t_switch = -np.log(x_switch) / (k * p)
    assert t[0] < t_switch < t[-1]
    before = t < t_switch
    x = np.where(
        before, np.exp(-k * p * t), x_switch * np.exp(-k * (t - t_switch))
    )
    dxdk = np.where(before, -p * t * x, -t * x)
    dxdp = np.where(
        before,
        -k * t * x,
        x * (dlogx_switch_dp * (1 - 1 / p) - k * t_switch / p),
    )

    measurements = np.full(len(t), 0.1)
    edata = ExpData(model)
    edata.set_measurements(measurements)
    edata.set_noise_scales(np.ones(len(t)))
    sllh_expected = -np.array(
        [np.sum((x - measurements) * dxdk), np.sum((x - measurements) * dxdp)]
    )
    assert np.all(np.abs(sllh_expected) > 0.05)

    for method in (SensitivityMethod.forward, SensitivityMethod.adjoint):
        solver = model.create_solver()
        solver.set_sensitivity_order(SensitivityOrder.first)
        solver.set_sensitivity_method(method)
        solver.set_relative_tolerance(1e-12)
        solver.set_absolute_tolerance(1e-12)
        rdata = run_simulation(model, solver, edata)
        assert rdata.status == AMICI_SUCCESS
        assert_allclose(
            rdata.sllh, sllh_expected, rtol=1e-6, err_msg=method.name
        )


@skip_on_valgrind
@pytest.mark.parametrize("from_trigger", [True, False])
def test_dae_simultaneous_events_adjoint(tempdir, from_trigger):
    """Adjoint sensitivities of a DAE with simultaneous events, triggered by
    an algebraic state, with assignments depending on an algebraic state
    whose equation switches at the same time."""
    module_name = f"dae_simultaneous_events_asa_{from_trigger}"
    antimony2amici(
        f"""
        model {module_name}
            const k = 0.5, p = 3
            x = 1; y = 3; z = 3; w = 0; v = 0
            x' = -k * z
            w' = 0
            v' = 0
            0 = y - p * x
            0 = z - piecewise(y, y > 0.9, 2 * y - 0.9)
            E1: at y < 0.9, priority=2, fromTrigger={str(from_trigger).lower()}:
                w = w + z
            E2: at y < 0.9, priority=1, fromTrigger={str(from_trigger).lower()}:
                v = v + 2 * z
        end
        """,
        model_name=module_name,
        output_dir=tempdir,
        observation_model=[
            MC("obs_x", formula="x"),
            MC("obs_w", formula="w"),
            MC("obs_v", formula="v"),
        ],
    )
    model = import_model_module(module_name, tempdir).get_model()
    t = np.linspace(0, 4, 9)
    model.set_timepoints(t)
    edata = ExpData(model)
    edata.set_measurements(np.full(len(t) * 3, 0.1))
    edata.set_noise_scales(np.ones(len(t) * 3))

    rdatas = {}
    for method in (SensitivityMethod.forward, SensitivityMethod.adjoint):
        solver = model.create_solver()
        solver.set_sensitivity_order(SensitivityOrder.first)
        solver.set_sensitivity_method(method)
        solver.set_relative_tolerance(1e-12)
        solver.set_absolute_tolerance(1e-12)
        rdatas[method] = run_simulation(model, solver, edata)
        assert rdatas[method].status == AMICI_SUCCESS

    rdata = rdatas[SensitivityMethod.forward]
    # both events fire once
    assert rdata.by_id("w")[-1] == pytest.approx(0.9)
    assert rdata.by_id("v")[-1] == pytest.approx(1.8)
    assert np.all(np.abs(rdata.sllh) > 0.05)
    assert_allclose(
        rdatas[SensitivityMethod.adjoint].sllh, rdata.sllh, rtol=1e-6
    )


@skip_on_valgrind
def test_dae_event_adjoint_sensitivities(tempdir):
    """Adjoint sensitivities of a DAE with an event that changes a
    differential state the algebraic state depends on."""
    ant_model = """
    model dae_event_asa
        const k = 0.5, p = 2, q = 0.1
        x = 1; y = 2
        x' = -k * x + q
        0 = y - p * x
        at (time > 1.3): x = x + 1
    end
    """
    module_name = "dae_event_asa"
    antimony2amici(
        ant_model,
        model_name=module_name,
        output_dir=tempdir,
        observation_model=[MC("obs_x", formula="x")],
    )
    model = import_model_module(module_name, tempdir).get_model()
    t = np.linspace(0.5, 4.5, 9)
    model.set_timepoints(t)
    edata = ExpData(model)
    edata.set_measurements(np.full(len(t), 1.0))
    edata.set_noise_scales(np.ones(len(t)))

    sllh = {}
    for method in (SensitivityMethod.forward, SensitivityMethod.adjoint):
        solver = model.create_solver()
        solver.set_sensitivity_order(SensitivityOrder.first)
        solver.set_sensitivity_method(method)
        solver.set_relative_tolerance(1e-12)
        solver.set_absolute_tolerance(1e-12)
        rdata = run_simulation(model, solver, edata)
        assert rdata.status == AMICI_SUCCESS, method
        sllh[method] = rdata.sllh
    assert np.all(np.abs(sllh[SensitivityMethod.forward][[0, 2]]) > 0.1)
    assert_allclose(
        sllh[SensitivityMethod.adjoint],
        sllh[SensitivityMethod.forward],
        rtol=1e-6,
        atol=1e-8,
    )


@skip_on_valgrind
def test_explicit_trigger_times_within_roundoff(tempdir):
    """Explicit trigger times within roundoff of each other or of the time
    the solver stopped at (gh-3393).

    The solver cannot step to a stop time within roundoff of the current
    time, so such events are handled at the current time.
    """
    model_name = "test_explicit_trigger_times_within_roundoff"
    antimony2amici(
        r"""
        p1 = 0.1
        p2 = 0.2
        p3 = 0.3
        s = 0
        x = 1
        x' = s - x
        E1: at time >= p1 + p2: s = s + 1
        E2: at time >= p3: s = s + 10
        """,
        observation_model=[MC("obs_x", formula="x")],
        model_name=model_name,
        output_dir=tempdir,
    )
    model = import_model_module(model_name, tempdir).get_model()
    # E1 triggers 1 ulp after 0.3
    assert 0.1 + 0.2 == np.nextafter(0.3, 1)

    def expected(t, t_e2):
        """Analytical s(t), x(t) and dx(t)/dp, for E1 at 0.3 and E2 at
        t_e2 >= 0.3."""
        s, x, t_prev, dx_dp = 0.0, 1.0, 0.0, np.zeros(3)
        for t_e, ds, dt_e_dp in ((0.3, 1, [1, 1, 0]), (t_e2, 10, [0, 0, 1])):
            if t_e > t:
                break
            x = s + (x - s) * np.exp(-(t_e - t_prev))
            t_prev, s = t_e, s + ds
            # the jump in x' = s - x shifts x by -ds * dt_e
            dx_dp -= ds * np.exp(-(t - t_e)) * np.array(dt_e_dp)
        return s, s + (x - s) * np.exp(-(t - t_prev)), dx_dp

    def check(timepoints, t_e2):
        model.set_timepoints(timepoints)
        edata = ExpData(model)
        # measure only at the last timepoint, where x is differentiable
        # w.r.t. the trigger times
        measurements = np.full(len(timepoints), np.nan)
        measurements[-1] = 5.0
        edata.set_measurements(measurements)
        edata.set_noise_scales(np.ones(len(timepoints)))
        s_expected, x_expected, dx_dp = zip(
            *(expected(t, t_e2) for t in timepoints), strict=True
        )
        sllh_expected = -(x_expected[-1] - 5.0) * dx_dp[-1]

        for sensi_method in (
            SensitivityMethod.none,
            SensitivityMethod.forward,
            SensitivityMethod.adjoint,
        ):
            solver = model.create_solver()
            solver.set_sensitivity_order(SensitivityOrder.first)
            solver.set_sensitivity_method(sensi_method)
            rdata = run_simulation(model, solver, edata)
            assert rdata.status == AMICI_SUCCESS, sensi_method
            assert_allclose(rdata.by_id("s"), s_expected, rtol=0, atol=0)
            assert_allclose(rdata.by_id("x"), x_expected, rtol=1e-6)
            if sensi_method != SensitivityMethod.none:
                assert_allclose(
                    rdata.sllh, sllh_expected, rtol=1e-5, err_msg=sensi_method
                )

    # (a) two trigger times 1 ulp apart
    model.set_free_parameter_by_id("p3", 0.3)
    check([0.0, 1.0], t_e2=0.3)

    # (b) trigger time 1 ulp after an output timepoint: E1 is handled
    #  before the output at t=0.3 is recorded
    model.set_free_parameter_by_id("p3", 0.7)
    check([0.0, 0.3, 1.0], t_e2=0.7)

    # (c) trigger time 1 ulp after the initial time
    model.set_t0(0.3)
    model.set_timepoints([0.3, 1.0])
    rdata = run_simulation(model, model.create_solver())
    assert rdata.status == AMICI_SUCCESS
    assert_allclose(rdata.by_id("s"), [1.0, 11.0], rtol=0, atol=0)
    model.set_t0(0.0)

    # (d) post-equilibration, two trigger times 1 ulp apart
    model.set_free_parameter_by_id("p3", 0.3)
    model.set_steady_state_computation_mode(
        SteadyStateComputationMode.integrationOnly
    )
    model.set_timepoints([0.0, np.inf])
    rdata = run_simulation(model, model.create_solver())
    assert rdata.status == AMICI_SUCCESS
    assert_allclose(rdata.by_id("s"), [0.0, 11.0], rtol=0, atol=0)
    assert_allclose(rdata.by_id("x")[-1], 11.0, rtol=1e-4)

    # (e) post-equilibration after a trigger time 1 ulp after the last
    #  finite output timepoint: E1 is handled once, during
    #  post-equilibration
    model.set_free_parameter_by_id("p3", 0.7)
    model.set_timepoints([0.0, 0.3, np.inf])
    rdata = run_simulation(model, model.create_solver())
    assert rdata.status == AMICI_SUCCESS
    assert_allclose(rdata.by_id("s"), [0.0, 0.0, 11.0], rtol=0, atol=0)
    assert_allclose(rdata.by_id("x")[-1], 11.0, rtol=1e-4)

    # (f) as (e), but with the trigger time 1 ulp after the initial time.
    #  After E1, x = s = 1 is a steady state, so post-equilibration ends
    #  before E2.
    model.set_t0(0.3)
    model.set_timepoints([0.3, np.inf])
    rdata = run_simulation(model, model.create_solver())
    assert rdata.status == AMICI_SUCCESS
    assert_allclose(rdata.by_id("s"), [0.0, 1.0], rtol=0, atol=0)
    assert_allclose(rdata.by_id("x"), [1.0, 1.0], rtol=0, atol=0)


@skip_on_valgrind
def test_explicit_trigger_times_handled_before_reached(tempdir):
    """Events whose explicit trigger times are within roundoff after the
    current time are handled before these times are reached (gh-3393).

    They fire even if their root function is already zero at the current
    time, and also if an event assignment moves another event's trigger time
    to such a time.
    """
    model_name = "test_explicit_trigger_times_handled_before_reached"
    antimony2amici(
        r"""
        pa = 1
        k = 1
        c = 1
        tc = 100
        s = 0
        x = 0
        x' = 1
        Ea: at time >= pa: s = s + 1
        Eb: at time / k >= c: s = s + 10, tc = c * k
        Ec: at time >= tc: s = s + 100
        """,
        model_name=model_name,
        output_dir=tempdir,
    )
    model = import_model_module(model_name, tempdir).get_model()
    # Eb triggers at c * k, 1 ulp after Ea, and its root function
    # `time / k - c` is zero at the trigger time of Ea. Eb moves the trigger
    # time of Ec to its own trigger time. The parameter values are set here,
    # as they would be rounded in the SBML model.
    pa, k, c = 1.7099052299459676, 3.3554874394253624, 0.5095847505955183
    assert c * k == np.nextafter(pa, np.inf)
    assert pa / k - c == 0.0
    model.set_free_parameter_by_id("pa", pa)
    model.set_free_parameter_by_id("k", k)
    model.set_free_parameter_by_id("c", c)

    model.set_timepoints([0.0, 6.0])
    rdata = run_simulation(model, model.create_solver())
    assert rdata.status == AMICI_SUCCESS
    assert_allclose(rdata.by_id("s"), [0.0, 111.0], rtol=0, atol=0)


@skip_on_valgrind
def test_secondary_event_heaviside_at_trigger(tempdir):
    """Assignments of events triggered by other events see the updated
    Heaviside variables (gh-3429)."""
    model_name = "test_secondary_event_heaviside_at_trigger"
    antimony2amici(
        r"""
        s = 0
        z1 = 0
        z2 = 0
        x = 1
        u = 0
        w := piecewise(1, s >= 1, 0)
        x' = -0.1 * x
        u' = w
        A: at time >= 1, priority=2: s = 1
        B1: at s >= 1, priority=1, fromTrigger=true: z1 = w
        B2: at s >= 1, priority=1, fromTrigger=false: z2 = w
        """,
        model_name=model_name,
        output_dir=tempdir,
    )
    model = import_model_module(model_name, tempdir).get_model()
    model.set_timepoints([0.0, 2.0])
    rdata = run_simulation(model, model.create_solver())
    assert rdata.status == AMICI_SUCCESS
    assert_allclose(rdata.by_id("s"), [0.0, 1.0], rtol=0, atol=0)
    assert_allclose(rdata.by_id("z1"), [0.0, 1.0], rtol=0, atol=0)
    assert_allclose(rdata.by_id("z2"), [0.0, 1.0], rtol=0, atol=0)


@skip_on_valgrind
def test_event_crosses_own_trigger_back(tempdir):
    """The Heaviside variables of roots that an event assignment crosses back
    in the same discontinuity are updated (gh-3425)."""
    model_name = "test_event_crosses_own_trigger_back"
    antimony2amici(
        r"""
        k1 = 0.7; k2 = 0.3; a = 1.5; thr = 0.5
        x = 1; y = 0.5
        x' = piecewise(-k1 * x, x > thr, -2 * k1 * x)
        y' = piecewise(-k2 * y, x > thr, -k2 * x * y)
        E1: at x < thr, fromTrigger=false: x = x * a
        """,
        model_name=model_name,
        output_dir=tempdir,
        observation_model=[
            MC("obs_x", formula="x"),
            MC("obs_y", formula="y"),
        ],
    )
    model = import_model_module(model_name, tempdir).get_model()
    k1, k2, a, thr = 0.7, 0.3, 1.5, 0.5
    t = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    model.set_timepoints(t)

    # `x` reaches `thr` at `t0`, `E1` moves it back above `thr`, from where it
    # decays to `thr` again within `period`. So, `x > thr` holds throughout,
    # and `x = a^m * exp(-k1 * t)` after `m` event occurrences.
    t0 = np.log(1 / thr) / k1
    period = np.log(a) / k1
    m = np.where(t < t0, 0, np.floor((t - t0) / period) + 1)
    assert np.sum(m >= 3) >= 2
    x = a**m * np.exp(-k1 * t)
    y = 0.5 * np.exp(-k2 * t)
    # The event times depend on `thr`, but the state jumps compensate the
    # shifted right-hand side discontinuities, so `dx/dthr = 0`.
    sx = {"k1": -t * x, "k2": 0 * t, "a": m * x / a, "thr": 0 * t}
    sy = {"k1": 0 * t, "k2": -t * y, "a": 0 * t, "thr": 0 * t}
    par_ids = list(model.get_free_parameter_ids())
    ix = model.get_state_ids().index("x")
    iy = model.get_state_ids().index("y")

    edata = ExpData(model)
    measurements = np.vstack([x, y]).T + 0.1
    edata.set_measurements(measurements.flatten())
    edata.set_noise_scales(np.ones_like(measurements).flatten())
    for sens_method in (SensitivityMethod.forward, SensitivityMethod.adjoint):
        solver = model.create_solver()
        solver.set_sensitivity_order(SensitivityOrder.first)
        solver.set_sensitivity_method(sens_method)
        solver.set_relative_tolerance(1e-12)
        solver.set_absolute_tolerance(1e-12)
        rdata = run_simulation(model, solver, edata=edata)
        assert rdata.status == AMICI_SUCCESS
        assert_allclose(rdata.by_id("x"), x, rtol=1e-8)
        assert_allclose(rdata.by_id("y"), y, rtol=1e-8)
        # the right-hand side at the final time point uses the `x > thr`
        # branch
        assert_allclose(rdata.xdot[ix], -k1 * x[-1], rtol=1e-8)

        expected_sllh = [
            np.sum(0.1 * (sx[par_id] + sy[par_id])) for par_id in par_ids
        ]
        assert_allclose(
            rdata.sllh,
            expected_sllh,
            rtol=1e-6,
            atol=1e-8,
            err_msg=str(sens_method),
        )
        if sens_method == SensitivityMethod.forward:
            for ip, par_id in enumerate(par_ids):
                assert_allclose(
                    rdata.sx[:, ip, ix], sx[par_id], rtol=1e-6, atol=1e-8
                )
                assert_allclose(
                    rdata.sx[:, ip, iy], sy[par_id], rtol=1e-6, atol=1e-8
                )


@skip_on_valgrind
def test_event_moves_own_trigger_time(tempdir):
    """An event that moves its own explicit trigger time fires again
    (gh-3425)."""
    model_name = "test_event_moves_own_trigger_time"
    antimony2amici(
        r"""
        s = 1; n = 0
        x = 1; x' = -x
        E: at time >= s: s = s + 1, n = n + 1
        """,
        model_name=model_name,
        output_dir=tempdir,
    )
    model = import_model_module(model_name, tempdir).get_model()
    assert model.ne_solver == 0
    model.set_timepoints([0.0, 1.5, 2.5, 3.5])
    rdata = run_simulation(model, model.create_solver())
    assert rdata.status == AMICI_SUCCESS
    assert_allclose(rdata.by_id("n"), [0.0, 1.0, 2.0, 3.0], rtol=0, atol=0)


@skip_on_valgrind
@pytest.mark.parametrize("x_new", ["xmax", "0.3"])
def test_event_assignment_onto_threshold(tempdir, x_new):
    """An event assignment that sets a just-crossed root function to zero,
    exactly or up to rounding errors, does not cross it back.

    `xmax` is computed at runtime, `0.1 + 0.2 = 0.30000000000000004`, so
    `x = 0.3` is one rounding error below it.
    """
    model_name = f"test_event_assignment_onto_threshold_{x_new != 'xmax'}"
    antimony2amici(
        rf"""
        r = 1; a = 0.1; b = 0.2
        xmax := a + b
        x = 0; x' = piecewise(r, x < xmax, 0)
        E: at x > xmax: x = {x_new}
        """,
        model_name=model_name,
        output_dir=tempdir,
    )
    model = import_model_module(model_name, tempdir).get_model()
    model.set_timepoints([0.0, 1.0, 2.0])
    rdata = run_simulation(model, model.create_solver())
    assert rdata.status == AMICI_SUCCESS
    x_expected = 0.1 + 0.2 if x_new == "xmax" else 0.3
    assert_allclose(rdata.by_id("x"), [0.0, x_expected, x_expected], rtol=1e-8)
