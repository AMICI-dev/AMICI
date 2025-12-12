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
from amici.importers.antimony import antimony2amici
from amici.sim.sundials import (
    AMICI_SUCCESS,
    ExpData,
    SensitivityMethod,
    SensitivityOrder,
    SteadyStateComputationMode,
    SteadyStateSensitivityMode,
    run_simulation,
)
from amici.sim.sundials.gradient_check import check_derivatives
from amici.testing import skip_on_valgrind
from amici.testing.models import create_amici_model, create_sbml_model
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
        # FIXME adjoint sensitivities w.r.t. decay_2 are slightly off
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

    # FIXME: For a few parameters of these models, adjoint sensitivities
    # are somewhat off. This needs to be investigated further.
    asa_xfail = amici_model.get_name() in (
        "events_plus_heavisides",
        "piecewise_plus_event_semi_complicated",
        "nested_events",
    )
    check_trajectories_with_adjoint_sensitivities(amici_model, asa_xfail)


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
    # watch out for too long path names on windows ...
    module_name = "tst_mltple_ea_w_cmprtmnt"
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
    # FIXME: sensitivities w.r.t. the bolus parameter are not correct
    model.set_parameter_list(
        [
            ip
            for ip, par in enumerate(model.get_free_parameter_ids())
            if par != "two"
        ]
    )

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
            # FIXME: adjoint sensitivities w.r.t. the bolus parameter `three`
            #  are wrong.
            #  maybe related to https://github.com/AMICI-dev/AMICI/issues/2805
            model.set_parameter_list(
                [
                    ip
                    for ip, par in enumerate(model.get_free_parameter_ids())
                    if par not in ["one", "three"]
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
    # FIXME: this is currently not supported
    #  (root-after-reinitialization when switching from pre-simulation#
    #   to main simulation)
    # edata = ExpData(rdata, 1, 0, 0)
    # edata.fixedParametersPresimulation = [1]
    # edata.t_presim = 10
    # rdata = runAmiciSimulation(model, solver, edata=edata)
    # assert rdata.status == AMICI_SUCCESS
    # assert list(rdata.by_id("target1")) == [1.0, 1.0, 1.0]
    # assert list(rdata.by_id("target2")) == [1.0, 1.0, 1.0]

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
