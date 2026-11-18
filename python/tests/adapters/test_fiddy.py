"""Tests for `amici.adapters.fiddy`."""

import sys
from pathlib import Path

import amici
import numpy as np
import pytest
from amici.adapters.fiddy import (
    output_labels_for_derivatives,
    run_simulation_to_function_and_derivative,
    simulate_petab_to_function_and_derivative,
)
from amici.importers.petab.v1 import import_petab_problem
from amici.sim.sundials import (
    ExpData,
    ParameterScaling,
    SensitivityOrder,
    SteadyStateSensitivityMode,
    run_simulation,
)
from fiddy import (
    FunctionEvaluationError,
    JoblibExecutor,
    SequentialExecutor,
    Type,
    check_gradient,
    check_jacobian,
    estimate_gradient,
)
from petab import v1


def lotka_volterra() -> tuple[v1.Problem, np.ndarray]:
    petab_problem = v1.Problem.from_yaml(
        str(
            Path(__file__).parents[1]
            / "petab_test_problems"
            / "lotka_volterra"
            / "petab"
            / "problem.yaml"
        )
    )
    point = np.array([2, 3], dtype=Type.SCALAR)
    return petab_problem, point


@pytest.fixture(scope="session")
def lotka_volterra_model_module():
    """Imports `lotka_volterra` model module."""
    petab_problem, _ = lotka_volterra()
    import_petab_problem(petab_problem)
    model_name = petab_problem.model.model_id
    return amici.import_model_module(
        model_name, amici.get_model_dir(model_name)
    )


def test_run_amici_simulation_to_function_and_derivative(
    lotka_volterra_model_module,
):
    petab_problem, point = lotka_volterra()
    timepoints = sorted(set(petab_problem.measurement_df.time))
    amici_model = lotka_volterra_model_module.get_model()
    amici_model.set_timepoints(timepoints)
    amici_solver = amici_model.create_solver()

    amici_solver.set_sensitivity_order(SensitivityOrder.first)

    parameter_ids = list(
        petab_problem.parameter_df[
            petab_problem.parameter_df.estimate == 1
        ].index
    )

    derivative_variables = ["x", "x0", "y", "sigmay"]
    function, derivative = run_simulation_to_function_and_derivative(
        free_parameter_ids=parameter_ids,
        amici_model=amici_model,
        amici_solver=amici_solver,
        derivative_variables=derivative_variables,
    )

    expected = derivative(point)
    output_labels = output_labels_for_derivatives(
        amici_model,
        derivative_variables=derivative_variables,
        timepoints=timepoints,
    )
    result = check_jacobian(
        function,
        point,
        expected,
        direction_labels=parameter_ids,
        output_labels=output_labels,
    )
    assert len(output_labels) == len(result.output_results)
    result.assert_success(always_print=True)


def test_run_simulation_to_function_and_derivative_raises_on_failure(
    lotka_volterra_model_module,
):
    """A failed simulation raises `FunctionEvaluationError`."""
    petab_problem, point = lotka_volterra()
    amici_model = lotka_volterra_model_module.get_model()
    amici_model.set_timepoints(sorted(set(petab_problem.measurement_df.time)))
    amici_solver = amici_model.create_solver()
    amici_solver.set_max_steps(1)

    function, derivative = run_simulation_to_function_and_derivative(
        free_parameter_ids=["alpha", "gamma"],
        amici_model=amici_model,
        amici_solver=amici_solver,
        derivative_variables=["llh"],
    )
    with pytest.raises(FunctionEvaluationError, match="status"):
        function(point)
    with pytest.raises(FunctionEvaluationError, match="status"):
        derivative(point)


def test_run_simulation_respects_a_customized_parameter_list(
    lotka_volterra_model_module,
):
    """Regression test for proper `plist` handling."""
    petab_problem, _ = lotka_volterra()
    timepoints = sorted(set(petab_problem.measurement_df.time))
    amici_model = lotka_volterra_model_module.get_model()
    amici_model.set_timepoints(timepoints)
    amici_solver = amici_model.create_solver()
    amici_solver.set_sensitivity_order(SensitivityOrder.first)

    free_parameter_ids = list(amici_model.get_free_parameter_ids())
    alpha_id = "alpha"
    sigma_id = "noiseParameter1_observable_prey"
    assert set(free_parameter_ids) == {alpha_id, "gamma", sigma_id}
    alpha_index = free_parameter_ids.index(alpha_id)
    sigma_index = free_parameter_ids.index(sigma_id)
    # Confirm the natural order actually puts alpha before sigma, so the
    # `[sigma_index, alpha_index]` plist below is really a reversal.
    assert alpha_index < sigma_index

    values = {alpha_id: 2.0, "gamma": 3.0, sigma_id: 1.0}
    full_point = np.array([values[pid] for pid in free_parameter_ids])

    # Baseline: default (identity) plist
    _, baseline_derivative = run_simulation_to_function_and_derivative(
        free_parameter_ids=free_parameter_ids,
        amici_model=amici_model,
        amici_solver=amici_solver,
        derivative_variables=["llh"],
    )
    baseline = baseline_derivative(full_point)["llh"]

    # A subset, reversed relative to the model's own natural order.
    amici_model.set_parameter_list([sigma_index, alpha_index])
    subset_ids = [alpha_id, sigma_id]
    subset_point = np.array([full_point[alpha_index], full_point[sigma_index]])
    _, subset_derivative = run_simulation_to_function_and_derivative(
        free_parameter_ids=subset_ids,
        amici_model=amici_model,
        amici_solver=amici_solver,
        derivative_variables=["llh"],
    )
    restricted = subset_derivative(subset_point)["llh"]

    np.testing.assert_allclose(restricted[0], baseline[alpha_index])
    np.testing.assert_allclose(restricted[1], baseline[sigma_index])


@pytest.mark.skipif(
    sys.platform == "win32",
    reason="Parallelization/pickling requires HDF5 support -- unavailable on Windows builds.",
)
def test_joblib_executor_agrees_with_sequential_executor(
    lotka_volterra_model_module,
):
    """Results from `SequentialExecutor` and `JoblibExecutor`
    must agree exactly.
    """
    petab_problem, point = lotka_volterra()
    timepoints = sorted(set(petab_problem.measurement_df.time))
    amici_model = lotka_volterra_model_module.get_model()
    amici_model.set_timepoints(timepoints)
    amici_solver = amici_model.create_solver()
    amici_solver.set_sensitivity_order(SensitivityOrder.first)

    parameter_ids = list(
        petab_problem.parameter_df[
            petab_problem.parameter_df.estimate == 1
        ].index
    )

    function, _ = run_simulation_to_function_and_derivative(
        free_parameter_ids=parameter_ids,
        amici_model=amici_model,
        amici_solver=amici_solver,
        derivative_variables=["x", "x0", "y", "sigmay"],
    )

    sequential = estimate_gradient(
        function, point, executor=SequentialExecutor()
    )
    parallel = estimate_gradient(
        function, point, executor=JoblibExecutor(n_jobs=4)
    )

    for s, p in zip(sequential, parallel, strict=True):
        assert s.value == p.value
        assert s.status == p.status


@pytest.mark.parametrize("scaled_parameters", (False, True))
def test_simulate_petab_to_function_and_derivative(
    scaled_parameters, lotka_volterra_model_module
):
    petab_problem, point = lotka_volterra()
    amici_model = lotka_volterra_model_module.get_model()
    amici_solver = amici_model.create_solver()

    if amici_model.get_name() == "simple":
        amici_model.set_steady_state_sensitivity_mode(
            SteadyStateSensitivityMode.integrationOnly
        )

    amici_solver.set_sensitivity_order(SensitivityOrder.first)

    if scaled_parameters:
        point = np.asarray(
            list(
                petab_problem.scale_parameters(
                    dict(
                        zip(
                            petab_problem.parameter_df.index,
                            point,
                            strict=True,
                        )
                    )
                ).values()
            )
        )

    function, derivative = simulate_petab_to_function_and_derivative(
        free_parameter_ids=petab_problem.parameter_df.index,
        petab_problem=petab_problem,
        amici_model=amici_model,
        solver=amici_solver,
        scaled_gradients=scaled_parameters,
        scaled_parameters=scaled_parameters,
    )

    expected = derivative(point)
    result = check_gradient(function, point, expected)
    result.assert_success(always_print=True)


def test_derivative_check_leaves_solver_and_model_unchanged(
    lotka_volterra_model_module,
):
    """The caller's solver and model are restored after each evaluation."""
    petab_problem, point = lotka_volterra()
    amici_model = lotka_volterra_model_module.get_model()
    amici_model.set_timepoints(sorted(set(petab_problem.measurement_df.time)))
    amici_solver = amici_model.create_solver()
    amici_solver.set_sensitivity_order(SensitivityOrder.first)

    parameter_ids = list(
        petab_problem.parameter_df[
            petab_problem.parameter_df.estimate == 1
        ].index
    )
    # evaluations at `point` must actually change the model's parameters
    amici_model.set_free_parameter_by_id(dict.fromkeys(parameter_ids, 1.0))
    assert not np.allclose(
        [amici_model.get_free_parameter_by_id(pid) for pid in parameter_ids],
        point,
    )
    original_parameters = amici_model.get_free_parameters()

    function, derivative = run_simulation_to_function_and_derivative(
        free_parameter_ids=parameter_ids,
        amici_model=amici_model,
        amici_solver=amici_solver,
        derivative_variables=["y"],
    )
    check_jacobian(function, point, derivative(point))
    assert amici_solver.get_sensitivity_order() == SensitivityOrder.first
    assert amici_model.get_free_parameters() == original_parameters

    function, derivative = simulate_petab_to_function_and_derivative(
        petab_problem=petab_problem,
        amici_model=amici_model,
        solver=amici_solver,
    )
    check_gradient(function, point, derivative(point))
    assert amici_solver.get_sensitivity_order() == SensitivityOrder.first
    assert amici_model.get_free_parameters() == original_parameters


def test_event_output_derivatives(tempdir):
    """The derivatives of event outputs are mapped to the right parameters.

    `E1` is time-triggered: AMICI's `sz` ignores the parameter dependence
    of state-dependent event times (#3361).
    """
    from amici.importers.antimony import antimony2amici
    from amici.importers.utils import MeasurementChannel as MC

    model_name = "test_fiddy_event_outputs"
    antimony2amici(
        """
        k = 0.5; q = 0.1; s = 0.3
        x = 1
        x' = -k * x + q
        E1: at time > 1.3: x = x + 1
        E2: at x > 5: x = 0
        """,
        model_name=model_name,
        output_dir=tempdir,
        observation_model=[
            MC("obs_x", formula="x", sigma="s"),
            MC("z1", formula="q * x", event_id="E1", sigma="2 * s"),
            # For events that did not occur, `sz` is 0, but `z` is the
            # output at the last time point: use a constant output.
            MC("z2", formula="1", event_id="E2", sigma="s * k"),
        ],
    )
    model = amici.import_model_module(model_name, tempdir).get_model()
    model.set_parameter_scale(ParameterScaling.none)
    model.set_timepoints([0.5, 2.0])
    model.set_n_max_event(1)
    edata = ExpData(model)
    edata.set_measurements(np.ones(2))
    edata.set_event_measurements(np.ones(model.nztrue))
    solver = model.create_solver()
    solver.set_relative_tolerance(1e-12)
    solver.set_absolute_tolerance(1e-15)
    solver.set_sensitivity_order(SensitivityOrder.first)

    rdata = run_simulation(model, solver, edata)
    # E1 triggers, E2 does not
    assert rdata.rz[0, 0] == 0
    assert rdata.rz[0, 1] != 0
    # nontrivial derivatives, w.r.t. different parameters
    assert np.count_nonzero(rdata.sz[0, :, 0]) == 2
    assert np.count_nonzero(rdata.srz[0, :, 1]) == 2
    assert np.count_nonzero(rdata.ssigmaz[0, :, 1]) == 2

    # not in the model's order, to check the reordering
    free_parameter_ids = ["s", "k", "q"]
    point = np.array(
        [model.get_free_parameter_by_id(p) for p in free_parameter_ids]
    )
    for variable in ("z", "rz", "sigmaz"):
        function, derivative = run_simulation_to_function_and_derivative(
            free_parameter_ids=free_parameter_ids,
            amici_model=model,
            amici_solver=solver,
            amici_edata=edata,
            derivative_variables=[variable],
        )
        expected = derivative(point)
        check_jacobian(
            function, point, expected, direction_labels=free_parameter_ids
        ).assert_success()
