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
from amici.sim.sundials import SensitivityOrder, SteadyStateSensitivityMode
from fiddy import (
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

    # `x_ss`/`llh`/`res` are excluded: this model has no steady state (a
    # pure oscillator, so `x_ss`/`sx_ss` are structurally undefined), and no
    # `amici_edata` is supplied here (this test is about plain-ReturnData
    # sensitivities, not PEtab-driven measurement fitting -- see
    # `test_simulate_petab_to_function_and_derivative` for the `llh`/`sllh`
    # case), so `llh`/`res` (which need measurements) are undefined too.
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
