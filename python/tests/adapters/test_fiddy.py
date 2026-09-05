"""Tests for `amici.adapters.fiddy`."""

from pathlib import Path

import numpy as np
import pytest
from amici.adapters.fiddy import (
    run_simulation_to_function_and_derivative,
    simulate_petab_to_function_and_derivative,
)
from amici.importers.petab.v1 import import_petab_problem
from amici.sim.sundials import SensitivityOrder, SteadyStateSensitivityMode
from fiddy import Type, check_gradient, check_jacobian
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


def test_run_amici_simulation_to_function_and_derivative():
    petab_problem, point = lotka_volterra()
    timepoints = sorted(set(petab_problem.measurement_df.time))
    amici_model = import_petab_problem(petab_problem)
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
    function, derivative = run_simulation_to_function_and_derivative(
        free_parameter_ids=parameter_ids,
        amici_model=amici_model,
        amici_solver=amici_solver,
        derivative_variables=["x", "x0", "y", "sigmay"],
    )

    expected = derivative(point)
    result = check_jacobian(function, point, expected)
    result.assert_success(always_print=True)


@pytest.mark.parametrize("scaled_parameters", (False, True))
def test_simulate_petab_to_function_and_derivative(scaled_parameters):
    petab_problem, point = lotka_volterra()
    amici_model = import_petab_problem(petab_problem)
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
