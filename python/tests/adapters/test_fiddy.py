"""Tests for `amici.adapters.fiddy`."""

from pathlib import Path

import amici
import amici.petab.petab_import
import amici.petab.simulations
import numpy as np
import pytest
from amici import SensitivityOrder, SteadyStateSensitivityMode
from amici.adapters.fiddy import (
    run_simulation_to_cached_functions,
    simulate_petab_to_cached_functions,
)
from fiddy import MethodId, Type, get_derivative
from fiddy.derivative_check import NumpyIsCloseDerivativeCheck
from fiddy.success import Consistency
from petab import v1

# Absolute and relative tolerances for finite difference gradient checks.
ATOL: float = 1e-3
RTOL: float = 1e-3


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


@pytest.mark.parametrize("problem_generator", [lotka_volterra])
def test_run_amici_simulation_to_functions(problem_generator):
    petab_problem, point = problem_generator()
    timepoints = sorted(set(petab_problem.measurement_df.time))
    amici_model = amici.petab.petab_import.import_petab_problem(petab_problem)
    amici_model.set_timepoints(timepoints)
    amici_solver = amici_model.create_solver()

    amici_solver.set_sensitivity_order(SensitivityOrder.first)

    parameter_ids = list(
        petab_problem.parameter_df[
            petab_problem.parameter_df.estimate == 1
        ].index
    )
    parameter_indices = [
        amici_model.get_parameter_ids().index(parameter_id)
        for parameter_id in parameter_ids
    ]

    (
        amici_function,
        amici_derivative,
        structures,
    ) = run_simulation_to_cached_functions(
        parameter_ids=parameter_ids,
        amici_model=amici_model,
        amici_solver=amici_solver,
    )

    expected_derivative = amici_derivative(point)[..., parameter_indices]

    derivative = get_derivative(
        function=amici_function,
        point=point,
        sizes=[1e-10, 1e-5],
        direction_ids=parameter_ids,
        method_ids=[MethodId.FORWARD, MethodId.BACKWARD, MethodId.CENTRAL],
        # analysis_classes=[],
        # analysis_classes=[
        #    lambda: TransformByDirectionScale(scales=parameter_scales),
        # ],
        success_checker=Consistency(atol=1e-2),
    )
    test_derivative = derivative.value

    # The test derivative is close to the expected derivative.
    assert np.isclose(
        test_derivative,
        expected_derivative,
        rtol=1e-1,
        atol=1e-1,
        equal_nan=True,
    ).all()

    # Same as above assert.
    check = NumpyIsCloseDerivativeCheck(
        derivative=derivative,
        expectation=expected_derivative,
        point=point,
    )
    result = check(rtol=1e-1, atol=1e-1, equal_nan=True)
    assert result.success


@pytest.mark.parametrize("problem_generator", [lotka_volterra])
@pytest.mark.parametrize("scaled_parameters", (False, True))
def test_simulate_petab_to_functions(problem_generator, scaled_parameters):
    petab_problem, point = problem_generator()
    amici_model = amici.petab.petab_import.import_petab_problem(petab_problem)
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

    amici_function, amici_derivative = simulate_petab_to_cached_functions(
        parameter_ids=petab_problem.parameter_df.index,
        petab_problem=petab_problem,
        amici_model=amici_model,
        solver=amici_solver,
        scaled_gradients=scaled_parameters,
        scaled_parameters=scaled_parameters,
    )

    expected_derivative = amici_derivative(point)

    parameter_ids = list(
        petab_problem.parameter_df[
            petab_problem.parameter_df.estimate == 1
        ].index
    )
    # parameter_scales = dict(
    #     petab_problem.parameter_df[
    #         petab_problem.parameter_df.estimate == 1
    #     ].parameterScale
    # )

    derivative = get_derivative(
        function=amici_function,
        point=point,
        sizes=[1e-10, 1e-5, 1e-3, 1e-1],
        direction_ids=parameter_ids,
        method_ids=[MethodId.FORWARD, MethodId.BACKWARD, MethodId.CENTRAL],
        success_checker=Consistency(),
    )

    check = NumpyIsCloseDerivativeCheck(
        derivative=derivative,
        expectation=expected_derivative,
        point=point,
    )
    result = check(rtol=1e-2)
    assert result.success
