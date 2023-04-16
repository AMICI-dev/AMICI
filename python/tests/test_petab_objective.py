"""Tests for petab_objective.py."""

from pathlib import Path

import amici
import amici.petab_import
import amici.petab_objective
import numpy as np
import petab
import pytest


# Absolute and relative tolerances for finite difference gradient checks.
ATOL: float = 1e-3
RTOL: float = 1e-3


@pytest.fixture
def lotka_volterra() -> petab.Problem:
    return petab.Problem.from_yaml(str(
        Path(__file__).parent
        / 'petab_test_problems'
        / 'lotka_volterra'
        / 'petab'
        / 'problem.yaml'
    ))


def test_simulate_petab_sensitivities(lotka_volterra):
    petab_problem = lotka_volterra
    amici_model = amici.petab_import.import_petab_problem(petab_problem)
    amici_solver = amici_model.getSolver()

    amici_solver.setSensitivityOrder(amici.SensitivityOrder_first)

    problem_parameters = dict(zip(
        petab_problem.x_ids,
        petab_problem.x_nominal,
    ))

    results = {}
    for scaled_parameters in [True, False]:
        for scaled_gradients in [True, False]:
            results[(scaled_parameters, scaled_gradients)] = (
                amici.petab_objective.check_grad_multi_eps(
                    petab_problem=petab_problem,
                    amici_model=amici_model,
                    amici_solver=amici_solver,
                    problem_parameters=(
                        petab_problem.scale_parameters(problem_parameters)
                        if scaled_parameters
                        else problem_parameters
                    ),
                    scaled_parameters=scaled_parameters,
                    scaled_gradients=scaled_gradients,
                )
            )

    # Each combination of `scaled_parameters` and `scaled_gradients` passes the
    # gradient check.
    for result_id, result in results.items():
        assert (result.rel_err < ATOL).all(), result_id
        assert (result.rel_err < RTOL).all(), result_id

    # `scaled_parameters` does not affect the output gradients, just
    # `scaled_gradients` in the gradient check (which affects `unscaled_gradients`
    # in `simulate_petab`)
    for scaled_gradients in [True, False]:
        assert np.isclose(
            results[(True, scaled_gradients)].grad,
            results[(False, scaled_gradients)].grad,
        ).all()

    # The gradient is transformed as expected.
    transformation_factor = np.array(petab_problem.x_nominal) * np.log(10)
    for scaled_parameters in [True, False]:
        assert np.isclose(
            results[(scaled_parameters, True)].grad,
            results[(scaled_parameters, False)].grad * transformation_factor,
        ).all()
