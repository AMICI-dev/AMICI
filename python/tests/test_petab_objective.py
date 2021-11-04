"""Tests for petab_objective.py."""

from pathlib import Path

import amici
import amici.petab_import
import amici.petab_objective
import petab
import pytest


# Absolute and relative tolerances for finite difference gradient checks.
ATOL: float = 1e-3
RTOL: float = 1e-3


@pytest.fixture
def lotka_volterra() -> petab.Problem:
    return petab.Problem.from_yaml(str(
        Path('petab_test_problems')
        / 'lotka_volterra'
        / 'petab'
        / 'problem.yaml'
    ))


def test_simulate_petab_sensitivities(lotka_volterra):
    petab_problem = lotka_volterra
    amici_model = amici.petab_import.import_petab_problem(petab_problem)
    amici_solver = amici_model.getSolver()

    amici_solver.setSensitivityOrder(amici.SensitivityOrder_first)

    result = amici.petab_objective.check_grad_multi_eps(
        petab_problem=petab_problem,
        amici_model=amici_model,
        amici_solver=amici_solver,
        detailed=True,
    )

    assert (result.abs_err < ATOL).all()
    assert (result.rel_err < RTOL).all()
