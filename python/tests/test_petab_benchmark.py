"""Tests for simulate_petab on PEtab benchmark problems."""

import copy  # FIXME remove
from functools import partial
from pathlib import Path

import amici
import amici.petab_import
import amici.petab_objective
import numpy as np
import pandas as pd
import petab
from petab.C import NOMINAL_VALUE, LOWER_BOUND, UPPER_BOUND, ESTIMATE
import pytest

import finite_difference_methods as fdm
from finite_difference_methods.extensions.amici import (
    simulate_petab_to_cached_functions,
)


# Absolute and relative tolerances for finite difference gradient checks.
ATOL: float = 1e-3
RTOL: float = 1e-3

benchmark_path = Path(__file__).parent.parent.parent.parent / "Benchmark-Models-PEtab" / "Benchmark-Models"
benchmark_yamls = [
    petab_path / (petab_path.stem + ".yaml")
    for petab_path in benchmark_path.glob("*") if petab_path.is_dir()
]
#benchmark_yamls = [benchmark_yamls[1]]  # FIXME remove


@pytest.mark.parametrize("petab_yaml", benchmark_yamls)
def test_benchmark_gradient(petab_yaml):
    petab_problem = petab.Problem.from_yaml(petab_yaml)

    # Only compute gradient for estimated parameters.
    index_estimated = petab_problem.parameter_df[ESTIMATE] == 1
    df = petab_problem.parameter_df.loc[index_estimated]
    parameter_ids = df.index
    # Set point to midpoint of bounds.
    # Hopefully no gradients are zero at this point...
    point = ((df[LOWER_BOUND] + df[UPPER_BOUND])/2).values

    # Setup AMICI objects.
    amici_model = amici.petab_import.import_petab_problem(petab_problem)
    amici_solver = amici_model.getSolver()
    amici_solver.setSensitivityOrder(amici.SensitivityOrder_first)

    function, gradient = simulate_petab_to_cached_functions(
        simulate_petab=amici.petab_objective.simulate_petab,
        parameter_ids=parameter_ids,
        petab_problem=petab_problem,
        amici_model=amici_model,
        solver=amici_solver,
        # TODO check if caching speeds up this test
        #      with either disk or ram caching
        cache=False,
    )

    if (gradient(point) == 0).any():
        raise ValueError(
            "The chosen point may be an issue: an expected gradient is zero."
        )

    gradient_check_partial = partial(
        fdm.gradient_check,
        function=function,
        point=point,
        gradient=gradient,
        sizes=[
            1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11,
            1e-12, 1e-13,
        ],
    )

    success, result_df = gradient_check_partial(
        fd_gradient_method='central',
    )

    # The gradients for all parameters are correct.
    assert success, result_df
