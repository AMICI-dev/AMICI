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

import fiddy
from fiddy.extensions.amici import (
    simulate_petab_to_cached_functions,
)
from fiddy.gradient_check import keep_lowest_error


# Absolute and relative tolerances for finite difference gradient checks.
ATOL: float = 1e-3
RTOL: float = 1e-3

benchmark_path = Path(__file__).parent.parent.parent.parent / "Benchmark-Models-PEtab" / "Benchmark-Models"
benchmark_yamls = [
    petab_path / (petab_path.stem + ".yaml")
    for petab_path in benchmark_path.glob("*") if petab_path.is_dir()
]
#benchmark_yamls = [benchmark_yamls[0]]  # FIXME remove

debug = True
if debug:
    debug_path = Path('debug')
    debug_full_path = debug_path / 'full'
    debug_minimal_path = debug_path / 'minimal'
    debug_full_path.mkdir(exist_ok=True, parents=True)
    debug_minimal_path.mkdir(exist_ok=True, parents=True)


@pytest.mark.parametrize("petab_yaml", benchmark_yamls)
def test_benchmark_gradient(petab_yaml):
    petab_problem = petab.Problem.from_yaml(petab_yaml)

    sizes = [
        #1e-1,
        #1e-2,
        1e-3,
        #1e-4,
        1e-5,
        #1e-6,
        1e-7,
        #1e-8,
        1e-9,
        #1e-10,
        #1e-11,
        #1e-12,
        #1e-13,
    ]


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
        fiddy.gradient_check,
        function=function,
        point=point,
        gradient=gradient,
        sizes=sizes,
    )

    success, result_df = gradient_check_partial(
        fd_gradient_method='central',
    )

    if debug:
        result_df.to_csv(debug_full_path / (petab_yaml.stem + ".tsv"), sep='\t')
        minimal_result_df = keep_lowest_error(result_df, inplace=False)
        minimal_result_df.to_csv(debug_minimal_path / (petab_yaml.stem + ".tsv"), sep='\t')

    # The gradients for all parameters are correct.
    assert success, result_df
