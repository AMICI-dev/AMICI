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
from fiddy import get_derivative, MethodId
from fiddy.success import Consistency
from fiddy.derivative_check import NumpyIsCloseDerivativeCheck
from fiddy.extensions.amici import (
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
#benchmark_yamls = [benchmark_yamls[0]]  # FIXME remove

debug = True
if debug:
    debug_path = Path('debug')
    debug_full_path = debug_path / 'full'
    debug_minimal_path = debug_path / 'minimal'
    debug_full_path.mkdir(exist_ok=True, parents=True)
    debug_minimal_path.mkdir(exist_ok=True, parents=True)


@pytest.mark.parametrize("scale", (False, True))
@pytest.mark.parametrize("petab_yaml", benchmark_yamls)
def test_benchmark_gradient(petab_yaml, scale):
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
        #1e-9,
        #1e-10,
        #1e-11,
        #1e-12,
        #1e-13,
    ]

    # Only compute gradient for estimated parameters.
    parameter_df_free = petab_problem.parameter_df.loc[petab_problem.x_free_ids]
    parameter_ids = list(parameter_df_free.index)
    parameter_scales = dict(parameter_df_free.parameterScale)
    # Set point to midpoint of bounds.
    # Hopefully no gradients are zero at this point...
    point = ((parameter_df_free[LOWER_BOUND] + parameter_df_free[UPPER_BOUND])/2).values

    # Setup AMICI objects.
    amici_model = amici.petab_import.import_petab_problem(petab_problem)
    amici_solver = amici_model.getSolver()
    amici_solver.setSensitivityOrder(amici.SensitivityOrder_first)
    amici_solver.setAbsoluteTolerance(1e-12)
    amici_solver.setRelativeTolerance(1e-12)

    #if amici_model.getName() == 'Bachmann_MSB2011':
    #    amici_solver.setMaxSteps(int(1e5))









    amici_function, amici_derivative = simulate_petab_to_cached_functions(
        parameter_ids=parameter_ids,
        petab_problem=petab_problem,
        amici_model=amici_model,
        solver=amici_solver,
        scaled_parameters=scale,
        scaled_gradients=scale,
        # TODO check if caching speeds up this test
        #      with either disk or ram caching
        cache=False,
    )

    expected_derivative = amici_derivative(point)
    if (expected_derivative == 0).any():
        raise ValueError(
            "The chosen point may be an issue: an expected gradient value is zero."
        )

    derivative = get_derivative(
        function=amici_function,
        point=point,
        sizes=sizes,
        direction_ids=parameter_ids,
        method_ids=[MethodId.FORWARD, MethodId.BACKWARD, MethodId.CENTRAL],
        #analysis_classes=[
        #    lambda: TransformByDirectionScale(scales=parameter_scales),
        #],
        success_checker=Consistency(),
    )
    test_value = derivative.value

    check = NumpyIsCloseDerivativeCheck(
        derivative=derivative,
        expectation=expected_derivative,
        point=point,
    )
    result = check(rtol=1e-2)
    assert result.success
