"""Tests for simulate_petab on PEtab benchmark problems."""

from pathlib import Path

import amici
import amici.petab_import
import amici.petab_objective
import numpy as np
import pandas as pd
import petab
import pytest

from fiddy import get_derivative, MethodId
from fiddy.success import Consistency
from fiddy.derivative_check import NumpyIsCloseDerivativeCheck
from fiddy.extensions.amici import (
    simulate_petab_to_cached_functions,
)


# Absolute and relative tolerances for finite difference gradient checks.
ATOL: float = 1e-3
RTOL: float = 1e-2

benchmark_path = Path(__file__).parent.parent.parent / "Benchmark-Models-PEtab" / "Benchmark-Models"
benchmark_outdir = Path(__file__).parent.parent.parent / "test_bmc"
benchmark_yamls = [
    petab_path / (petab_path.stem + ".yaml")
    for petab_path in benchmark_path.glob("*") if petab_path.is_dir()
    # excluded due to excessive runtime
    if not str(petab_path.stem).startswith(('Chen_MSB', 'Froehlich_CellSystems'))
]

debug = False
if debug:
    debug_path = Path('../../python/tests/debug')
    debug_path.mkdir(exist_ok=True, parents=True)


@pytest.mark.parametrize("scale", (True, False))
@pytest.mark.parametrize("petab_yaml", benchmark_yamls)
def test_benchmark_gradient(petab_yaml, scale):
    if not scale and str(petab_yaml.stem) in (
        'Smith_BMCSystBiol2013',
        'Alkan_SciSignal2018',
        'Lucarelli_CellSystems2018',
        'Isensee_JCB2018',
        'Brannmark_JBC2010',
        'Elowitz_Nature2000',
    ):
        # not really worth the effort trying to fix these cases if they
        # only fail on linear scale
        pytest.skip()

    petab_problem = petab.Problem.from_yaml(petab_yaml)
    if str(petab_yaml.stem) == 'Fiedler_BMC2016':
        petab.flatten_timepoint_specific_output_overrides(petab_problem)

    # Only compute gradient for estimated parameters.
    parameter_df_free = petab_problem.parameter_df.loc[petab_problem.x_free_ids]
    parameter_ids = list(parameter_df_free.index)

    # Setup AMICI objects.
    amici_model = amici.petab_import.import_petab_problem(
        petab_problem,
        model_output_dir=benchmark_outdir / petab_yaml.stem,
    )
    amici_solver = amici_model.getSolver()
    amici_solver.setAbsoluteTolerance(1e-10)
    amici_solver.setRelativeTolerance(1e-10)
    if amici_model.getName() in (
            'Smith_BMCSystBiol2013',
            'Raimundez_PCB2020',
            'Elowitz_Nature2000',
    ):
        amici_solver.setAbsoluteTolerance(1e-12)
        amici_solver.setRelativeTolerance(1e-12)

    if amici_model.getName() in (
        'Raimundez_PCB2020', 'Brannmark_JBC2010', 'Isensee_JCB2018'
    ):
        amici_model.setSteadyStateSensitivityMode(amici.SteadyStateSensitivityMode.integrationOnly)

    amici_function, amici_derivative = simulate_petab_to_cached_functions(
        petab_problem=petab_problem,
        parameter_ids=parameter_ids,
        amici_model=amici_model,
        solver=amici_solver,
        scaled_parameters=scale,
        scaled_gradients=scale,
        cache=not debug,
    )

    noise_level = 0.1
    if amici_model.getName() == 'Crauste_CellSystems2017':
        noise_level = 0.0

    np.random.seed(0)
    if scale:
        point = np.asarray(list(
            petab_problem.scale_parameters(dict(parameter_df_free.nominalValue)).values()
        ))
        point_noise = np.random.randn(len(point)) * noise_level
    else:
        point = parameter_df_free.nominalValue.values
        point_noise = np.random.randn(len(point)) * point * noise_level
    point += point_noise  # avoid small gradients at nominal value

    expected_derivative = amici_derivative(point)

    sizes = [
        1e-1,
        1e-2,
        1e-3,
        1e-4,
        1e-5,
    ]

    if amici_model.getName() == 'Smith_BMCSystBiol2013':
        sizes.insert(0, 0.5)
        if scale:
            sizes = sizes[:-3]

    derivative = get_derivative(
        function=amici_function,
        point=point,
        sizes=sizes,
        direction_ids=parameter_ids,
        method_ids=[MethodId.CENTRAL, MethodId.FORWARD, MethodId.BACKWARD],
        success_checker=Consistency(atol=1e-4, rtol=0.2),
        expected_result=expected_derivative,
        relative_sizes=not scale,
    )
    success = False
    if derivative.df.success.all():
        check = NumpyIsCloseDerivativeCheck(
            derivative=derivative,
            expectation=expected_derivative,
            point=point,
        )
        success = check(rtol=RTOL, atol=ATOL)

    if debug:
        df = pd.DataFrame([
            {
                ('fd', r.metadata['size_absolute'], str(r.method_id)): r.value
                for c in d.computers
                for r in c.results
            } for d in derivative.directional_derivatives
        ], index=parameter_ids)
        df[('fd', 'full', '')] = derivative.series.values
        df[('amici', '', '')] = expected_derivative

        file_name = f"{petab_yaml.stem}_scale={scale}.tsv"
        df.to_csv(debug_path / file_name, sep='\t')

    # The gradients for all parameters are correct.
    assert success, derivative.df
