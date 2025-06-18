import copy
from pathlib import Path

import pytest
import petab.v1 as petab
from petab.v1.lint import measurement_table_has_timepoint_specific_mappings

import benchmark_models_petab
from amici.petab.petab_import import import_petab_problem

from test_petab_benchmark import problems

script_dir = Path(__file__).parent
repo_root = script_dir.parent.parent
benchmark_outdir = repo_root / "test_bmc"


@pytest.fixture(scope="session", params=problems, ids=problems)
def benchmark_problem(request):
    """Fixture providing model and PEtab problem for a benchmark model."""
    problem_id = request.param
    petab_problem = benchmark_models_petab.get_problem(problem_id)
    flat_petab_problem = copy.deepcopy(petab_problem)
    if measurement_table_has_timepoint_specific_mappings(
        petab_problem.measurement_df,
    ):
        petab.flatten_timepoint_specific_output_overrides(flat_petab_problem)

    amici_model = import_petab_problem(
        flat_petab_problem,
        model_output_dir=benchmark_outdir / problem_id,
    )
    return problem_id, flat_petab_problem, petab_problem, amici_model
