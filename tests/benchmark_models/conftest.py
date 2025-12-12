import copy
import sys
from pathlib import Path

import benchmark_models_petab
import petab.v1 as petab
import pytest
from amici import get_model_root_dir
from amici.importers.petab.v1 import import_petab_problem
from petab.v1.lint import measurement_table_has_timepoint_specific_mappings

script_dir = Path(__file__).parent.resolve()
if str(script_dir) not in sys.path:
    sys.path.insert(0, str(script_dir))

from test_petab_benchmark import problems

benchmark_outdir = get_model_root_dir() / "test_bmc"


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
        output_dir=benchmark_outdir / problem_id,
    )
    return problem_id, flat_petab_problem, petab_problem, amici_model
