"""Diagnostic canary for AMICI-dev/AMICI#3078.

Evaluates the plain (no-sensitivity) simulated log-likelihood of the
Weber_BMC2015 PEtab benchmark model at a fixed point, multiple times within
this one process, and prints each result at full precision.

This model/point/seed reproduces the `a32`-direction finite-difference
instability behind #3078's intermittent CI failures. This script does not
itself assert anything -- it is a data-collection aid, run unconditionally
(pass or fail) alongside the real benchmark tests, to build up a dataset of
(runner hardware, `AMICI_EXTRACT_CSE` setting, computed value) correlations
without needing to catch an actual test failure. Assumes the Weber_BMC2015
model has already been built by the main test run (uses the default
output directory and does not force recompilation).
"""

import os

import benchmark_models_petab
import numpy as np
import petab.v1 as petab
from amici import get_model_root_dir
from amici.adapters.fiddy import simulate_petab_to_cached_functions
from amici.importers.petab.v1 import import_petab_problem
from amici.sim.sundials import SensitivityMethod, SteadyStateSensitivityMode
from petab.v1.lint import measurement_table_has_timepoint_specific_mappings

N_REPEATS = 3
PROBLEM_ID = "Weber_BMC2015"
RNG_SEED = 2
NOISE_LEVEL = 0.05


def main():
    petab_problem = benchmark_models_petab.get_problem(PROBLEM_ID)
    if measurement_table_has_timepoint_specific_mappings(
        petab_problem.measurement_df
    ):
        petab.flatten_timepoint_specific_output_overrides(petab_problem)

    parameter_ids = petab_problem.x_free_ids
    output_dir = get_model_root_dir() / "test_bmc" / PROBLEM_ID

    amici_model = import_petab_problem(petab_problem, output_dir=output_dir)
    amici_solver = amici_model.create_solver()
    amici_solver.set_absolute_tolerance(1e-12)
    amici_solver.set_relative_tolerance(1e-12)
    amici_solver.set_max_steps(2 * 10**5)
    amici_solver.set_sensitivity_method(SensitivityMethod.forward)
    amici_model.set_steady_state_sensitivity_mode(
        SteadyStateSensitivityMode.integrateIfNewtonFails
    )

    amici_function, _ = simulate_petab_to_cached_functions(
        petab_problem=petab_problem,
        free_parameter_ids=parameter_ids,
        amici_model=amici_model,
        solver=amici_solver,
        scaled_parameters=False,
        scaled_gradients=False,
        cache=False,
        num_threads=os.cpu_count(),
    )

    np.random.seed(RNG_SEED)
    point = petab_problem.x_nominal_free
    point = point + np.random.randn(len(point)) * point * NOISE_LEVEL

    print(f"canary: AMICI_EXTRACT_CSE={os.environ.get('AMICI_EXTRACT_CSE')!r}")
    for i in range(N_REPEATS):
        llh = float(amici_function(point))
        print(f"canary: run={i} llh={llh!r}")


if __name__ == "__main__":
    main()
