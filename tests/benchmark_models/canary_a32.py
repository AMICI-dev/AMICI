"""Diagnostic canary for AMICI-dev/AMICI#3078.

Evaluates the plain (no-sensitivity) simulated log-likelihood of the
Weber_BMC2015 PEtab benchmark model at a fixed point, and the finite
-difference central estimate for the `a32` direction at the two step sizes
known to produce catastrophically-cancelled (and hence anomaly-prone)
estimates -- the actual quantity that flips between "fine" and failing the
real test's tolerance. Each is repeated multiple times within this one
process and printed at full precision (`repr()` and `.hex()`, for
unambiguous bit-level comparison across runs).

This script does not itself assert anything -- it is a data-collection aid,
run unconditionally (pass or fail) alongside the real benchmark tests, to
build up a dataset of (runner hardware, `AMICI_EXTRACT_CSE` setting,
computed values) correlations without needing to catch an actual test
failure. Repeats within one process additionally distinguish per-process
variability (would point to a runtime race) from per-job/VM-only
variability (consistent with genuine hardware/microcode floating-point
differences). Assumes the Weber_BMC2015 model has already been built by the
main test run (uses the default output directory and does not force
recompilation).
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
DIRECTION = "a32"
# the two step sizes (out of the test's default 8) that produce wildly
# inaccurate finite differences for `a32`, due to catastrophic cancellation
ANOMALOUS_STEP_SIZES = (1e-4, 1e-5)


def fmt(x: float) -> str:
    return f"{x!r} ({x.hex()})"


def main():
    petab_problem = benchmark_models_petab.get_problem(PROBLEM_ID)
    if measurement_table_has_timepoint_specific_mappings(
        petab_problem.measurement_df
    ):
        petab.flatten_timepoint_specific_output_overrides(petab_problem)

    parameter_ids = list(petab_problem.x_free_ids)
    a32_idx = parameter_ids.index(DIRECTION)
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
    h0 = float(point[a32_idx])

    print(f"canary: AMICI_EXTRACT_CSE={os.environ.get('AMICI_EXTRACT_CSE')!r}")

    # plain nominal-point log-likelihood, repeated within-process
    for i in range(N_REPEATS):
        llh = float(amici_function(point))
        print(f"canary: nominal run={i} llh={fmt(llh)}")

    # central-difference estimate for `a32`, at the two anomalous step
    # sizes, repeated within-process
    for size in ANOMALOUS_STEP_SIZES:
        h = size * h0
        p_plus = point.copy()
        p_plus[a32_idx] = h0 + h
        p_minus = point.copy()
        p_minus[a32_idx] = h0 - h
        for i in range(N_REPEATS):
            f_plus = float(amici_function(p_plus))
            f_minus = float(amici_function(p_minus))
            central = (f_plus - f_minus) / (2 * h)
            print(
                f"canary: fd size={size} run={i} "
                f"f+={fmt(f_plus)} f-={fmt(f_minus)} central={fmt(central)}"
            )


if __name__ == "__main__":
    main()
