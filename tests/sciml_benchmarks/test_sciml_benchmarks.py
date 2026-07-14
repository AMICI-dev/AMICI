import os
from contextlib import contextmanager
from pathlib import Path

import jax
import numpy as np
import pytest
from amici.importers.petab import *
from amici.sim.jax import run_simulations
from petab import v2
from yaml import safe_load


@contextmanager
def change_directory(destination):
    # Save the current working directory
    original_directory = os.getcwd()
    try:
        # Change to the new directory
        os.chdir(destination)
        yield
    finally:
        # Change back to the original directory
        os.chdir(original_directory)


jax.config.update("jax_enable_x64", True)

cases_dir = Path(__file__).parent / "testsuite" / "Benchmark-Models"


@pytest.mark.parametrize(
    "test", sorted([d.name for d in cases_dir.iterdir() if d.is_dir()])
)
def test_benchmark(test):
    test_dir = cases_dir / test

    for split in ("train", "validation"):
        with open(test_dir / f"{split}_{test}.yaml") as f:
            petab_yaml = safe_load(f)
        with open(test_dir / f"expected_{test}.yaml") as f:
            solutions = safe_load(f)

        with change_directory(test_dir):
            petab_problem = v2.Problem.from_yaml(petab_yaml)

            pi = PetabImporter(
                petab_problem=petab_problem,
                module_name="hybrid" + test,
                compile_=True,
                jax=jax,
                validate=False,
            )

            jax_problem = pi.create_simulator(
                force_import=True,
            )

        # llh
        llh, _ = run_simulations(jax_problem)
        solution_split = "training" if split == "train" else "validation"
        np.testing.assert_allclose(
            llh,
            solutions["llh"][solution_split],
            atol=1e-3,
            rtol=1e-3,
        )
