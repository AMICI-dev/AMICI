#!/usr/bin/env python3
"""
Run SBML Test Suite and verify simulation results
[https://github.com/sbmlteam/sbml-test-suite/releases]

Usage:
    pytest tests.sbml.testSBMLSuite -n CORES --dist=loadgroup --cases=SELECTION
        CORES can be an integer or `auto` for all available cores.
        `--dist=loadgroup` is required whenever `-n` is used: each case's
        simulation and sensitivity checks share one compiled model
        (`compiled_case` fixture) and must run on the same xdist worker.
        SELECTION can be e.g.: `1`, `1,3`, `-3,4,6-7`, or `100-` to select
        specific test cases. If `--cases` is omitted, all cases are run.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import amici
import diffrax
import jax
import jax.numpy as jnp
import numpy as np
import optimistix
import pandas as pd
import pytest
from amici.sim.jax.petab import (
    DEFAULT_CONTROLLER_SETTINGS,
    DEFAULT_ROOT_FINDER_SETTINGS,
)
from amici.sim.sundials import (
    AMICI_SUCCESS,
    Model,
    SensitivityMethod,
    SensitivityOrder,
    Solver,
    run_simulation,
)
from amici.sim.sundials.gradient_check import check_derivatives
from utils import (
    apply_settings,
    find_model_file,
    read_settings_file,
    verify_results,
    write_result_file,
)

# test cases for which sensitivities are to be checked
#  key: case ID; value: epsilon for finite differences
_SENSITIVITY_CHECK_CASES = {
    # parameter-dependent conservation laws
    "00783": 1.5e-2,
    # initial events
    "00995": 1e-3,
}


@pytest.fixture(scope="session")
def compiled_case(test_id, sbml_semantic_cases_dir):
    """Compile and configure a case's model, shared by its simulation and
    sensitivity test nodes.

    Model import/compilation is expensive and must happen at most once per
    `test_id` -- see the `xdist_group`/`scope="session"` comment in
    `pytest_generate_tests` (conftest.py) for why that requires both
    consuming test nodes to run on the same xdist worker.
    """
    current_test_path = sbml_semantic_cases_dir / test_id
    model_dir = Path(__file__).parent / "SBMLTestModels" / test_id
    try:
        model, solver, wrapper = compile_model(
            current_test_path,
            test_id,
            model_dir,
            generate_sensitivity_code=True,
        )
    except amici.importers.sbml.SBMLException as err:
        # `pytest.skip` from inside a fixture correctly propagates as
        # SKIPPED to every dependent test, not as a fixture ERROR.
        pytest.skip(str(err))

    settings = read_settings_file(current_test_path, test_id)
    atol, rtol = apply_settings(settings, solver, model, test_id)
    solver.set_sensitivity_order(SensitivityOrder.first)
    solver.set_sensitivity_method(SensitivityMethod.forward)
    if test_id == "00885":
        # 00885: root-after-reinitialization with FSA with default settings
        solver.set_absolute_tolerance(1e-16)
        solver.set_relative_tolerance(1e-15)

    yield model, solver, wrapper, settings, atol, rtol, current_test_path

    shutil.rmtree(model_dir, ignore_errors=True)


def _check_simulation_status(rdata, test_id: str) -> None:
    """Skip/fail consistently for a known-bad vs. a genuinely unexpected
    base simulation failure."""
    if rdata["status"] != AMICI_SUCCESS:
        if test_id in ("00748", "00374", "00369"):
            pytest.skip("Simulation Failed expectedly")
        raise RuntimeError("Simulation failed unexpectedly")


def test_sbml_testsuite_case(test_id, compiled_case, result_path):
    model, solver, wrapper, settings, atol, rtol, current_test_path = (
        compiled_case
    )

    # parse expected results
    results_file = current_test_path / f"{test_id}-results.csv"
    results = pd.read_csv(results_file, delimiter=",")
    results.rename(
        columns={c: c.replace(" ", "") for c in results.columns},
        inplace=True,
    )

    # simulate model
    rdata = run_simulation(model, solver)
    _check_simulation_status(rdata, test_id)

    # verify
    simulated = verify_results(
        settings, rdata, results, wrapper, model, atol, rtol
    )

    # record results
    write_result_file(simulated, test_id, result_path)


def test_sbml_testsuite_case_sensitivity(test_id, compiled_case):
    """Check the model's sensitivities for the selected cases in
    `_SENSITIVITY_CHECK_CASES`, via finite differences and against JAX
    autodiff.

    Reported independently from `test_sbml_testsuite_case`: whether AMICI
    can correctly *simulate* an SBML feature and whether its *sensitivities*
    for that feature are correct are different questions, and a case whose
    simulation is right but whose sensitivities aren't (yet) checked
    shouldn't count against basic SBML simulation support.
    """
    model, solver, wrapper, settings, atol, rtol, current_test_path = (
        compiled_case
    )

    epsilon = _SENSITIVITY_CHECK_CASES.get(test_id)
    if epsilon is None:
        pytest.skip("No sensitivity check defined for this case")

    rdata = run_simulation(model, solver)
    _check_simulation_status(rdata, test_id)

    check_derivatives(model, solver=solver, epsilon=epsilon)
    jax_sensitivity_check(
        current_test_path,
        test_id,
        model,
        rdata,
        atol,
        rtol,
    )


def compile_model(
    sbml_dir: Path,
    test_id: str,
    model_dir: Path,
    generate_sensitivity_code: bool = False,
) -> tuple[Model, Solver, amici.SbmlImporter]:
    """Import the given test model to AMICI"""
    model_dir.mkdir(parents=True, exist_ok=True)

    sbml_file = find_model_file(sbml_dir, test_id)
    sbml_importer = amici.SbmlImporter(sbml_file)

    model_name = f"SBMLTest{test_id}"
    sbml_importer.sbml2amici(
        model_name,
        output_dir=model_dir,
        generate_sensitivity_code=generate_sensitivity_code,
    )

    # settings
    model_module = amici.import_model_module(model_name, model_dir)

    model = model_module.get_model()
    solver = model.create_solver()

    return model, solver, sbml_importer


def compile_model_jax(sbml_dir: Path, test_id: str, model_dir: Path):
    """Import the given test model as JAX model"""
    model_dir.mkdir(parents=True, exist_ok=True)
    sbml_file = find_model_file(sbml_dir, test_id)
    sbml_importer = amici.SbmlImporter(sbml_file)
    model_name = f"SBMLTest{test_id}_jax"
    sbml_importer.sbml2jax(model_name, output_dir=model_dir)
    model_module = amici.import_model_module(model_dir.name, model_dir.parent)
    jax_model = model_module.Model()
    return jax_model, sbml_importer


def jax_sensitivity_check(
    sbml_dir: Path,
    test_id: str,
    amici_model: Model,
    rdata: dict,
    atol: float,
    rtol: float,
):
    """Compare AMICI forward sensitivities against JAX autodiff"""
    model_dir = Path(__file__).parent / "SBMLTestModelsJaxGrad" / test_id
    try:
        jax_model, _ = compile_model_jax(sbml_dir, test_id, model_dir)
    except NotImplementedError as err:
        if "The JAX backend does not support" in str(err):
            pytest.skip(str(err))
        raise

    try:
        ts = rdata["ts"]
        p = jax_model.parameters
        ts_jnp = jnp.asarray(ts, dtype=float)
        zeros = jnp.zeros_like(ts_jnp)
        tol_factor = 1e2
        if int(test_id) in (
            191,
            192,
            193,
            194,
            198,
            199,
            201,
            270,
            272,
            273,
            274,
            276,
            277,
            279,
            1148,
            1159,
            1160,
            1161,
            1395,
        ):
            tol_factor = 1e4

        solver = diffrax.Kvaerno5()
        controller = diffrax.PIDController(
            rtol=rtol / tol_factor,
            atol=atol / tol_factor,
            pcoeff=DEFAULT_CONTROLLER_SETTINGS["pcoeff"],
            icoeff=DEFAULT_CONTROLLER_SETTINGS["icoeff"],
            dcoeff=DEFAULT_CONTROLLER_SETTINGS["dcoeff"],
        )
        root_finder = optimistix.Newton(**DEFAULT_ROOT_FINDER_SETTINGS)

        def simulate(pars):
            x, _ = jax_model.simulate_condition(
                pars,
                ts_jnp,
                jnp.array([]),
                zeros,
                jnp.zeros_like(ts_jnp, dtype=int),
                jnp.zeros_like(ts_jnp, dtype=int),
                jnp.zeros((ts_jnp.shape[0], 0)),
                jnp.zeros((ts_jnp.shape[0], 0)),
                solver,
                controller,
                root_finder,
                diffrax.DirectAdjoint(),
                diffrax.SteadyStateEvent(),
                2**10,
                ret=amici.sim.jax.ReturnValue.x,
            )
            return x

        x = simulate(p)
        sx = jax.jacfwd(simulate)(p)
        par_idx = [
            jax_model.parameter_ids.index(pid)
            for pid in amici_model.get_free_parameter_ids()
        ]
        sx = jnp.transpose(sx[:, :, par_idx], (0, 2, 1))

        if rdata["sx"] is None:
            solver_amici = amici_model.create_solver()
            solver_amici.set_sensitivity_order(SensitivityOrder.first)
            solver_amici.set_sensitivity_method(SensitivityMethod.forward)
            rdata = run_simulation(amici_model, solver_amici)

        np.testing.assert_allclose(x, rdata["x"], rtol=rtol, atol=atol)
        np.testing.assert_allclose(sx, rdata["sx"], rtol=rtol, atol=atol)
    finally:
        shutil.rmtree(model_dir, ignore_errors=True)
