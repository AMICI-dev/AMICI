#!/usr/bin/env python3
"""Run SBML Test Suite with JAX and verify simulation results."""

import shutil
from pathlib import Path

import amici
import jax
import jax.numpy as jnp
import numpy as np
import optimistix
import pytest
from amici.sim.sundials import AMICI_SUCCESS
from utils import (
    CaseSettings,
    compile_model_jax,
    read_expected_results,
    simulate_jax,
    verify_results,
    write_result_file,
)

jax.config.update("jax_enable_x64", True)


@pytest.fixture(scope="session")
def result_path_jax() -> Path:
    return Path(__file__).parent / "amici-semantic-results-jax"


def run_jax_simulation(model, ts, atol, rtol, tol_factor=1e2):
    p = model.parameters
    ts_jnp = jnp.asarray(ts, dtype=float)
    x, stats = simulate_jax(
        model,
        p,
        ts_jnp,
        atol,
        rtol,
        tol_factor,
        root_finder=optimistix.Newton(atol=atol, rtol=rtol),
    )
    y = jax.vmap(
        lambda t, x_solver, x_rdata, hs: model._y(
            t,
            x_solver,
            p,
            model._tcl(x_rdata, p),
            hs,
            jnp.zeros(len(model.observable_ids)),
        )
    )(ts_jnp, stats["x"], x, stats["hs"])
    w = jax.vmap(
        lambda t, x_solver, x_rdata, hs: model._w(
            t, x_solver, p, model._tcl(x_rdata, p), hs
        )
    )(ts_jnp, stats["x"], x, stats["hs"])

    class RData(dict):
        __getattr__ = dict.__getitem__

    return RData(
        ts=np.asarray(ts_jnp).copy(),
        y=np.asarray(y).copy(),
        w=np.asarray(w).copy(),
        status=AMICI_SUCCESS,
    )


def test_sbml_testsuite_case_jax(
    test_id, result_path_jax, sbml_semantic_cases_dir
):
    model_dir = Path(__file__).parent / "SBMLTestModelsJax" / test_id
    try:
        case_dir = sbml_semantic_cases_dir / test_id
        model, importer = compile_model_jax(case_dir, test_id, model_dir)
        settings = CaseSettings.from_file(case_dir, test_id)

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
            356,
            357,
            752,
            1148,
            1159,
            1160,
            1161,
            1395,
        ):
            tol_factor = 1e4
        rdata = run_jax_simulation(
            model,
            settings.timepoints,
            settings.atol,
            settings.rtol,
            tol_factor,
        )
        simulated = verify_results(
            settings,
            rdata,
            read_expected_results(case_dir, test_id),
            importer,
            parameters=dict(
                zip(
                    model.parameter_ids,
                    map(float, model.parameters),
                    strict=True,
                )
            ),
            expression_ids=list(model.expression_ids),
        )
        write_result_file(simulated, test_id, result_path_jax)
    except amici.importers.sbml.SBMLException as err:
        pytest.skip(str(err))
    except NotImplementedError as err:
        if "The JAX backend does not support" in str(err):
            pytest.skip(str(err))
        raise err
    finally:
        shutil.rmtree(model_dir, ignore_errors=True)
