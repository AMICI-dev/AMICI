#!/usr/bin/env python3
"""Run SBML Test Suite with JAX and verify simulation results."""

import copy
import os
import shutil
import sys
from pathlib import Path

import amici
import jax
import jax.numpy as jnp
import diffrax
import numpy as np
import pandas as pd
import pytest
from amici.jax.petab import DEFAULT_CONTROLLER_SETTINGS


from tests.sbml.testSBMLSuite import (
    verify_results,
    write_result_file,
    find_model_file,
    read_settings_file,
)
from tests.sbml.conftest import format_test_id

jax.config.update("jax_enable_x64", True)


@pytest.fixture(scope="session")
def result_path_jax() -> Path:
    return Path(__file__).parent / "amici-semantic-results-jax"


@pytest.fixture(scope="function", autouse=True)
def sbml_test_dir():
    old_cwd = os.getcwd()
    old_path = copy.copy(sys.path)
    yield
    os.chdir(old_cwd)
    sys.path = old_path


class DummyModel:
    def __init__(self, jax_model: object, importer: amici.SbmlImporter):
        self.jax_model = jax_model
        self.importer = importer

    def getParameterIds(self):
        return list(self.jax_model.parameter_ids)

    def getParameterById(self, pid: str):
        return float(
            self.jax_model.parameters[
                list(self.jax_model.parameter_ids).index(pid)
            ]
        )

    def getExpressionIds(self):
        return list(self.jax_model.expression_ids)


def compile_model_jax(sbml_dir: Path, test_id: str, model_dir: Path):
    model_dir.mkdir(parents=True, exist_ok=True)
    sbml_file = find_model_file(sbml_dir, test_id)
    sbml_importer = amici.SbmlImporter(sbml_file)
    model_name = f"SBMLTest{test_id}_jax"
    sbml_importer.sbml2jax(model_name, output_dir=model_dir)
    model_module = amici.import_model_module(model_dir.name, model_dir.parent)
    jax_model = model_module.Model()
    return jax_model, sbml_importer


def run_jax_simulation(model, importer, ts, atol, rtol, tol_factor=1e2):
    p = model.parameters
    ts_jnp = jnp.asarray(ts, dtype=float)
    zeros = jnp.zeros_like(ts_jnp)
    solver = diffrax.Kvaerno5()
    controller = diffrax.PIDController(
        rtol=rtol / tol_factor,
        atol=atol / tol_factor,
        pcoeff=DEFAULT_CONTROLLER_SETTINGS["pcoeff"],
        icoeff=DEFAULT_CONTROLLER_SETTINGS["icoeff"],
        dcoeff=DEFAULT_CONTROLLER_SETTINGS["dcoeff"],
    )
    x, stats = model.simulate_condition(
        p,
        ts_jnp,
        jnp.array([]),
        zeros,
        jnp.zeros_like(ts_jnp, dtype=int),
        jnp.zeros_like(ts_jnp, dtype=int),
        jnp.zeros((ts_jnp.shape[0], 0)),
        jnp.zeros((ts_jnp.shape[0], 0)),
        solver,
        controller,
        diffrax.DirectAdjoint(),
        diffrax.SteadyStateEvent(),
        2**10,
        ret=amici.jax.ReturnValue.x,
    )
    y = jax.vmap(
        lambda t, x_solver, x_rdata: model._y(
            t,
            x_solver,
            p,
            model._tcl(x_rdata, p),
            jnp.zeros(len(model.observable_ids)),
        )
    )(ts_jnp, stats["x"], x)
    w = jax.vmap(
        lambda t, x_solver, x_rdata: model._w(
            t, x_solver, p, model._tcl(x_rdata, p)
        )
    )(ts_jnp, stats["x"], x)

    class RData(dict):
        __getattr__ = dict.__getitem__

    return RData(
        ts=np.asarray(ts_jnp).copy(),
        y=np.asarray(y).copy(),
        w=np.asarray(w).copy(),
        status=amici.AMICI_SUCCESS,
    )


def test_sbml_testsuite_case_jax(
    test_number, result_path_jax, sbml_semantic_cases_dir
):
    test_id = format_test_id(test_number)
    model_dir = Path(__file__).parent / "SBMLTestModelsJax" / test_id
    try:
        current_test_path = sbml_semantic_cases_dir / test_id
        results_file = current_test_path / f"{test_id}-results.csv"
        results = pd.read_csv(results_file, delimiter=",")
        results.rename(
            columns={c: c.replace(" ", "") for c in results.columns},
            inplace=True,
        )

        model, wrapper = compile_model_jax(
            current_test_path, test_id, model_dir
        )
        settings = read_settings_file(current_test_path, test_id)
        ts = np.linspace(
            float(settings["start"] or 0),
            float(settings["start"] or 0) + float(settings["duration"] or 0),
            int(settings["steps"] or 0) + 1,
        )
        atol = float(settings["absolute"])
        rtol = float(settings["relative"])

        dummy = DummyModel(model, wrapper)

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
        rdata = run_jax_simulation(model, wrapper, ts, atol, rtol, tol_factor)
        simulated = verify_results(
            settings, rdata, results, wrapper, dummy, atol, rtol
        )
        write_result_file(simulated, test_id, result_path_jax)
    except amici.sbml_import.SBMLException as err:
        pytest.skip(str(err))
    except NotImplementedError as err:
        if "The JAX backend does not support" in str(err):
            pytest.skip(str(err))
        raise err
    finally:
        shutil.rmtree(model_dir, ignore_errors=True)
