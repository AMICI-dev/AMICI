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
from amici.constants import SymbolId


def _steady_state_event(rtol=None, atol=None, norm=None):
    """Replicate :func:`diffrax.steady_state_event` without using it."""

    def cond_fn(t, y, args, *, terms, solver, stepsize_controller, **kwargs):
        del kwargs
        _rtol = rtol if rtol is not None else stepsize_controller.rtol
        _atol = atol if atol is not None else stepsize_controller.atol
        _norm = norm if norm is not None else stepsize_controller.norm
        vf = solver.func(terms, t, y, args)
        return _norm(vf) < _atol + _rtol * _norm(y)

    return cond_fn


from tests.sbml.tests.testSBMLSuite import (
    verify_results,
    write_result_file,
    find_model_file,
    read_settings_file,
)
from tests.sbml.tests.conftest import format_test_id


@pytest.fixture(scope="session")
def result_path_jax() -> Path:
    return Path(__file__).parents[1] / "amici-semantic-results-jax"


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
        return list(map(str, self.jax_model.parameter_ids))

    def getParameterById(self, pid: str):
        par = self.importer.sbml.getParameter(pid)
        return par.getValue() if par else np.nan

    def getExpressionIds(self):
        exprs = list(map(str, self.importer.flux_ids))
        exprs += [
            str(k)
            for k in self.importer.symbols.get(SymbolId.EXPRESSION, {}).keys()
        ]
        return exprs


def compile_model_jax(sbml_dir: Path, test_id: str, model_dir: Path):
    model_dir.mkdir(parents=True, exist_ok=True)
    sbml_file = find_model_file(sbml_dir, test_id)
    sbml_importer = amici.SbmlImporter(sbml_file)
    model_name = f"SBMLTest{test_id}_jax"
    sbml_importer.sbml2jax(model_name, output_dir=model_dir)
    model_module = amici.import_model_module(model_dir.name, model_dir.parent)
    jax_model = model_module.Model()
    return jax_model, sbml_importer


def run_jax_simulation(model, importer, ts, atol, rtol):
    p = jnp.array(
        [
            importer.sbml.getParameter(pid).getValue()
            for pid in model.parameter_ids
        ]
    )
    ts_jnp = jnp.asarray(ts, dtype=float)
    zeros = jnp.zeros_like(ts_jnp)
    solver = diffrax.Kvaerno5()
    controller = diffrax.PIDController(rtol=rtol / 1e4, atol=atol / 1e4)
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
        _steady_state_event(),
        2**8,
        ret=amici.jax.ReturnValue.x,
    )
    tcl = model._tcl(x[0], p)
    y = jax.vmap(
        lambda t, xs: model._y(
            t, xs, p, tcl, jnp.zeros(len(model.observable_ids))
        )
    )(ts_jnp, stats["x"])
    w = jax.vmap(lambda t, xs: model._w(t, xs, p, tcl))(ts_jnp, stats["x"])

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

        rdata = run_jax_simulation(model, wrapper, ts, atol, rtol)
        dummy = DummyModel(model, wrapper)
        simulated = verify_results(
            settings, rdata, results, wrapper, dummy, atol, rtol
        )
        write_result_file(simulated, test_id, result_path_jax)
    except amici.sbml_import.SBMLException as err:
        pytest.skip(str(err))
    finally:
        shutil.rmtree(model_dir, ignore_errors=True)
