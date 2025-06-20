#!/usr/bin/env python3
"""
Run SBML Test Suite and verify simulation results
[https://github.com/sbmlteam/sbml-test-suite/releases]

Usage:
    pytest tests.sbml.testSBMLSuite -n CORES --cases=SELECTION
        CORES can be an integer or `auto` for all available cores.
        SELECTION can be e.g.: `1`, `1,3`, `-3,4,6-7`, or `100-` to select
        specific test cases. If `--cases` is omitted, all cases are run.
"""

import shutil
from pathlib import Path

import amici
import pandas as pd
import pytest
from amici.gradient_check import check_derivatives

from utils import (
    verify_results,
    write_result_file,
    find_model_file,
    read_settings_file,
    apply_settings,
)


@pytest.fixture(scope="session")
def result_path() -> Path:
    return Path(__file__).parent / "amici-semantic-results"


def test_sbml_testsuite_case(test_id, result_path, sbml_semantic_cases_dir):
    model_dir = None

    # test cases for which sensitivities are to be checked
    #  key: case ID; value: epsilon for finite differences
    sensitivity_check_cases = {
        # parameter-dependent conservation laws
        "00783": 1.5e-2,
        # initial events
        "00995": 1e-3,
    }

    try:
        current_test_path = sbml_semantic_cases_dir / test_id

        # parse expected results
        results_file = current_test_path / f"{test_id}-results.csv"
        results = pd.read_csv(results_file, delimiter=",")
        results.rename(
            columns={c: c.replace(" ", "") for c in results.columns},
            inplace=True,
        )

        # setup model
        model_dir = Path(__file__).parent / "SBMLTestModels" / test_id
        model, solver, wrapper = compile_model(
            current_test_path,
            test_id,
            model_dir,
            generate_sensitivity_code=test_id in sensitivity_check_cases,
        )
        settings = read_settings_file(current_test_path, test_id)

        atol, rtol = apply_settings(settings, solver, model, test_id)

        # simulate model
        rdata = amici.runAmiciSimulation(model, solver)
        if rdata["status"] != amici.AMICI_SUCCESS:
            if test_id in ("00748", "00374", "00369"):
                pytest.skip("Simulation Failed expectedly")
            else:
                raise RuntimeError("Simulation failed unexpectedly")

        # verify
        simulated = verify_results(
            settings, rdata, results, wrapper, model, atol, rtol
        )

        # record results
        write_result_file(simulated, test_id, result_path)

        # check sensitivities for selected models
        if epsilon := sensitivity_check_cases.get(test_id):
            solver.setSensitivityOrder(amici.SensitivityOrder.first)
            solver.setSensitivityMethod(amici.SensitivityMethod.forward)
            check_derivatives(model, solver, epsilon=epsilon)

    except amici.sbml_import.SBMLException as err:
        pytest.skip(str(err))
    finally:
        if model_dir:
            shutil.rmtree(model_dir, ignore_errors=True)


def compile_model(
    sbml_dir: Path,
    test_id: str,
    model_dir: Path,
    generate_sensitivity_code: bool = False,
):
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

    model = model_module.getModel()
    solver = model.getSolver()

    return model, solver, sbml_importer
