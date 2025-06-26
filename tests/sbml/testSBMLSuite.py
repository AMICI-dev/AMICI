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
from fiddy import MethodId, get_derivative
from fiddy.extensions.amici import (
    reshape,
    run_amici_simulation_to_cached_functions,
)
from fiddy.success import Consistency
from utils import (
    verify_results,
    write_result_file,
    find_model_file,
    read_settings_file,
    apply_settings,
)
import libsbml
import numpy as np


@pytest.fixture(scope="session")
def result_path() -> Path:
    return Path(__file__).parent / "amici-semantic-results"


def test_sbml_testsuite_case(test_id, result_path, sbml_semantic_cases_dir):
    model_dir = None

    # test cases for which sensitivities are to be checked
    #  key: case ID; value: epsilon for finite differences
    # sensitivity_check_cases = {
    #     # parameter-dependent conservation laws
    #     "00783": 1.5e-2,
    #     # initial events
    #     "00995": 1e-3,
    # }

    try:
        current_test_path = sbml_semantic_cases_dir / test_id

        # parse expected results
        results_file = current_test_path / f"{test_id}-results.csv"
        results = pd.read_csv(results_file, delimiter=",")
        results.rename(
            columns={c: c.replace(" ", "") for c in results.columns},
            inplace=True,
        )

        # TODO remove after https://github.com/AMICI-dev/AMICI/pull/2101
        #   and https://github.com/AMICI-dev/AMICI/issues/2106
        # Don't attempt to generate sensitivity code for models with events+algebraic rules, which will fail
        sbml_file = find_model_file(current_test_path, test_id)
        sbml_document = libsbml.SBMLReader().readSBMLFromFile(str(sbml_file))
        sbml_model = sbml_document.getModel()
        has_events = sbml_model.getNumEvents() > 0
        has_algebraic_rules = any(
            rule.getTypeCode() == libsbml.SBML_ALGEBRAIC_RULE
            for rule in sbml_model.getListOfRules()
        )
        generate_sensitivity_code = not (has_events and has_algebraic_rules)
        # TODO https://github.com/AMICI-dev/AMICI/issues/2109
        generate_sensitivity_code &= test_id not in {"01240"}
        # ^^^^^^^^

        # setup model
        model_dir = Path(__file__).parent / "SBMLTestModels" / test_id
        model, solver, wrapper = compile_model(
            current_test_path,
            test_id,
            model_dir,
            generate_sensitivity_code=generate_sensitivity_code,
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

        # verify simulation results
        simulated = verify_results(
            settings, rdata, results, wrapper, model, atol, rtol
        )

        # record results
        write_result_file(simulated, test_id, result_path)

        # test sensitivities
        if not model.getParameters():
            pytest.skip("No parameters -> no sensitivities to check")

        # TODO see https://github.com/AMICI-dev/AMICI/pull/2101
        if not generate_sensitivity_code:
            pytest.skip("Sensitivity analysis is known to fail.")
        if any(id_ == 0 for id_ in model.idlist):
            pytest.skip("Sensitivity analysis for DAE is known to fail.")

        solver.setSensitivityOrder(amici.SensitivityOrder.first)
        solver.setSensitivityMethod(amici.SensitivityMethod.forward)
        # currently only checking "x"/"sx" for FSA
        (
            amici_function_f,
            amici_derivative_f,
            structures_f,
        ) = run_amici_simulation_to_cached_functions(
            amici_model=model,
            amici_solver=solver,
            derivative_variables=["x"],
            cache=False,
        )
        rdata_f = amici.runAmiciSimulation(model, solver)

        # solver.setSensitivityMethod(amici.SensitivityMethod.adjoint)
        # (
        #    amici_function_a,
        #    amici_derivative_a,
        #    structures_a,
        # ) = run_amici_simulation_to_cached_functions(
        #    amici_model=model,
        #    amici_solver=solver,
        #    derivative_variables=["x"],
        #    cache=False,
        # )
        # rdata_a = amici.runAmiciSimulation(model, solver)

        point = np.asarray(model.getParameters())

        derivative = get_derivative(
            # can use `_f` or `_a` here, should be no difference
            function=amici_function_f,
            point=point,
            sizes=[1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1],
            direction_ids=model.getParameterIds(),
            method_ids=[MethodId.FORWARD, MethodId.BACKWARD, MethodId.CENTRAL],
            relative_sizes=True,
            success_checker=Consistency(rtol=1e-2, atol=1e-4),
        )

        derivative_fd = reshape(
            derivative.value.flat,
            structures_f["derivative"],
            sensitivities=True,
        )["x"]
        derivative_fsa = rdata_f.sx
        # derivative_asa = rdata_a.sllh  # currently None, define some objective?

        # could alternatively use a `fiddy.DerivativeCheck` class
        if not np.isclose(
            derivative_fd, derivative_fsa, rtol=5e-2, atol=5e-2
        ).all():
            raise ValueError("Gradients were not validated.")

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
