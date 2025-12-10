"""Tests for SBML events, including piecewise expressions."""

import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from amici.importers.sbml import SbmlImporter
from amici.sim.sundials import (
    AMICI_SUCCESS,
    AmiciModel,
    ExpData,
    SensitivityMethod,
    SensitivityOrder,
    import_model_module,
    run_simulation,
)
from amici.sim.sundials.gradient_check import _check_close
from numpy.testing import assert_allclose


def create_amici_model(sbml_model, model_name, **kwargs) -> AmiciModel:
    """
    Import an sbml file and create an AMICI model from it
    """
    sbml_test_models_output_dir = Path("amici_models")
    sbml_test_models_output_dir.mkdir(parents=True, exist_ok=True)

    sbml_importer = SbmlImporter(sbml_model)
    # try not to exceed the stupid maximum path length on windows 💩
    output_dir = (
        sbml_test_models_output_dir / model_name
        if sys.platform != "win32"
        else tempfile.mkdtemp()
    )

    sbml_importer.sbml2amici(
        model_name=model_name, output_dir=output_dir, **kwargs
    )

    model_module = import_model_module(model_name, output_dir)
    return model_module.get_model()


def check_trajectories_without_sensitivities(
    amici_model: AmiciModel,
    result_expected_x: np.ndarray,
):
    """
    Check whether the AMICI simulation matches a known solution
    (ideally an analytically calculated one).
    """
    solver = amici_model.create_solver()
    solver.set_absolute_tolerance(1e-15)
    solver.set_relative_tolerance(1e-12)
    rdata = run_simulation(amici_model, solver=solver)
    _check_close(
        rdata["x"], result_expected_x, field="x", rtol=5e-9, atol=1e-13
    )


def check_trajectories_with_forward_sensitivities(
    amici_model: AmiciModel,
    result_expected_x: np.ndarray,
    result_expected_sx: np.ndarray,
):
    """
    Check whether the forward sensitivities of the AMICI simulation match
    a known solution (ideally an analytically calculated one).
    """
    solver = amici_model.create_solver()
    solver.set_sensitivity_order(SensitivityOrder.first)
    solver.set_sensitivity_method(SensitivityMethod.forward)
    solver.set_absolute_tolerance(1e-15)
    solver.set_relative_tolerance(1e-13)
    solver.set_absolute_tolerance_fsa(1e-15)
    solver.set_relative_tolerance_fsa(1e-13)
    rdata = run_simulation(amici_model, solver=solver)
    _check_close(
        rdata["x"], result_expected_x, field="x", rtol=1e-10, atol=1e-12
    )
    _check_close(
        rdata["sx"], result_expected_sx, field="sx", rtol=1e-7, atol=1e-9
    )


def check_trajectories_with_adjoint_sensitivities(
    amici_model: AmiciModel, asa_xfail: bool = False
):
    """
    Check whether adjoint sensitivities match forward sensitivities and finite
    differences.

    :param amici_model: AMICI model to test
    :param asa_xfail: If True, deviations between adjoint and forward
        sensitivities will not raise an AssertionError.
    """
    # First compute dummy experimental data to use adjoints
    solver = amici_model.create_solver()
    rdata = run_simulation(amici_model, solver=solver)
    assert rdata.status == AMICI_SUCCESS
    rng_seed = 42
    edata = ExpData(rdata, 1.0, 1.0, rng_seed)

    # FSA
    solver.set_sensitivity_order(SensitivityOrder.first)
    solver.set_sensitivity_method(SensitivityMethod.forward)
    solver.set_absolute_tolerance(1e-15)
    solver.set_relative_tolerance(1e-13)
    rdata_fsa = run_simulation(amici_model, solver=solver, edata=edata)
    assert rdata_fsa.status == AMICI_SUCCESS

    # ASA
    solver.set_sensitivity_method(SensitivityMethod.adjoint)
    solver.set_absolute_tolerance(1e-16)
    solver.set_relative_tolerance(1e-14)
    solver.set_absolute_tolerance_b(1e-16)
    solver.set_relative_tolerance_b(1e-15)
    solver.set_absolute_tolerance_quadratures(1e-14)
    solver.set_relative_tolerance_quadratures(1e-8)
    rdata_asa = run_simulation(amici_model, solver=solver, edata=edata)
    assert rdata_asa.status == AMICI_SUCCESS

    assert_allclose(rdata_fsa.x, rdata_asa.x, atol=1e-14, rtol=1e-10)
    assert_allclose(rdata_fsa.llh, rdata_asa.llh, atol=1e-14, rtol=1e-10)
    df = pd.DataFrame(
        {
            "fsa": rdata_fsa["sllh"],
            "asa": rdata_asa["sllh"],
            "fd": np.nan,
        },
        index=list(amici_model.get_free_parameter_ids()),
    )
    df["abs_diff"] = df["fsa"] - df["asa"]
    df["rel_diff"] = df["abs_diff"] / df["fsa"]

    # Also test against finite differences
    parameters = amici_model.get_unscaled_parameters()
    solver.set_sensitivity_order(SensitivityOrder.none)
    sllh_fd = []
    eps = 1e-5
    for i_par, par in enumerate(parameters):
        tmp_par = np.array(parameters[:])
        tmp_par[i_par] += eps
        amici_model.set_free_parameters(tmp_par)
        rdata_p = run_simulation(amici_model, solver=solver, edata=edata)
        tmp_par = np.array(parameters[:])
        tmp_par[i_par] -= eps
        amici_model.set_free_parameters(tmp_par)
        rdata_m = run_simulation(amici_model, solver=solver, edata=edata)
        sllh_fd.append((rdata_p["llh"] - rdata_m["llh"]) / (2 * eps))
    df["fd"] = sllh_fd
    df["asa_matches_fsa"] = np.isclose(
        df["asa"], df["fsa"], rtol=1e-8, atol=1e-12
    )
    print()
    with pd.option_context(
        "display.max_rows",
        None,
        "display.max_columns",
        None,
        "display.width",
        None,
    ):
        print(df)

    if asa_xfail:
        if df["asa_matches_fsa"].all():
            # Note that this is machine-dependent...
            print("Incorrectly marked as xfail?")
        print(
            "FIXME: Ignoring differences between adjoint and forward sensitivities."
        )
        return

    assert_allclose(
        sllh_fd,
        rdata_fsa["sllh"],
        rtol=1e-5,
        atol=1e-8,
        err_msg="Finite differences and forward sensitivities do not match.",
    )

    assert_allclose(
        rdata_fsa["sllh"],
        rdata_asa["sllh"],
        rtol=1e-8,
        atol=1e-12,
        err_msg="Adjoint and forward sensitivities do not match.",
    )
    assert_allclose(
        sllh_fd,
        rdata_asa["sllh"],
        rtol=1e-5,
        atol=1e-6,
        err_msg="Finite differences and adjoint sensitivities do not match.",
    )
