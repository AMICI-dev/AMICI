#!/usr/bin/env python3
"""
Run SBML Test Suite and verify simulation results
[https://github.com/sbmlteam/sbml-test-suite/releases]

Usage:
    pytest tests.sbml.testSBMLSuite -n CORES --dist=loadgroup --cases=SELECTION
        CORES can be an integer or `auto` for all available cores.
        `--dist=loadgroup` is required whenever `-n` is used: each case's
        simulation and sensitivity checks share one compiled model module
        (`compiled_case` fixture) and must run on the same xdist worker.
        SELECTION can be e.g.: `1`, `1,3`, `-3,4,6-7`, or `100-` to select
        specific test cases. If `--cases` is omitted, all cases are run.
"""

from __future__ import annotations

import logging
import shutil
from collections.abc import Callable
from pathlib import Path

import amici
import diffrax
import jax
import jax.numpy as jnp
import libsbml
import numpy as np
import optimistix
import pandas as pd
import pytest
from amici.adapters.fiddy import run_simulation_to_function_and_derivative
from amici.sim.jax.petab import (
    DEFAULT_CONTROLLER_SETTINGS,
    DEFAULT_ROOT_FINDER_SETTINGS,
)
from amici.sim.sundials import (
    AMICI_SUCCESS,
    ExpData,
    Model,
    SensitivityMethod,
    SensitivityOrder,
    Solver,
    run_simulation,
)
from fiddy import FunctionEvaluationError, check_gradient, check_jacobian
from utils import (
    apply_settings,
    find_model_file,
    read_settings_file,
    verify_results,
    write_result_file,
)

# test cases for which the separate, autodiff-based JAX sensitivity
# cross-check is additionally run (see `jax_sensitivity_check`)
_JAX_CHECK_CASES = {
    # parameter-dependent conservation laws
    "00783",
    # initial events
    "00995",
}


@pytest.fixture(scope="session")
def compiled_case(test_id, sbml_semantic_cases_dir):
    """Compile a case's model once, reuse simulation and sensitivity tests.

    For use with pytest-xdist, see conftest.py.
    """
    current_test_path = sbml_semantic_cases_dir / test_id
    model_dir = Path(__file__).parent / "SBMLTestModels" / test_id
    try:
        model_module, sbml_importer = compile_model(
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

    yield model_module, sbml_importer, settings, current_test_path

    shutil.rmtree(model_dir, ignore_errors=True)


def _fresh_model_and_solver(
    model_module, settings: dict, test_id: str
) -> tuple[Model, Solver, float, float]:
    """Build a fresh `Model`/`Solver` and set up for test."""
    model = model_module.get_model()
    solver = model.create_solver()
    atol, rtol = apply_settings(settings, solver, model, test_id)
    solver.set_sensitivity_order(SensitivityOrder.first)
    solver.set_sensitivity_method(SensitivityMethod.forward)
    if test_id == "00885":
        # 00885: root-after-reinitialization with FSA with default settings
        solver.set_absolute_tolerance(1e-16)
        solver.set_relative_tolerance(1e-15)
    return model, solver, atol, rtol


def _check_simulation_status(rdata, test_id: str) -> None:
    """Skip/fail consistently for a known-bad vs. a genuinely unexpected
    base simulation failure."""
    if rdata["status"] != AMICI_SUCCESS:
        if test_id in ("00748", "00374", "00369"):
            pytest.skip("Simulation Failed expectedly")
        raise RuntimeError("Simulation failed unexpectedly")


def test_sbml_testsuite_case(test_id, compiled_case, result_path):
    model_module, sbml_importer, settings, current_test_path = compiled_case
    model, solver, atol, rtol = _fresh_model_and_solver(
        model_module, settings, test_id
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
        settings, rdata, results, sbml_importer, model, atol, rtol
    )

    # record results
    write_result_file(simulated, test_id, result_path)


# FIXME: Skip list - to be investigated further
#  test_id -> adjoint_only (whether forward is unaffected)
_OTHER_KNOWN_SENSITIVITY_CHECK_ISSUES = {
    "00048": True,
    "00066": True,
    "00208": True,
    "00589": True,
    "00879": False,
    "01530": False,
    "01104": True,
    "01107": True,
    "01148": True,
}


def _sensitivity_preflight_checks(
    model: Model, sbml_importer, test_id: str, uses_adjoint: bool
):
    """Skip if a sensitivity check wouldn't be meaningful/known-correct for
    this SBML feature.

    :param uses_adjoint: Whether the caller's check involves adjoint
        sensitivities.
    :return: The current libsbml model.
    """
    if not model.get_free_parameter_ids():
        pytest.skip("No free parameters to differentiate w.r.t.")

    sbml_model = sbml_importer.sbml_model
    if any(
        rule.getTypeCode() == libsbml.SBML_ALGEBRAIC_RULE
        for rule in sbml_model.getListOfRules()
    ):
        pytest.skip(
            "Sensitivities for AlgebraicRule models are known to "
            "be wrong -- see "
            "https://github.com/AMICI-dev/AMICI/issues/3250"
        )
    if uses_adjoint and model.nx_rdata == 0:
        pytest.skip(
            "Adjoint sensitivities for zero-state models are known to crash."
        )
    if test_id in _OTHER_KNOWN_SENSITIVITY_CHECK_ISSUES:
        adjoint_only = _OTHER_KNOWN_SENSITIVITY_CHECK_ISSUES[test_id]
        if uses_adjoint or not adjoint_only:
            pytest.skip("Known sensitivity-check issue, not yet investigated.")
    return sbml_model


def test_sbml_testsuite_case_sensitivity_forward(test_id, compiled_case):
    """Finite-difference-check the model's forward sensitivities."""
    model_module, sbml_importer, settings, current_test_path = compiled_case
    model, solver, atol, rtol = _fresh_model_and_solver(
        model_module, settings, test_id
    )
    sbml_model = _sensitivity_preflight_checks(
        model, sbml_importer, test_id, uses_adjoint=False
    )

    # Test whether base-simulation succeeds. If not, we can skip fail right
    # away with a clearer message than the FD check's own failure would give.
    rdata = run_simulation(model, solver)
    _check_simulation_status(rdata, test_id)

    def check(sensi_solver, point, bounds, retried):
        sensi_solver.set_sensitivity_method(SensitivityMethod.forward)
        function, derivative = run_simulation_to_function_and_derivative(
            amici_model=model,
            amici_solver=sensi_solver,
            derivative_variables=["x", "x0", "y", "sigmay"],
        )
        expected = derivative(point)
        # A bare parameter-only model has nothing for `x`/`x0`/`y`/`sigmay` to
        # report -- `derivative`/`function` then both return an empty
        # dict
        if not expected:
            return None
        result = check_jacobian(function, point, expected, bounds=bounds)
        result.assert_success(always_print=True)
        _assert_check_confirms_something(result, has_events, retried)

    has_events = sbml_model.getNumEvents() > 0
    _run_sensitivity_check(model, settings, test_id, has_events, check)

    # additionally cross-check against JAX autodiff for a couple of
    # historically tricky cases
    if test_id in _JAX_CHECK_CASES:
        jax_sensitivity_check(
            current_test_path,
            test_id,
            model,
            rdata,
            atol,
            rtol,
        )


@pytest.mark.filterwarnings(
    # https://github.com/AMICI-dev/AMICI/issues/18
    "ignore:Adjoint sensitivity analysis for models with discontinuous "
    "right hand sides .*:UserWarning",
)
def test_sbml_testsuite_case_sensitivity_adjoint(test_id, compiled_case):
    """Finite-difference-check the model's adjoint sensitivities."""
    model_module, sbml_importer, settings, current_test_path = compiled_case
    model, solver, atol, rtol = _fresh_model_and_solver(
        model_module, settings, test_id
    )
    sbml_model = _sensitivity_preflight_checks(
        model, sbml_importer, test_id, uses_adjoint=True
    )

    # generate synthetic measurements
    rdata = run_simulation(model, solver)
    _check_simulation_status(rdata, test_id)
    # `amici.ExpData(rdata, sigma_y, sigma_z, seed)`
    edata = ExpData(rdata, 1.0, 1.0, 42)

    def check(sensi_solver, point, bounds, retried):
        sensi_solver.set_sensitivity_method(SensitivityMethod.adjoint)
        sensi_solver.set_absolute_tolerance_b(
            sensi_solver.get_absolute_tolerance()
        )
        sensi_solver.set_relative_tolerance_b(
            sensi_solver.get_relative_tolerance()
        )

        def _run(atol_quad, rtol_quad):
            sensi_solver.set_absolute_tolerance_quadratures(atol_quad)
            sensi_solver.set_relative_tolerance_quadratures(rtol_quad)
            function, derivative = run_simulation_to_function_and_derivative(
                amici_model=model,
                amici_solver=sensi_solver,
                amici_edata=edata,
                derivative_variables=["llh"],
            )
            return function, derivative(point)

        # Try the tightest quadrature tolerance first: loosening it can
        # silently corrupt an otherwise-successful backward integration's
        # result (confirmed on case 00945: keeping this fixed tight gives
        # ~0.1% error at the same state/backward tolerance that gives ~7%
        # error if quadrature tolerance also scales alongside it). Only if
        # that fails outright (NaN) set the quadrature
        # tolerance to the already-escalated backward tolerance -- some
        # cases' backward integration genuinely cannot converge *at all*
        # without that (confirmed on cases 00754/00755/00756).
        function, expected = _run(1e-16, 1e-15)
        if expected and np.any(np.isnan(expected["llh"])):
            function, expected = _run(
                sensi_solver.get_absolute_tolerance_b(),
                sensi_solver.get_relative_tolerance_b(),
            )
        if not expected:
            return None
        # Raise on simulation failure (NaN)
        if np.any(np.isnan(expected["llh"])):
            raise FunctionEvaluationError("Simulation failed.")
        result = check_gradient(
            function, point, expected["llh"], bounds=bounds
        )
        result.assert_success(always_print=True)
        _assert_check_confirms_something(result, has_events, retried)

    has_events = sbml_model.getNumEvents() > 0
    _run_sensitivity_check(
        model, settings, test_id, has_events, check, jitter_seed=43
    )


@pytest.mark.filterwarnings(
    # https://github.com/AMICI-dev/AMICI/issues/18
    "ignore:Adjoint sensitivity analysis for models with discontinuous "
    "right hand sides .*:UserWarning",
)
def test_sbml_testsuite_case_sensitivity_consistency(test_id, compiled_case):
    """Compare forward vs. adjoint sensitivities (`sllh`) against *each other*."""
    model_module, sbml_importer, settings, current_test_path = compiled_case
    model, solver, atol, rtol = _fresh_model_and_solver(
        model_module, settings, test_id
    )
    sbml_model = _sensitivity_preflight_checks(
        model, sbml_importer, test_id, uses_adjoint=True
    )

    rdata = run_simulation(model, solver)
    _check_simulation_status(rdata, test_id)
    edata = ExpData(rdata, 1.0, 1.0, 42)

    def check(sensi_solver, point, bounds, retried):
        # FSA
        sensi_solver.set_sensitivity_method(SensitivityMethod.forward)
        _, fsa_derivative = run_simulation_to_function_and_derivative(
            amici_model=model,
            amici_solver=sensi_solver,
            amici_edata=edata,
            derivative_variables=["llh"],
        )
        fsa_expected = fsa_derivative(point)
        if not fsa_expected:
            return None

        # ASA
        asa_solver = model.create_solver()
        apply_settings(settings, asa_solver, model, test_id)
        asa_solver.set_absolute_tolerance(
            sensi_solver.get_absolute_tolerance()
        )
        asa_solver.set_relative_tolerance(
            sensi_solver.get_relative_tolerance()
        )
        asa_solver.set_absolute_tolerance_b(
            sensi_solver.get_absolute_tolerance()
        )
        asa_solver.set_relative_tolerance_b(
            sensi_solver.get_relative_tolerance()
        )
        asa_solver.set_sensitivity_order(SensitivityOrder.first)
        asa_solver.set_sensitivity_method(SensitivityMethod.adjoint)

        def _run_asa(atol_quad, rtol_quad):
            asa_solver.set_absolute_tolerance_quadratures(atol_quad)
            asa_solver.set_relative_tolerance_quadratures(rtol_quad)
            _, asa_derivative = run_simulation_to_function_and_derivative(
                amici_model=model,
                amici_solver=asa_solver,
                amici_edata=edata,
                derivative_variables=["llh"],
            )
            return asa_derivative(point)

        # Same two-tier quadrature-tolerance strategy (tightest first,
        # falling back to matching the already-escalated backward
        # tolerance only if that's NaN) as
        # `test_sbml_testsuite_case_sensitivity_adjoint`'s `check` -- see
        # its comment for why.
        asa_expected = _run_asa(1e-16, 1e-15)
        if asa_expected and np.any(np.isnan(asa_expected["llh"])):
            asa_expected = _run_asa(
                asa_solver.get_absolute_tolerance_b(),
                asa_solver.get_relative_tolerance_b(),
            )
        if not asa_expected:
            return None

        fsa_sllh = fsa_expected["llh"]
        asa_sllh = asa_expected["llh"]
        if np.any(np.isnan(fsa_sllh)) or np.any(np.isnan(asa_sllh)):
            raise FunctionEvaluationError(
                "AMICI's forward or adjoint sensitivity computation "
                "returned NaN (likely a near-tangent event crossing)."
            )
        # Tight by default; deliberately loosened once we've already had
        # to retry at a much looser solver tolerance to get a finite
        # result at all -- that widening alone causes up to ~1% residual
        # disagreement between otherwise-correct forward and adjoint
        # sensitivities (case 00753)
        rtol = 1e-6 if not retried else 1e-2
        np.testing.assert_allclose(
            asa_sllh,
            fsa_sllh,
            rtol=rtol,
            atol=1e-8,
            err_msg="Forward and adjoint sensitivities disagree",
        )

    has_events = sbml_model.getNumEvents() > 0
    _run_sensitivity_check(
        model, settings, test_id, has_events, check, jitter_seed=44
    )


def compile_model(
    sbml_dir: Path,
    test_id: str,
    model_dir: Path,
    generate_sensitivity_code: bool = False,
):
    """Import the given test model."""
    model_dir.mkdir(parents=True, exist_ok=True)

    sbml_file = find_model_file(sbml_dir, test_id)
    sbml_importer = amici.SbmlImporter(sbml_file)

    model_name = f"SBMLTest{test_id}"
    sbml_importer.sbml2amici(
        model_name,
        output_dir=model_dir,
        generate_sensitivity_code=generate_sensitivity_code,
    )

    model_module = amici.import_model_module(model_name, model_dir)

    return model_module, sbml_importer


# Case 01395-specific: `v1_h`..`v15_h` are Hill coefficients with nominal
# value exactly 1, acting on species (`p1`, `p2`, ...) that start at
# exactly 0. The forward-sensitivity RHS (not the state RHS itself)
# contains a `species**(h-1)` term; reducing `h` even infinitesimally
# below 1 makes that exponent negative, and `0**(negative)` is `+inf`.
_HILL_EXPONENT_LOWER_BOUND_OVERRIDES = {
    "01395": tuple(f"v{i}_h" for i in range(1, 16)),
}


def _derive_generic_bounds(
    point: np.ndarray, param_ids: tuple[str, ...], test_id: str
) -> tuple[np.ndarray, np.ndarray]:
    """Derive a generic, per-parameter valid domain from nominal values
    alone.

    Without bounds, fiddy might try to evaluate infeasible points.
    A bound of nominal * (1e-6, 1000) on the same-signed side (never
    including 0 itself -- needed for case 00313, singular at 0) fixes
    69/76 previously-failing cases outright. Parameters nominally at 0
    fall back to (0, 1000): every free parameter in this suite is
    semantically non-negative.

    :param point: Nominal free-parameter values.
    :param param_ids: `point`'s parameter IDs, same order -- only used to
        apply `_HILL_EXPONENT_LOWER_BOUND_OVERRIDES`.
    :param test_id: The SBML semantic test suite case ID -- ditto.
    :return: A `(lower, upper)` bounds tuple, same shape as `point`.
    """
    lower = np.empty_like(point)
    upper = np.empty_like(point)
    positive = point > 0
    negative = point < 0
    zero = ~positive & ~negative
    lower[positive] = np.maximum(point[positive] * 1e-6, 1e-9)
    upper[positive] = np.maximum(point[positive] * 1000, 1000.0)
    lower[negative] = np.minimum(point[negative] * 1000, -1000.0)
    upper[negative] = np.minimum(point[negative] * 1e-6, -1e-9)
    lower[zero] = 0.0
    upper[zero] = 1000.0
    for param_id in _HILL_EXPONENT_LOWER_BOUND_OVERRIDES.get(test_id, ()):
        idx = param_ids.index(param_id)
        lower[idx] = point[idx]
    return lower, upper


def _assert_check_confirms_something(
    result, has_events: bool, retried: bool
) -> None:
    """`check_jacobian`/`check_gradient`'s `success` is `True` as long as no
    direction is confidently *wrong* -- a check where every direction came
    back "inconclusive" (noise-dominated/discontinuity-suspected) would
    still report success, having actually confirmed nothing. Require at
    least one direction, in at least one output component, to have been
    confirmed converged, so a silent coverage regression fails loudly instead
    of passing vacuously.

    Exception: for an event-triggered model, if every non-passing direction was
    flagged `"discontinuity_suspected"`, skip instead of failing.
    Same for `"noise_dominated"`, but only once we've already had to retry at
    a much looser tolerance to get a finite result at all (`retried`).
    """
    direction_results = (
        [d for o in result.output_results for d in o.direction_results]
        if hasattr(result, "output_results")
        else result.direction_results
    )
    if any(r.outcome == "passed" for r in direction_results):
        return
    acceptable_statuses = {"discontinuity_suspected"}
    if retried:
        acceptable_statuses.add("noise_dominated")
    if has_events and all(
        r.estimate.status in acceptable_statuses for r in direction_results
    ):
        pytest.skip(
            "Every checked direction was inconclusive "
            f"({sorted({r.estimate.status for r in direction_results})})."
        )
    raise AssertionError(
        "check reported success, but every direction was inconclusive -- "
        "nothing was actually confirmed correct."
    )


def _run_sensitivity_check(
    model: Model,
    settings: dict,
    test_id: str,
    has_events: bool,
    check: Callable[
        [Solver, np.ndarray, tuple[np.ndarray, np.ndarray]], object
    ],
    jitter_seed: int | None = None,
) -> None:
    """Retry `check` (a sensitivity FD check via fiddy)
    at looser solver tolerance (event-triggered models only) if the
    simulation itself fails to produce a finite value.

    An FD-perturbed parameter point can turn a clean event trigger crossing
    into a near-tangent one (confirmed on cases 00375/00754), triggering
    AMICI's "root after reinitialization" error. To avoid this error,
    we retry the check at a looser solver tolerance.

    :param model: The AMICI model.
    :param settings: This case's parsed `{test_id}-settings.txt`.
    :param test_id: The SBML semantic test suite case ID.
    :param has_events: Whether the model has any SBML events -- gates the
        tolerance-retry loop; non-event models get exactly one attempt.
    :param check: Called as `check(sensi_solver, point, bounds, retried)`
        (`retried` is `True` once a looser-than-default tolerance was
        needed to get this far) -- responsible for configuring
        `sensi_solver`'s sensitivity method, running the check, and
        asserting success itself. Returns `None` if there's nothing to check
        for this model -- in which case `_run_sensitivity_check` returns
        immediately without retrying; any other return value is ignored.
    :param jitter_seed: If given, `point` is perturbed by fixed-seed 5%
        relative Gaussian noise (then clipped back into `bounds`) before
        checking -- the same "avoid small gradients at nominal value"
        mitigation `tests/benchmark_models/test_petab_benchmark.py`'s
        gradient check already uses.
    """
    # Forward only ever needs up to 100x (e.g. 00375/00754);
    # adjoint's backward+quadrature integration can need
    # much looser tolerance still to get through some event crossings at
    # all -- confirmed on cases 00026/00041/00074/00745/00746/00747/
    # 00789/00845/00945 (all fail with CVodeB's error-test repeatedly
    # failing/`|h| = hmin` up to 1e4x, several needing up to 1e8x).
    tolerance_multipliers = (1, 10, 100, 1e4, 1e6, 1e8) if has_events else (1,)
    param_ids = model.get_free_parameter_ids()
    point = np.asarray(
        [
            model.get_free_parameter_by_id(parameter_id)
            for parameter_id in param_ids
        ]
    )
    bounds = _derive_generic_bounds(point, param_ids, test_id)
    if jitter_seed is not None:
        rng = np.random.default_rng(jitter_seed)
        point = point + rng.standard_normal(len(point)) * point * 0.05
        point = np.clip(point, *bounds)
    error = None
    for multiplier in tolerance_multipliers:
        # Use a fresh solver every attempt
        # avoids failures for zero-state models
        # (this calls `CVodeSensSStolerances` without
        # having gone through a matching `CVodeSensReInit` for the updated
        # parameters first, which CVODES rejects with `CV_ILL_INPUT`
        # "CVODE routine CVodeSensSStolerances failed with error code -40").
        sensi_solver = model.create_solver()
        apply_settings(settings, sensi_solver, model, test_id)
        # start with tightest tolerances, then loosen if needed
        sensi_solver.set_absolute_tolerance(1e-16)
        sensi_solver.set_relative_tolerance(1e-15)
        if multiplier != 1:
            sensi_solver.set_absolute_tolerance(
                sensi_solver.get_absolute_tolerance() * multiplier
            )
            sensi_solver.set_relative_tolerance(
                sensi_solver.get_relative_tolerance() * multiplier
            )
        sensi_solver.set_sensitivity_order(SensitivityOrder.first)

        # Skip, not fail on "root after reinitialization" errors during FD
        # checks that occur at any tolerance
        # (e.g., for 00369/00754/00755/00756/00883/00885).
        root_after_reinit_messages = []
        log_handler = logging.Handler()
        log_handler.emit = lambda record: (
            root_after_reinit_messages.append(record.getMessage())
            if "root after reinitialization" in record.getMessage()
            else None
        )
        amici_logger = logging.getLogger("amici.sim.sundials._swig_wrappers")
        amici_logger.addHandler(log_handler)
        try:
            if check(sensi_solver, point, bounds, multiplier != 1) is None:
                return
        except FunctionEvaluationError as err:
            if root_after_reinit_messages:
                pytest.skip(
                    "AMICI/CVODES cannot integrate through a near-tangent "
                    f"event crossing: {root_after_reinit_messages[0]}"
                )
            error = err
            continue
        finally:
            amici_logger.removeHandler(log_handler)
        return
    raise error


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
