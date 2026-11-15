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

import shutil
from collections.abc import Callable
from pathlib import Path

import amici
import diffrax
import jax
import jax.numpy as jnp
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
    return model, solver, atol, rtol


def _check_simulation_status(rdata) -> None:
    """Fail on an unsuccessful base simulation."""
    if rdata["status"] != AMICI_SUCCESS:
        raise RuntimeError("Simulation failed unexpectedly")


def test_sbml_testsuite_case(test_id, compiled_case, result_path):
    model_module, sbml_importer, settings, current_test_path = compiled_case
    model, solver, atol, rtol = _fresh_model_and_solver(
        model_module, settings, test_id
    )

    # parse expected results
    results_file = current_test_path / f"{test_id}-results.csv"
    results = pd.read_csv(results_file, delimiter=",")
    results = results.rename(
        columns={c: c.replace(" ", "") for c in results.columns},
    )

    # simulate model
    rdata = run_simulation(model, solver)
    _check_simulation_status(rdata)

    # verify
    simulated = verify_results(
        settings, rdata, results, sbml_importer, model, atol, rtol
    )

    # record results
    write_result_file(simulated, test_id, result_path)


# Known sensitivity check failures: test_id -> (affected checks, reason)
_KNOWN_SENSITIVITY_CHECK_ISSUES = {
    "01530": (
        {"forward", "adjoint"},
        "The free parameters are unused and nominally at their lower "
        "bound 0, so fiddy can't evaluate any direction.",
    ),
}


def _sensitivity_preflight_checks(
    model: Model, sbml_importer, test_id: str, check: str
):
    """Skip if a sensitivity check wouldn't be meaningful for this case.

    :param check: The calling check: "forward", "adjoint", or "consistency".
    :return: The current libsbml model.
    """
    if not model.get_free_parameter_ids():
        pytest.skip("No free parameters to differentiate w.r.t.")
    if not np.all(np.isfinite(model.get_free_parameters())):
        pytest.skip("Non-finite nominal parameter values.")
    if check == "forward":
        if model.nx_rdata == 0 and model.ny == 0:
            pytest.skip("No states or observables.")
    else:
        if model.ny == 0:
            pytest.skip("No observables, the likelihood is constant.")
        if model.nx_rdata == 0:
            pytest.skip(
                "Adjoint sensitivity analysis fails for models without "
                "state variables (#3546)."
            )
    if test_id in _KNOWN_SENSITIVITY_CHECK_ISSUES:
        checks, reason = _KNOWN_SENSITIVITY_CHECK_ISSUES[test_id]
        if check in checks:
            pytest.skip(f"Known sensitivity check issue: {reason}")
    return sbml_importer.sbml_model


def test_sbml_testsuite_case_sensitivity_forward(test_id, compiled_case):
    """Finite-difference-check the model's forward sensitivities."""
    model_module, sbml_importer, settings, current_test_path = compiled_case
    model, solver, atol, rtol = _fresh_model_and_solver(
        model_module, settings, test_id
    )
    sbml_model = _sensitivity_preflight_checks(
        model, sbml_importer, test_id, check="forward"
    )

    # Test whether base-simulation succeeds. If not, we can fail right
    # away with a clearer message than the FD check's own failure would give.
    rdata = run_simulation(model, solver)
    _check_simulation_status(rdata)

    def check(sensi_solver, point, bounds, retried):
        sensi_solver.set_sensitivity_method(SensitivityMethod.forward)
        function, derivative = run_simulation_to_function_and_derivative(
            amici_model=model,
            amici_solver=sensi_solver,
            derivative_variables=["x", "x0", "y", "sigmay"],
        )
        expected = derivative(point)
        result = check_jacobian(function, point, expected, bounds=bounds)
        result.assert_success(always_print=True)
        _assert_check_confirms_something(result, has_events)

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
        model, sbml_importer, test_id, check="adjoint"
    )

    # generate synthetic measurements
    rdata = run_simulation(model, solver)
    _check_simulation_status(rdata)
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
        # that fails outright (failure or NaN) set the quadrature
        # tolerance to the already-escalated backward tolerance -- some
        # cases' backward integration genuinely cannot converge *at all*
        # without that (confirmed on cases 00754/00755/00756).
        try:
            function, expected = _run(1e-16, 1e-15)
        except FunctionEvaluationError:
            expected = None
        if expected is None or np.any(np.isnan(expected["llh"])):
            function, expected = _run(
                sensi_solver.get_absolute_tolerance_b(),
                sensi_solver.get_relative_tolerance_b(),
            )
        if np.any(np.isnan(expected["llh"])):
            raise FunctionEvaluationError("Gradient is NaN.")
        result = check_gradient(
            function, point, expected["llh"], bounds=bounds
        )
        result.assert_success(always_print=True)
        _assert_check_confirms_something(result, has_events)

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
        model, sbml_importer, test_id, check="consistency"
    )

    rdata = run_simulation(model, solver)
    _check_simulation_status(rdata)
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
        asa_solver.set_max_steps_backward_problem(
            sensi_solver.get_max_steps_backward_problem()
        )

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
        # tolerance only if that fails) as
        # `test_sbml_testsuite_case_sensitivity_adjoint`'s `check` -- see
        # its comment for why.
        try:
            asa_expected = _run_asa(1e-16, 1e-15)
        except FunctionEvaluationError:
            asa_expected = None
        if asa_expected is None or np.any(np.isnan(asa_expected["llh"])):
            asa_expected = _run_asa(
                asa_solver.get_absolute_tolerance_b(),
                asa_solver.get_relative_tolerance_b(),
            )

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


# Case-specific bounds overriding `_derive_generic_bounds`:
#  test_id -> {parameter_id: (lower, upper)}, as factors of the nominal
#  value (assumed positive); `None` keeps the generic bound.
_BOUND_OVERRIDES = {
    # `v1_h`..`v15_h` are Hill coefficients with nominal value exactly 1,
    # acting on species (`p1`, `p2`, ...) that start at exactly 0. The
    # forward-sensitivity RHS (not the state RHS itself) contains a
    # `species**(h-1)` term; reducing `h` even infinitesimally below 1
    # makes that exponent negative, and `0**(negative)` is `+inf`.
    "01395": {f"v{i}_h": (1.0, None) for i in range(1, 16)},
    # The species grows like `exp(c * t**2)` with
    # `c ~ extentconv**3 / timeconv**6`: overflow at nominal * (0.5, 2).
    "01148": {
        f"{prefix}{name}": (0.9, 1.1)
        for prefix in ("", "sub1__", "sub1__sub1__")
        for name in ("timeconv", "extentconv")
    },
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
        apply `_BOUND_OVERRIDES`.
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
    for param_id, (lower_factor, upper_factor) in _BOUND_OVERRIDES.get(
        test_id, {}
    ).items():
        idx = param_ids.index(param_id)
        if lower_factor is not None:
            lower[idx] = point[idx] * lower_factor
        if upper_factor is not None:
            upper[idx] = point[idx] * upper_factor
    return lower, upper


#: Maximum relative difference between finite difference and sensitivity for
#: an inconclusive direction to count as confirmed.
_AGREEMENT_RTOL = 1e-6


def _agrees_closely(direction_result) -> bool:
    """Check whether a direction's finite difference and sensitivity agree.

    Exactly matching zeros don't count.
    """
    scale = max(abs(direction_result.test), abs(direction_result.expectation))
    return (
        abs(direction_result.test - direction_result.expectation)
        < _AGREEMENT_RTOL * scale
    )


def _assert_check_confirms_something(result, has_events: bool) -> None:
    """`check_jacobian`/`check_gradient`'s `success` is `True` as long as no
    direction is confidently *wrong* -- a check where every direction came
    back "inconclusive" (noise-dominated/discontinuity-suspected) would
    still report success, having actually confirmed nothing. Require at
    least one direction, in at least one output component, to have been
    confirmed converged, so a silent coverage regression fails loudly instead
    of passing vacuously.

    An inconclusive direction also counts as confirmed if the finite
    difference and the sensitivity agree to within a relative difference of
    `_AGREEMENT_RTOL`. fiddy's classification of well-resolved directions as
    converged is sensitive to round-off errors, e.g., on different CPUs.

    Exception: for an event-triggered model, if every non-passing direction was
    flagged `"discontinuity_suspected"` or `"noise_dominated"`, skip instead
    of failing.
    """
    direction_results = (
        [d for o in result.output_results for d in o.direction_results]
        if hasattr(result, "output_results")
        else result.direction_results
    )
    if any(
        r.outcome == "passed"
        or (r.outcome == "inconclusive" and _agrees_closely(r))
        for r in direction_results
    ):
        return
    acceptable_statuses = {"discontinuity_suspected", "noise_dominated"}
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
        [Solver, np.ndarray, tuple[np.ndarray, np.ndarray], bool], None
    ],
    jitter_seed: int | None = None,
) -> None:
    """Retry `check` (a sensitivity FD check via fiddy)
    at looser solver tolerance (event-triggered models and DAEs only) if
    the simulation itself fails.

    :param model: The AMICI model.
    :param settings: This case's parsed `{test_id}-settings.txt`.
    :param test_id: The SBML semantic test suite case ID.
    :param has_events: Whether the model has any SBML events -- gates the
        tolerance-retry loop, together with the model being a DAE;
        other models get exactly one attempt.
    :param check: Called as `check(sensi_solver, point, bounds, retried)`
        (`retried` is `True` once a looser tolerance than the first attempted
        one was needed to get this far) -- responsible for configuring
        `sensi_solver`'s sensitivity method, running the check, and
        asserting success itself.
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
    # IDAS fails at tolerances close to machine precision (forward
    # sensitivities silently become NaN, backward integration fails), but
    # always succeeded at 100x in the semantic test suite.
    is_dae = any(i == 0.0 for i in model.get_id_list())
    if is_dae:
        tolerance_multipliers = (100, 1e4, 1e6, 1e8)
    elif has_events:
        tolerance_multipliers = (1, 10, 100, 1e4, 1e6, 1e8)
    else:
        tolerance_multipliers = (1,)
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
        # The default (100 * forward max steps) lets failing backward
        # integrations run for minutes.
        sensi_solver.set_max_steps_backward_problem(int(1e6))

        try:
            check(
                sensi_solver,
                point,
                bounds,
                multiplier != tolerance_multipliers[0],
            )
        except FunctionEvaluationError as err:
            error = err
            continue
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
