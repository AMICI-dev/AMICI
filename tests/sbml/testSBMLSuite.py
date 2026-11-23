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
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType

import amici
import diffrax
import jax
import jax.numpy as jnp
import numpy as np
import optimistix
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
    ReturnData,
    SensitivityMethod,
    SensitivityOrder,
    Solver,
    run_simulation,
)
from fiddy import FunctionEvaluationError, check_gradient, check_jacobian
from utils import (
    CaseSettings,
    apply_settings,
    find_model_file,
    read_expected_results,
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


@dataclass
class CompiledCase:
    """A test case with its compiled model."""

    model_module: ModuleType
    importer: amici.SbmlImporter
    settings: CaseSettings
    #: Directory with the test case files
    case_dir: Path


@pytest.fixture(scope="session")
def compiled_case(test_id, sbml_semantic_cases_dir):
    """Compile a case's model once, reuse simulation and sensitivity tests.

    For use with pytest-xdist, see conftest.py.
    """
    case_dir = sbml_semantic_cases_dir / test_id
    model_dir = Path(__file__).parent / "SBMLTestModels" / test_id
    try:
        model_module, sbml_importer = compile_model(
            case_dir,
            test_id,
            model_dir,
            generate_sensitivity_code=True,
        )
    except amici.importers.sbml.SBMLException as err:
        # `pytest.skip` from inside a fixture correctly propagates as
        # SKIPPED to every dependent test, not as a fixture ERROR.
        pytest.skip(str(err))

    yield CompiledCase(
        model_module=model_module,
        importer=sbml_importer,
        settings=CaseSettings.from_file(case_dir, test_id),
        case_dir=case_dir,
    )

    shutil.rmtree(model_dir, ignore_errors=True)


def _fresh_model_and_solver(
    compiled_case: CompiledCase, test_id: str
) -> tuple[Model, Solver]:
    """Build a fresh `Model`/`Solver` and set up for test."""
    model = compiled_case.model_module.get_model()
    if test_id in _NON_NEGATIVE_STATES:
        model.set_all_states_non_negative()
    solver = model.create_solver()
    apply_settings(compiled_case.settings, solver, model, test_id)
    solver.set_sensitivity_order(SensitivityOrder.first)
    solver.set_sensitivity_method(SensitivityMethod.forward)
    return model, solver


def _check_simulation_status(rdata) -> None:
    """Fail on an unsuccessful base simulation."""
    if rdata["status"] != AMICI_SUCCESS:
        raise RuntimeError("Simulation failed unexpectedly")


def test_sbml_testsuite_case(test_id, compiled_case, result_path):
    model, solver = _fresh_model_and_solver(compiled_case, test_id)

    rdata = run_simulation(model, solver)
    _check_simulation_status(rdata)

    simulated = verify_results(
        compiled_case.settings,
        rdata,
        read_expected_results(compiled_case.case_dir, test_id),
        compiled_case.importer,
        parameters=dict(
            zip(
                model.get_free_parameter_ids(),
                model.get_free_parameters(),
                strict=True,
            )
        ),
        expression_ids=model.get_expression_ids(),
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


def _sensitivity_preflight_checks(model: Model, test_id: str, check: str):
    """Skip if a sensitivity check wouldn't be meaningful for this case.

    :param check: The calling check: "forward", "adjoint", or "consistency".
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
    if test_id in _KNOWN_SENSITIVITY_CHECK_ISSUES:
        checks, reason = _KNOWN_SENSITIVITY_CHECK_ISSUES[test_id]
        if check in checks:
            pytest.skip(f"Known sensitivity check issue: {reason}")


@dataclass
class _SensitivityCase:
    """A test case set up for a sensitivity check."""

    model: Model
    settings: CaseSettings
    #: Result of the nominal simulation with forward sensitivities
    rdata: ReturnData
    has_events: bool
    #: Synthetic measurements generated from `rdata`
    edata: ExpData


def _setup_sensitivity_case(
    test_id: str, compiled_case: CompiledCase, check: str
) -> _SensitivityCase:
    """Set up a sensitivity check and run the nominal simulation.

    Skips if the check isn't meaningful for this case and fails if the
    nominal simulation fails.
    """
    model, solver = _fresh_model_and_solver(compiled_case, test_id)
    _sensitivity_preflight_checks(model, test_id, check)
    rdata = run_simulation(model, solver)
    _check_simulation_status(rdata)
    return _SensitivityCase(
        model=model,
        settings=compiled_case.settings,
        rdata=rdata,
        has_events=compiled_case.importer.sbml_model.getNumEvents() > 0,
        # `ExpData(rdata, sigma_y, sigma_z, seed)`
        edata=ExpData(rdata, 1.0, 1.0, 42),
    )


def _llh_gradient(model: Model, solver: Solver, edata: ExpData, point):
    """Return the log-likelihood function and its gradient at `point`."""
    function, derivative = run_simulation_to_function_and_derivative(
        amici_model=model,
        amici_solver=solver,
        amici_edata=edata,
        derivative_variables=["llh"],
    )
    gradient = derivative(point)["llh"]
    if np.any(np.isnan(gradient)):
        raise FunctionEvaluationError("Gradient is NaN.")
    return function, gradient


def test_sbml_testsuite_case_sensitivity_forward(test_id, compiled_case):
    """Finite-difference-check the model's forward sensitivities."""
    case = _setup_sensitivity_case(test_id, compiled_case, check="forward")

    def check(sensi_solver, point, bounds):
        sensi_solver.set_sensitivity_method(SensitivityMethod.forward)
        function, derivative = run_simulation_to_function_and_derivative(
            amici_model=case.model,
            amici_solver=sensi_solver,
            # `sigmay` is constant without an observation model
            derivative_variables=["x", "x0", "y"],
        )
        expected = derivative(point)
        result = check_jacobian(function, point, expected, bounds=bounds)
        result.assert_success(always_print=True)
        _assert_check_confirms_something(result, case.has_events)

    _run_sensitivity_check(case, test_id, check)

    # additionally cross-check against JAX autodiff for a couple of
    # historically tricky cases
    if test_id in _JAX_CHECK_CASES:
        jax_sensitivity_check(
            compiled_case.case_dir,
            test_id,
            case.model,
            case.rdata,
            case.settings.atol,
            case.settings.rtol,
        )


@pytest.mark.filterwarnings(
    # https://github.com/AMICI-dev/AMICI/issues/18
    "ignore:Adjoint sensitivity analysis for models with discontinuous "
    "right hand sides .*:UserWarning",
)
def test_sbml_testsuite_case_sensitivity_adjoint(test_id, compiled_case):
    """Finite-difference-check the model's adjoint sensitivities."""
    case = _setup_sensitivity_case(test_id, compiled_case, check="adjoint")

    def check(sensi_solver, point, bounds):
        sensi_solver.set_sensitivity_method(SensitivityMethod.adjoint)
        function, gradient = _llh_gradient(
            case.model, sensi_solver, case.edata, point
        )
        result = check_gradient(function, point, gradient, bounds=bounds)
        result.assert_success(always_print=True)
        _assert_check_confirms_something(result, case.has_events)

    _run_sensitivity_check(case, test_id, check, jitter_seed=43)


@pytest.mark.filterwarnings(
    # https://github.com/AMICI-dev/AMICI/issues/18
    "ignore:Adjoint sensitivity analysis for models with discontinuous "
    "right hand sides .*:UserWarning",
)
def test_sbml_testsuite_case_sensitivity_consistency(test_id, compiled_case):
    """Compare forward and adjoint sensitivities (`sllh`)."""
    case = _setup_sensitivity_case(test_id, compiled_case, check="consistency")

    def check(sensi_solver, point, bounds):
        asa_solver = sensi_solver.clone()
        sensi_solver.set_sensitivity_method(SensitivityMethod.forward)
        asa_solver.set_sensitivity_method(SensitivityMethod.adjoint)
        _, fsa_sllh = _llh_gradient(
            case.model, sensi_solver, case.edata, point
        )
        _, asa_sllh = _llh_gradient(case.model, asa_solver, case.edata, point)
        np.testing.assert_allclose(
            asa_sllh,
            fsa_sllh,
            rtol=1e-6,
            atol=1e-8,
            err_msg="Forward and adjoint sensitivities disagree",
        )

    _run_sensitivity_check(case, test_id, check, jitter_seed=44)


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
    # `v*_h` are Hill coefficients (nominal values 1, 2, or 4), acting on
    # species (`p1`, `p2`, ...) that start at exactly 0. The
    # forward-sensitivity RHS (not the state RHS itself) contains a
    # `species**(h-1)` term; reducing `h` even infinitesimally below 1
    # makes that exponent negative, and `0**(negative)` is `+inf`.
    # Therefore, keep `h >= 1`.
    "01395": {
        f"v{i}_h": (1.0 / nominal, None)
        for nominal, indices in (
            (1, (1, 3, 15)),
            (2, (11, 12, 13)),
            (4, (2, 4, 5, 6, 7, 8, 9, 10, 14)),
        )
        for i in indices
    },
    # The species grows like `exp(c * t**2)` with
    # `c ~ extentconv**3 / timeconv**6`: overflow at nominal * (0.5, 2).
    "01148": {
        f"{prefix}{name}": (0.9, 1.1)
        for prefix in ("", "sub1__", "sub1__sub1__")
        for name in ("timeconv", "extentconv")
    },
}

# Cases whose states are kept non-negative: slightly negative values due to
# integration errors would make non-integer powers of them NaN, so that the
# simulation fails at isolated parameter values.
#  - 01395: `(p / Kd)**h` with Hill coefficients `h` (see above).
_NON_NEGATIVE_STATES = {"01395"}


def _derive_generic_bounds(
    point: np.ndarray, param_ids: tuple[str, ...], test_id: str
) -> tuple[np.ndarray, np.ndarray]:
    """Derive per-parameter bounds for finite differences.

    Without bounds, fiddy might evaluate infeasible points. Bounds are
    nominal * (1e-6, 1000) on the same-signed side, excluding 0 (e.g., case
    00313 is singular at 0). Parameters nominally at 0 get (0, 1000): the
    free parameters in this suite are non-negative.

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
    case: _SensitivityCase,
    test_id: str,
    check: Callable[[Solver, np.ndarray, tuple[np.ndarray, np.ndarray]], None],
    jitter_seed: int | None = None,
) -> None:
    """Run a sensitivity check at tight solver tolerances.

    If a simulation fails, retry once at 100x looser tolerances.

    :param case: The test case.
    :param test_id: The SBML semantic test suite case ID.
    :param check: Called as `check(sensi_solver, point, bounds)`. Sets the
        sensitivity method, runs the check, and asserts success.
    :param jitter_seed: If given, `point` is perturbed by fixed-seed 5%
        relative Gaussian noise (then clipped to `bounds`) to avoid
        special cases at the nominal parameters.
    """
    model = case.model
    # IDAS fails at tolerances close to machine precision (forward
    # sensitivities silently become NaN, backward integration fails).
    is_dae = any(i == 0.0 for i in model.get_id_list())
    first_multiplier = 100 if is_dae else 1
    param_ids = model.get_free_parameter_ids()
    point = np.asarray(model.get_free_parameters())
    bounds = _derive_generic_bounds(point, param_ids, test_id)
    if jitter_seed is not None:
        rng = np.random.default_rng(jitter_seed)
        point = point + rng.standard_normal(len(point)) * point * 0.05
        point = np.clip(point, *bounds)
    error = None
    for multiplier in (first_multiplier, first_multiplier * 100):
        # A fresh solver per attempt; changing the tolerances of a used
        # solver fails for models without states (CV_ILL_INPUT).
        sensi_solver = model.create_solver()
        apply_settings(case.settings, sensi_solver, model, test_id)
        sensi_solver.set_absolute_tolerance(1e-16 * multiplier)
        sensi_solver.set_relative_tolerance(1e-15 * multiplier)
        sensi_solver.set_absolute_tolerance_b(1e-16 * multiplier)
        sensi_solver.set_relative_tolerance_b(1e-15 * multiplier)
        # Looser quadrature tolerances may silently corrupt the gradient
        sensi_solver.set_absolute_tolerance_quadratures(1e-16)
        sensi_solver.set_relative_tolerance_quadratures(1e-15)
        sensi_solver.set_sensitivity_order(SensitivityOrder.first)
        # The default (100 * forward max steps) lets failing backward
        # integrations run for minutes.
        sensi_solver.set_max_steps_backward_problem(int(1e6))
        try:
            check(sensi_solver, point, bounds)
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

        np.testing.assert_allclose(x, rdata["x"], rtol=rtol, atol=atol)
        np.testing.assert_allclose(sx, rdata["sx"], rtol=rtol, atol=atol)
    finally:
        shutil.rmtree(model_dir, ignore_errors=True)
