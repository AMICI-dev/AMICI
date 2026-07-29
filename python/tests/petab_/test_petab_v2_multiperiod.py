"""Tests for the JAX PEtab v2 simulator's native N-period chaining.

These exercise PEtab v2 experiments with more than two periods (pre-
equilibration + several sequential dosing/condition-switching periods),
which the JAX backend now simulates by chaining one ODE integration per
period directly, instead of collapsing them into SBML events at import
time (see ``amici.sim.jax.petab.JAXProblem``).
"""

import jax
import numpy as np
import pytest
from petab.v2 import C, Problem, ProblemConfig
from petab.v2.models.sbml_model import SbmlModel

from amici.importers.petab import PetabImporter
from amici.sim.jax import ReturnValue, run_simulations

jax.config.update("jax_enable_x64", True)


def _linear_decay_problem() -> Problem:
    """A single-species linear-decay model (dx/dt = -k*x)."""
    problem = Problem()
    problem.config = ProblemConfig()
    problem.model = SbmlModel.from_antimony("xx = 1; xx' = -kk*xx;")
    problem.add_observable("obs1", "xx", noise_formula="1")
    return problem


def _import_jax(problem: Problem, module_name: str, tmp_path):
    pi = PetabImporter(
        petab_problem=problem,
        module_name=module_name,
        output_dir=tmp_path / module_name,
        compile_=True,
        jax=True,
        verbose=False,
    )
    return pi.create_simulator(force_import=True)


def test_three_period_chain_matches_analytical_solution(tmp_path):
    """A pre-equilibration followed by three sequential dosing periods,
    with measurements split across the 2nd and 3rd periods, must match a
    closed-form (segment-wise exponential decay) reference solution."""
    problem = _linear_decay_problem()
    problem.add_condition("cond_preeq", kk=0.5)
    # period 1: t in [0, 2), k=0.3, dosed to xx=3.0 at t=0
    problem.add_condition("cond_p1", kk=0.3, xx=3.0)
    # period 2: t in [2, 4), k=0.6, no reinit (state carries over from p1)
    problem.add_condition("cond_p2", kk=0.6)
    # period 3: t >= 4, k=0.2, dosed again to xx=1.5
    problem.add_condition("cond_p3", kk=0.2, xx=1.5)
    problem.add_experiment(
        "exp1",
        C.TIME_PREEQUILIBRATION,
        "cond_preeq",
        0.0,
        "cond_p1",
        2.0,
        "cond_p2",
        4.0,
        "cond_p3",
    )
    measurement_times = (0.5, 1.5, 2.5, 3.5, 4.5, 5.5)
    for t in measurement_times:
        problem.add_measurement(
            "obs1", time=t, measurement=0.0, experiment_id="exp1"
        )

    jax_problem = _import_jax(problem, "test_three_period_chain", tmp_path)
    assert jax_problem._max_periods == 3

    x, _ = run_simulations(jax_problem, ret=ReturnValue.x)
    ts_mask = np.asarray(jax_problem._ts_masks)[0].reshape(-1)
    actual = np.asarray(x)[0].reshape(-1)[ts_mask]

    x2 = 3.0 * np.exp(-0.3 * 2.0)

    def analytical(t):
        if t < 2.0:
            return 3.0 * np.exp(-0.3 * t)
        if t < 4.0:
            return x2 * np.exp(-0.6 * (t - 2.0))
        return 1.5 * np.exp(-0.2 * (t - 4.0))

    expected = np.array([analytical(t) for t in measurement_times])
    np.testing.assert_allclose(actual, expected, rtol=1e-4)

    llh, _ = run_simulations(jax_problem, ret=ReturnValue.llh)
    assert np.isfinite(llh)


def test_gradient_through_multiperiod_chain_matches_analytical_derivative(
    tmp_path,
):
    """Gradients of the log-likelihood w.r.t. an estimated parameter that
    is referenced from a non-first period's condition table must flow
    correctly through the chained periods.

    The reference gradient is derived analytically rather than by finite
    differences: within a segment, ``xx(t) = xx0 * exp(-k*(t - t0))``, and
    with a unit-sigma Gaussian noise model the per-measurement
    log-likelihood contribution is ``-0.5*log(2*pi) - 0.5*(m-y)^2``, so
    ``d(llh)/dk = sum_i (m_i - y_i) * dy_i/dk``. Period 3 reinitialises
    ``xx`` to a literal value, so it -- and its measurements -- no longer
    depend on ``k_free`` at all, giving a zero gradient contribution.
    """
    problem = _linear_decay_problem()
    problem.add_parameter(
        "k_free", nominal_value=0.3, estimate=True, lb=0.01, ub=2
    )
    problem.add_condition("cond_preeq", kk=0.5)
    problem.add_condition("cond_p1", kk="k_free", xx=3.0)
    problem.add_condition("cond_p2", kk=0.6)
    problem.add_condition("cond_p3", kk=0.2, xx=1.5)
    problem.add_experiment(
        "exp1",
        C.TIME_PREEQUILIBRATION,
        "cond_preeq",
        0.0,
        "cond_p1",
        2.0,
        "cond_p2",
        4.0,
        "cond_p3",
    )
    measurement_times = (0.5, 1.5, 2.5, 3.5, 4.5, 5.5)
    measurement_value = 1.0
    for t in measurement_times:
        problem.add_measurement(
            "obs1",
            time=t,
            measurement=measurement_value,
            experiment_id="exp1",
        )

    jax_problem = _import_jax(problem, "test_multiperiod_gradient", tmp_path)

    def llh_fn(p):
        return run_simulations(jax_problem.update_parameters(p))[0]

    p0 = jax_problem.parameters
    grad = jax.grad(llh_fn)(p0)

    k = float(p0[0])
    x_p1_end = 3.0 * np.exp(-k * 2.0)  # xx at t=2, end of period 1
    dx_p1_end_dk = -2.0 * x_p1_end

    def y_and_dy_dk(t: float) -> tuple[float, float]:
        if t < 2.0:
            y = 3.0 * np.exp(-k * t)
            return y, -t * y
        if t < 4.0:
            y = x_p1_end * np.exp(-0.6 * (t - 2.0))
            return y, dx_p1_end_dk * np.exp(-0.6 * (t - 2.0))
        # period 3 reinitialises xx to a literal, independent of k_free
        return 1.5 * np.exp(-0.2 * (t - 4.0)), 0.0

    expected_grad = sum(
        (measurement_value - y) * dy_dk
        for y, dy_dk in map(y_and_dy_dk, measurement_times)
    )
    np.testing.assert_allclose(float(grad[0]), expected_grad, rtol=1e-4)


def _threshold_piecewise_decay_problem() -> Problem:
    """A single-species model whose decay rate depends on whether the
    species concentration is above or below a threshold, via a
    ``piecewise`` rate law. The JAX backend compiles this to a
    root-finding/heaviside event, rather than a state-assignment event."""
    problem = Problem()
    problem.config = ProblemConfig()
    problem.model = SbmlModel.from_antimony(
        "xx = 1; kfast = 1.0; "
        "xx' = piecewise(-kfast*xx, xx > 2, -0.1*xx);"
    )
    problem.add_observable("obs1", "xx", noise_formula="1")
    return problem


def test_event_heaviside_state_reevaluated_after_period_reinit(
    tmp_path,
):
    """A period-boundary reinitialisation that crosses the threshold of a
    ``piecewise`` rate law must select the branch matching the
    *reinitialised* state at the new period's t0, not whatever heaviside
    state the previous period happened to end with.

    ``JAXModel._handle_t0_event`` re-evaluates the trigger condition
    against the actual incoming state at every period boundary (and after
    preequilibration), using the previous heaviside state only as the
    pre-transition reference for detecting a crossing -- it does not carry
    it over unconditionally.
    """
    problem = _threshold_piecewise_decay_problem()
    # period 1: xx starts at 5 (above the threshold of 2) and decays past
    # it, ending (at t=1) below the threshold.
    problem.add_condition("cond_p1", xx=5.0)
    # period 2 reinitialises xx to 3, back above the threshold: the
    # trigger must be re-evaluated at t0 of period 2 rather than reusing
    # period 1's ending ("below threshold") heaviside state.
    problem.add_condition("cond_p2", xx=3.0)
    problem.add_experiment("exp1", 0.0, "cond_p1", 1.0, "cond_p2")
    measurement_times = (0.3, 0.8, 1.3, 1.8)
    for t in measurement_times:
        problem.add_measurement(
            "obs1", time=t, measurement=0.0, experiment_id="exp1"
        )

    jax_problem = _import_jax(
        problem, "test_event_reinit_heaviside_reevaluated", tmp_path
    )
    assert jax_problem._max_periods == 2

    x, _ = run_simulations(jax_problem, ret=ReturnValue.x)
    ts_mask = np.asarray(jax_problem._ts_masks)[0].reshape(-1)
    actual = np.asarray(x)[0].reshape(-1)[ts_mask]

    # period 1 crosses the threshold via root-finding during integration
    # (heaviside starts "above" since xx=5 > 2 at t=0).
    t_cross1 = np.log(5.0 / 2.0)

    def period1(t):
        if t < t_cross1:
            return 5.0 * np.exp(-1.0 * t)
        return 2.0 * np.exp(-0.1 * (t - t_cross1))

    # period 2 must also start "above" (heaviside re-evaluated against the
    # reinitialised xx=3.0 at period 2's t0), decaying fast until it
    # crosses the threshold again, then slow.
    t_cross2 = np.log(3.0 / 2.0)

    def period2(t_local):
        if t_local < t_cross2:
            return 3.0 * np.exp(-1.0 * t_local)
        return 2.0 * np.exp(-0.1 * (t_local - t_cross2))

    expected = np.array(
        [
            period1(0.3),
            period1(0.8),
            period2(0.3),
            period2(0.8),
        ]
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-4)


@pytest.mark.parametrize("jax_flag", [True])
def test_petab_importer_skips_event_conversion_for_jax(tmp_path, jax_flag):
    """The JAX backend must not run ExperimentsToSbmlConverter (the
    experiments-to-events conversion), even for experiments with more than
    two periods, since it now chains periods natively."""
    problem = _linear_decay_problem()
    problem.add_condition("cond_preeq", kk=0.5)
    problem.add_condition("cond_p1", kk=0.3, xx=3.0)
    problem.add_condition("cond_p2", kk=0.6)
    problem.add_condition("cond_p3", kk=0.2, xx=1.5)
    problem.add_experiment(
        "exp1",
        C.TIME_PREEQUILIBRATION,
        "cond_preeq",
        0.0,
        "cond_p1",
        2.0,
        "cond_p2",
        4.0,
        "cond_p3",
    )
    problem.add_measurement("obs1", time=0.5, measurement=0.0, experiment_id="exp1")

    pi = PetabImporter(
        petab_problem=problem,
        module_name="test_no_event_conversion",
        output_dir=tmp_path / "test_no_event_conversion",
        compile_=True,
        jax=jax_flag,
        verbose=False,
    )
    # experiment periods are untouched (still 4: preeq + 3 dosing periods)
    assert len(pi.petab_problem.experiments[0].periods) == 4
    assert pi._unconverted_problem is None
