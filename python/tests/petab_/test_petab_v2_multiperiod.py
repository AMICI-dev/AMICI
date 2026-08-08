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


def _two_species_decay_problem() -> Problem:
    """Two independently decaying species, ``AA' = -kA*AA``, ``BB' = -kB*BB``."""
    problem = Problem()
    problem.config = ProblemConfig()
    problem.model = SbmlModel.from_antimony(
        "AA = 1; BB = 2; kA = 0.5; kB = 0.25; AA' = -kA*AA; BB' = -kB*BB;"
    )
    problem.add_observable("obs_a", "AA", noise_formula="1")
    problem.add_observable("obs_b", "BB", noise_formula="1")
    return problem


def _state_from_state_problem(
    reinit_expr: str = "AA + scale*BB",
) -> Problem:
    """A two-period problem whose second period initialises ``AA`` from the
    *simulated* value of both ``AA`` and ``BB`` at the period boundary."""
    problem = _two_species_decay_problem()
    problem.add_parameter(
        "scale", nominal_value=5.0, estimate=True, lb=0.1, ub=10
    )
    problem.add_condition("c_p1", AA=1.0, BB=2.0)
    # the previously-unsupported case: a state's (re)initial value is an
    # expression in other states' simulated values, mixed with an estimated
    # PEtab parameter that is not a model parameter
    problem.add_condition("c_p2", AA=reinit_expr)
    problem.add_experiment("exp1", 0.0, "c_p1", 1.0, "c_p2")
    return problem


#: measurement times of :func:`_state_from_state_problem`, one in the first
#: period and two in the second
_STATE_FROM_STATE_TIMES = (0.5, 1.5, 2.5)


def _state_from_state_analytical(
    scale: float, t: float
) -> tuple[float, float]:
    """``(AA(t), BB(t))`` of :func:`_state_from_state_problem`, segment-wise.

    ``AA`` and ``BB`` decay independently at rates 0.5 and 0.25 from
    ``AA(0)=1``/``BB(0)=2``. At ``t=1`` the second period's condition replaces
    ``AA`` by ``AA(1) + scale*BB(1)``, leaving ``BB`` untouched.
    """
    if t < 1.0:
        return np.exp(-0.5 * t), 2.0 * np.exp(-0.25 * t)
    aa1 = np.exp(-0.5)
    bb1 = 2.0 * np.exp(-0.25)
    return (
        (aa1 + scale * bb1) * np.exp(-0.5 * (t - 1.0)),
        bb1 * np.exp(-0.25 * (t - 1.0)),
    )


def test_period_condition_initialises_state_from_another_state(tmp_path):
    """A non-first period's condition may set a state from other states'
    *simulated* values (``AA = AA + scale*BB``).

    This is only possible because the reinitialisation expression is
    generated into the model file (``JAXModel._x_reinit``) and evaluated at
    the period boundary against the state actually reached there. Resolving
    it in Python ahead of the simulation -- as the PEtab layer used to --
    raised ``NotImplementedError``, since no simulated state exists yet at
    that point.
    """
    problem = _state_from_state_problem()
    for t in _STATE_FROM_STATE_TIMES:
        for obs in ("obs_a", "obs_b"):
            problem.add_measurement(
                obs, time=t, measurement=0.0, experiment_id="exp1"
            )

    jax_problem = _import_jax(problem, "test_reinit_state_ref", tmp_path)
    assert jax_problem._max_periods == 2

    # the reinitialisation expression really is in the generated model file,
    # not resolved by the PEtab layer
    assert jax_problem.model.reinitialisation_targets == (
        ("c_p1", "AA", "1.00000000000000"),
        ("c_p1", "BB", "2.00000000000000"),
        ("c_p2", "AA", "AA + BB*scale"),
    )
    assert jax_problem.model.reinitialisation_parameter_ids == ("scale",)
    assert "AA + BB*scale" in jax_problem.model.jax_py_file.read_text()

    x, _ = run_simulations(jax_problem, ret=ReturnValue.x)
    ts_mask = np.asarray(jax_problem._ts_masks)[0].reshape(-1)
    actual = np.asarray(x)[0].reshape(-1, len(jax_problem.model.state_ids))
    actual = actual[ts_mask]

    scale = float(
        jax_problem.parameters[jax_problem.parameter_ids.index("scale")]
    )
    i_aa = list(jax_problem.model.state_ids).index("AA")
    i_bb = list(jax_problem.model.state_ids).index("BB")
    # two measurements (obs_a, obs_b) share each timepoint
    expected = np.array(
        [
            _state_from_state_analytical(scale, t)
            for t in _STATE_FROM_STATE_TIMES
            for _ in range(2)
        ]
    )
    np.testing.assert_allclose(actual[:, i_aa], expected[:, 0], rtol=1e-6)
    np.testing.assert_allclose(actual[:, i_bb], expected[:, 1], rtol=1e-6)


def test_gradient_through_state_referencing_reinitialisation(tmp_path):
    """The reinitialisation ``AA = AA + scale*BB`` must stay differentiable
    w.r.t. the estimated parameter ``scale`` appearing in it.

    ``scale`` is a PEtab parameter that is *not* a model parameter (it only
    ever appears in a condition-table ``targetValue``), so it reaches the
    generated ``_x_reinit`` through the model's
    ``reinitialisation_parameter_ids`` channel. With a unit-sigma Gaussian
    noise model, ``d(llh)/d(scale) = sum_i (m_i - y_i) * dy_i/d(scale)``, and
    only the second period's ``obs_a`` measurements depend on ``scale``, with
    ``dAA/d(scale) = BB(1) * exp(-0.5*(t-1))``.
    """
    problem = _state_from_state_problem()
    measurement_value = 1.0
    for t in _STATE_FROM_STATE_TIMES:
        problem.add_measurement(
            "obs_a",
            time=t,
            measurement=measurement_value,
            experiment_id="exp1",
        )

    jax_problem = _import_jax(problem, "test_reinit_state_ref_grad", tmp_path)
    i_scale = jax_problem.parameter_ids.index("scale")

    def llh_fn(p):
        return run_simulations(jax_problem.update_parameters(p))[0]

    p0 = jax_problem.parameters
    grad = jax.grad(llh_fn)(p0)

    scale = float(p0[i_scale])
    bb1 = 2.0 * np.exp(-0.25)
    expected_grad = 0.0
    for t in _STATE_FROM_STATE_TIMES:
        if t < 1.0:
            continue  # first period is independent of `scale`
        y = _state_from_state_analytical(scale, t)[0]
        dy_dscale = bb1 * np.exp(-0.5 * (t - 1.0))
        expected_grad += (measurement_value - y) * dy_dscale
    np.testing.assert_allclose(float(grad[i_scale]), expected_grad, rtol=1e-5)

    # cross-check against central finite differences
    eps = 1e-6
    fd = (
        float(llh_fn(p0.at[i_scale].add(eps)))
        - float(llh_fn(p0.at[i_scale].add(-eps)))
    ) / (2 * eps)
    np.testing.assert_allclose(float(grad[i_scale]), fd, rtol=1e-5, atol=1e-6)


def test_stale_model_condition_table_mismatch_is_detected(tmp_path):
    """Editing a ``targetValue`` and reusing a cached model directory must
    fail loudly.

    Reinitialisation values are generated into the model file, so a model
    generated from a different condition table would otherwise silently keep
    simulating the value it was generated with.
    """
    problem = _state_from_state_problem()
    for t in _STATE_FROM_STATE_TIMES:
        problem.add_measurement(
            "obs_a", time=t, measurement=0.0, experiment_id="exp1"
        )
    jax_problem = _import_jax(problem, "test_reinit_stale", tmp_path)

    from amici.sim.jax.petab import JAXProblem

    # same (already generated) model, but the condition table now says
    # something else
    edited = _state_from_state_problem("AA + 2*scale*BB")
    for t in _STATE_FROM_STATE_TIMES:
        edited.add_measurement(
            "obs_a", time=t, measurement=0.0, experiment_id="exp1"
        )

    with pytest.raises(ValueError, match="must be regenerated"):
        JAXProblem(jax_problem.model, edited)


def test_reinitialisation_referencing_assignment_rule_is_rejected(tmp_path):
    """A ``targetValue`` referencing an assignment-rule target cannot be
    represented and must be rejected at code generation time.

    Reinitialisation happens at a period boundary *before* the period's
    expressions (``w``) and conservation laws (``tcl``) are computed -- they
    are computed *from* the reinitialised state -- so the generated
    ``_x_reinit`` has no way to evaluate one.
    """
    problem = Problem()
    problem.config = ProblemConfig()
    problem.model = SbmlModel.from_antimony(
        "AA = 1; kA = 0.5; scaled := 3*kA; AA' = -kA*AA;"
    )
    problem.add_observable("obs_a", "AA", noise_formula="1")
    problem.add_condition("c_p1", AA=1.0)
    problem.add_condition("c_p2", AA="scaled")
    problem.add_experiment("exp1", 0.0, "c_p1", 1.0, "c_p2")
    problem.add_measurement(
        "obs_a", time=0.5, measurement=0.0, experiment_id="exp1"
    )

    with pytest.raises(NotImplementedError, match="assignment rule"):
        _import_jax(problem, "test_reinit_w_ref", tmp_path)


def _threshold_piecewise_decay_problem() -> Problem:
    """A single-species model whose decay rate depends on whether the
    species concentration is above or below a threshold, via a
    ``piecewise`` rate law. The JAX backend compiles this to a
    root-finding/heaviside event, rather than a state-assignment event."""
    problem = Problem()
    problem.config = ProblemConfig()
    problem.model = SbmlModel.from_antimony(
        "xx = 1; kfast = 1.0; xx' = piecewise(-kfast*xx, xx > 2, -0.1*xx);"
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
    problem.add_measurement(
        "obs1", time=0.5, measurement=0.0, experiment_id="exp1"
    )

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


def test_posteq_runs_at_own_last_period_of_a_short_experiment(tmp_path):
    """A ragged problem: post-equilibration must be applied per experiment.

    ``exp_long`` has two dynamic periods and so sets the problem-wide
    period count; ``exp_short`` has one, plus a ``time=inf`` measurement.
    The steady state of ``xx' = -kk*xx`` is ``0``, so a post-equilibrated
    value is analytically distinguishable from the end-of-dynamics value.
    """
    problem = _linear_decay_problem()
    problem.add_condition("c_a", kk=0.5)
    problem.add_condition("c_b", kk=0.9)

    # two dynamic periods -> problem-wide maximum is 2
    problem.add_experiment("exp_long", 0.0, "c_a", 2.0, "c_b")
    for t in (0.5, 2.5):
        problem.add_measurement(
            "obs1", time=t, measurement=0.0, experiment_id="exp_long"
        )

    # a single dynamic period, i.e. fewer than the problem-wide maximum,
    # carrying a post-equilibration measurement
    problem.add_experiment("exp_short", 0.0, "c_a")
    problem.add_measurement(
        "obs1", time=1.0, measurement=0.0, experiment_id="exp_short"
    )
    problem.add_measurement(
        "obs1", time=np.inf, measurement=0.0, experiment_id="exp_short"
    )

    jax_problem = _import_jax(problem, "test_ragged_posteq", tmp_path)
    assert jax_problem._max_periods == 2

    x, _ = run_simulations(jax_problem, ret=ReturnValue.x)
    exp_ids = [e.id for e in jax_problem._petab_problem.experiments]
    i_short = exp_ids.index("exp_short")
    mask = np.asarray(jax_problem._ts_masks)[i_short].reshape(-1)
    actual = np.asarray(x)[i_short].reshape(-1)[mask]

    # dynamic measurement at t=1.0, then the post-equilibrated steady state
    np.testing.assert_allclose(actual[0], np.exp(-0.5 * 1.0), rtol=1e-4)
    np.testing.assert_allclose(actual[-1], 0.0, atol=1e-6)


def test_simulation_experiments_can_select_a_subset(tmp_path):
    """``simulation_experiments`` must accept a proper subset.

    The measurement arrays are built once for every experiment of the
    problem, while the parameter arrays are built only for the experiments
    being simulated. Passing both to the same ``vmap`` only lined up when
    the "subset" happened to be all of them; any real subset raised
    ``vmap got inconsistent sizes for array axes to be mapped``.
    """
    problem = _linear_decay_problem()
    problem.add_condition("c1", kk=0.5)
    problem.add_condition("c2", kk=0.9)
    for exp_id, cond_id in (("e1", "c1"), ("e2", "c2")):
        problem.add_experiment(exp_id, 0.0, cond_id)
        for t in (0.5, 1.5):
            problem.add_measurement(
                "obs1", time=t, measurement=0.0, experiment_id=exp_id
            )

    jax_problem = _import_jax(problem, "test_experiment_subset", tmp_path)

    llh_all, _ = run_simulations(jax_problem, ret=ReturnValue.llh)
    llh_e1, _ = run_simulations(
        jax_problem, simulation_experiments=["e1"], ret=ReturnValue.llh
    )
    llh_e2, _ = run_simulations(
        jax_problem, simulation_experiments=["e2"], ret=ReturnValue.llh
    )

    # experiments contribute independently, so the parts must sum to the
    # whole -- a subset that merely runs but mixes up which experiment's
    # measurements it scores would not satisfy this
    np.testing.assert_allclose(
        float(llh_e1) + float(llh_e2), float(llh_all), rtol=1e-8
    )
    assert float(llh_e1) != float(llh_e2)

    # selection is by id, not by position
    llh_reversed, _ = run_simulations(
        jax_problem, simulation_experiments=["e2", "e1"], ret=ReturnValue.llh
    )
    np.testing.assert_allclose(float(llh_reversed), float(llh_all), rtol=1e-8)


def test_run_simulations_under_filter_jit(tmp_path):
    """``run_simulations`` must work inside ``eqx.filter_jit``.

    ``JAXProblem`` is an ``eqx.Module``, so its array fields -- including
    the measurement masks -- become tracers under ``filter_jit``. Reading
    concrete values out of them (rather than just their shapes) raises
    ``TracerArrayConversionError``, which plain (non-jitted) calls never
    surface. This is the pattern the JAX PEtab example notebook uses to
    build an optimisation step.
    """
    import equinox as eqx

    problem = _linear_decay_problem()
    problem.add_parameter(
        "k_free", nominal_value=0.4, estimate=True, lb=0.01, ub=10
    )
    # an estimated parameter drives the first period, so the gradient is
    # non-trivial and must reach `problem.parameters`
    problem.add_condition("c_p1", kk="k_free")
    problem.add_condition("c_p2", kk=0.8)
    problem.add_experiment("e1", 0.0, "c_p1", 2.0, "c_p2")
    for t in (0.5, 1.5, 2.5, 3.5):
        problem.add_measurement(
            "obs1", time=t, measurement=0.0, experiment_id="e1"
        )
    # a post-equilibration row, so the post-equilibration bookkeeping is
    # exercised too
    problem.add_measurement(
        "obs1", time=np.inf, measurement=0.0, experiment_id="e1"
    )

    jax_problem = _import_jax(problem, "test_filter_jit", tmp_path)

    @eqx.filter_jit
    def llh_and_grad(p):
        return eqx.filter_value_and_grad(
            lambda q: run_simulations(q, ret=ReturnValue.llh)[0]
        )(p)

    llh, grad = llh_and_grad(jax_problem)
    assert np.isfinite(float(llh))
    # the gradient must actually reach the estimated parameters
    assert np.any(np.asarray(grad.parameters) != 0.0)


def test_measurement_row_indices_stay_integral(tmp_path):
    """Measurement row indices must not be widened to float.

    The per-period arrays are assembled by padding a dynamic and a
    post-equilibration part and concatenating them. An absent
    post-equilibration part is an empty array, which numpy types as
    float64 regardless of what the data is -- so concatenating it against
    the integer row indices silently produced a float index, which no
    longer compares equal to the measurement table's integer index.
    """
    problem = _linear_decay_problem()
    problem.add_condition("c_p1", kk=0.4)
    problem.add_condition("c_p2", kk=0.8)
    problem.add_experiment("e1", 0.0, "c_p1", 2.0, "c_p2")
    for t in (0.5, 1.5, 2.5):
        problem.add_measurement(
            "obs1", time=t, measurement=0.0, experiment_id="e1"
        )

    jax_problem = _import_jax(problem, "test_index_dtype", tmp_path)

    assert jax_problem._petab_measurement_indices.dtype.kind in "iu", (
        "measurement row indices must stay integral, got "
        f"{jax_problem._petab_measurement_indices.dtype}"
    )
