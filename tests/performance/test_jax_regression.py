"""JAX regression test suite.

Each test exercises one (model, operation) pair and records:
- llh + solver step counts (deterministic, primary regression signal)
- wall-clock execution time (secondary, soft signal)

Results are written to jax_regression_results.json when --results-dir is
passed.
"""

import diffrax
import equinox as eqx
import jax
import jax.numpy as jnp
import pytest

pytest.importorskip("jax")


from tests.performance._utils import measure_exec_time

# ── Tier 1 test cases (synthetic models) ──────────────────────────────────

TIER1_FWD_CASES = [
    "LinearDecay",
    "Robertson",
    "LotkaVolterra",
    "ConservationLaw",
    "SingleEvent",
    "MultiEvent",
]

# ── Helper: build simulate_condition_unjitted kwargs for each model ─────────


def _sim_kwargs(model, solver_kwargs) -> dict:
    """Return keyword arguments for simulate_condition_unjitted."""
    import tests.performance.synthetic_models.conservation_law as cl
    import tests.performance.synthetic_models.linear_decay as ld
    import tests.performance.synthetic_models.lotka_volterra as lv
    import tests.performance.synthetic_models.multi_event as me
    import tests.performance.synthetic_models.robertson as rob
    import tests.performance.synthetic_models.single_event as se

    _MAP = {
        "LinearDecay": (
            ld.TS_DYN,
            ld.MY,
            ld.IYS,
            ld.IY_TRAFOS,
            ld.OPS,
            ld.NPS,
        ),
        "Robertson": (
            rob.TS_DYN,
            rob.MY,
            rob.IYS,
            rob.IY_TRAFOS,
            rob.OPS,
            rob.NPS,
        ),
        "LotkaVolterra": (
            lv.TS_DYN,
            lv.MY,
            lv.IYS,
            lv.IY_TRAFOS,
            lv.OPS,
            lv.NPS,
        ),
        "ConservationLaw": (
            cl.TS_DYN,
            cl.MY,
            cl.IYS,
            cl.IY_TRAFOS,
            cl.OPS,
            cl.NPS,
        ),
        "SingleEvent": (
            se.TS_DYN,
            se.MY,
            se.IYS,
            se.IY_TRAFOS,
            se.OPS,
            se.NPS,
        ),
        "MultiEvent": (me.TS_DYN, me.MY, me.IYS, me.IY_TRAFOS, me.OPS, me.NPS),
    }

    model_name = type(model).__name__
    ts_dyn, my, iys, iy_trafos, ops, nps = _MAP[model_name]

    # `simulate_condition`/`simulate_condition_unjitted` chain one ODE
    # integration per experiment period; add a leading period axis of
    # size 1 since these synthetic models are all single-period.
    return dict(
        ts_dyn=ts_dyn[None, :],
        ts_posteq=jnp.zeros((1, 0)),
        my=my[None, :],
        iys=iys[None, :],
        iy_trafos=iy_trafos[None, :],
        ops=ops[None, :, :],
        nps=nps[None, :, :],
        **solver_kwargs,
    )


def _extract_stats(stats: dict) -> dict:
    """Pull step counts out of the stats dict returned by simulate_condition."""
    out = {}
    for key in ("stats_dyn", "stats_posteq"):
        s = stats.get(key)
        # `stats_dyn` is a list with one entry per experiment period (see
        # JAXModel._simulate_period); these synthetic models are all
        # single-period, so use the (only) entry.
        if isinstance(s, list):
            s = next((entry for entry in reversed(s) if entry is not None), None)
        if s is None:
            out[key] = None
        else:
            out[key] = {
                k: int(v)
                for k, v in s.items()
                if k
                in ("num_accepted_steps", "num_rejected_steps", "num_steps")
            }
    return out


# ── Forward simulation ─────────────────────────────────────────────


@pytest.mark.parametrize("model_id", TIER1_FWD_CASES)
def test_tier1_fwd_sim(
    model_id, tier1_models, solver_kwargs, results_collector
):
    model = tier1_models[model_id]
    # add a leading period axis of size 1 (single-period simulation)
    p = model.parameters[None, :]
    kwargs = _sim_kwargs(model, solver_kwargs)

    # Deterministic run (unjitted, for exact step counts)
    llh, stats = model.simulate_condition_unjitted(p, **kwargs)

    # Timing (JIT-compiled path)
    sim_fn = model.simulate_condition
    t_first, t_exec = measure_exec_time(sim_fn, p, **kwargs)

    results_collector.add(
        model_id,
        "fwd_sim",
        {
            "llh": float(llh),
            **_extract_stats(stats),
            "t_first_s": t_first,
            "t_exec_s": t_exec,
        },
    )


# ── Adjoint (gradient) ────────────────────────────────────────────


@pytest.mark.parametrize("model_id", TIER1_FWD_CASES)
def test_tier1_adj(model_id, tier1_models, solver_kwargs, results_collector):
    model = tier1_models[model_id]
    # add a leading period axis of size 1 (single-period simulation)
    p = model.parameters[None, :]
    kwargs = _sim_kwargs(model, solver_kwargs)

    def _fn(p):
        return model.simulate_condition(p, **kwargs)

    t_first, t_exec = measure_exec_time(
        eqx.filter_value_and_grad(_fn, has_aux=True), p
    )

    results_collector.add(
        model_id,
        "adj",
        {"t_first_s": t_first, "t_exec_s": t_exec},
    )


# ── Forward sensitivities (small models only) ─────────────────────


@pytest.mark.parametrize(
    "model_id", ["LinearDecay", "ConservationLaw", "SingleEvent"]
)
def test_tier1_fwd_sens(
    model_id, tier1_models, solver_kwargs, results_collector
):
    model = tier1_models[model_id]
    # add a leading period axis of size 1 (single-period simulation)
    p = model.parameters[None, :]
    # Forward sensitivity requires DirectAdjoint
    kwargs = {
        **_sim_kwargs(model, solver_kwargs),
        "adjoint": diffrax.DirectAdjoint(),
    }

    def _fn(p):
        return model.simulate_condition(p, **kwargs)

    t_first, t_exec = measure_exec_time(jax.jacfwd(_fn, has_aux=True), p)

    results_collector.add(
        model_id,
        "fwd_sens",
        {"t_first_s": t_first, "t_exec_s": t_exec},
    )
