"""Simulation subroutines."""

# ruff: noqa: F821 F722

import os
from collections.abc import Callable

import diffrax
import equinox.internal as eqxi
import jax
import jax.numpy as jnp
import jax.tree_util as jtu
import jaxtyping as jt
from optimistix import AbstractRootFinder

STARTING_STATS = {
    "max_steps": 0,
    "num_accepted_steps": 0,
    "num_rejected_steps": 0,
    "num_steps": 0,
}


def eq(
    p: jt.Float[jt.Array, "np"],
    tcl: jt.Float[jt.Array, "ncl"],
    h0: jt.Float[jt.Array, "ne"],
    x0: jt.Float[jt.Array, "nxs"],
    h_mask: jt.Bool[jt.Array, "ne"],
    solver: diffrax.AbstractSolver,
    controller: diffrax.AbstractStepSizeController,
    root_finder: AbstractRootFinder,
    steady_state_event: Callable[..., diffrax._custom_types.BoolScalarLike],
    term: diffrax.ODETerm,
    root_cond_fns: list[Callable],
    root_cond_fn: Callable,
    delta_x: Callable,
    known_discs: jt.Float[jt.Array, "*nediscs"],
    max_steps: jnp.int_,
) -> tuple[jt.Float[jt.Array, "nxs"], jt.Float[jt.Array, "ne"], dict]:
    """
    Simulate the ODE system until steady state.

    :param p:
        parameters
    :param tcl:
        total values for conservation laws
    :param h:
        heaviside variables
    :param x0:
        initial state vector
    :param solver:
        ODE solver
    :param controller:
        step size controller
    :param root_finder:
        root finder for discontinuities
    :param steady_state_event:
        event function for steady state
    :param term:
        ODE term
    :param root_cond_fns:
        list of individual root condition functions for discontinuities
    :param root_cond_fn:
        root condition function for all discontinuities
    :param delta_x:
        function to compute state changes at events
    :param known_discs:
        known discontinuities, used to clip the step size controller
    :param max_steps:
        maximum number of steps
    :return:
        steady state solution, heaviside variables, and statistics
    """
    # if there are no events, we can avoid expensive looping and just run a single segment
    if not root_cond_fns:
        sol, _, stats = _run_segment(
            0.0,
            jnp.inf,
            x0,
            p,
            tcl,
            h0,
            solver,
            controller,
            root_finder,
            max_steps,
            diffrax.DirectAdjoint(),
            [steady_state_event],
            [None],
            diffrax.SaveAt(t1=True),
            term,
            dict(**STARTING_STATS),
        )
        y1 = jnp.where(
            diffrax.is_event(sol.result),
            sol.ys[-1],
            jnp.inf * jnp.ones_like(sol.ys[-1]),
        )
        return y1, h0, stats

    def cond_fn(carry):
        _, y0, _, event_index, _ = carry
        return jnp.logical_and(
            event_index != 0,  # has not reached steady state yet
            jnp.isfinite(
                y0
            ).all(),  # y0 is finite, also used to track integration failure
        )

    def body_fn(carry):
        t_start, y0, h, _, _ = carry

        # skip the solve if already at steady state
        # TODO: remove once https://github.com/patrick-kidger/diffrax/issues/720 is resolved?
        at_steady_state = steady_state_event(
            t_start, 
            y0, 
            (p,tcl,h), 
            terms=term, 
            solver=solver, 
            stepsize_controller=controller
        )

        def at_steady_state_fn(carry):
            return carry

        def not_at_steady_state_fn(carry):
            t_start, y0, h, event_index, stats = carry

            sol, event_index, stats = _run_segment(
                t_start,
                jnp.inf,
                y0,
                p,
                tcl,
                h,
                solver,
                controller,
                root_finder,
                max_steps,  # TODO: figure out how to pass `max_steps - stats['num_steps']` here
                diffrax.DirectAdjoint(),
                [steady_state_event] + root_cond_fns,
                [None] + [True] * len(root_cond_fns),
                diffrax.SaveAt(t1=True),
                term,
                stats,
            )
            y0_next = jnp.where(
                jnp.logical_or(
                    diffrax.is_successful(sol.result),
                    diffrax.is_event(sol.result),
                ),
                sol.ys[-1],
                jnp.inf * jnp.ones_like(sol.ys[-1]),
            )
            t0_next = jnp.where(jnp.isfinite(sol.ts), sol.ts, -jnp.inf).max()

            y0_next, h_next, stats = _handle_event(
                t0_next,
                y0_next,
                p,
                tcl,
                h,
                root_finder,
                term,
                root_cond_fn,
                delta_x,
                h_mask,
                stats,
            )

            return (
                t0_next,
                y0_next,
                h_next,
                event_index,
                dict(**STARTING_STATS),
            )
        
        return jax.lax.cond(at_steady_state, at_steady_state_fn, not_at_steady_state_fn, carry)

    # run the loop until no event is triggered (which will also be the case if we run out of steps)
    _, y1, h, _, stats = eqxi.while_loop(
        cond_fn,
        body_fn,
        (0.0, x0, h0, -1, dict(**STARTING_STATS)),
        kind="bounded",
        max_steps=2**6,
    )

    return y1, h, stats


def solve(
    p: jt.Float[jt.Array, "np"],
    t0: jnp.float_,
    ts: jt.Float[jt.Array, "nt_dyn"],
    tcl: jt.Float[jt.Array, "ncl"],
    h: jt.Float[jt.Array, "ne"],
    x0: jt.Float[jt.Array, "nxs"],
    h_mask: jt.Bool[jt.Array, "ne"],
    solver: diffrax.AbstractSolver,
    controller: diffrax.AbstractStepSizeController,
    root_finder: AbstractRootFinder,
    max_steps: jnp.int_,
    adjoint: diffrax.AbstractAdjoint,
    term: diffrax.ODETerm,
    root_cond_fns: list[Callable],
    root_cond_fn: Callable,
    delta_x: Callable,
    known_discs: jt.Float[jt.Array, "*nediscs"],
    observable_ids: list[str],
) -> tuple[jt.Float[jt.Array, "nt nxs"], jt.Float[jt.Array, "nt ne"], dict]:
    """
    Simulate the ODE system for the specified timepoints.

    :param p:
        parameters
    :param t0:
        initial time point
    :param ts:
        time points at which solutions are evaluated
    :param tcl:
        total values for conservation laws
    :param x0:
        initial state vector
    :param solver:
        ODE solver
    :param controller:
        step size controller
    :param max_steps:
        maximum number of steps
    :param adjoint:
        adjoint method
        :param term:
        ODE term
    :param root_cond_fns:
        list of individual root condition functions for discontinuities
    :param root_cond_fn:
        root condition function for all discontinuities
    :param delta_x:
        function to compute state changes at events
    :param known_discs:
        known discontinuities, used to clip the step size controller
    :param observable_ids:
        list of observable IDs
    :return:
        solution+heaviside variables at time points ts and statistics
    """
    # if there are no events, we can avoid expensive looping and just run a single segment
    if not root_cond_fns:
        # no events, we can just run a single segment
        sol, _, stats = _run_segment(
            t0,
            ts[-1],
            x0,
            p,
            tcl,
            h,
            solver,
            controller,
            root_finder,
            max_steps,
            adjoint,
            [],
            [],
            diffrax.SaveAt(ts=ts),
            term,
            dict(**STARTING_STATS),
        )
        return sol.ys, jnp.repeat(h[None, :], sol.ys.shape[0]), stats

    def cond_fn(carry):
        _, t_start, y0, _, _, stats = carry
        # check if we have reached the end of the time points
        return jnp.logical_and(
            t_start < ts[-1],  # final time point not reached
            jnp.isfinite(
                y0
            ).all(),  # y0 is finite, also used to track integration failure
        )

    def body_fn(carry):
        ys, t_start, y0, hs, h, stats = carry
        sol, _, stats = _run_segment(
            t_start,
            ts[-1],
            y0,
            p,
            tcl,
            h,
            solver,
            controller,
            root_finder,
            max_steps,  # TODO: figure out how to pass `max_steps - stats['num_steps']` here
            adjoint,
            root_cond_fns,
            [True] * len(root_cond_fns),
            diffrax.SaveAt(
                subs=[
                    diffrax.SubSaveAt(
                        ts=jnp.where(ts >= t_start, ts, t_start)
                    ),  # datapoints
                    diffrax.SubSaveAt(t1=True),  # events
                ]
            ),
            term,
            stats,
        )
        # update the solution for all timepoints in the simulated segment
        was_simulated = jnp.isin(ts, sol.ts[0])
        ys = jnp.where(was_simulated[:, None], sol.ys[0], ys)
        hs = jnp.where(was_simulated[:, None], h[None, :], hs)

        t0_next = sol.ts[1][
            -1
        ]  # next start time is the end of the current segment
        y0_next = sol.ys[1][
            -1
        ]  # next initial state is the last state of the current segment

        y0_next, h_next, stats = _handle_event(
            t0_next,
            y0_next,
            p,
            tcl,
            h,
            root_finder,
            term,
            root_cond_fn,
            delta_x,
            h_mask,
            stats,
        )

        after_event = sol.ts[1] < ts
        hs = jnp.where(after_event[:, None], h_next[None, :], hs)

        return ys, t0_next, y0_next, hs, h_next, stats

    # run the loop until we have reached the end of the time points
    ys, t0_next, y0_next, hs, _, stats = eqxi.while_loop(
        cond_fn,
        body_fn,
        (
            jnp.zeros((ts.shape[0], x0.shape[0]), dtype=x0.dtype) + x0,
            t0,
            x0,
            jnp.zeros((ts.shape[0], h.shape[0]), dtype=h.dtype),
            h,
            dict(**STARTING_STATS),
        ),
        kind="bounded",
        max_steps=2**6,
    )

    mask = ts == t0_next
    n_obs = len(observable_ids)
    y0_obs = y0_next[None, -n_obs:]
    updated_last = jnp.where(mask[:, None], y0_obs, ys[:, -n_obs:])
    ys = ys.at[:, -n_obs:].set(updated_last)

    return ys, hs, stats


def _run_segment(
    t_start: float,
    t_end: float,
    y0: jt.Float[jt.Array, "nxs"],
    p: jt.Float[jt.Array, "np"],
    tcl: jt.Float[jt.Array, "ncl"],
    h: jt.Float[jt.Array, "ne"],
    solver: diffrax.AbstractSolver,
    controller: diffrax.AbstractStepSizeController,
    root_finder: AbstractRootFinder,
    max_steps: jnp.int_,
    adjoint: diffrax.AbstractAdjoint,
    cond_fns: list[Callable],
    cond_dirs: list[None | bool],
    saveat: diffrax.SaveAt,
    term: diffrax.ODETerm,
    stats: dict,
) -> tuple[diffrax.Solution, int, dict]:
    """Solve a single integration segment and return triggered event index, start time for the next segment,
    aggregated statistics

    The returned index corresponds to the event in ``cond_fns`` that was
    triggered during the integration. ``None`` indicates that the solver
    reached ``t_end`` without any event firing.
    """
    # combine all discontinuity conditions into a single diffrax.Event
    event = (
        diffrax.Event(
            cond_fn=cond_fns,
            root_finder=root_finder,
            direction=cond_dirs,
        )
        if cond_fns
        else None
    )

    sol = diffrax.diffeqsolve(
        term,
        solver,
        args=(p, tcl, h),
        t0=t_start,
        t1=t_end,
        dt0=None,
        y0=y0,
        stepsize_controller=controller,
        max_steps=max_steps,
        adjoint=adjoint,
        saveat=saveat,
        event=event,
        throw=False,
    )
    if os.getenv("JAX_DEBUG") == "1":
        jax.debug.print(
            "Segment: t0 = {}, t1 = {}, y0 = {}, y1 = {}, numsteps = {}, numrejected = {}, result = {}",
            t_start,
            t_end,
            y0,
            sol.ys[-1],
            sol.stats["num_steps"],
            sol.stats["num_rejected_steps"],
            sol.result,
        )

    # extract the event index
    if sol.event_mask is None:
        # no events were triggered
        event_index = -1
    else:
        event_index = jnp.where(
            diffrax.is_event(sol.result),
            jnp.argmax(jnp.array(sol.event_mask)),
            -1,
        )

    # aggregate statistics
    stats = jtu.tree_map(lambda x, y: x + y, stats, sol.stats)

    return sol, event_index, stats


def _handle_event(
    t0_next: float,
    y0_next: jt.Float[jt.Array, "nxs"],
    p: jt.Float[jt.Array, "np"],
    tcl: jt.Float[jt.Array, "ncl"],
    h: jt.Float[jt.Array, "ne"],
    root_finder: AbstractRootFinder,
    term: diffrax.ODETerm,
    root_cond_fn: Callable,
    delta_x: Callable,
    h_mask: jt.Bool[jt.Array, "ne"],
    stats: dict,
):
    args = (p, tcl, h)
    rootvals = root_cond_fn(t0_next, y0_next, args)
    roots_found = jnp.isclose(
        rootvals, 0.0, atol=root_finder.atol, rtol=root_finder.rtol
    )
    droot_dt = (
        # ∂root_cond_fn/∂t
        jax.jacfwd(root_cond_fn, argnums=0)(t0_next, y0_next, args)
        +
        # ∂root_cond_fn/∂y * ∂y/∂t
        jax.jacfwd(root_cond_fn, argnums=1)(t0_next, y0_next, args)
        @ term.vector_field(t0_next, y0_next, args)
    )
    roots_dir = jnp.sign(droot_dt)  # direction of the root condition function

    y0_next, h_next = _apply_event_assignments(
        roots_found,
        roots_dir,
        y0_next,
        p,
        tcl,
        h,
        delta_x,
    )

    h_next = jnp.where(h_mask, h_next, h)

    if os.getenv("JAX_DEBUG") == "1":
        jax.debug.print(
            "rootvals: {}, roots_found: {}, roots_dir: {}, h: {}, h_next: {}",
            rootvals,
            roots_found,
            roots_dir,
            h,
            h_next,
        )

    y0_next = _check_cascading_events(
        t0_next,
        y0_next,
        rootvals,
        p,
        tcl,
        h_next,
        root_finder,
        term,
        root_cond_fn,
        delta_x,
        h_mask,
    )

    return y0_next, h_next, stats

def _check_cascading_events(
    t0_next: float,
    y0_next: jt.Float[jt.Array, "nxs"],
    rootval_prev: jt.Float[jt.Array, "nroots"],
    p: jt.Float[jt.Array, "np"],
    tcl: jt.Float[jt.Array, "ncl"],
    h: jt.Float[jt.Array, "ne"],
    root_finder: AbstractRootFinder,
    term: diffrax.ODETerm,
    root_cond_fn: Callable,
    delta_x: Callable,
    h_mask: jt.Bool[jt.Array, "ne"],
):
    args = (p, tcl, h)
    rootvals = root_cond_fn(t0_next, y0_next, args)
    root_vals_changed_sign = jnp.sign(rootvals) != jnp.sign(rootval_prev)
    roots_dir = jnp.sign(jnp.sign(rootvals) - jnp.sign(rootval_prev))

    y0_next, _ = _apply_event_assignments(
        root_vals_changed_sign,
        roots_dir,
        y0_next,
        p,
        tcl,
        h,
        delta_x,
    )

    return y0_next

def _apply_event_assignments(
    roots_found,
    roots_dir,
    y0_next,
    p,
    tcl,
    h,
    delta_x,
):
    h_next = jnp.where(
        roots_found,
        jnp.logical_not(h),
        h,
    )  # update heaviside variables based on the root condition function

    mask = jnp.array(
        [
            (roots_found & (roots_dir > 0.0) & (h == 0.0))
            for _ in range(y0_next.shape[0])
        ]
    ).T
    delx = delta_x(y0_next, p, tcl)
    if y0_next.size:
        delx = delx.reshape(delx.size // y0_next.shape[0], y0_next.shape[0],)
    y0_up = jnp.where(mask, delx, 0.0)
    y0_next = y0_next + jnp.sum(y0_up, axis=0)

    return y0_next, h_next
