"""Model simulation using JAX."""

# ruff: noqa: F821 F722

from abc import abstractmethod
from pathlib import Path
import enum
from dataclasses import field

import diffrax
import equinox as eqx
import equinox.internal as eqxi
import jax.numpy as jnp
import jax
import jax.tree_util as jtu
from optimistix import AbstractRootFinder
import jaxtyping as jt

from collections.abc import Callable


class ReturnValue(enum.Enum):
    llh = "log-likelihood"
    nllhs = "pointwise negative log-likelihood"
    x0 = "full initial state vector"
    x0_solver = "reduced initial state vector"
    x = "full state vector"
    x_solver = "reduced state vector"
    y = "observables"
    sigmay = "standard deviations of the observables"
    tcl = "total values for conservation laws"
    res = "residuals"
    chi2 = "sum(((observed - simulated) / sigma ) ** 2)"


STARTING_STATS = {
    "max_steps": 0,
    "num_accepted_steps": 0,
    "num_rejected_steps": 0,
    "num_steps": 0,
}


class JAXModel(eqx.Module):
    """
    JAXModel provides an abstract base class for a JAX-based implementation of an AMICI model. The class implements
    routines for simulation and evaluation of derived quantities, model specific implementations need to be provided by
    classes inheriting from JAXModel.

    :ivar api_version:
        API version of the derived class. Needs to match the API version of the base class (MODEL_API_VERSION).
    :ivar MODEL_API_VERSION:
        API version of the base class.
    :ivar jax_py_file:
        Path to the JAX model file.
    """

    MODEL_API_VERSION = "0.0.4"
    api_version: str
    jax_py_file: Path
    parameters: jnp.ndarray = field(default_factory=lambda: jnp.array([]))

    def __init__(self):
        if self.api_version != self.MODEL_API_VERSION:
            raise ValueError(
                "JAXModel API version mismatch, please regenerate the model class."
            )
        super().__init__()

    @abstractmethod
    def _xdot(
        self,
        t: jnp.float_,
        x: jt.Float[jt.Array, "nxs"],
        args: tuple[
            jt.Float[jt.Array, "np"],
            jt.Float[jt.Array, "ncl"],
            jt.Float[jt.Array, "ne"],
        ],
    ) -> jt.Float[jt.Array, "nxs"]:
        """
        Right-hand side of the ODE system.

        :param t: time point
        :param x: state vector
        :param args: tuple of parameters, total values for conservation laws and heaviside variables
        :return:
            Temporal derivative of the state vector x at time point t.
        """
        ...

    @abstractmethod
    def _w(
        self,
        t: jt.Float[jt.Array, ""],
        x: jt.Float[jt.Array, "nxs"],
        p: jt.Float[jt.Array, "np"],
        tcl: jt.Float[jt.Array, "ncl"],
        h: jt.Float[jt.Array, "ne"],
    ) -> jt.Float[jt.Array, "nw"]:
        """
        Compute the expressions, i.e. derived quantities that are used in other parts of the model.

        :param t: time point
        :param x: state vector
        :param p: parameters
        :param tcl: total values for conservation laws
        :param h: heaviside variables
        :return:
            Expression values.
        """
        ...

    @abstractmethod
    def _x0(
        self, t: jnp.float_, p: jt.Float[jt.Array, "np"]
    ) -> jt.Float[jt.Array, "nx"]:
        """
        Compute the initial state vector.

        :param t: initial time point
        :param p: parameters
        :return:
            Initial state vector.
        """
        ...

    @abstractmethod
    def _x_solver(
        self, x: jt.Float[jt.Array, "nx"]
    ) -> jt.Float[jt.Array, "nxs"]:
        """
        Transform the full state vector to the reduced state vector for ODE solving.

        :param x:
            full state vector
        :return:
            reduced state vector
        """
        ...

    @abstractmethod
    def _x_rdata(
        self, x: jt.Float[jt.Array, "nxs"], tcl: jt.Float[jt.Array, "ncl"]
    ) -> jt.Float[jt.Array, "nx"]:
        """
        Compute the full state vector from the reduced state vector and conservation laws.

        :param x:
            reduced state vector
        :param tcl:
            total values for conservation laws
        :return:
            full state vector
        """
        ...

    @abstractmethod
    def _tcl(
        self, x: jt.Float[jt.Array, "nx"], p: jt.Float[jt.Array, "np"]
    ) -> jt.Float[jt.Array, "ncl"]:
        """
        Compute the total values for conservation laws.

        :param x:
            state vector
        :param p:
            parameters
        :return:
            total values for conservation laws
        """
        ...

    @abstractmethod
    def _y(
        self,
        t: jt.Float[jt.Scalar, ""],
        x: jt.Float[jt.Array, "nxs"],
        p: jt.Float[jt.Array, "np"],
        tcl: jt.Float[jt.Array, "ncl"],
        h: jt.Float[jt.Array, "ne"],
        op: jt.Float[jt.Array, "ny"],
    ) -> jt.Float[jt.Array, "ny"]:
        """
        Compute the observables.

        :param t:
            time point
        :param x:
            state vector
        :param p:
            parameters
        :param tcl:
            total values for conservation laws
        :param h:
            heaviside variables
        :param op:
            observable parameters
        :return:
            observables
        """
        ...

    @abstractmethod
    def _sigmay(
        self,
        y: jt.Float[jt.Array, "ny"],
        p: jt.Float[jt.Array, "np"],
        np: jt.Float[jt.Array, "ny"],
    ) -> jt.Float[jt.Array, "ny"]:
        """
        Compute the standard deviations of the observables.

        :param y:
            observables
        :param p:
            parameters
        :param np:
            noise parameters
        :return:
            standard deviations of the observables
        """
        ...

    @abstractmethod
    def _nllh(
        self,
        t: jt.Float[jt.Scalar, ""],
        x: jt.Float[jt.Array, "nxs"],
        p: jt.Float[jt.Array, "np"],
        tcl: jt.Float[jt.Array, "ncl"],
        h: jt.Float[jt.Array, "ne"],
        my: jt.Float[jt.Array, ""],
        iy: jt.Int[jt.Array, ""],
        op: jt.Float[jt.Array, "ny"],
        np: jt.Float[jt.Array, "ny"],
    ) -> jt.Float[jt.Scalar, ""]:
        """
        Compute the negative log-likelihood of the observable for the specified observable index.

        :param t:
            time point
        :param x:
            state vector
        :param p:
            parameters
        :param tcl:
            total values for conservation laws
        :param h:
            heaviside variables
        :param my:
            observed data
        :param iy:
            observable index
        :param op:
            observable parameters
        :param np:
            noise parameters
        :return:
            log-likelihood of the observable
        """
        ...

    @abstractmethod
    def _known_discs(
        self, p: jt.Float[jt.Array, "np"]
    ) -> jt.Float[jt.Array, "ndiscs"]:
        """
        Compute the known discontinuity points of the ODE system.

        :param p:
            parameters
        :return:
            known discontinuity points in the ODE system
        """
        ...

    @abstractmethod
    def _root_cond_fns(
        self,
    ) -> tuple[
        Callable[[float, jt.Float[jt.Array, "nxs"], tuple], jt.Float], ...
    ]:
        """Return condition functions for implicit discontinuities.

        These functions are passed to :class:`diffrax.Event` and must evaluate
        to zero when a discontinuity is triggered.

        :param p:
            model parameters
        :return:
            tuple of callable root functions
        """
        ...

    @abstractmethod
    def _root_cond_fn(
        self,
        t: jt.Float[jt.Scalar, ""],
        y: jt.Float[jt.Array, "nxs"],
        args: tuple[
            jt.Float[jt.Array, "np"],
            jt.Float[jt.Array, "ncl"],
            jt.Float[jt.Array, "ne"],
        ],
    ) -> jt.Float[jt.Array, "ne"]:
        """
        Root condition function for implicit discontinuities.

        :param t:
            time point
        :param y:
            state vector
        :param args:
            tuple of parameters, total values for conservation laws and heaviside variables
        :return:
            root condition values
        """
        ...

    @property
    @abstractmethod
    def state_ids(self) -> list[str]:
        """
        Get the state ids of the model.

        :return:
            State ids
        """
        ...

    @property
    @abstractmethod
    def observable_ids(self) -> list[str]:
        """
        Get the observable ids of the model.

        :return:
            Observable ids
        """
        ...

    @property
    @abstractmethod
    def parameter_ids(self) -> list[str]:
        """
        Get the parameter ids of the model.

        :return:
            Parameter ids
        """
        ...

    @property
    @abstractmethod
    def expression_ids(self) -> list[str]:
        """
        Get the expression ids of the model.

        :return:
            Expression ids
        """
        ...

    def _eq(
        self,
        p: jt.Float[jt.Array, "np"],
        tcl: jt.Float[jt.Array, "ncl"],
        h: jt.Float[jt.Array, "ne"],
        x0: jt.Float[jt.Array, "nxs"],
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        root_finder: AbstractRootFinder,
        steady_state_event: Callable[
            ..., diffrax._custom_types.BoolScalarLike
        ],
        max_steps: jnp.int_,
    ) -> tuple[jt.Float[jt.Array, "1 nxs"], dict]:
        """
        Solve the steady state equation.

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
        :param max_steps:
            maximum number of steps
        :return:
        """
        # if there are no events, we can avoid expensive looping and just run a single segment
        if not self._root_cond_fns():
            sol, _, stats = self._run_segment(
                0.0,
                jnp.inf,
                x0,
                p,
                tcl,
                h,
                solver,
                controller,
                root_finder,
                max_steps,
                diffrax.DirectAdjoint(),
                [steady_state_event],
                diffrax.SaveAt(t1=True),
                dict(**STARTING_STATS),
            )
            return sol.ys[-1][None, :], stats

        def body_fn(carry):
            t_start, y0, event_index, stats = carry
            sol, event_index, stats = self._run_segment(
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
                [steady_state_event] + list(self._root_cond_fns()),
                diffrax.SaveAt(t1=True),
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
            t_start_next = jnp.where(
                jnp.isfinite(sol.ts), sol.ts, -jnp.inf
            ).max()
            return t_start_next, y0_next, event_index, dict(**STARTING_STATS)

        def cond_fn(carry):
            _, y0, event_index, _ = carry
            return jnp.logical_and(
                event_index != 0,  # has not reached steady state yet
                jnp.isfinite(
                    y0
                ).all(),  # y0 is finite, used to indicate integration failure
            )

        # run the loop until no event is triggered (which will also be the case if we run out of steps)
        _, ys, _, stats = eqxi.while_loop(
            cond_fn,
            body_fn,
            (0.0, x0, -1, dict(**STARTING_STATS)),
            kind="bounded",
            max_steps=2**6,
        )

        return ys, stats

    def _solve(
        self,
        p: jt.Float[jt.Array, "np"],
        ts: jt.Float[jt.Array, "nt_dyn"],
        tcl: jt.Float[jt.Array, "ncl"],
        h: jt.Float[jt.Array, "ne"],
        x0: jt.Float[jt.Array, "nxs"],
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        root_finder: AbstractRootFinder,
        max_steps: jnp.int_,
        adjoint: diffrax.AbstractAdjoint,
    ) -> tuple[
        jt.Float[jt.Array, "nt nxs"], jt.Float[jt.Array, "nt ne"], dict
    ]:
        """
        Solve the ODE system.

        :param p:
            parameters
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
        :return:
            solution+heaviside variables at time points ts and statistics
        """
        # if there are no events, we can avoid expensive looping and just run a single segment
        if not self._root_cond_fns():
            # no events, we can just run a single segment
            sol, _, stats = self._run_segment(
                0.0,
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
                diffrax.SaveAt(ts=ts),
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
                ).all(),  # y0 is finite, used to indicate integration failure
            )

        def body_fn(carry):
            ys, t_start, y0, hs, h, stats = carry
            sol, idx, stats = self._run_segment(
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
                list(self._root_cond_fns()),
                diffrax.SaveAt(
                    subs=[
                        diffrax.SubSaveAt(
                            ts=jnp.where(ts >= t_start, ts, t_start)
                        ),  # datapoints
                        diffrax.SubSaveAt(t1=True),  # events
                    ]
                ),
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

            args = (p, tcl, h)
            roots_found = jnp.isclose(
                self._root_cond_fn(t0_next, y0_next, args), 0.0
            )
            droot_dt = (
                # ∂root_cond_fn/∂t
                jax.jacfwd(self._root_cond_fn, argnums=0)(
                    t0_next, y0_next, args
                )
                +
                # ∂root_cond_fn/∂y * ∂y/∂t
                jax.jacfwd(self._root_cond_fn, argnums=1)(
                    t0_next, y0_next, args
                )
                @ self._xdot(t0_next, y0_next, args)
            )
            roots_dir = jnp.sign(
                droot_dt
            )  # direction of the root condition function

            h_next = (
                h
                + jnp.where(
                    roots_found,
                    roots_dir,
                    jnp.zeros_like(h),
                )
            )  # update heaviside variables based on the root condition function
            was_event = jnp.isin(ts, sol.ts[1])
            hs = jnp.where(was_event[:, None], h_next[None, :], hs)
            jax.debug.print(
                "roots_found: {}, roots_dir: {}, h: {}, h_next: {}",
                roots_found,
                roots_dir,
                h,
                h_next,
            )

            ts_next = jnp.where(
                ts > t0_next, ts, ts[-1]
            ).min()  # timepoint of next datapoint, don't step over that
            # while we are trying to step out of the event
            y0_next, t0_next, stats = self._step_out_of_event(
                t0_next,
                ts_next,
                y0_next,
                p,
                tcl,
                h_next,
                solver,
                controller,
                adjoint,
                stats,
            )
            return ys, t0_next, y0_next, hs, h_next, stats

        # run the loop until we have reached the end of the time points
        ys, _, _, hs, _, stats = eqxi.while_loop(
            cond_fn,
            body_fn,
            (
                jnp.zeros((ts.shape[0], x0.shape[0]), dtype=x0.dtype) + x0,
                0.0,
                x0,
                jnp.zeros((ts.shape[0], h.shape[0]), dtype=h.dtype),
                h,
                dict(**STARTING_STATS),
            ),
            kind="bounded",
            max_steps=2**6,
        )

        return ys, hs, stats

    def _run_segment(
        self,
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
        saveat: diffrax.SaveAt,
        stats: dict,
    ) -> tuple[diffrax.Solution, int, dict]:
        """Solve a single integration segment and return triggered event index, start time for the next segment,
        aggregated statistics

        The returned index corresponds to the event in ``cond_fns`` that was
        triggered during the integration. ``None`` indicates that the solver
        reached ``t_end`` without any event firing.
        """

        # combine all discontinuity conditions into a single diffrax.Event
        event = diffrax.Event(cond_fns, root_finder) if cond_fns else None

        # manage events with explicit discontinuities
        controller = (
            diffrax.ClipStepSizeController(
                controller,
                jump_ts=self._known_discs(p),
            )
            if self._known_discs(p).size
            else controller
        )

        sol = diffrax.diffeqsolve(
            diffrax.ODETerm(self._xdot),
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

    def _step_out_of_event(
        self,
        t_start: float,
        t_end: float,
        y0: jt.Float[jt.Array, "nxs"],
        p: jt.Float[jt.Array, "np"],
        tcl: jt.Float[jt.Array, "ncl"],
        h: jt.Float[jt.Array, "ne"],
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        adjoint: diffrax.AbstractAdjoint,
        stats: dict,
    ):
        # take a single step, without event, to step out of the event
        # TODO: loop until we are out of the event
        sol = diffrax.diffeqsolve(
            diffrax.ODETerm(self._xdot),
            solver,
            args=(p, tcl, h),
            t0=t_start,
            t1=t_end,
            dt0=None,
            y0=y0,
            stepsize_controller=controller,
            max_steps=10,
            adjoint=adjoint,
            saveat=diffrax.SaveAt(
                subs=[
                    diffrax.SubSaveAt(t1=True),  # we manage to reach t_end
                    diffrax.SubSaveAt(
                        steps=True
                    ),  # we want to save the state at t_end
                ]
            ),
            throw=False,
        )
        jax.debug.print(
            "StepOutOfEvent: t0 = {}, t1 = {}, y0 = {}, y1 = {}, numsteps = {}, num_rejected = {}, result = {}",
            t_start,
            t_end,
            y0,
            sol.ys[0][-1],
            sol.stats["num_steps"],
            sol.stats["num_rejected_steps"],
            sol.result,
        )
        jax.debug.print("ys: {}", sol.ys)
        # aggregate stats
        stats = jtu.tree_map(lambda x, y: x + y, stats, sol.stats)
        return sol.ys[0][-1], sol.ts[0][-1], stats

    def _initialise_heaviside_variables(
        self,
        t0: jt.Float[jt.Scalar, ""],
        x_solver: jt.Float[jt.Array, "nxs"],
        p: jt.Float[jt.Array, "np"],
        tcl: jt.Float[jt.Array, "ncl"],
    ) -> jt.Float[jt.Array, "ne"]:
        """
        Initialise the heaviside variables.

        :param t0:
            initial time point
        :param x_solver:
            reduced initial state vector
        :param p:
            parameters
        :param tcl:
            total values for conservation laws
        :return:
            heaviside variables
        """
        h0 = jnp.zeros((len(self._root_cond_fns()),))  # dummy values
        roots_found = self._root_cond_fn(t0, x_solver, (p, tcl, h0))
        return jnp.where(
            roots_found >= 0.0, jnp.ones_like(h0), jnp.zeros_like(h0)
        )

    def _x_rdatas(
        self, x: jt.Float[jt.Array, "nt nxs"], tcl: jt.Float[jt.Array, "ncl"]
    ) -> jt.Float[jt.Array, "nt nx"]:
        """
        Compute the full state vector from the reduced state vector and conservation laws.

        :param x:
            reduced state vector
        :param tcl:
            total values for conservation laws
        :return:
            full state vector
        """
        return jax.vmap(self._x_rdata, in_axes=(0, None))(x, tcl)

    def _nllhs(
        self,
        ts: jt.Float[jt.Array, "nt nx"],
        xs: jt.Float[jt.Array, "nt nxs"],
        p: jt.Float[jt.Array, "np"],
        tcl: jt.Float[jt.Array, "ncl"],
        h: jt.Float[jt.Array, "nt ne"],
        mys: jt.Float[jt.Array, "nt"],
        iys: jt.Int[jt.Array, "nt"],
        ops: jt.Float[jt.Array, "nt *nop"],
        nps: jt.Float[jt.Array, "nt *nnp"],
    ) -> jt.Float[jt.Array, "nt"]:
        """
        Compute the negative log-likelihood for each observable.

        :param ts:
            time points
        :param xs:
            state vectors
        :param p:
            parameters
        :param tcl:
            total values for conservation laws
        :param h:
            heaviside variables
        :param mys:
            observed data
        :param iys:
            observable indices
        :param ops:
            observable parameters
        :param nps:
            noise parameters
        :return:
            negative log-likelihoods of the observables
        """
        return jax.vmap(self._nllh, in_axes=(0, 0, None, None, 0, 0, 0, 0, 0))(
            ts, xs, p, tcl, h, mys, iys, ops, nps
        )

    def _ys(
        self,
        ts: jt.Float[jt.Array, "nt"],
        xs: jt.Float[jt.Array, "nt nxs"],
        p: jt.Float[jt.Array, "np"],
        tcl: jt.Float[jt.Array, "ncl"],
        h: jt.Float[jt.Array, "nt ne"],
        iys: jt.Float[jt.Array, "nt"],
        ops: jt.Float[jt.Array, "nt *nop"],
    ) -> jt.Int[jt.Array, "nt"]:
        """
        Compute the observables.

        :param ts:
            time points
        :param xs:
            state vectors
        :param p:
            parameters
        :param tcl:
            total values for conservation laws
        :param h:
            heaviside variables
        :param iys:
            observable indices
        :param ops:
            observables parameters
        :return:
            observables
        """
        return jax.vmap(
            lambda t, x, p, tcl, iy, op: self._y(t, x, p, tcl, h, op)
            .at[iy]
            .get(),
            in_axes=(0, 0, None, None, 0, 0, 0),
        )(ts, xs, p, tcl, h, iys, ops)

    def _sigmays(
        self,
        ts: jt.Float[jt.Array, "nt"],
        xs: jt.Float[jt.Array, "nt nxs"],
        p: jt.Float[jt.Array, "np"],
        tcl: jt.Float[jt.Array, "ncl"],
        h: jt.Float[jt.Array, "nt ne"],
        iys: jt.Int[jt.Array, "nt"],
        ops: jt.Float[jt.Array, "nt *nop"],
        nps: jt.Float[jt.Array, "nt *nnp"],
    ):
        """
        Compute the standard deviations of the observables.

        :param ts:
            time points
        :param xs:
            state vectors
        :param p:
            parameters
        :param tcl:
            total values for conservation laws
        :param h:
            heaviside variables
        :param iys:
            observable indices
        :param ops:
            observable parameters
        :param nps:
            noise parameters
        :return:
            standard deviations of the observables
        """
        return jax.vmap(
            lambda t, x, p, tcl, iy, op, np: self._sigmay(
                self._y(t, x, p, tcl, h, op), p, np
            )
            .at[iy]
            .get(),
            in_axes=(0, 0, None, None, 0, 0, 0, 0),
        )(ts, xs, p, tcl, h, iys, ops, nps)

    @eqx.filter_jit
    def simulate_condition(
        self,
        p: jt.Float[jt.Array, "np"] | None,
        ts_dyn: jt.Float[jt.Array, "nt_dyn"],
        ts_posteq: jt.Float[jt.Array, "nt_posteq"],
        my: jt.Float[jt.Array, "nt"],
        iys: jt.Int[jt.Array, "nt"],
        iy_trafos: jt.Int[jt.Array, "nt"],
        ops: jt.Float[jt.Array, "nt *nop"],
        nps: jt.Float[jt.Array, "nt *nnp"],
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        root_finder: AbstractRootFinder,
        adjoint: diffrax.AbstractAdjoint,
        steady_state_event: Callable[
            ..., diffrax._custom_types.BoolScalarLike
        ],
        max_steps: int | jnp.int_,
        x_preeq: jt.Float[jt.Array, "*nx"] = jnp.array([]),
        mask_reinit: jt.Bool[jt.Array, "*nx"] = jnp.array([]),
        x_reinit: jt.Float[jt.Array, "*nx"] = jnp.array([]),
        ts_mask: jt.Bool[jt.Array, "nt"] = jnp.array([]),
        ret: ReturnValue = ReturnValue.llh,
    ) -> tuple[jt.Float[jt.Array, "nt *nx"] | jnp.float_, dict]:
        r"""
        Simulate a condition.

        :param p:
            parameters for simulation ordered according to ids in :ivar parameter_ids:. If ``None``,
            the values stored in :attr:`parameters` are used.
        :param ts_dyn:
            time points for dynamic simulation. Sorted in monotonically increasing order but duplicate time points are
            allowed to facilitate the evaluation of multiple observables at specific time points.
        :param ts_posteq:
            time points for post-equilibration. Usually valued \Infty, but needs to be shaped according to
            the number of observables that are evaluated after post-equilibration.
        :param my:
            observed data
        :param iys:
            indices of the observables according to ordering in :ivar observable_ids:
        :param iy_trafos:
            indices of transformations for observables
        :param ops:
            observable parameters
        :param nps:
            noise parameters
        :param solver:
            ODE solver
        :param controller:
            step size controller
        :param adjoint:
            adjoint method. Recommended values are `diffrax.DirectAdjoint()` for jax.jacfwd (with vector-valued
            outputs) and  `diffrax.RecursiveCheckpointAdjoint()` for jax.grad (for scalar-valued outputs).
        :param steady_state_event:
            event function for steady state. See :func:`diffrax.steady_state_event` for details.
        :param max_steps:
            maximum number of solver steps
        :param x_preeq:
            initial state vector for pre-equilibration. If not provided, the initial state vector is computed using
            :meth:`_x0`.
        :param mask_reinit:
            mask for re-initialization. If `True`, the corresponding state variable is re-initialized.
        :param x_reinit:
            re-initialized state vector. If not provided, the state vector is not re-initialized.
        :param ts_mask:
            mask to remove (padded) time points. If `True`, the corresponding time point is used for the evaluation of
            the output. Only applied if ret is ReturnValue.llh, ReturnValue.nllhs, ReturnValue.res, or ReturnValue.chi2.
        :param ret:
            which output to return. See :class:`ReturnValue` for available options.
        :return:
            output according to `ret` and general results/statistics
        """
        t0 = 0.0
        if p is None:
            p = self.parameters

        if x_preeq.shape[0]:
            x = x_preeq
        else:
            x = self._x0(t0, p)

        if not ts_mask.shape[0]:
            ts_mask = jnp.ones_like(my, dtype=jnp.bool_)

        # Re-initialization
        if x_reinit.shape[0]:
            x = jnp.where(mask_reinit, x_reinit, x)
        x_solver = self._x_solver(x)
        tcl = self._tcl(x, p)
        h = self._initialise_heaviside_variables(t0, x_solver, p, tcl)

        # Dynamic simulation
        if ts_dyn.shape[0]:
            x_dyn, h_dyn, stats_dyn = self._solve(
                p,
                ts_dyn,
                tcl,
                h,
                x_solver,
                solver,
                controller,
                root_finder,
                max_steps,
                adjoint,
            )
            x_solver = x_dyn[-1, :]
        else:
            x_dyn = jnp.repeat(x_solver[None, :], ts_dyn.shape[0], axis=0)
            h_dyn = jnp.repeat(h[None, :], ts_dyn.shape[0], axis=0)
            stats_dyn = None

        # Post-equilibration
        if ts_posteq.shape[0]:
            x_solver, h, stats_posteq = self._eq(
                p,
                tcl,
                h,
                x_solver,
                solver,
                controller,
                root_finder,
                steady_state_event,
                max_steps,
            )
        else:
            stats_posteq = None

        x_posteq = jnp.repeat(x_solver[None, :], ts_posteq.shape[0], axis=0)
        h_posteq = jnp.repeat(h[None, :], ts_posteq.shape[0], axis=0)

        ts = jnp.concatenate((ts_dyn, ts_posteq), axis=0)
        if len(self._root_cond_fns()):
            hs = jnp.concatenate((h_dyn, h_posteq), axis=0)
        else:
            hs = jnp.zeros((ts.shape[0], h.shape[0]))
        x = jnp.concatenate((x_dyn, x_posteq), axis=0)

        nllhs = self._nllhs(ts, x, p, tcl, hs, my, iys, ops, nps)
        nllhs = jnp.where(ts_mask, nllhs, 0.0)
        llh = -jnp.sum(nllhs)

        stats = dict(
            ts=ts,
            x=x,
            hs=hs,
            llh=llh,
            stats_dyn=stats_dyn,
            stats_posteq=stats_posteq,
        )
        if ret == ReturnValue.llh:
            output = llh
        elif ret == ReturnValue.nllhs:
            output = nllhs
        elif ret == ReturnValue.x:
            output = self._x_rdatas(x, tcl)
        elif ret == ReturnValue.x_solver:
            output = x
        elif ret == ReturnValue.y:
            output = self._ys(ts, x, p, tcl, hs, iys, ops)
        elif ret == ReturnValue.sigmay:
            output = self._sigmays(ts, x, p, tcl, hs, iys, ops, nps)
        elif ret == ReturnValue.x0:
            output = self._x_rdata(x[0, :], tcl)
        elif ret == ReturnValue.x0_solver:
            output = x[0, :]
        elif ret == ReturnValue.tcl:
            output = tcl
        elif ret in (ReturnValue.res, ReturnValue.chi2):
            obs_trafo = jax.vmap(
                lambda y, iy_trafo: jnp.array(
                    # needs to follow order in amici.jax.petab.SCALE_TO_INT
                    [y, safe_log(y), safe_log(y) / jnp.log(10)]
                )
                .at[iy_trafo]
                .get(),
            )
            ys_obj = obs_trafo(
                self._ys(ts, x, p, tcl, hs, iys, ops), iy_trafos
            )
            m_obj = obs_trafo(my, iy_trafos)
            if ret == ReturnValue.chi2:
                sigma_obj = self._sigmays(ts, x, p, tcl, hs, iys, ops, nps)
                chi2 = jnp.square((ys_obj - m_obj) / sigma_obj)
                chi2 = jnp.where(ts_mask, chi2, 0.0)
                output = jnp.sum(chi2)
            else:
                output = ys_obj - m_obj
                output = jnp.where(ts_mask, output, 0.0)
        else:
            raise NotImplementedError(f"Return value {ret} not implemented.")

        return output, stats

    @eqx.filter_jit
    def preequilibrate_condition(
        self,
        p: jt.Float[jt.Array, "np"] | None,
        x_reinit: jt.Float[jt.Array, "*nx"],
        mask_reinit: jt.Bool[jt.Array, "*nx"],
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        root_finder: AbstractRootFinder,
        steady_state_event: Callable[
            ..., diffrax._custom_types.BoolScalarLike
        ],
        max_steps: int | jnp.int_,
    ) -> tuple[jt.Float[jt.Array, "nx"], dict]:
        r"""
        Simulate a condition.

        :param p:
            parameters for simulation ordered according to ids in :ivar parameter_ids:. If ``None``,
            the values stored in :attr:`parameters` are used.
        :param x_reinit:
            re-initialized state vector. If not provided, the state vector is not re-initialized.
        :param mask_reinit:
            mask for re-initialization. If `True`, the corresponding state variable is re-initialized.
        :param solver:
            ODE solver
        :param controller:
            step size controller
        :param max_steps:
            maximum number of solver steps
        :return:
            pre-equilibrated state variables and statistics
        """
        # Pre-equilibration
        t0 = 0.0
        if p is None:
            p = self.parameters

        x0 = self._x0(t0, p)
        if x_reinit.shape[0]:
            x0 = jnp.where(mask_reinit, x_reinit, x0)
        tcl = self._tcl(x0, p)
        h = self.initialise_heaviside_variables(t0, self._x_solver(x0), p, tcl)
        current_x = self._x_solver(x0)
        current_x, _, stats_preeq = self._eq(
            p,
            tcl,
            h,
            current_x,
            solver,
            controller,
            root_finder,
            steady_state_event,
            max_steps,
        )

        return self._x_rdata(current_x, tcl), dict(stats_preeq=stats_preeq)


def safe_log(x: jnp.float_) -> jnp.float_:
    """
    Safe logarithm that returns `jnp.log(jnp.finfo(jnp.float_).eps)` for x <= 0.

    :param x:
        input
    :return:
        logarithm of x
    """
    # see https://docs.kidger.site/equinox/api/debug/, need double jnp.where to guard
    # against nans in forward & backward passes
    safe_x = jnp.where(
        x > jnp.finfo(jnp.float_).eps, x, jnp.finfo(jnp.float_).eps
    )
    return jnp.where(
        x > 0, jnp.log(safe_x), jnp.log(jnp.finfo(jnp.float_).eps)
    )


def safe_div(x: jnp.float_, y: jnp.float_) -> jnp.float_:
    """
    Safe division that returns `x/jnp.finfo(jnp.float_).eps` for `y == 0`.

    :param x:
        numerator
    :param y:
        denominator
    :return:
        x / y
    """
    # see https://docs.kidger.site/equinox/api/debug/, need double jnp.where to guard
    # against nans in forward & backward passes
    safe_y = jnp.where(y != 0, y, jnp.finfo(jnp.float_).eps)
    return jnp.where(y != 0, x / safe_y, x / jnp.finfo(jnp.float_).eps)
