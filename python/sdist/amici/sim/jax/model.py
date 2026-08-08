"""Model simulation using JAX."""

# ruff: noqa: F821 F722

import enum
import os
import warnings
from abc import abstractmethod
from collections.abc import Callable
from dataclasses import field
from pathlib import Path

import diffrax
import equinox as eqx
import jax
import jax.numpy as jnp
import jaxtyping as jt
from optimistix import AbstractRootFinder

from amici.sim.jax._simulation import _apply_event_assignments, eq, solve


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
    nns: dict
    parameters: jnp.ndarray = field(default_factory=lambda: jnp.array([]))
    _array_inputs: dict = field(default_factory=dict)
    _array_input_index: jnp.int32 = field(default_factory=lambda: jnp.int32(0))

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

    @abstractmethod
    def _delta_x(
        self, y: jt.Float[jt.Array, "nxs"]
    ) -> jt.Float[jt.Array, "nxs"]:
        """
        Compute the state vector changes at discontinuities.

        :param y:
            state vector
        :return:
            changes in the state vector at discontinuities
        """
        ...

    @property
    @abstractmethod
    def n_events(self) -> int:
        """
        Get the number of events (that require implicit tracking)

        :return:
            number
        """

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

    def _root_cond_fn_event(
        self, ie: int, t: float, y: jt.Float[jt.Array, "nxs"], args: tuple, **_
    ):
        """
        Root condition function for a specific event index.

        :param ie:
            event index
        :param t:
            time point
        :param y:
            state vector
        :param args:
            tuple of arguments required for _root_cond_fn
        :return:
            mask of root condition value for the specified event index
        """
        __, __, h = args
        rval = self._root_cond_fn(t, y, args, **_)
        # only allow root triggers where trigger function is negative (heaviside == 0)
        masked_rval = jnp.where(h == 0.0, rval, 1.0)
        return masked_rval.at[ie].get()

    def _root_cond_fns(
        self,
    ) -> list[Callable[[float, jt.Float[jt.Array, "nxs"], tuple], jt.Float]]:
        """Return condition functions for implicit discontinuities.

        These functions are passed to :class:`diffrax.Event` and must evaluate
        to zero when a discontinuity is triggered.

        :param p:
            model parameters
        :return:
            iterable of callable root functions
        """
        return [
            eqx.Partial(self._root_cond_fn_event, ie)
            for ie in range(self.n_events)
        ]

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
        h0 = self.event_initial_values.astype(float)
        if os.getenv("JAX_DEBUG") == "1":
            jax.debug.print(
                "h0: {}",
                h0,
            )
        roots_found = self._root_cond_fn(t0, x_solver, (p, tcl, h0))
        return jnp.where(
            jnp.logical_and(roots_found >= 0.0, h0 == 1.0),
            jnp.ones_like(h0),
            jnp.zeros_like(h0),
        )

    def _x_rdatas(
        self,
        x: jt.Float[jt.Array, "nt nxs"],
        tcl: jt.Float[jt.Array, "nt ncl"],
    ) -> jt.Float[jt.Array, "nt nx"]:
        """
        Compute the full state vector from the reduced state vector and conservation laws.

        :param x:
            reduced state vector
        :param tcl:
            total values for conservation laws, per time point (conservation
            laws may change across period boundaries as parameters change)
        :return:
            full state vector
        """
        return jax.vmap(self._x_rdata, in_axes=(0, 0))(x, tcl)

    def _nllhs(
        self,
        ts: jt.Float[jt.Array, "nt nx"],
        xs: jt.Float[jt.Array, "nt nxs"],
        p: jt.Float[jt.Array, "nt np"],
        tcl: jt.Float[jt.Array, "nt ncl"],
        hs: jt.Float[jt.Array, "nt ne"],
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
            parameters, per time point (parameters may differ across period
            boundaries)
        :param tcl:
            total values for conservation laws, per time point
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
        return jax.vmap(self._nllh, in_axes=(0, 0, 0, 0, 0, 0, 0, 0, 0))(
            ts, xs, p, tcl, hs, mys, iys, ops, nps
        )

    def _ys(
        self,
        ts: jt.Float[jt.Array, "nt"],
        xs: jt.Float[jt.Array, "nt nxs"],
        p: jt.Float[jt.Array, "nt np"],
        tcl: jt.Float[jt.Array, "nt ncl"],
        hs: jt.Float[jt.Array, "nt ne"],
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
            parameters, per time point (parameters may differ across period
            boundaries)
        :param tcl:
            total values for conservation laws, per time point
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
            lambda t, x, p, tcl, h, iy, op: (
                self._y(t, x, p, tcl, h, op).at[iy].get()
            ),
            in_axes=(0, 0, 0, 0, 0, 0, 0),
        )(ts, xs, p, tcl, hs, iys, ops)

    def _sigmays(
        self,
        ts: jt.Float[jt.Array, "nt"],
        xs: jt.Float[jt.Array, "nt nxs"],
        p: jt.Float[jt.Array, "nt np"],
        tcl: jt.Float[jt.Array, "nt ncl"],
        hs: jt.Float[jt.Array, "nt ne"],
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
            parameters, per time point (parameters may differ across period
            boundaries)
        :param tcl:
            total values for conservation laws, per time point
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
            lambda t, x, p, tcl, h, iy, op, np: (
                self._sigmay(self._y(t, x, p, tcl, h, op), p, np).at[iy].get()
            ),
            in_axes=(0, 0, 0, 0, 0, 0, 0, 0),
        )(ts, xs, p, tcl, hs, iys, ops, nps)

    def _simulate_period(
        self,
        p: jt.Float[jt.Array, "np"],
        t0: jnp.float_,
        ts_dyn: jt.Float[jt.Array, "nt_dyn"],
        ts_posteq: jt.Float[jt.Array, "nt_posteq"],
        tcl: jt.Float[jt.Array, "ncl"],
        h: jt.Bool[jt.Array, "ne"],
        h_mask: jt.Bool[jt.Array, "ne"],
        x_solver: jt.Float[jt.Array, "nxs"],
        has_posteq_slot: bool,
        do_posteq: jt.Bool[jt.Scalar, ""],  # noqa: F722
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        root_finder: AbstractRootFinder,
        adjoint: diffrax.AbstractAdjoint,
        steady_state_event: Callable[
            ..., diffrax._custom_types.BoolScalarLike
        ],
        max_steps: int | jnp.int_,
    ):
        """
        Simulate a single experiment period, starting from ``x_solver``/``h``
        at time ``t0`` with parameters ``p``/``tcl``.

        Only an experiment's own final period is post-equilibrated (see
        ``has_posteq_slot``/``do_posteq``); earlier periods only integrate
        up to (and including) the synthetic hand-off time point appended to
        ``ts_dyn`` by :meth:`amici.sim.jax.petab.JAXProblem._get_measurements`.

        :return:
            Tuple of (state trajectory, heaviside trajectory, ending reduced
            state, ending heaviside state, dynamic simulation statistics,
            post-equilibration statistics).
        """
        if ts_dyn.shape[0]:
            x_dyn, h_dyn, stats_dyn = solve(
                p,
                t0,
                ts_dyn,
                tcl,
                h,
                x_solver,
                h_mask,
                solver,
                controller,
                root_finder,
                max_steps,
                adjoint,
                diffrax.ODETerm(self._xdot),
                self._root_cond_fns(),
                self._root_cond_fn,
                self._delta_x,
                self._known_discs(p),
                self.observable_ids,
            )
            x_solver = x_dyn[-1, :]
            h = h_dyn[-1, :]
        else:
            x_dyn = jnp.repeat(x_solver[None, :], ts_dyn.shape[0], axis=0)
            h_dyn = jnp.repeat(h[None, :], ts_dyn.shape[0], axis=0)
            stats_dyn = None

        # Post-equilibration. Which period is an experiment's last is a
        # per-experiment property, but this loop is shared across the
        # vmapped experiment axis, so the two halves of the decision are
        # split: ``has_posteq_slot`` is static and decides whether the
        # steady-state solve is traced here at all (no experiment carries
        # post-equilibration rows in this slot -> no solve), while
        # ``do_posteq`` is a per-experiment traced flag deciding whether
        # the result is actually adopted. The blend matters: ``eq``
        # replaces the state that is handed to the next period, so
        # equilibrating unconditionally would clobber the trajectory of
        # any experiment whose chain continues past this slot.
        if has_posteq_slot and ts_posteq.shape[0]:
            x_eq, h_eq, stats_posteq = eq(
                p,
                tcl,
                h,
                x_solver,
                h_mask,
                solver,
                controller,
                root_finder,
                steady_state_event,
                diffrax.ODETerm(self._xdot),
                self._root_cond_fns(),
                self._root_cond_fn,
                self._delta_x,
                self._known_discs(p),
                max_steps,
            )
            x_solver = jnp.where(do_posteq, x_eq, x_solver)
            if h.shape[0]:
                h = jnp.where(do_posteq, h_eq, h)
        else:
            stats_posteq = None

        x_posteq = jnp.repeat(x_solver[None, :], ts_posteq.shape[0], axis=0)
        h_posteq = jnp.repeat(h[None, :], ts_posteq.shape[0], axis=0)

        ts = jnp.concatenate((ts_dyn, ts_posteq), axis=0)
        if h.shape[0]:
            hs = jnp.concatenate((h_dyn, h_posteq), axis=0)
        else:
            hs = jnp.zeros((ts.shape[0], h.shape[0]))
        xs = jnp.concatenate((x_dyn, x_posteq), axis=0)

        return ts, xs, hs, x_solver, h, stats_dyn, stats_posteq

    def simulate_experiment_unjitted(
        self,
        p: jt.Float[jt.Array, "P np"],
        ts_dyn: jt.Float[jt.Array, "P nt_dyn"],
        ts_posteq: jt.Float[jt.Array, "P nt_posteq"],
        my: jt.Float[jt.Array, "P nt"],
        iys: jt.Int[jt.Array, "P nt"],
        iy_trafos: jt.Int[jt.Array, "P nt"],
        ops: jt.Float[jt.Array, "P nt *nop"],
        nps: jt.Float[jt.Array, "P nt *nnp"],
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        root_finder: AbstractRootFinder,
        adjoint: diffrax.AbstractAdjoint,
        steady_state_event: Callable[
            ..., diffrax._custom_types.BoolScalarLike
        ],
        max_steps: int | jnp.int_,
        x_preeq: jt.Float[jt.Array, "*nx"] | None = None,
        h_preeq: jt.Float[jt.Array, "*ne"] | None = None,
        mask_reinit: jt.Bool[jt.Array, "P *nx"] | None = None,
        x_reinit: jt.Float[jt.Array, "P *nx"] | None = None,
        init_override: jt.Float[jt.Array, "*nx"] | None = None,
        init_override_mask: jt.Bool[jt.Array, "*nx"] | None = None,
        ts_mask: jt.Bool[jt.Array, "P nt"] | None = None,
        h_mask: jt.Bool[jt.Array, "ne"] | None = None,
        t_zero: jt.Float[jt.Array, "P"] | None = None,
        posteq_mask: jt.Bool[jt.Array, "P"] | None = None,
        posteq_slots: tuple[bool, ...] | None = None,
        ret: ReturnValue = ReturnValue.llh,
    ) -> tuple[jt.Float[jt.Array, "*nt"], dict]:
        """
        Unjitted version of simulate_experiment.

        Chains one ODE integration per experiment period (the leading axis,
        of static size ``P``, of ``p``/``ts_dyn``/``ts_posteq``/``my``/
        ``iys``/``iy_trafos``/``ops``/``nps``/``mask_reinit``/``x_reinit``/
        ``ts_mask``/``t_zero``), carrying the ODE state and heaviside/event
        state from the end of one period into the start of the next, in
        lieu of encoding period transitions as model events. ``P == 1``
        reduces to a single, non-chained simulation.

        See :meth:`simulate_experiment` for full documentation.
        """
        n_periods = p.shape[0]

        # Normalize omitted optional arrays here, at call time, rather than
        # via eager `jnp.array(...)`-valued default arguments: a default
        # constructed once at function-definition time freezes to
        # float32 if `jax_enable_x64` is enabled only after this module is
        # first imported, silently diverging in dtype from every other
        # (call-time-constructed) array flowing through the same call.
        if x_preeq is None:
            x_preeq = jnp.array([])
        if h_preeq is None:
            h_preeq = jnp.array([])
        if mask_reinit is None:
            mask_reinit = jnp.array([])
        if x_reinit is None:
            x_reinit = jnp.array([])
        if init_override is None:
            init_override = jnp.array([])
        if init_override_mask is None:
            init_override_mask = jnp.array([])
        if ts_mask is None:
            ts_mask = jnp.array([])
        if h_mask is None:
            h_mask = jnp.array([])
        if t_zero is None:
            t_zero = jnp.zeros(n_periods)
        if posteq_mask is None:
            # default: only the last period post-equilibrates, which is
            # what a single-experiment (or uniform-period-count) caller
            # wants and reproduces the previous `is_final` behaviour
            posteq_mask = jnp.arange(n_periods) == n_periods - 1
        if posteq_slots is None:
            posteq_slots = tuple(
                i == n_periods - 1 for i in range(n_periods)
            )

        if not h_mask.shape[0]:
            h_mask = jnp.ones(self.n_events, dtype=jnp.bool_)

        if not ts_mask.shape[0]:
            ts_mask = jnp.ones_like(my, dtype=jnp.bool_)

        has_reinit = x_reinit.shape[-1] > 0

        t0_0 = t_zero[0]
        if x_preeq.shape[0]:
            x = x_preeq
        elif init_override.shape[0]:
            x_def = self._x0(t0_0, p[0])
            x = jnp.where(init_override_mask, init_override, x_def)
        else:
            x = self._x0(t0_0, p[0])

        h = h_preeq
        x_solver = None
        tcl_prev = None

        ts_list = []
        x_list = []
        h_list = []
        tcl_list = []
        p_list = []
        stats_dyn_list = []
        stats_posteq_final = None

        for i in range(n_periods):
            p_i = p[i]
            t0_i = t_zero[i]

            if i == 0:
                x_i_full = x
            else:
                # carry reduced state from the end of the previous period
                # back to full state space
                x_i_full = self._x_rdata(x_solver, tcl_prev)

            if has_reinit:
                x_i_full = jnp.where(mask_reinit[i], x_reinit[i], x_i_full)

            x_solver = self._x_solver(x_i_full)
            tcl_i = self._tcl(x_i_full, p_i)

            x_solver, _, h, _ = self._handle_t0_event(
                t0_i,
                x_solver,
                p_i,
                tcl_i,
                root_finder,
                self._root_cond_fn,
                self._delta_x,
                h_mask,
                h,
                {},
            )

            is_final = i == n_periods - 1
            ts_i, xs_i, hs_i, x_solver, h, stats_dyn_i, stats_posteq_i = (
                self._simulate_period(
                    p_i,
                    t0_i,
                    ts_dyn[i],
                    ts_posteq[i],
                    tcl_i,
                    h,
                    h_mask,
                    x_solver,
                    posteq_slots[i],
                    posteq_mask[i],
                    solver,
                    controller,
                    root_finder,
                    adjoint,
                    steady_state_event,
                    max_steps,
                )
            )
            if is_final:
                stats_posteq_final = stats_posteq_i

            ts_list.append(ts_i)
            x_list.append(xs_i)
            h_list.append(hs_i)
            tcl_list.append(jnp.repeat(tcl_i[None, :], ts_i.shape[0], axis=0))
            p_list.append(jnp.repeat(p_i[None, :], ts_i.shape[0], axis=0))
            stats_dyn_list.append(stats_dyn_i)

            tcl_prev = tcl_i

        ts = jnp.concatenate(ts_list, axis=0)
        x = jnp.concatenate(x_list, axis=0)
        hs = jnp.concatenate(h_list, axis=0)
        tcls = jnp.concatenate(tcl_list, axis=0)
        ps = jnp.concatenate(p_list, axis=0)

        my = my.reshape(-1)
        iys = iys.reshape(-1)
        iy_trafos = iy_trafos.reshape(-1)
        # avoid `-1` in reshape: it errors on a `math.prod(...) == 0`
        # trailing shape (e.g. no observable/noise parameter overrides)
        ops = ops.reshape(ops.shape[0] * ops.shape[1], *ops.shape[2:])
        nps = nps.reshape(nps.shape[0] * nps.shape[1], *nps.shape[2:])
        ts_mask = ts_mask.reshape(-1)

        nllhs = self._nllhs(ts, x, ps, tcls, hs, my, iys, ops, nps)
        nllhs = jnp.where(ts_mask, nllhs, 0.0)
        llh = -jnp.sum(nllhs)

        stats = dict(
            ts=ts,
            x=x,
            hs=hs,
            llh=llh,
            stats_dyn=stats_dyn_list,
            stats_posteq=stats_posteq_final,
        )
        if ret == ReturnValue.llh:
            output = llh
        elif ret == ReturnValue.nllhs:
            output = nllhs
        elif ret == ReturnValue.x:
            output = self._x_rdatas(x, tcls)
        elif ret == ReturnValue.x_solver:
            output = x
        elif ret == ReturnValue.y:
            output = self._ys(ts, x, ps, tcls, hs, iys, ops)
        elif ret == ReturnValue.sigmay:
            output = self._sigmays(ts, x, ps, tcls, hs, iys, ops, nps)
        elif ret == ReturnValue.x0:
            output = self._x_rdata(x[0, :], tcls[0])
        elif ret == ReturnValue.x0_solver:
            output = x[0, :]
        elif ret == ReturnValue.tcl:
            output = tcls[0]
        elif ret in (ReturnValue.res, ReturnValue.chi2):
            obs_trafo = jax.vmap(
                lambda y, iy_trafo: (
                    jnp.array(
                        # needs to follow order in amici.jax.petab.SCALE_TO_INT
                        [y, safe_log(y), safe_log(y) / jnp.log(10)]
                    )
                    .at[iy_trafo]
                    .get()
                ),
            )
            ys_obj = obs_trafo(
                self._ys(ts, x, ps, tcls, hs, iys, ops), iy_trafos
            )
            m_obj = obs_trafo(my, iy_trafos)
            if ret == ReturnValue.chi2:
                sigma_obj = self._sigmays(ts, x, ps, tcls, hs, iys, ops, nps)
                chi2 = jnp.square((m_obj - ys_obj) / sigma_obj)
                chi2 = jnp.where(ts_mask, chi2, 0.0)
                output = jnp.sum(chi2)
            else:
                output = m_obj - ys_obj
                output = jnp.where(ts_mask, output, 0.0)
        else:
            raise NotImplementedError(f"Return value {ret} not implemented.")

        return output, stats

    def simulate_condition_unjitted(self, *args, **kwargs):
        """Deprecated alias for :meth:`simulate_experiment_unjitted`."""
        warnings.warn(
            "`simulate_condition_unjitted` has been renamed to "
            "`simulate_experiment_unjitted` and will be removed in a "
            "future release.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.simulate_experiment_unjitted(*args, **kwargs)

    @eqx.filter_jit
    def simulate_experiment(
        self,
        p: jt.Float[jt.Array, "P np"],
        ts_dyn: jt.Float[jt.Array, "P nt_dyn"],
        ts_posteq: jt.Float[jt.Array, "P nt_posteq"],
        my: jt.Float[jt.Array, "P nt"],
        iys: jt.Int[jt.Array, "P nt"],
        iy_trafos: jt.Int[jt.Array, "P nt"],
        ops: jt.Float[jt.Array, "P nt *nop"],
        nps: jt.Float[jt.Array, "P nt *nnp"],
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        root_finder: AbstractRootFinder,
        adjoint: diffrax.AbstractAdjoint,
        steady_state_event: Callable[
            ..., diffrax._custom_types.BoolScalarLike
        ],
        max_steps: int | jnp.int_,
        x_preeq: jt.Float[jt.Array, "*nx"] | None = None,
        h_preeq: jt.Bool[jt.Array, "*ne"] | None = None,
        mask_reinit: jt.Bool[jt.Array, "P *nx"] | None = None,
        x_reinit: jt.Float[jt.Array, "P *nx"] | None = None,
        init_override: jt.Float[jt.Array, "*nx"] | None = None,
        init_override_mask: jt.Bool[jt.Array, "*nx"] | None = None,
        ts_mask: jt.Bool[jt.Array, "P nt"] | None = None,
        h_mask: jt.Bool[jt.Array, "ne"] | None = None,
        t_zero: jt.Float[jt.Array, "P"] | None = None,
        posteq_mask: jt.Bool[jt.Array, "P"] | None = None,
        posteq_slots: tuple[bool, ...] | None = None,
        ret: ReturnValue = ReturnValue.llh,
    ) -> tuple[jt.Float[jt.Array, "*nt"], dict]:
        r"""
        Simulate an experiment (JIT-compiled version).

        This is the JIT-compiled version for optimal performance. For runtime type checking
        with beartype, use :meth:`simulate_experiment_unjitted` instead.

        Chains one ODE integration per experiment period (the leading axis,
        of static size ``P``, of ``p``/``ts_dyn``/``ts_posteq``/``my``/
        ``iys``/``iy_trafos``/``ops``/``nps``/``mask_reinit``/``x_reinit``/
        ``ts_mask``/``t_zero``); ``P == 1`` reduces to a single, non-chained
        simulation.

        :param p:
            parameters for simulation ordered according to ids in :ivar parameter_ids:, one row per
            experiment period.
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
        :param init_override:
            override model input e.g. with neural net outputs. If not provided, the inputs are not overridden.
        :param init_override_mask:
            mask for input override. If `True`, the corresponding input is replaced with the corresponding value from `init_override`.
        :param ts_mask:
            mask to remove (padded) time points. If `True`, the corresponding time point is used for the evaluation of
            the output. Only applied if ret is ReturnValue.llh, ReturnValue.nllhs, ReturnValue.res, or ReturnValue.chi2.
        :param h_mask:
            mask for heaviside variables. If `True`, the corresponding heaviside variable is updated during simulation, otherwise it
            it marked as 1.0.
        :param ret:
            which output to return. See :class:`ReturnValue` for available options.
        :return:
            output according to `ret` and general results/statistics
        """
        return self.simulate_experiment_unjitted(
            p,
            ts_dyn,
            ts_posteq,
            my,
            iys,
            iy_trafos,
            ops,
            nps,
            solver,
            controller,
            root_finder,
            adjoint,
            steady_state_event,
            max_steps,
            x_preeq,
            h_preeq,
            mask_reinit,
            x_reinit,
            init_override,
            init_override_mask,
            ts_mask,
            h_mask,
            t_zero,
            posteq_mask,
            posteq_slots,
            ret,
        )

    def simulate_condition(self, *args, **kwargs):
        """Deprecated alias for :meth:`simulate_experiment`."""
        warnings.warn(
            "`simulate_condition` has been renamed to `simulate_experiment` "
            "and will be removed in a future release.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.simulate_experiment(*args, **kwargs)

    @eqx.filter_jit
    def preequilibrate_condition(
        self,
        p: jt.Float[jt.Array, "np"] | None,
        x_reinit: jt.Float[jt.Array, "*nx"],
        mask_reinit: jt.Bool[jt.Array, "*nx"],
        h_mask: jt.Bool[jt.Array, "ne"],
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
        :param h_mask:
            mask for heaviside variables. If `True`, the corresponding heaviside variable is updated during simulation, otherwise it
            it marked as 1.0.
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

        if not h_mask.shape[0]:
            h_mask = jnp.ones(self.n_events, dtype=jnp.bool_)

        x0 = self._x0(t0, p)
        if x_reinit.shape[0]:
            x0 = jnp.where(mask_reinit, x_reinit, x0)
        tcl = self._tcl(x0, p)

        current_x = self._x_solver(x0)

        current_x, _, h, _ = self._handle_t0_event(
            t0,
            self._x_solver(x0),
            p,
            tcl,
            root_finder,
            self._root_cond_fn,
            self._delta_x,
            h_mask,
            jnp.array([]),
            {},
        )

        current_x, h, stats_preeq = eq(
            p,
            tcl,
            h,
            current_x,
            h_mask,
            solver,
            controller,
            root_finder,
            steady_state_event,
            diffrax.ODETerm(self._xdot),
            self._root_cond_fns(),
            self._root_cond_fn,
            self._delta_x,
            self._known_discs(p),
            max_steps,
        )

        return self._x_rdata(current_x, tcl), dict(stats_preeq=stats_preeq), h

    def _handle_t0_event(
        self,
        t0_next: float,
        y0_next: jt.Float[jt.Array, "nxs"],
        p: jt.Float[jt.Array, "np"],
        tcl: jt.Float[jt.Array, "ncl"],
        root_finder: AbstractRootFinder,
        root_cond_fn: Callable,
        delta_x: Callable,
        h_mask: jt.Bool[jt.Array, "ne"],
        h_prev: jt.Bool[jt.Array, "ne"],
        stats: dict,
    ):
        rf0 = self.event_initial_values - 0.5

        if h_prev.shape[0]:
            # `h_prev` is the heaviside state at wherever `y0_next` came
            # from (a preceding preequilibration, or the end of the
            # previous experiment period). It is not necessarily the
            # trigger state *at* `(t0_next, y0_next)`: a period boundary
            # (or the reinitialisation applied after preequilibration) may
            # have crossed an event's trigger threshold, or the dynamic
            # phase may use different parameters than preequilibration did
            # (e.g. a stimulus onset time that is inactive during
            # preequilibration), without the ODE integrator ever seeing it.
            # So the trigger condition is always re-evaluated below against
            # the actual incoming state, exactly as for a genuine t=0;
            # `h_prev` only supplies the pre-transition reference value.
            h = jnp.where(h_mask, h_prev, jnp.ones_like(h_prev))
            rf0 = jnp.where(h > 0.5, 0.5, -0.5)
        else:
            h = jnp.where(h_mask, jnp.heaviside(rf0, 0.0), jnp.ones_like(rf0))
        args = (p, tcl, h)
        rfx = root_cond_fn(t0_next, y0_next, args)
        roots_dir = jnp.sign(rfx - rf0)
        roots_found = jnp.sign(rfx) != jnp.sign(rf0)

        y0_next, h_next = _apply_event_assignments(
            roots_found,
            roots_dir,
            y0_next,
            p,
            tcl,
            h,
            delta_x,
        )

        roots_zero = jnp.isclose(
            rfx, 0.0, atol=root_finder.atol, rtol=root_finder.rtol
        )
        droot_dt = (
            # ∂root_cond_fn/∂t
            jax.jacfwd(root_cond_fn, argnums=0)(t0_next, y0_next, args)
            +
            # ∂root_cond_fn/∂y * ∂y/∂t
            jax.jacfwd(root_cond_fn, argnums=1)(t0_next, y0_next, args)
            @ self._xdot(t0_next, y0_next, args)
        )
        h_next = jnp.where(
            roots_zero,
            droot_dt >= 0.0,
            h_next,
        )

        if os.getenv("JAX_DEBUG") == "1":
            jax.debug.print(
                "h: {}, rf0: {}, rfx: {}, roots_found: {}, roots_dir: {}, h_next: {}",
                h,
                rf0,
                rfx,
                roots_found,
                roots_dir,
                h_next,
            )

        return y0_next, t0_next, h_next, stats


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
