"""Model simulation using JAX."""

# ruff: noqa: F821 F722

from abc import abstractmethod

import diffrax
import equinox as eqx
import jax.numpy as jnp
import jax
import jaxtyping as jt


class JAXModel(eqx.Module):
    """
    JAXModel provides an abstract base class for a JAX-based implementation of an AMICI model. The class implements
    routines for simulation and evaluation of derived quantities, model specific implementations need to be provided by
    classes inheriting from JAXModel.
    """

    MODEL_API_VERSION = "0.0.1"
    api_version: str

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
        args: tuple[jt.Float[jt.Array, "np"], jt.Float[jt.Array, "ncl"]],
    ) -> jt.Float[jt.Array, "nxs"]:
        """
        Right-hand side of the ODE system.

        :param t: time point
        :param x: state vector
        :param args: tuple of parameters and total values for conservation laws
        :return:
            Temporal derivative of the state vector x at time point t.
        """
        ...

    @abstractmethod
    def _w(
        self,
        t: jt.Float[jt.Array, ""],
        x: jt.Float[jt.Array, "nxs"],
        pk: jt.Float[jt.Array, "np"],
        tcl: jt.Float[jt.Array, "ncl"],
    ) -> jt.Float[jt.Array, "nw"]:
        """
        Compute the expressions, i.e. derived quantities that are used in other parts of the model.

        :param t: time point
        :param x: state vector
        :param pk: parameters
        :param tcl: total values for conservation laws
        :return:
            Expression values.
        """
        ...

    @abstractmethod
    def _x0(self, pk: jt.Float[jt.Array, "np"]) -> jt.Float[jt.Array, "nx"]:
        """
        Compute the initial state vector.

        :param pk: parameters
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
        self, x: jt.Float[jt.Array, "nx"], pk: jt.Float[jt.Array, "np"]
    ) -> jt.Float[jt.Array, "ncl"]:
        """
        Compute the total values for conservation laws.

        :param x:
            state vector
        :param pk:
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
        pk: jt.Float[jt.Array, "np"],
        tcl: jt.Float[jt.Array, "ncl"],
    ) -> jt.Float[jt.Array, "ny"]:
        """
        Compute the observables.

        :param t:
            time point
        :param x:
            state vector
        :param pk:
            parameters
        :param tcl:
            total values for conservation laws
        :return:
            observables
        """
        ...

    @abstractmethod
    def _sigmay(
        self, y: jt.Float[jt.Array, "ny"], pk: jt.Float[jt.Array, "np"]
    ) -> jt.Float[jt.Array, "ny"]:
        """
        Compute the standard deviations of the observables.

        :param y:
            observables
        :param pk:
            parameters
        :return:
            standard deviations of the observables
        """
        ...

    @abstractmethod
    def _nllh(
        self,
        t: jt.Float[jt.Scalar, ""],
        x: jt.Float[jt.Array, "nxs"],
        pk: jt.Float[jt.Array, "np"],
        tcl: jt.Float[jt.Array, "ncl"],
        my: jt.Float[jt.Array, ""],
        iy: jt.Int[jt.Array, ""],
    ) -> jt.Float[jt.Scalar, ""]:
        """
        Compute the negative log-likelihood of the observable for the specified observable index.

        :param t:
            time point
        :param x:
            state vector
        :param pk:
            parameters
        :param tcl:
            total values for conservation laws
        :param my:
            observed data
        :param iy:
            observable index
        :return:
            log-likelihood of the observable
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

    def _eq(
        self,
        p: jt.Float[jt.Array, "np"],
        tcl: jt.Float[jt.Array, "ncl"],
        x0: jt.Float[jt.Array, "nxs"],
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        max_steps: jnp.int_,
    ) -> tuple[jt.Float[jt.Array, "1 nxs"], dict]:
        """
        Solve the steady state equation.

        :param p:
            parameters
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
        :return:
        """
        sol = diffrax.diffeqsolve(
            diffrax.ODETerm(self._xdot),
            solver,
            args=(p, tcl),
            t0=0.0,
            t1=jnp.inf,
            dt0=None,
            y0=x0,
            stepsize_controller=controller,
            max_steps=max_steps,
            adjoint=diffrax.DirectAdjoint(),
            event=diffrax.Event(cond_fn=diffrax.steady_state_event()),
            throw=False,
        )
        return sol.ys[-1, :], sol.stats

    def _solve(
        self,
        p: jt.Float[jt.Array, "np"],
        ts: jt.Float[jt.Array, "nt_dyn"],
        tcl: jt.Float[jt.Array, "ncl"],
        x0: jt.Float[jt.Array, "nxs"],
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        max_steps: jnp.int_,
        adjoint: diffrax.AbstractAdjoint,
    ) -> tuple[jt.Float[jt.Array, "nt nxs"], dict]:
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
            solution at time points ts and statistics
        """
        sol = diffrax.diffeqsolve(
            diffrax.ODETerm(self._xdot),
            solver,
            args=(p, tcl),
            t0=0.0,
            t1=ts[-1],
            dt0=None,
            y0=x0,
            stepsize_controller=controller,
            max_steps=max_steps,
            adjoint=adjoint,
            saveat=diffrax.SaveAt(ts=ts),
            throw=False,
        )
        return sol.ys, sol.stats

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
        mys: jt.Float[jt.Array, "nt"],
        iys: jt.Int[jt.Array, "nt"],
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
        :param mys:
            observed data
        :param iys:
            observable indices
        :return:
            negative log-likelihoods of the observables
        """
        return jax.vmap(self._nllh, in_axes=(0, 0, None, None, 0, 0))(
            ts, xs, p, tcl, mys, iys
        )

    def _ys(
        self,
        ts: jt.Float[jt.Array, "nt"],
        xs: jt.Float[jt.Array, "nt nxs"],
        p: jt.Float[jt.Array, "np"],
        tcl: jt.Float[jt.Array, "ncl"],
        iys: jt.Float[jt.Array, "nt"],
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
        :param iys:
            observable indices
        :return:
            observables
        """
        return jax.vmap(
            lambda t, x, p, tcl, iy: self._y(t, x, p, tcl).at[iy].get(),
            in_axes=(0, 0, None, None, 0),
        )(ts, xs, p, tcl, iys)

    def _sigmays(
        self,
        ts: jt.Float[jt.Array, "nt"],
        xs: jt.Float[jt.Array, "nt nxs"],
        p: jt.Float[jt.Array, "np"],
        tcl: jt.Float[jt.Array, "ncl"],
        iys: jt.Int[jt.Array, "nt"],
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
        :param iys:
            observable indices
        :return:
            standard deviations of the observables
        """
        return jax.vmap(
            lambda t, x, p, tcl, iy: self._sigmay(self._y(t, x, p, tcl), p)
            .at[iy]
            .get(),
            in_axes=(0, 0, None, None, 0),
        )(ts, xs, p, tcl, iys)

    @eqx.filter_jit
    def simulate_condition(
        self,
        p: jt.Float[jt.Array, "np"],
        p_preeq: jt.Float[jt.Array, "*np"],
        ts_preeq: jt.Float[jt.Array, "nt_preeq"],
        ts_dyn: jt.Float[jt.Array, "nt_dyn"],
        ts_posteq: jt.Float[jt.Array, "nt_posteq"],
        my: jt.Float[jt.Array, "nt"],
        iys: jt.Int[jt.Array, "nt"],
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        adjoint: diffrax.AbstractAdjoint,
        max_steps: int | jnp.int_,
        ret: str = "llh",
    ) -> tuple[jt.Float[jt.Array, "nt *nx"] | jnp.float_, dict]:
        r"""
        Simulate a condition.

        :param p:
            parameters for simulation ordered according to ids in :ivar parameter_ids:
        :param p_preeq:
            parameters for pre-equilibration ordered according to ids in :ivar parameter_ids:. May be empty to
            disable pre-equilibration.
        :param ts_preeq:
            time points for pre-equilibration. Usually valued 0.0, but needs to be shaped according to
            the number of observables that are evaluated after pre-equilibration.
        :param ts_dyn:
            time points for dynamic simulation. Usually valued > 0.0 and sorted in monotonically increasing order.
            Duplicate time points are allowed to facilitate the evaluation of multiple observables at specific time
            points.
        :param ts_posteq:
            time points for post-equilibration. Usually valued \Infty, but needs to be shaped according to
            the number of observables that are evaluated after post-equilibration.
        :param my:
            observed data
        :param iys:
            indices of the observables according to ordering in :ivar observable_ids:
        :param solver:
            ODE solver
        :param controller:
            step size controller
        :param adjoint:
            adjoint method. Recommended values are `diffrax.DirectAdjoint()` for jax.jacfwd (with vector-valued
            outputs) and  `diffrax.RecursiveCheckpointAdjoint()` for jax.grad (for scalar-valued outputs).
        :param max_steps:
            maximum number of solver steps
        :param ret:
            which output to return. Valid values are
                - `llh`: log-likelihood (default)
                - `nllhs`: negative log-likelihood at each time point
                - `x0`: full initial state vector (after pre-equilibration)
                - `x0_solver`: reduced initial state vector (after pre-equilibration)
                - `x`: full state vector
                - `x_solver`: reduced state vector
                - `y`: observables
                - `sigmay`: standard deviations of the observables
                - `tcl`: total values for conservation laws (at final timepoint)
                - `res`: residuals (observed - simulated)
        :return:
            output according to `ret` and statistics
        """
        # Pre-equilibration
        if p_preeq.shape[0] > 0:
            x0 = self._x0(p_preeq)
            tcl = self._tcl(x0, p_preeq)
            current_x = self._x_solver(x0)
            current_x, stats_preeq = self._eq(
                p_preeq, tcl, current_x, solver, controller, max_steps
            )
            # update tcl with new parameters
            tcl = self._tcl(self._x_rdata(current_x, tcl), p)
        else:
            x0 = self._x0(p)
            current_x = self._x_solver(x0)
            stats_preeq = None

            tcl = self._tcl(x0, p)
        x_preq = jnp.repeat(
            current_x.reshape(1, -1), ts_preeq.shape[0], axis=0
        )

        # Dynamic simulation
        if ts_dyn.shape[0] > 0:
            x_dyn, stats_dyn = self._solve(
                p,
                ts_dyn,
                tcl,
                current_x,
                solver,
                controller,
                max_steps,
                adjoint,
            )
            current_x = x_dyn[-1, :]
        else:
            x_dyn = jnp.repeat(
                current_x.reshape(1, -1), ts_dyn.shape[0], axis=0
            )
            stats_dyn = None

        # Post-equilibration
        if ts_posteq.shape[0] > 0:
            current_x, stats_posteq = self._eq(
                p, tcl, current_x, solver, controller, max_steps
            )
        else:
            stats_posteq = None

        x_posteq = jnp.repeat(
            current_x.reshape(1, -1), ts_posteq.shape[0], axis=0
        )

        ts = jnp.concatenate((ts_preeq, ts_dyn, ts_posteq), axis=0)
        x = jnp.concatenate((x_preq, x_dyn, x_posteq), axis=0)

        nllhs = self._nllhs(ts, x, p, tcl, my, iys)
        llh = -jnp.sum(nllhs)
        return {
            "llh": llh,
            "nllhs": nllhs,
            "x": self._x_rdatas(x, tcl),
            "x_solver": x,
            "y": self._ys(ts, x, p, tcl, iys),
            "sigmay": self._sigmays(ts, x, p, tcl, iys),
            "x0": self._x_rdata(x[0, :], tcl),
            "x0_solver": x[0, :],
            "tcl": tcl,
            "res": self._ys(ts, x, p, tcl, iys) - my,
        }[ret], dict(
            ts=ts,
            x=x,
            stats_preeq=stats_preeq,
            stats_dyn=stats_dyn,
            stats_posteq=stats_posteq,
        )
