from abc import abstractmethod

import diffrax
import equinox as eqx
import jax.numpy as jnp
import jax

# always use 64-bit precision. No-brainer on CPUs and GPUs don't make sense for stiff systems.
jax.config.update("jax_enable_x64", True)


class JAXModel(eqx.Module):
    """
    JAXModel provides an abstract base class for a JAX-based implementation of an AMICI model. Models inheriting from
    JAXModel must provide model specific implementations of abstract methods.
    """

    @abstractmethod
    def xdot(
        self,
        t: jnp.float_,
        x: jnp.ndarray,
        args: tuple[jnp.ndarray, jnp.ndarray],
    ) -> jnp.ndarray:
        """
        Right-hand side of the ODE system.

        :param t: time point
        :param x: state vector
        :param args: tuple of parameters, fixed parameters and total values for conservation laws
        :return:
            Derivative of the state vector at time point, same data structure as x.
        """
        ...

    @staticmethod
    @abstractmethod
    def _w(
        t: jnp.float_, x: jnp.ndarray, pk: jnp.ndarray, tcl: jnp.ndarray
    ) -> jnp.ndarray:
        """
        Compute the expressions (algebraic variables) of the model.

        :param t: time point
        :param x: state vector
        :param pk: parameters
        :param tcl: total values for conservation laws
        :return:
            Expression values.
        """
        ...

    @staticmethod
    @abstractmethod
    def x0(pk: jnp.ndarray) -> jnp.ndarray:
        """
        Compute the initial state vector.

        :param pk: parameters
        """
        ...

    @staticmethod
    @abstractmethod
    def x_solver(x: jnp.ndarray) -> jnp.ndarray:
        """
        Transform the full state vector to the reduced state vector for ODE solving.

        :param x:
            full state vector
        :return:
            reduced state vector
        """
        ...

    @staticmethod
    @abstractmethod
    def x_rdata(x: jnp.ndarray, tcl: jnp.ndarray) -> jnp.ndarray:
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

    @staticmethod
    @abstractmethod
    def tcl(x: jnp.ndarray, pk: jnp.ndarray) -> jnp.ndarray:
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
    def y(
        self, t: jnp.float_, x: jnp.ndarray, pk: jnp.ndarray, tcl: jnp.ndarray
    ) -> jnp.ndarray:
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
            observable
        """
        ...

    @staticmethod
    @abstractmethod
    def sigmay(y: jnp.ndarray, pk: jnp.ndarray) -> jnp.ndarray:
        """
        Compute the standard deviations of the observables.
        :param y:
            observable for the specified observable id
        :param pk:
            parameters
        :return:
            standard deviations of the observables
        """
        ...

    @abstractmethod
    def llh(
        self,
        t: jnp.float_,
        x: jnp.ndarray,
        pk: jnp.ndarray,
        tcl: jnp.ndarray,
        iy: int,
    ) -> jnp.float_:
        """
        Compute the log-likelihood of the observable for the specified observable id.
        :param t:
            time point
        :param x:
            state vector
        :param pk:
            parameters
        :param tcl:
            total values for conservation laws
        :param iy:
            observable id
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

    def _eq(self, p, tcl, x0, solver, controller, max_steps):
        sol = diffrax.diffeqsolve(
            diffrax.ODETerm(self.xdot),
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
        )
        return sol.ys[-1, :], sol.stats

    def _solve(self, p, ts, tcl, x0, solver, controller, max_steps, adjoint):
        sol = diffrax.diffeqsolve(
            diffrax.ODETerm(self.xdot),
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

    def _x_rdata(self, x, tcl):
        return jax.vmap(self.x_rdata, in_axes=(0, None))(x, tcl)

    def _outputs(self, ts, x, p, tcl, my, iys) -> jnp.float_:
        return jax.vmap(self.llh, in_axes=(0, 0, None, None, 0, 0))(
            ts, x, p, tcl, my, iys
        )

    def _y(self, ts, xs, p, tcl, iys):
        return jax.vmap(
            lambda t, x, p, tcl, iy: self.y(t, x, p, tcl).at[iy].get(),
            in_axes=(0, 0, None, None, 0),
        )(ts, xs, p, tcl, iys)

    def _sigmay(self, ts, xs, p, tcl, iys):
        return jax.vmap(
            lambda t, x, p, tcl, iy: self.sigmay(self.y(t, x, p, tcl), p)
            .at[iy]
            .get(),
            in_axes=(0, 0, None, None, 0),
        )(ts, xs, p, tcl, iys)

    # @eqx.filter_jit
    def simulate_condition(
        self,
        p: jnp.ndarray,
        p_preeq: jnp.ndarray,
        ts_preeq: jnp.ndarray,
        ts_dyn: jnp.ndarray,
        ts_posteq: jnp.ndarray,
        my: jnp.ndarray,
        iys: jnp.ndarray,
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        adjoint: diffrax.AbstractAdjoint,
        max_steps: int,
        ret: str = "llh",
    ):
        # Pre-equilibration
        if p_preeq.shape[0] > 0:
            x0 = self.x0(p_preeq)
            tcl = self.tcl(x0, p_preeq)
            current_x = self.x_solver(x0)
            current_x, stats_preeq = self._eq(
                p_preeq, tcl, current_x, solver, controller, max_steps
            )
            # update tcl with new parameters
            tcl = self.tcl(self.x_rdata(current_x, tcl), p)
        else:
            x0 = self.x0(p)
            current_x = self.x_solver(x0)
            stats_preeq = None

            tcl = self.tcl(x0, p)
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

        llhs = self._outputs(ts, x, p, tcl, my, iys)
        llh = -jnp.sum(llhs)
        return {
            "llh": llh,
            "llhs": llhs,
            "x": self._x_rdata(x, tcl),
            "x_solver": x,
            "y": self._y(ts, x, p, tcl, iys),
            "sigmay": self._sigmay(ts, x, p, tcl, iys),
            "x0": self.x_rdata(x_preq[-1, :], tcl),
            "x0_solver": x_preq[-1, :],
            "tcl": tcl,
            "res": self._y(ts, x, p, tcl, iys) - my,
        }[ret], dict(
            stats_preeq=stats_preeq,
            stats_dyn=stats_dyn,
            stats_posteq=stats_posteq,
        )
