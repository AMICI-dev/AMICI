from abc import abstractmethod

import diffrax
import equinox as eqx
import jax.numpy as jnp
import numpy as np
import jax

# always use 64-bit precision. No-brainer on CPUs and GPUs don't make sense for stiff systems.
jax.config.update("jax_enable_x64", True)


class JAXModel(eqx.Module):
    """
    JAXModel provides an abstract base class for a JAX-based implementation of an AMICI model. Models inheriting from
    JAXModel must provide model specific implementations of abstract methods.
    """

    @staticmethod
    @abstractmethod
    def xdot(
        t: jnp.float_, x: jnp.ndarray, args: tuple[jnp.ndarray, jnp.ndarray]
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

    def _preeq(self, p, solver, controller, max_steps):
        """
        Pre-equilibration of the model.
        :param p:
            parameters
        :return:
            Initial state vector
        """
        x0 = self.x_solver(self.x0(p))
        tcl = self.tcl(x0, p)
        return self._eq(p, tcl, x0, solver, controller, max_steps)

    def _posteq(self, p, x, tcl, solver, controller, max_steps):
        return self._eq(p, tcl, x, solver, controller, max_steps)

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
            event=diffrax.Event(cond_fn=diffrax.steady_state_event()),
        )
        return sol.ys[-1, :]

    def _solve(self, ts, p, x0, solver, controller, max_steps):
        tcl = self.tcl(x0, p)
        sol = diffrax.diffeqsolve(
            diffrax.ODETerm(self.xdot),
            solver,
            args=(p, tcl),
            t0=0.0,
            t1=ts[-1],
            dt0=None,
            y0=self.x_solver(x0),
            stepsize_controller=controller,
            max_steps=max_steps,
            adjoint=diffrax.RecursiveCheckpointAdjoint(),
            saveat=diffrax.SaveAt(ts=ts),
            throw=False,
        )
        return sol.ys, tcl, sol.stats

    def _x_rdata(self, x, tcl):
        return jax.vmap(self.x_rdata, in_axes=(0, None))(x, tcl)

    def _outputs(self, ts, x, p, tcl, my, iys) -> jnp.float_:
        return jax.vmap(self.llh, in_axes=(0, 0, None, None, 0, 0))(
            ts, x, p, tcl, my, iys
        )

    # @eqx.filter_jit
    def simulate_condition(
        self,
        ts: np.ndarray,
        ts_dyn: np.ndarray,
        my: np.ndarray,
        iys: np.ndarray,
        p: jnp.ndarray,
        p_preeq: jnp.ndarray,
        dynamic: bool,
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        max_steps: int,
    ):
        # Pre-equilibration
        if p_preeq.shape[0] > 0:
            x0 = self._preeq(p_preeq, solver, controller, max_steps)
        else:
            x0 = self.x0(p)

        # Dynamic simulation
        if dynamic:
            x, tcl, stats = self._solve(
                ts_dyn, p, x0, solver, controller, max_steps
            )
        else:
            x = jnp.repeat(
                self.x_solver(x0).reshape(1, -1),
                len(ts_dyn),
                axis=0,
            )
            tcl = self.tcl(x0, p)
            stats = None

        # Post-equilibration
        if len(ts) > len(ts_dyn):
            if len(ts_dyn) > 0:
                x_final = x[-1, :]
            else:
                x_final = self.x_solver(x0)
            x_posteq = self._posteq(
                p, x_final, tcl, solver, controller, max_steps
            )
            x_posteq = jnp.repeat(
                x_posteq.reshape(1, -1),
                len(ts) - len(ts_dyn),
                axis=0,
            )
            if len(ts_dyn) > 0:
                x = jnp.concatenate((x, x_posteq), axis=0)
            else:
                x = x_posteq

        outputs = self._outputs(ts, x, p, tcl, my, iys)
        llh = -jnp.sum(outputs[:, 0])
        obs = outputs[:, 1]
        sigmay = outputs[:, 2]
        x_rdata = jnp.stack(self._x_rdata(x, tcl), axis=1)
        return llh, dict(llh=llh, x=x_rdata, y=obs, sigmay=sigmay, stats=stats)
