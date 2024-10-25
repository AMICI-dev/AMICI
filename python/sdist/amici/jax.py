from abc import abstractmethod
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor

import diffrax
import equinox as eqx
import jax.numpy as jnp
import numpy as np
import jax
from collections.abc import Iterable

import amici

jax.config.update("jax_enable_x64", True)


class JAXModel(eqx.Module):
    _unscale_funs = {
        amici.ParameterScaling.none: lambda x: x,
        amici.ParameterScaling.ln: lambda x: jnp.exp(x),
        amici.ParameterScaling.log10: lambda x: jnp.power(10, x),
    }
    solver: diffrax.AbstractSolver
    controller: diffrax.AbstractStepSizeController
    atol: float
    rtol: float
    pcoeff: float
    icoeff: float
    dcoeff: float
    maxsteps: int
    term: diffrax.ODETerm

    def __init__(self):
        self.solver = diffrax.Kvaerno5()
        self.atol: float = 1e-8
        self.rtol: float = 1e-8
        self.pcoeff: float = 0.4
        self.icoeff: float = 0.3
        self.dcoeff: float = 0.0
        self.maxsteps: int = 2**14
        self.controller = diffrax.PIDController(
            rtol=self.rtol,
            atol=self.atol,
            pcoeff=self.pcoeff,
            icoeff=self.icoeff,
            dcoeff=self.dcoeff,
        )
        self.term = diffrax.ODETerm(self.xdot)

    @staticmethod
    @abstractmethod
    def xdot(t, x, args): ...

    @staticmethod
    @abstractmethod
    def _w(t, x, p, k, tcl): ...

    @staticmethod
    @abstractmethod
    def x0(p, k): ...

    @staticmethod
    @abstractmethod
    def x_solver(x): ...

    @staticmethod
    @abstractmethod
    def x_rdata(x, tcl): ...

    @staticmethod
    @abstractmethod
    def tcl(x, p, k): ...

    @staticmethod
    @abstractmethod
    def y(t, x, p, k, tcl): ...

    @staticmethod
    @abstractmethod
    def sigmay(y, p, k): ...

    @staticmethod
    @abstractmethod
    def Jy(y, my, sigmay): ...

    def unscale_p(self, p, pscale):
        return jax.vmap(
            lambda p_i, pscale_i: jnp.stack(
                (p_i, jnp.exp(p_i), jnp.power(10, p_i))
            )
            .at[pscale_i]
            .get()
        )(p, pscale)

    def _preeq(self, p, k):
        x0 = self.x_solver(self.x0(p, k))
        tcl = self.tcl(x0, p, k)
        return self._eq(p, k, tcl, x0)

    def _posteq(self, p, k, x, tcl):
        return self._eq(p, k, tcl, x)

    def _eq(self, p, k, tcl, x0):
        sol = diffrax.diffeqsolve(
            self.term,
            self.solver,
            args=(p, k, tcl),
            t0=0.0,
            t1=jnp.inf,
            dt0=None,
            y0=x0,
            stepsize_controller=self.controller,
            max_steps=self.maxsteps,
            event=diffrax.Event(cond_fn=diffrax.steady_state_event()),
        )
        return sol.ys

    def _solve(self, ts, p, k, x0, checkpointed):
        tcl = self.tcl(x0, p, k)
        sol = diffrax.diffeqsolve(
            self.term,
            self.solver,
            args=(p, k, tcl),
            t0=0.0,
            t1=ts[-1],
            dt0=None,
            y0=self.x_solver(x0),
            stepsize_controller=self.controller,
            max_steps=self.maxsteps,
            adjoint=diffrax.RecursiveCheckpointAdjoint()
            if checkpointed
            else diffrax.DirectAdjoint(),
            saveat=diffrax.SaveAt(ts=ts),
            throw=False,
        )
        return sol.ys, tcl, sol.stats

    def _obs(self, ts, x, p, k, tcl):
        return jax.vmap(self.y, in_axes=(0, 0, None, None, None))(
            ts, x, p, k, tcl
        )

    def _sigmay(self, obs, p, k):
        return jax.vmap(self.sigmay, in_axes=(0, None, None))(obs, p, k)

    def _x_rdata(self, x, tcl):
        return jax.vmap(self.x_rdata, in_axes=(0, None))(x, tcl)

    def _loss(self, obs: jnp.ndarray, sigmay: jnp.ndarray, my: np.ndarray):
        loss_fun = jax.vmap(self.Jy, in_axes=(0, 0, 0))
        return -jnp.sum(loss_fun(obs, my, sigmay))

    def _run(
        self,
        ts: np.ndarray,
        ts_dyn: np.ndarray,
        p: np.ndarray,
        k: jnp.ndarray,
        k_preeq: jnp.ndarray,
        my: jnp.ndarray,
        pscale: np.ndarray,
        checkpointed=True,
        dynamic=True,
    ):
        ps = self.unscale_p(p, pscale)

        # Pre-equilibration
        if k_preeq.shape[0] > 0:
            x0 = self._preeq(ps, k_preeq)
        else:
            x0 = self.x0(ps, k)

        # Dynamic simulation
        if dynamic == "true":
            x, tcl, stats = self._solve(
                ts_dyn, ps, k, x0, checkpointed=checkpointed
            )
        else:
            x = tuple(
                jnp.array([x0_i] * len(ts_dyn)) for x0_i in self.x_solver(x0)
            )
            tcl = self.tcl(x0, ps, k)
            stats = None

        # Post-equilibration
        if len(ts) > len(ts_dyn):
            if len(ts_dyn) > 0:
                x_final = tuple(x_i[-1] for x_i in x)
            else:
                x_final = self.x_solver(x0)
            x_posteq = self._posteq(ps, k, x_final, tcl)
            x_posteq = tuple(
                jnp.array([x0_i] * (len(ts) - len(ts_dyn)))
                for x0_i in x_posteq
            )
            if len(ts_dyn) > 0:
                x = tuple(
                    jnp.concatenate((x_i, x_posteq_i), axis=0)
                    for x_i, x_posteq_i in zip(x, x_posteq)
                )
            else:
                x = x_posteq

        obs = self._obs(ts, x, ps, k, tcl)
        my_r = my.reshape((len(ts), -1))
        sigmay = self._sigmay(obs, ps, k)
        llh = self._loss(obs, sigmay, my_r)
        x_rdata = self._x_rdata(x, tcl)
        return llh, (x_rdata, obs, stats)

    @eqx.filter_jit
    def run(
        self,
        ts: np.ndarray,
        ts_dyn: np.ndarray,
        p: jnp.ndarray,
        k: np.ndarray,
        k_preeq: np.ndarray,
        my: np.ndarray,
        pscale: np.ndarray,
        dynamic=True,
    ):
        return self._run(
            ts, ts_dyn, p, k, k_preeq, my, pscale, dynamic=dynamic
        )

    @eqx.filter_jit
    def srun(
        self,
        ts: np.ndarray,
        ts_dyn: np.ndarray,
        p: jnp.ndarray,
        k: np.ndarray,
        k_preeq: np.ndarray,
        my: np.ndarray,
        pscale: np.ndarray,
        dynamic=True,
    ):
        (llh, (x, obs, stats)), sllh = (
            jax.value_and_grad(self._run, 2, True)
        )(ts, ts_dyn, p, k, k_preeq, my, pscale, dynamic=dynamic)
        return llh, sllh, (x, obs, stats)

    @eqx.filter_jit
    def s2run(
        self,
        ts: np.ndarray,
        ts_dyn: np.ndarray,
        p: jnp.ndarray,
        k: np.ndarray,
        k_preeq: np.ndarray,
        my: np.ndarray,
        pscale: np.ndarray,
        dynamic=True,
    ):
        (llh, (x, obs, stats)), sllh = (
            jax.value_and_grad(self._run, 2, True)
        )(ts, ts_dyn, p, k, k_preeq, my, pscale, dynamic=dynamic)

        s2llh = jax.hessian(self._run, 2, True)(
            ts, ts_dyn, p, k, k_preeq, my, pscale, dynamic=dynamic
        )

        return llh, sllh, s2llh, (x, obs, stats)

    def run_simulation(
        self, edata: amici.ExpData, sensitivity_order: amici.SensitivityOrder
    ):
        ts = np.asarray(edata.getTimepoints())
        p = jnp.asarray(edata.parameters)
        k = np.asarray(edata.fixedParameters)
        k_preeq = np.asarray(edata.fixedParametersPreequilibration)
        my = np.asarray(edata.getObservedData())
        pscale = np.asarray(edata.pscale)
        ts_dyn = ts[np.isfinite(ts)]
        dynamic = "true" if len(ts_dyn) and np.max(ts_dyn) > 0 else "false"

        rdata_kwargs = dict()

        if sensitivity_order == amici.SensitivityOrder.none:
            (
                rdata_kwargs["llh"],
                (rdata_kwargs["x"], rdata_kwargs["y"], rdata_kwargs["stats"]),
            ) = self.run(
                ts, ts_dyn, p, k, k_preeq, my, pscale, dynamic=dynamic
            )
        elif sensitivity_order == amici.SensitivityOrder.first:
            (
                rdata_kwargs["llh"],
                rdata_kwargs["sllh"],
                (rdata_kwargs["x"], rdata_kwargs["y"], rdata_kwargs["stats"]),
            ) = self.srun(
                ts, ts_dyn, p, k, k_preeq, my, pscale, dynamic=dynamic
            )
        elif sensitivity_order == amici.SensitivityOrder.second:
            (
                rdata_kwargs["llh"],
                rdata_kwargs["sllh"],
                rdata_kwargs["s2llh"],
                (rdata_kwargs["x"], rdata_kwargs["y"], rdata_kwargs["stats"]),
            ) = self.s2run(
                ts, ts_dyn, p, k, k_preeq, my, pscale, dynamic=dynamic
            )

        for field in rdata_kwargs.keys():
            if field == "llh":
                rdata_kwargs[field] = np.float64(rdata_kwargs[field])
            elif field not in ["sllh", "s2llh"]:
                rdata_kwargs[field] = np.asarray(rdata_kwargs[field]).T
                if rdata_kwargs[field].ndim == 1:
                    rdata_kwargs[field] = np.expand_dims(
                        rdata_kwargs[field], 1
                    )

        return ReturnDataJAX(**rdata_kwargs)

    def run_simulations(
        self,
        edatas: Iterable[amici.ExpData],
        sensitivity_order: amici.SensitivityOrder = amici.SensitivityOrder.none,
        num_threads: int = 1,
    ):
        fun = eqx.Partial(
            self.run_simulation, sensitivity_order=sensitivity_order
        )
        if num_threads > 1:
            with ThreadPoolExecutor(max_workers=num_threads) as pool:
                results = pool.map(fun, edatas)
        else:
            results = map(fun, edatas)
        return list(results)


@dataclass
class ReturnDataJAX(dict):
    x: np.array = None
    sx: np.array = None
    y: np.array = None
    sy: np.array = None
    sigmay: np.array = None
    ssigmay: np.array = None
    llh: np.array = None
    sllh: np.array = None
    stats: dict = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__dict__ = self
