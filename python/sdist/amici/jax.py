from abc import abstractmethod
from dataclasses import dataclass

import diffrax
import jax.numpy as jnp
import numpy as np
import jax
from functools import partial

import amici


class JAXModel(object):
    _unscale_funs = {
        amici.ParameterScaling.none: lambda x: x,
        amici.ParameterScaling.ln: lambda x: jnp.exp(x),
        amici.ParameterScaling.log10: lambda x: jnp.power(10, x)
    }

    @abstractmethod
    def xdot(self, t, x, args):
        ...

    @abstractmethod
    def _w(self, x, p, k, tcl):
        ...

    @abstractmethod
    def x0(self, p, k):
        ...

    @abstractmethod
    def y(self, x, p, k, tcl):
        ...

    @abstractmethod
    def sigmay(self, y, p, k):
        ...

    @abstractmethod
    def Jy(self, y, my, sigmay):
        ...

    def unscale_p(self, p, pscale):
        return jnp.stack([
            self._unscale_funs[pscale_i](p_i)
            for p_i, pscale_i in zip(p, pscale)
        ])

    def get_solver(self):
        return JAXSolver(model=self)


class JAXSolver(object):
    def __init__(self, model: JAXModel):
        self.model: JAXModel = model
        self.solver: diffrax.AbstractSolver = diffrax.Tsit5()
        self.atol: float = 1e-8
        self.rtol: float = 1e-8
        self.sensi_mode: amici.SensitivityMethod = \
            amici.SensitivityMethod.adjoint
        self.sensi_order: amici.SensitivityOrder = \
            amici.SensitivityOrder.none

    def _solve(self, ts, p, k):
        y0 = self.model.x0(p, k)
        tcl = 0
        sol = diffrax.diffeqsolve(
            diffrax.ODETerm(self.model.xdot),
            self.solver,
            args=(p, k, tcl),
            t0=ts[0],
            t1=ts[-1],
            dt0=ts[1] - ts[0],
            y0=y0,
            stepsize_controller=diffrax.PIDController(
                rtol=self.rtol,
                atol=self.atol
            ),
            saveat=diffrax.SaveAt(ts=ts)
        )
        return sol.ys

    def _obs(self, x, p, k, tcl):
        y = jnp.apply_along_axis(
            lambda x: self.model.y(x, p, k, tcl),
            axis=1,
            arr=x
        )
        return y

    def _sigmay(self, obs, p, k):
        sigmay = jnp.apply_along_axis(
            lambda y: self.model.sigmay(y, p, k),
            axis=1,
            arr=obs
        )
        return sigmay

    def _loss(self, obs: jnp.ndarray, sigmay: jnp.ndarray, my: np.ndarray):
        llh = - jnp.sum(jnp.stack(
            [self.model.Jy(obs[i, :], my[i, :], sigmay[i, :])
             for i in range(my.shape[0])]
        ))
        return llh

    @partial(jax.jit, static_argnames=('self', 'ts', 'k', 'my', 'pscale'))
    def run(self,
            ts: tuple,
            p: jnp.ndarray,
            k: tuple,
            my: tuple,
            pscale: tuple):
        ps = self.model.unscale_p(p, pscale)
        x = self._solve(ts, ps, k)
        tcl = 0
        obs = self._obs(x, ps, k, tcl)
        my_r = np.asarray(my).reshape(obs.shape)
        sigmay = self._sigmay(obs, ps, k)
        llh = self._loss(obs, sigmay, my_r)
        return llh, (x, obs)

    @partial(jax.jit, static_argnames=('self', 'ts', 'k', 'my', 'pscale'))
    def srun(self,
             ts: tuple,
             p: jnp.ndarray,
             k: tuple,
             my: tuple,
             pscale: tuple):
        (llh, (x, obs)), sllh = (jax.value_and_grad(self.run, 1, True))(
            ts, p, k, my, pscale
        )
        return llh, sllh, (x, obs)

    @partial(jax.jit, static_argnames=('self', 'ts', 'k', 'my', 'pscale'))
    def s2run(self,
              ts: tuple,
              p: jnp.ndarray,
              k: tuple,
              my: tuple,
              pscale: tuple):
        (llh, (x, obs)), sllh = (jax.value_and_grad(self.run, 1, True))(
            ts, p, k, my, pscale
        )
        s2llh, (x, obs) = jax.jacfwd(jax.grad(self.run, 1, True), 1, True)(
            ts, p, k, my, pscale
        )
        return llh, sllh, s2llh, (x, obs)


def runAmiciSimulationJAX(model: JAXModel,
                          solver: JAXSolver,
                          edata: amici.ExpData):
    ts = tuple(edata.getTimepoints())
    p = jnp.asarray(edata.parameters)
    k = tuple(edata.fixedParameters)
    my = tuple(edata.getObservedData())
    pscale = tuple(edata.pscale)

    rdata_kwargs = dict()

    if solver.sensi_order == amici.SensitivityOrder.none:
        rdata_kwargs['llh'], (rdata_kwargs['x'], rdata_kwargs['y']) = \
            solver.run(ts, p, k, my, pscale)
    elif solver.sensi_order == amici.SensitivityOrder.first:
        rdata_kwargs['llh'], rdata_kwargs['sllh'], (
            rdata_kwargs['x'], rdata_kwargs['y']
        ) = solver.srun(ts, p, k, my, pscale)
    elif solver.sensi_order == amici.SensitivityOrder.second:
        rdata_kwargs['llh'], rdata_kwargs['sllh'], rdata_kwargs['s2llh'], (
            rdata_kwargs['x'], rdata_kwargs['y']
        ) = solver.s2run(ts, p, k, my, pscale)

    return ReturnDataJAX(**rdata_kwargs)


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

    def __init__(self, *args, **kwargs):
        super(ReturnDataJAX, self).__init__(*args, **kwargs)
        self.__dict__ = self
