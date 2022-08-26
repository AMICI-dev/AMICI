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
    def x_solver(self, x):
        ...

    @abstractmethod
    def x_rdata(self, x, tcl):
        ...

    @abstractmethod
    def tcl(self, x, p, k):
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
        x0 = self.model.x0(p, k)
        tcl = self.model.tcl(x0, p, k)
        sol = diffrax.diffeqsolve(
            diffrax.ODETerm(self.model.xdot),
            self.solver,
            args=(p, k, tcl),
            t0=ts[0],
            t1=ts[-1],
            dt0=ts[1] - ts[0],
            y0=self.model.x_solver(x0),
            stepsize_controller=diffrax.PIDController(
                rtol=self.rtol,
                atol=self.atol
            ),
            saveat=diffrax.SaveAt(ts=ts)
        )
        return sol.ys, tcl

    def _obs(self, x, p, k, tcl):
        return jax.vmap(self.model.y, in_axes=(0, None, None, None))(
            x, p, k, tcl
        )

    def _sigmay(self, obs, p, k):
        return jax.vmap(self.model.sigmay, in_axes=(0, None, None))(obs, p, k)

    def _x_rdata(self, x, tcl):
        return jax.vmap(self.model.x_rdata, in_axes=(0, None))(x, tcl)

    def _loss(self, obs: jnp.ndarray, sigmay: jnp.ndarray, my: np.ndarray):
        loss_fun = jax.vmap(self.model.Jy, in_axes=(0, 0, 0))
        return - jnp.sum(loss_fun(obs, my, sigmay))

    def _run(self,
             ts: tuple,
             p: jnp.ndarray,
             k: tuple,
             my: tuple,
             pscale: tuple):
        ps = self.model.unscale_p(p, pscale)
        x, tcl = self._solve(ts, ps, k)
        obs = self._obs(x, ps, k, tcl)
        my_r = np.asarray(my).reshape(obs.shape)
        sigmay = self._sigmay(obs, ps, k)
        llh = self._loss(obs, sigmay, my_r)
        x_rdata = self._x_rdata(x, tcl)
        return llh, (x_rdata, obs)

    @partial(jax.jit, static_argnames=('self', 'ts', 'k', 'my', 'pscale'))
    def run(self,
            ts: tuple,
            p: jnp.ndarray,
            k: tuple,
            my: tuple,
            pscale: tuple):
        return self._run(ts, p, k, my, pscale)

    @partial(jax.jit, static_argnames=('self', 'ts', 'k', 'my', 'pscale'))
    def srun(self,
             ts: tuple,
             p: jnp.ndarray,
             k: tuple,
             my: tuple,
             pscale: tuple):
        (llh, (x, obs)), sllh = (jax.value_and_grad(self._run, 1, True))(
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
        (llh, (x, obs)), sllh = (jax.value_and_grad(self._run, 1, True))(
            ts, p, k, my, pscale
        )
        s2llh, (x, obs) = jax.jacfwd(jax.grad(self._run, 1, True), 1, True)(
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

    for field in rdata_kwargs.keys():
        if field == 'llh':
            rdata_kwargs[field] = np.float(rdata_kwargs[field])
        elif field not in ['sllh', 's2llh']:
            rdata_kwargs[field] = np.asarray(rdata_kwargs[field]).T
            if rdata_kwargs[field].ndim == 1:
                rdata_kwargs[field] = np.expand_dims(rdata_kwargs[field], 1)

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
