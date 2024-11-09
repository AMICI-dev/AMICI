from abc import abstractmethod
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor
from numbers import Number

import diffrax
import equinox as eqx
import jax.numpy as jnp
import numpy as np
import pandas as pd
import jax
import petab.v1 as petab

import amici
from amici.petab.parameter_mapping import (
    ParameterMapping,
    ParameterMappingForCondition,
)
from amici.petab.conditions import (
    _get_timepoints_with_replicates,
    _get_measurements_and_sigmas,
)

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

    @property
    @abstractmethod
    def state_ids(self): ...

    @property
    @abstractmethod
    def observable_ids(self): ...

    @property
    @abstractmethod
    def parameter_ids(self): ...

    @property
    @abstractmethod
    def fixed_parameter_ids(self): ...

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
        p: jnp.ndarray,
        k: np.ndarray,
        k_preeq: np.ndarray,
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
            ts,
            ts_dyn,
            p,
            k,
            k_preeq,
            my,
            pscale,
            checkpointed=False,
            dynamic=dynamic,
        )

        return llh, sllh, s2llh, (x, obs, stats)

    def run_simulation(
        self,
        parameter_mapping: ParameterMappingForCondition = None,
        measurements: pd.DataFrame = None,
        parameters: pd.DataFrame = None,
        sensitivity_order: amici.SensitivityOrder = amici.SensitivityOrder.none,
    ):
        cond_id, measurements_df = measurements
        ts = _get_timepoints_with_replicates(measurements_df)
        p = jnp.array(
            [
                pval
                if isinstance(
                    pval := parameter_mapping.map_sim_var[par], Number
                )
                else petab.scale(
                    parameters.loc[pval, petab.NOMINAL_VALUE],
                    parameters.loc[pval, petab.PARAMETER_SCALE],
                )
                for par in self.parameter_ids
            ]
        )
        pscale = jnp.array(
            [
                0 if s == petab.LIN else 1 if s == petab.LOG else 2
                for s in parameter_mapping.scale_map_sim_var.values()
            ]
        )
        k_sim = np.array(
            [
                parameter_mapping.map_sim_fix[k]
                for k in self.fixed_parameter_ids
            ]
        )
        k_preeq = np.array(
            [
                parameter_mapping.map_preeq_fix[k]
                for k in self.fixed_parameter_ids
                if k in parameter_mapping.map_preeq_fix
            ]
        )
        my = _get_measurements_and_sigmas(
            measurements_df, ts, self.observable_ids
        )[0].flatten()
        ts = np.array(ts)
        ts_dyn = ts[np.isfinite(ts)]
        dynamic = "true" if len(ts_dyn) and np.max(ts_dyn) > 0 else "false"

        rdata_kwargs = dict()

        if sensitivity_order == amici.SensitivityOrder.none:
            (
                rdata_kwargs["llh"],
                (rdata_kwargs["x"], rdata_kwargs["y"], rdata_kwargs["stats"]),
            ) = self.run(
                ts, ts_dyn, p, k_sim, k_preeq, my, pscale, dynamic=dynamic
            )
        elif sensitivity_order == amici.SensitivityOrder.first:
            (
                rdata_kwargs["llh"],
                rdata_kwargs["sllh"],
                (rdata_kwargs["x"], rdata_kwargs["y"], rdata_kwargs["stats"]),
            ) = self.srun(
                ts, ts_dyn, p, k_sim, k_preeq, my, pscale, dynamic=dynamic
            )
        elif sensitivity_order == amici.SensitivityOrder.second:
            (
                rdata_kwargs["llh"],
                rdata_kwargs["sllh"],
                rdata_kwargs["s2llh"],
                (rdata_kwargs["x"], rdata_kwargs["y"], rdata_kwargs["stats"]),
            ) = self.s2run(
                ts, ts_dyn, p, k_sim, k_preeq, my, pscale, dynamic=dynamic
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
        sensitivity_order: amici.SensitivityOrder = amici.SensitivityOrder.none,
        num_threads: int = 1,
        parameter_mappings: ParameterMapping = None,
        parameters: pd.DataFrame = None,
        simulation_conditions: pd.DataFrame = None,
        measurements: pd.DataFrame = None,
    ):
        fun = eqx.Partial(
            self.run_simulation,
            sensitivity_order=sensitivity_order,
            parameters=parameters,
        )
        gb = (
            [
                petab.PREEQUILIBRATION_CONDITION_ID,
                petab.SIMULATION_CONDITION_ID,
            ]
            if petab.PREEQUILIBRATION_CONDITION_ID in measurements.columns
            and petab.PREEQUILIBRATION_CONDITION_ID in simulation_conditions
            else petab.SIMULATION_CONDITION_ID
        )

        per_condition_measurements = measurements.groupby(gb)

        order_conditions = [
            tuple(c) if isinstance(c, np.ndarray) else c
            for c in simulation_conditions[gb].values
        ]

        sorted_mappings = [
            parameter_mappings[order_conditions.index(condition)]
            for condition in per_condition_measurements.groups.keys()
        ]

        if num_threads > 1:
            with ThreadPoolExecutor(max_workers=num_threads) as pool:
                results = pool.map(
                    fun, sorted_mappings, per_condition_measurements
                )
        else:
            results = map(fun, sorted_mappings, per_condition_measurements)
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
