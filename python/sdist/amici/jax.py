from abc import abstractmethod
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor
from numbers import Number

import diffrax
import equinox as eqx
import jax.numpy as jnp
import numpy as np
import jax
import pandas as pd
import petab.v1 as petab

import amici
from amici.petab.parameter_mapping import (
    ParameterMappingForCondition,
    create_parameter_mapping,
)
from amici.petab.conditions import (
    _get_timepoints_with_replicates,
    _get_measurements_and_sigmas,
)

jax.config.update("jax_enable_x64", True)


class JAXModel(eqx.Module):
    solver: diffrax.AbstractSolver
    controller: diffrax.AbstractStepSizeController
    maxsteps: int
    parameters: jnp.ndarray
    parameter_mappings: dict[tuple[str], ParameterMappingForCondition] | None
    measurements: dict[tuple[str], pd.DataFrame] | None
    petab_problem: petab.Problem | None

    def __init__(self):
        self.solver = diffrax.Kvaerno5()
        self.maxsteps: int = 2**14
        self.controller = diffrax.PIDController(
            rtol=1e-8,
            atol=1e-8,
            pcoeff=0.4,
            icoeff=0.3,
            dcoeff=0.0,
        )
        self.petab_problem = None
        self.parameter_mappings = None
        self.measurements = None
        self.parameters = jnp.array([])

    def _set_parameter_mappings(
        self, simulation_conditions: pd.DataFrame
    ) -> "JAXModel":
        mappings = create_parameter_mapping(
            petab_problem=self.petab_problem,
            simulation_conditions=simulation_conditions,
            scaled_parameters=False,
            amici_model=self,
        )

        parameter_mappings = {
            tuple(simulation_condition.values): mapping
            for (_, simulation_condition), mapping in zip(
                simulation_conditions.iterrows(), mappings
            )
        }

        is_leaf = (  # noqa: E731
            lambda x: x is None if self.parameter_mappings is None else None
        )
        return eqx.tree_at(
            lambda x: x.parameter_mappings,
            self,
            parameter_mappings,
            is_leaf=is_leaf,
        )

    def _set_measurements(
        self, simulation_conditions: pd.DataFrame
    ) -> "JAXModel":
        measurements = dict()
        for _, simulation_condition in simulation_conditions.iterrows():
            measurements_df = self.petab_problem.measurement_df
            for k, v in simulation_condition.items():
                measurements_df = measurements_df.query(f"{k} == '{v}'")

            ts = _get_timepoints_with_replicates(measurements_df)
            my = _get_measurements_and_sigmas(
                measurements_df, ts, self.observable_ids
            )[0].flatten()
            measurements[tuple(simulation_condition)] = np.array(ts), my
        is_leaf = (  # noqa: E731
            lambda x: x is None if self.measurements is None else None
        )
        return eqx.tree_at(
            lambda x: x.measurements,
            self,
            measurements,
            is_leaf=is_leaf,
        )

    def _set_nominal_parameter_values(self) -> "JAXModel":
        nominal_values = jnp.array(
            [
                petab.scale(
                    self.petab_problem.parameter_df.loc[
                        pval, petab.NOMINAL_VALUE
                    ],
                    self.petab_problem.parameter_df.loc[
                        pval, petab.PARAMETER_SCALE
                    ],
                )
                for pval in self.petab_parameter_ids()
            ]
        )
        return eqx.tree_at(lambda x: x.parameters, self, nominal_values)

    def _set_petab_problem(self, petab_problem: petab.Problem) -> "JAXModel":
        is_leaf = lambda x: x is None if self.petab_problem is None else None  # noqa: E731
        return eqx.tree_at(
            lambda x: x.petab_problem,
            self,
            petab_problem,
            is_leaf=is_leaf,
        )

    def set_petab_problem(self, petab_problem: petab.Problem) -> "JAXModel":
        """
        Set the PEtab problem for the model and updates parameters to the nominal values.
        :param petab_problem:
            Petab problem to set.
        :return: JAXModel instance
        """

        model = self._set_petab_problem(petab_problem)
        simulation_conditions = (
            petab_problem.get_simulation_conditions_from_measurement_df()
        )
        model = model._set_parameter_mappings(simulation_conditions)
        model = model._set_measurements(simulation_conditions)
        return model._set_nominal_parameter_values()

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

    def getParameterIds(self) -> list[str]:  # noqa: N802
        """
        Get the parameter ids of the model. Adds compatibility with AmiciModel, added to enable generation of
        parameter mappings via :func:`amici.petab.create_parameter_mapping`.
        :return:
        """
        return self.parameter_ids

    def getFixedParameterIds(self) -> list[str]:  # noqa: N802
        """
        Get the fixed parameter ids of the model. Adds compatibility with AmiciModel, added to enable generation of
        parameter mappings via :func:`amici.petab.create_parameter_mapping`.
        :return:
        """
        return self.fixed_parameter_ids

    def petab_parameter_ids(self) -> list[str]:
        return self.petab_problem.parameter_df[
            self.petab_problem.parameter_df[petab.ESTIMATE] == 1
        ].index.tolist()

    def get_petab_parameter_by_name(self, name: str) -> jnp.float_:
        return self.parameters[self.petab_parameter_ids().index(name)]

    def _unscale_p(self, p, pscale):
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
            diffrax.ODETerm(self.xdot),
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
            diffrax.ODETerm(self.xdot),
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

    def run_condition(
        self,
        ts: jnp.ndarray,
        ts_dyn: jnp.ndarray,
        p: jnp.ndarray,
        k: jnp.ndarray,
        k_preeq: jnp.ndarray,
        my: jnp.ndarray,
        pscale: jnp.ndarray,
        checkpointed=True,
        dynamic="true",
    ):
        ps = self._unscale_p(p, pscale)

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

        obs = jnp.stack(self._obs(ts, x, ps, k, tcl), axis=1)
        my_r = my.reshape((len(ts), -1))
        sigmay = self._sigmay(obs, ps, k)
        llh = self._loss(obs, sigmay, my_r)
        x_rdata = jnp.stack(self._x_rdata(x, tcl), axis=1)
        return llh, (x_rdata, obs, stats)

    @eqx.filter_jit
    def _fun(
        self,
        ts: jnp.ndarray,
        ts_dyn: jnp.ndarray,
        p: jnp.ndarray,
        k: jnp.ndarray,
        k_preeq: jnp.ndarray,
        my: jnp.ndarray,
        pscale: jnp.ndarray,
        dynamic="true",
    ):
        return self.run_condition(
            ts, ts_dyn, p, k, k_preeq, my, pscale, dynamic=dynamic
        )

    @eqx.filter_jit
    def _grad(
        self,
        ts: jnp.ndarray,
        ts_dyn: jnp.ndarray,
        p: jnp.ndarray,
        k: jnp.ndarray,
        k_preeq: jnp.ndarray,
        my: jnp.ndarray,
        pscale: jnp.ndarray,
        dynamic="true",
    ):
        (llh, (x, obs, stats)), sllh = (
            jax.value_and_grad(self.run_condition, 2, True)
        )(ts, ts_dyn, p, k, k_preeq, my, pscale, dynamic=dynamic)
        return llh, sllh, (x, obs, stats)

    @eqx.filter_jit
    def _hessian(
        self,
        ts: jnp.ndarray,
        ts_dyn: jnp.ndarray,
        p: jnp.ndarray,
        k: jnp.ndarray,
        k_preeq: jnp.ndarray,
        my: jnp.ndarray,
        pscale: jnp.ndarray,
        dynamic="true",
    ):
        (llh, (x, obs, stats)), sllh = (
            jax.value_and_grad(self.run_condition, 2, True)
        )(ts, ts_dyn, p, k, k_preeq, my, pscale, dynamic=dynamic)

        s2llh = jax.hessian(self.run_condition, 2, True)(
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
        simulation_condition: tuple[str],
        sensitivity_order: amici.SensitivityOrder = amici.SensitivityOrder.none,
    ):
        parameter_mapping = self.parameter_mappings[simulation_condition]
        ts, my = self.measurements[simulation_condition]
        p = jnp.array(
            [
                pval
                if isinstance(
                    pval := parameter_mapping.map_sim_var[par], Number
                )
                else self.get_petab_parameter_by_name(pval)
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

        ts_dyn = ts[np.isfinite(ts)]
        dynamic = "true" if len(ts_dyn) and np.max(ts_dyn) > 0 else "false"

        rdata_kwargs = dict(
            simulation_condition=simulation_condition,
        )

        if sensitivity_order == amici.SensitivityOrder.none:
            (
                rdata_kwargs["llh"],
                (rdata_kwargs["x"], rdata_kwargs["y"], rdata_kwargs["stats"]),
            ) = self._fun(
                ts, ts_dyn, p, k_sim, k_preeq, my, pscale, dynamic=dynamic
            )
        elif sensitivity_order == amici.SensitivityOrder.first:
            (
                rdata_kwargs["llh"],
                rdata_kwargs["sllh"],
                (rdata_kwargs["x"], rdata_kwargs["y"], rdata_kwargs["stats"]),
            ) = self._grad(
                ts, ts_dyn, p, k_sim, k_preeq, my, pscale, dynamic=dynamic
            )
        elif sensitivity_order == amici.SensitivityOrder.second:
            (
                rdata_kwargs["llh"],
                rdata_kwargs["sllh"],
                rdata_kwargs["s2llh"],
                (rdata_kwargs["x"], rdata_kwargs["y"], rdata_kwargs["stats"]),
            ) = self._hessian(
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
        simulation_conditions: tuple[tuple[str]] = None,
    ):
        fun = eqx.Partial(
            self.run_simulation,
            sensitivity_order=sensitivity_order,
        )

        if num_threads > 1:
            with ThreadPoolExecutor(max_workers=num_threads) as pool:
                results = pool.map(fun, simulation_conditions)
        else:
            results = map(fun, simulation_conditions)
        return list(results)


@dataclass
class ReturnDataJAX(dict):
    simulation_condition: tuple[str] = None
    x: np.array = None
    y: np.array = None
    sigmay: np.array = None
    llh: np.array = None
    sllh: np.array = None
    s2llh: np.array = None
    stats: dict = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__dict__ = self
