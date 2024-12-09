"""PEtab wrappers for JAX models.""" ""
import shutil
from numbers import Number
from collections.abc import Iterable
from pathlib import Path


import diffrax
import equinox as eqx
import jaxtyping as jt
import jax.lax
import jax.numpy as jnp
import numpy as np
import pandas as pd
import petab.v1 as petab

from amici import _module_from_path
from amici.petab.parameter_mapping import (
    ParameterMappingForCondition,
    create_parameter_mapping,
)
from amici.jax.model import JAXModel, ReturnValue

DEFAULT_CONTROLLER_SETTINGS = {
    "atol": 1e-8,
    "rtol": 1e-8,
    "pcoeff": 0.4,
    "icoeff": 0.3,
    "dcoeff": 0.0,
}

SCALE_TO_INT = {
    petab.LIN: 0,
    petab.LOG: 1,
    petab.LOG10: 2,
}


def jax_unscale(
    parameter: jnp.float_,
    scale_str: str,
) -> jnp.float_:
    """Unscale parameter according to ``scale_str``.

    Arguments:
        parameter:
            Parameter to be unscaled.
        scale_str:
            One of ``petab.LIN``, ``petab.LOG``, ``petab.LOG10``.

    Returns:
        The unscaled parameter.
    """
    if scale_str == petab.LIN or not scale_str:
        return parameter
    if scale_str == petab.LOG:
        return jnp.exp(parameter)
    if scale_str == petab.LOG10:
        return jnp.power(10, parameter)
    raise ValueError(f"Invalid parameter scaling: {scale_str}")


class JAXProblem(eqx.Module):
    """
    PEtab problem wrapper for JAX models.

    :ivar parameters:
        Values for the model parameters. Do not change dimensions, values may be changed during, e.g. model training.
    :ivar model:
        JAXModel instance to use for simulation.
    :ivar _parameter_mappings:
        :class:`ParameterMappingForCondition` instances for each simulation condition.
    :ivar _measurements:
        Subset measurement dataframes for each simulation condition.
    :ivar _petab_problem:
        PEtab problem to simulate.
    """

    parameters: jnp.ndarray
    model: JAXModel
    _parameter_mappings: dict[str, ParameterMappingForCondition]
    _measurements: dict[
        tuple[str, ...],
        tuple[
            np.ndarray,
            np.ndarray,
            np.ndarray,
            np.ndarray,
            np.ndarray,
            np.ndarray,
        ],
    ]
    _petab_measurement_indices: dict[tuple[str, ...], tuple[int, ...]]
    _petab_problem: petab.Problem

    def __init__(self, model: JAXModel, petab_problem: petab.Problem):
        """
        Initialize a JAXProblem instance with a model and a PEtab problem.

        :param model:
            JAXModel instance to use for simulation.
        :param petab_problem:
            PEtab problem to simulate.
        """
        self.model = model
        scs = petab_problem.get_simulation_conditions_from_measurement_df()
        self._petab_problem = petab_problem
        self._parameter_mappings = self._get_parameter_mappings(scs)
        self._measurements, self._petab_measurement_indices = (
            self._get_measurements(scs)
        )
        self.parameters = self._get_nominal_parameter_values()

    def save(self, directory: Path):
        """
        Save the problem to a directory.

        :param directory:
            Directory to save the problem to.
        """
        self._petab_problem.to_files(
            prefix_path=directory,
            model_file="model",
            condition_file="conditions.tsv",
            measurement_file="measurements.tsv",
            parameter_file="parameters.tsv",
            observable_file="observables.tsv",
            yaml_file="problem.yaml",
        )
        shutil.copy(self.model.jax_py_file, directory / "jax_py_file.py")
        with open(directory / "parameters.pkl", "wb") as f:
            eqx.tree_serialise_leaves(f, self)

    @classmethod
    def load(cls, directory: Path):
        """
        Load a problem from a directory.

        :param directory:
            Directory to load the problem from.

        :return:
            Loaded problem instance.
        """
        petab_problem = petab.Problem.from_yaml(
            directory / "problem.yaml",
        )
        model = _module_from_path("jax", directory / "jax_py_file.py").Model()
        problem = cls(model, petab_problem)
        with open(directory / "parameters.pkl", "rb") as f:
            return eqx.tree_deserialise_leaves(f, problem)

    def _get_parameter_mappings(
        self, simulation_conditions: pd.DataFrame
    ) -> dict[str, ParameterMappingForCondition]:
        """
        Create parameter mappings for the provided simulation conditions.

        :param simulation_conditions:
            Simulation conditions to create parameter mappings for. Same format as returned by
            :meth:`petab.Problem.get_simulation_conditions_from_measurement_df`.
        :return:
            Dictionary mapping simulation conditions to parameter mappings.
        """
        scs = list(set(simulation_conditions.values.flatten()))
        mappings = create_parameter_mapping(
            petab_problem=self._petab_problem,
            simulation_conditions=[
                {petab.SIMULATION_CONDITION_ID: sc} for sc in scs
            ],
            scaled_parameters=False,
        )
        for mapping in mappings:
            for sim_var, value in mapping.map_sim_var.items():
                if isinstance(value, Number) and not np.isfinite(value):
                    mapping.map_sim_var[sim_var] = 1.0
        return dict(zip(scs, mappings, strict=True))

    def _get_measurements(
        self, simulation_conditions: pd.DataFrame
    ) -> tuple[
        dict[
            tuple[str, ...],
            tuple[
                np.ndarray,
                np.ndarray,
                np.ndarray,
                np.ndarray,
                np.ndarray,
                np.ndarray,
            ],
        ],
        dict[tuple[str, ...], tuple[int, ...]],
    ]:
        """
        Get measurements for the model based on the provided simulation conditions.

        :param simulation_conditions:
            Simulation conditions to create parameter mappings for. Same format as returned by
            :meth:`petab.Problem.get_simulation_conditions_from_measurement_df`.
        :return:
            Dictionary mapping simulation conditions to measurements (tuple of pre-equilibrium, dynamic,
            post-equilibrium time points; measurements and observable indices).
        """
        measurements = dict()
        indices = dict()
        for _, simulation_condition in simulation_conditions.iterrows():
            query = " & ".join(
                [f"{k} == '{v}'" for k, v in simulation_condition.items()]
            )
            m = self._petab_problem.measurement_df.query(query).sort_values(
                by=petab.TIME
            )

            ts = m[petab.TIME]
            ts_preeq = ts[np.isfinite(ts) & (ts == 0)]
            ts_dyn = ts[np.isfinite(ts) & (ts > 0)]
            ts_posteq = ts[np.logical_not(np.isfinite(ts))]
            index = pd.concat([ts_preeq, ts_dyn, ts_posteq]).index
            ts_preeq = ts_preeq.values
            ts_dyn = ts_dyn.values
            ts_posteq = ts_posteq.values
            my = m[petab.MEASUREMENT].values
            iys = np.array(
                [
                    self.model.observable_ids.index(oid)
                    for oid in m[petab.OBSERVABLE_ID].values
                ]
            )
            if (
                petab.OBSERVABLE_TRANSFORMATION
                in self._petab_problem.observable_df
            ):
                iy_trafos = np.array(
                    [
                        SCALE_TO_INT[
                            self._petab_problem.observable_df.loc[
                                oid, petab.OBSERVABLE_TRANSFORMATION
                            ]
                        ]
                        for oid in m[petab.OBSERVABLE_ID].values
                    ]
                )
            else:
                iy_trafos = np.zeros_like(iys)

            measurements[tuple(simulation_condition)] = (
                ts_preeq,
                ts_dyn,
                ts_posteq,
                my,
                iys,
                iy_trafos,
            )
            indices[tuple(simulation_condition)] = tuple(index.tolist())
        return measurements, indices

    def get_all_simulation_conditions(self) -> tuple[tuple[str, ...], ...]:
        simulation_conditions = (
            self._petab_problem.get_simulation_conditions_from_measurement_df()
        )
        return tuple(tuple(row) for _, row in simulation_conditions.iterrows())

    def _get_nominal_parameter_values(self) -> jt.Float[jt.Array, "np"]:
        """
        Get the nominal parameter values for the model based on the nominal values in the PEtab problem.

        :return:
            jax array with nominal parameter values
        """
        return jnp.array(
            [
                petab.scale(
                    self._petab_problem.parameter_df.loc[
                        pval, petab.NOMINAL_VALUE
                    ],
                    self._petab_problem.parameter_df.loc[
                        pval, petab.PARAMETER_SCALE
                    ],
                )
                for pval in self.parameter_ids
            ]
        )

    @property
    def parameter_ids(self) -> list[str]:
        """
        Parameter ids that are estimated in the PEtab problem. Same ordering as values in :attr:`parameters`.

        :return:
            PEtab parameter ids
        """
        return self._petab_problem.parameter_df[
            self._petab_problem.parameter_df[petab.ESTIMATE] == 1
        ].index.tolist()

    def get_petab_parameter_by_id(self, name: str) -> jnp.float_:
        """
        Get the value of a PEtab parameter by name.

        :param name:
            PEtab parameter id, as returned by :attr:`parameter_ids`.
        :return:
            Value of the parameter
        """
        return self.parameters[self.parameter_ids.index(name)]

    def _unscale(
        self, p: jt.Float[jt.Array, "np"], scales: tuple[str, ...]
    ) -> jt.Float[jt.Array, "np"]:
        """
        Unscaling of parameters.

        :param p:
            Parameter values
        :param scales:
            Parameter scalings
        :return:
            Unscaled parameter values
        """
        return jnp.array(
            [jax_unscale(pval, scale) for pval, scale in zip(p, scales)]
        )

    def load_parameters(
        self, simulation_condition: str
    ) -> jt.Float[jt.Array, "np"]:
        """
        Load parameters for a simulation condition.

        :param simulation_condition:
            Simulation condition to load parameters for.
        :return:
            Parameters for the simulation condition.
        """
        mapping = self._parameter_mappings[simulation_condition]
        p = jnp.array(
            [
                pval
                if isinstance(pval := mapping.map_sim_var[pname], Number)
                else self.get_petab_parameter_by_id(pval)
                for pname in self.model.parameter_ids
            ]
        )
        pscale = tuple(
            [
                mapping.scale_map_sim_var[pname]
                for pname in self.model.parameter_ids
            ]
        )
        return self._unscale(p, pscale)

    def _state_needs_reinitialisation(
        self,
        simulation_condition: str,
        state_id: str,
    ) -> bool:
        """
        Check if a state needs reinitialisation for a simulation condition.

        :param simulation_condition:
            simulation condition to check reinitialisation for
        :param state_id:
            state id to check reinitialisation for
        :return:
            True if state needs reinitialisation, False otherwise
        """
        if state_id not in self._petab_problem.condition_df:
            return False
        xval = self._petab_problem.condition_df.loc[
            simulation_condition, state_id
        ]
        if isinstance(xval, Number) and np.isnan(xval):
            return False
        return True

    def _state_reinitialisation_value(
        self,
        simulation_condition: str,
        state_id: str,
        p: jt.Float[jt.Array, "np"],
    ) -> jt.Float[jt.Scalar, ""] | float:  # noqa: F722
        """
        Get the reinitialisation value for a state.

        :param simulation_condition:
            simulation condition to get reinitialisation value for
        :param state_id:
            state id to get reinitialisation value for
        :param p:
            parameters for the simulation condition
        :return:
            reinitialisation value for the state
        """
        if state_id not in self._petab_problem.condition_df:
            # no reinitialisation, return dummy value
            return 0.0
        xval = self._petab_problem.condition_df.loc[
            simulation_condition, state_id
        ]
        if isinstance(xval, Number) and np.isnan(xval):
            # no reinitialisation, return dummy value
            return 0.0
        if isinstance(xval, Number):
            # numerical value, return as is
            return xval
        if xval in self.model.parameter_ids:
            # model parameter, return value
            return p[self.model.parameter_ids.index(xval)]
        if xval in self.parameter_ids:
            # estimated PEtab parameter, return unscaled value
            return jax_unscale(
                self.get_petab_parameter_by_id(xval),
                self._petab_problem.parameter_df.loc[
                    xval, petab.PARAMETER_SCALE
                ],
            )
        # only remaining option is nominal value for PEtab parameter
        # that is not estimated, return nominal value
        return self._petab_problem.parameter_df.loc[xval, petab.NOMINAL_VALUE]

    def load_reinitialisation(
        self,
        simulation_condition: str,
        p: jt.Float[jt.Array, "np"],
    ) -> tuple[jt.Bool[jt.Array, "nx"], jt.Float[jt.Array, "nx"]]:  # noqa: F821
        """
        Load reinitialisation values and mask for the state vector for a simulation condition.

        :param simulation_condition:
            Simulation condition to load reinitialisation for.
        :param p:
            Parameters for the simulation condition.
        :return:
            Tuple of reinitialisation masm and value for states.
        """
        if not any(
            x_id in self._petab_problem.condition_df
            for x_id in self.model.state_ids
        ):
            return jnp.array([]), jnp.array([])

        mask = jnp.array(
            [
                self._state_needs_reinitialisation(simulation_condition, x_id)
                for x_id in self.model.state_ids
            ]
        )
        reinit_x = jnp.array(
            [
                self._state_reinitialisation_value(
                    simulation_condition, x_id, p
                )
                for x_id in self.model.state_ids
            ]
        )
        return mask, reinit_x

    def update_parameters(self, p: jt.Float[jt.Array, "np"]) -> "JAXProblem":
        """
        Update parameters for the model.

        :param p:
            New problem instance with updated parameters.
        """
        return eqx.tree_at(lambda p: p.parameters, self, p)

    def run_simulation(
        self,
        simulation_condition: tuple[str, ...],
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        max_steps: jnp.int_,
        x_preeq: jt.Float[jt.Array, "*nx"] = jnp.array([]),  # noqa: F821, F722
        ret: ReturnValue = ReturnValue.llh,
    ) -> tuple[jnp.float_, dict]:
        """
        Run a simulation for a given simulation condition.

        :param simulation_condition:
            Simulation condition to run simulation for.
        :param solver:
            ODE solver to use for simulation
        :param controller:
            Step size controller to use for simulation
        :param max_steps:
            Maximum number of steps to take during simulation
        :param x_preeq:
            Pre-equilibration state if available
        :param ret:
            which output to return. See :class:`ReturnValue` for available options.
        :return:
            Tuple of output value and simulation statistics
        """
        ts_preeq, ts_dyn, ts_posteq, my, iys, iy_trafos = self._measurements[
            simulation_condition
        ]
        p = self.load_parameters(simulation_condition[0])
        mask_reinit, x_reinit = self.load_reinitialisation(
            simulation_condition[0], p
        )
        return self.model.simulate_condition(
            p=eqx.debug.backward_nan(p),
            ts_init=jax.lax.stop_gradient(jnp.array(ts_preeq)),
            ts_dyn=jax.lax.stop_gradient(jnp.array(ts_dyn)),
            ts_posteq=jax.lax.stop_gradient(jnp.array(ts_posteq)),
            my=jax.lax.stop_gradient(jnp.array(my)),
            iys=jax.lax.stop_gradient(jnp.array(iys)),
            iy_trafos=jax.lax.stop_gradient(jnp.array(iy_trafos)),
            x_preeq=x_preeq,
            mask_reinit=jax.lax.stop_gradient(mask_reinit),
            x_reinit=x_reinit,
            solver=solver,
            controller=controller,
            max_steps=max_steps,
            adjoint=diffrax.RecursiveCheckpointAdjoint()
            if ret in (ReturnValue.llh, ReturnValue.chi2)
            else diffrax.DirectAdjoint(),
            ret=ret,
        )

    def run_preequilibration(
        self,
        simulation_condition: str,
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        max_steps: jnp.int_,
    ) -> tuple[jt.Float[jt.Array, "nx"], dict]:  # noqa: F821
        """
        Run a pre-equilibration simulation for a given simulation condition.

        :param simulation_condition:
            Simulation condition to run simulation for.
        :param solver:
            ODE solver to use for simulation
        :param controller:
            Step size controller to use for simulation
        :param max_steps:
            Maximum number of steps to take during simulation
        :return:
            Pre-equilibration state
        """
        p = self.load_parameters(simulation_condition)
        mask_reinit, x_reinit = self.load_reinitialisation(
            simulation_condition, p
        )
        return self.model.preequilibrate_condition(
            p=eqx.debug.backward_nan(p),
            mask_reinit=mask_reinit,
            x_reinit=x_reinit,
            solver=solver,
            controller=controller,
            max_steps=max_steps,
        )


def run_simulations(
    problem: JAXProblem,
    simulation_conditions: Iterable[tuple[str, ...]] | None = None,
    solver: diffrax.AbstractSolver = diffrax.Kvaerno5(),
    controller: diffrax.AbstractStepSizeController = diffrax.PIDController(
        **DEFAULT_CONTROLLER_SETTINGS
    ),
    max_steps: int = 2**10,
    ret: ReturnValue | str = ReturnValue.llh,
):
    """
    Run simulations for a problem.

    :param problem:
        Problem to run simulations for.
    :param simulation_conditions:
        Simulation conditions to run simulations for.
    :param solver:
        ODE solver to use for simulation.
    :param controller:
        Step size controller to use for simulation.
    :param max_steps:
        Maximum number of steps to take during simulation.
    :param ret:
        which output to return. See :class:`ReturnValue` for available options.
    :return:
        Overall output value and condition specific results and statistics.
    """
    if isinstance(ret, str):
        ret = ReturnValue[ret]

    if simulation_conditions is None:
        simulation_conditions = problem.get_all_simulation_conditions()

    preeqs = {
        sc: problem.run_preequilibration(sc, solver, controller, max_steps)
        # only run preequilibration once per condition
        for sc in {sc[1] for sc in simulation_conditions if len(sc) > 1}
    }

    results = {
        sc: problem.run_simulation(
            sc,
            solver,
            controller,
            max_steps,
            preeqs.get(sc[1])[0] if len(sc) > 1 else jnp.array([]),
            ret=ret,
        )
        for sc in simulation_conditions
    }
    stats = {
        sc: res[1] | preeqs[sc[1]][1] if len(sc) > 1 else res[1]
        for sc, res in results.items()
    }
    if ret in (ReturnValue.llh, ReturnValue.chi2):
        output = sum(r for r, _ in results.values())
    else:
        output = {sc: res[0] for sc, res in results.items()}

    return output, stats


def petab_simulate(
    problem: JAXProblem,
    solver: diffrax.AbstractSolver = diffrax.Kvaerno5(),
    controller: diffrax.AbstractStepSizeController = diffrax.PIDController(
        **DEFAULT_CONTROLLER_SETTINGS
    ),
    max_steps: int = 2**10,
):
    """
    Run simulations for a problem and return the results as a petab simulation dataframe.

    :param problem:
        Problem to run simulations for.
    :param solver:
        ODE solver to use for simulation.
    :param controller:
        Step size controller to use for simulation.
    :param max_steps:
        Maximum number of steps to take during simulation.
    :return:
        petab simulation dataframe.
    """
    y, r = run_simulations(
        problem,
        solver=solver,
        controller=controller,
        max_steps=max_steps,
        ret=ReturnValue.y,
    )
    dfs = []
    for sc, ys in y.items():
        obs = [
            problem.model.observable_ids[io]
            for io in problem._measurements[sc][4]
        ]
        t = jnp.concat(problem._measurements[sc][:2])
        df_sc = pd.DataFrame(
            {
                petab.SIMULATION: ys,
                petab.TIME: t,
                petab.OBSERVABLE_ID: obs,
                petab.SIMULATION_CONDITION_ID: [sc[0]] * len(t),
            },
            index=problem._petab_measurement_indices[sc],
        )
        if (
            petab.OBSERVABLE_PARAMETERS
            in problem._petab_problem.measurement_df
        ):
            df_sc[petab.OBSERVABLE_PARAMETERS] = (
                problem._petab_problem.measurement_df.query(
                    f"{petab.SIMULATION_CONDITION_ID} == '{sc[0]}'"
                )[petab.OBSERVABLE_PARAMETERS]
            )
        if petab.NOISE_PARAMETERS in problem._petab_problem.measurement_df:
            df_sc[petab.NOISE_PARAMETERS] = (
                problem._petab_problem.measurement_df.query(
                    f"{petab.SIMULATION_CONDITION_ID} == '{sc[0]}'"
                )[petab.NOISE_PARAMETERS]
            )
        if (
            petab.PREEQUILIBRATION_CONDITION_ID
            in problem._petab_problem.measurement_df
        ):
            df_sc[petab.PREEQUILIBRATION_CONDITION_ID] = sc[1]
        dfs.append(df_sc)
    return pd.concat(dfs).sort_index()
