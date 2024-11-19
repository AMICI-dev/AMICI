"""PEtab wrappers for JAX models.""" ""

from numbers import Number
from collections.abc import Iterable

import diffrax
import equinox as eqx
import jaxtyping as jt
import jax.lax
import jax.numpy as jnp
import numpy as np
import pandas as pd
import petab.v1 as petab

from amici.petab.parameter_mapping import (
    ParameterMappingForCondition,
    create_parameter_mapping,
)
from amici.jax.model import JAXModel


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
        tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    ]
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
        self._measurements = self._get_measurements(scs)
        self.parameters = self._get_nominal_parameter_values()

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
    ) -> dict[
        tuple[str],
        tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray],
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
        for _, simulation_condition in simulation_conditions.iterrows():
            query = " & ".join(
                [f"{k} == '{v}'" for k, v in simulation_condition.items()]
            )
            m = self._petab_problem.measurement_df.query(query).sort_values(
                by=petab.TIME
            )

            ts = m[petab.TIME].values
            ts_preeq = ts[np.isfinite(ts) & (ts == 0)]
            ts_dyn = ts[np.isfinite(ts) & (ts > 0)]
            ts_posteq = ts[np.logical_not(np.isfinite(ts))]
            my = m[petab.MEASUREMENT].values
            iys = np.array(
                [
                    self.model.observable_ids.index(oid)
                    for oid in m[petab.OBSERVABLE_ID].values
                ]
            )

            measurements[tuple(simulation_condition)] = (
                ts_preeq,
                ts_dyn,
                ts_posteq,
                my,
                iys,
            )
        return measurements

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
    ) -> tuple[jnp.float_, dict]:
        """
        Run a simulation for a given simulation condition.

        :param simulation_condition:
            Tuple of simulation conditions to run the simulation for. can be a single string (simulation only) or a
            tuple of strings (pre-equilibration followed by simulation).
        :param solver:
            ODE solver to use for simulation
        :param controller:
            Step size controller to use for simulation
        :param max_steps:
            Maximum number of steps to take during simulation
        :return:
            Tuple of log-likelihood and simulation statistics
        """
        ts_preeq, ts_dyn, ts_posteq, my, iys = self._measurements[
            simulation_condition
        ]
        p = self.load_parameters(simulation_condition[0])
        p_preeq = (
            self.load_parameters(simulation_condition[1])
            if len(simulation_condition) > 1
            else jnp.array([])
        )
        return self.model.simulate_condition(
            p=p,
            p_preeq=p_preeq,
            ts_preeq=jax.lax.stop_gradient(jnp.array(ts_preeq)),
            ts_dyn=jax.lax.stop_gradient(jnp.array(ts_dyn)),
            ts_posteq=jax.lax.stop_gradient(jnp.array(ts_posteq)),
            my=jax.lax.stop_gradient(jnp.array(my)),
            iys=jax.lax.stop_gradient(jnp.array(iys)),
            solver=solver,
            controller=controller,
            max_steps=max_steps,
            adjoint=diffrax.RecursiveCheckpointAdjoint(),
        )


def run_simulations(
    problem: JAXProblem,
    simulation_conditions: Iterable[tuple] | None = None,
    solver: diffrax.AbstractSolver = diffrax.Kvaerno5(),
    controller: diffrax.AbstractStepSizeController = diffrax.PIDController(
        rtol=1e-8,
        atol=1e-8,
        pcoeff=0.4,
        icoeff=0.3,
        dcoeff=0.0,
    ),
    max_steps: int = 2**10,
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
    :return:
        Overall negative log-likelihood and condition specific results and statistics.
    """
    if simulation_conditions is None:
        simulation_conditions = problem.get_all_simulation_conditions()

    results = {
        sc: problem.run_simulation(sc, solver, controller, max_steps)
        for sc in simulation_conditions
    }
    return sum(llh for llh, _ in results.values()), results
