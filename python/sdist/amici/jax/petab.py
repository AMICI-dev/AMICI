"""
JAX
----
This module provides functions and classes to enable the use of JAX-based ODE solvers (currently diffrax) to simulate
 AMICI processed models. The API in this module is experimental. Expect substantial changes and do not use in production
 code.

Loading this module will automatically enable 64-bit precision for JAX.
"""

from numbers import Number
from collections.abc import Iterable

import diffrax
import equinox as eqx
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
            One of ``'lin'`` (synonymous with ``''``), ``'log'``, ``'log10'``.

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
    :ivar solver:
        Diffrax solver to use for model simulation
    :ivar controller:
        Step-size controller to use for model simulation
    :ivar max_steps:
        Maximum number of steps to take during a simulation
    :ivar parameters:
        Values for the model parameters. Only populated after setting the PEtab problem via :meth:`set_petab_problem`.
        Do not change dimensions, values may be changed during, e.g. model training.
    :ivar parameter_mappings:
        :class:`ParameterMappingForCondition` instances for each simulation condition. Only populated after setting the
        PEtab problem via :meth:`set_petab_problem`. Do not set manually unless you know what you are doing.
    :ivar measurements:
        Subset measurement dataframes for each simulation condition. Only populated after setting the PEtab problem
        via :meth:`set_petab_problem`. Do not set manually unless you know what you are doing.
    :ivar petab_problem:
        PEtab problem to simulate. Set via :meth:`set_petab_problem`.
    """

    parameters: jnp.ndarray
    model: JAXModel
    parameter_mappings: dict[tuple[str], ParameterMappingForCondition] = (
        eqx.field(static=True)
    )
    measurements: dict[
        tuple[str],
        tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, str],
    ] = eqx.field(static=True)
    petab_problem: petab.Problem

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
        self.petab_problem = petab_problem
        self.parameter_mappings = self._get_parameter_mappings(scs)
        self.measurements = self._get_measurements(scs)
        self.parameters = self._get_nominal_parameter_values()

    def _get_parameter_mappings(self, simulation_conditions: pd.DataFrame):
        scs = list(set(simulation_conditions.values.flatten()))
        mappings = create_parameter_mapping(
            petab_problem=self.petab_problem,
            simulation_conditions=[
                {petab.SIMULATION_CONDITION_ID: sc} for sc in scs
            ],
            scaled_parameters=False,
        )
        for mapping in mappings:
            for sim_var, value in mapping.map_sim_var.items():
                if isinstance(value, Number) and not np.isfinite(value):
                    mapping.map_sim_var[sim_var] = 1.0
        return dict(zip(scs, mappings))

    def _get_measurements(self, simulation_conditions: pd.DataFrame):
        """
        Set measurements for the model based on the provided simulation conditions.
        :param simulation_conditions:
            Simulation conditions to create parameter mappings for. Same format as returned by
            :meth:`petab.Problem.get_simulation_conditions_from_measurement_df`.
        :return:
            JAXModel instance with measurements set.
        """
        measurements = dict()
        for _, simulation_condition in simulation_conditions.iterrows():
            measurements_df = self.petab_problem.measurement_df
            for k, v in simulation_condition.items():
                measurements_df = measurements_df.query(f"{k} == '{v}'")

            measurements_df.sort_values(by=petab.TIME, inplace=True)

            ts = measurements_df[petab.TIME].values
            ts_dyn = [t for t in ts if np.isfinite(t)]
            my = measurements_df[petab.MEASUREMENT].values
            iys = np.array(
                [
                    self.model.observable_ids.index(oid)
                    for oid in measurements_df[petab.OBSERVABLE_ID].values
                ]
            )

            # using strings here prevents tracing in jax
            dynamic = ts_dyn and max(ts_dyn) > 0
            measurements[tuple(simulation_condition)] = (
                np.array(ts),
                np.array(ts_dyn),
                my,
                iys,
                dynamic,
            )
        return measurements

    def _get_nominal_parameter_values(self) -> jnp.ndarray:
        """
        Set the nominal parameter values for the model based on the nominal values in the PEtab problem.
        :return:
            JAXModel instance with parameter values set to the nominal values.
        """
        if self.petab_problem is None:
            raise ValueError(
                "PEtab problem not set, cannot set nominal values."
            )
        return jnp.array(
            [
                petab.scale(
                    self.petab_problem.parameter_df.loc[
                        pval, petab.NOMINAL_VALUE
                    ],
                    self.petab_problem.parameter_df.loc[
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
        return self.petab_problem.parameter_df[
            self.petab_problem.parameter_df[petab.ESTIMATE] == 1
        ].index.tolist()

    def get_petab_parameter_by_id(self, name: str) -> jnp.float_:
        """
        Get the value of a PEtab parameter by name.
        :param name:
            PEtab parameter id
        :return:
            Value of the parameter
        """
        return self.parameters[self.parameter_ids.index(name)]

    def _unscale_p(
        self, p: jnp.ndarray, pscale: tuple[str, ...]
    ) -> jnp.ndarray:
        """
        Unscaling of parameters.

        :param p:
            Parameter values
        :param pscale:
            Parameter scaling
        :return:
            Unscaled parameter values
        """
        return jnp.array(
            [jax_unscale(pval, scale) for pval, scale in zip(p, pscale)]
        )

    def load_parameters(self, simulation_condition) -> jnp.ndarray:
        mapping = self.parameter_mappings[simulation_condition]
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
        return self._unscale_p(p, pscale)

    def run_simulation(
        self,
        simulation_condition: tuple[str],
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        max_steps: int,
    ):
        ts, ts_dyn, my, iys, dynamic = self.measurements[simulation_condition]
        p = self.load_parameters(simulation_condition[0])
        p_preeq = (
            self.load_parameters(simulation_condition[1])
            if len(simulation_condition) > 1
            else jnp.array([])
        )
        return self.model.simulate_condition(
            ts,
            ts_dyn,
            my,
            iys,
            p,
            p_preeq,
            dynamic,
            solver,
            controller,
            max_steps,
        )


def run_simulations(
    problem: JAXProblem,
    simulation_conditions: Iterable[tuple] = None,
    solver: diffrax.AbstractSolver = diffrax.Kvaerno5(),
    controller: diffrax.AbstractStepSizeController = diffrax.PIDController(
        rtol=1e-8,
        atol=1e-8,
        pcoeff=0.4,
        icoeff=0.3,
        dcoeff=0.0,
    ),
    max_steps: int = 2**14,
):
    results = {
        sc: problem.run_simulation(sc, solver, controller, max_steps)
        for sc in simulation_conditions
    }
    return sum(llh for llh, _ in results.values()), results
