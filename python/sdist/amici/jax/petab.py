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
    _inputs: dict[str, dict[str, np.ndarray]]
    _petab_problem: petab.Problem

    def __init__(self, model: JAXModel, petab_problem: petab.Problem):
        """
        Initialize a JAXProblem instance with a model and a PEtab problem.

        :param model:
            JAXModel instance to use for simulation.
        :param petab_problem:
            PEtab problem to simulate.
        """
        scs = petab_problem.get_simulation_conditions_from_measurement_df()
        self._petab_problem = petab_problem
        self.parameters, self.model = self._get_nominal_parameter_values(model)
        self._parameter_mappings = self._get_parameter_mappings(scs)
        self._measurements = self._get_measurements(scs)
        self._inputs = self._get_inputs()

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

    def _get_nominal_parameter_values(
        self, model: JAXModel
    ) -> tuple[jt.Float[jt.Array, "np"], JAXModel]:
        """
        Get the nominal parameter values for the model based on the nominal values in the PEtab problem.
        Also set nominal values in the model (where applicable).

        :return:
            jax array with nominal parameter values and model with nominal parameter values set.
        """
        # initialize everything with zeros
        model_pars = {
            net_id: {
                layer_id: {
                    attribute: jnp.zeros_like(getattr(layer, attribute))
                    for attribute in ["weight", "bias"]
                    if hasattr(layer, attribute)
                }
                for layer_id, layer in nn.layers.items()
            }
            for net_id, nn in model.nns.items()
        }
        # extract nominal values from petab problem
        for pname, row in self._petab_problem.parameter_df.iterrows():
            if (net := pname.split("_")[0]) in model.nns:
                nn = model_pars[net]
                layer = nn[pname.split("_")[1]]
                attribute = pname.split("_")[2]
                index = tuple(np.array(pname.split("_")[3:]).astype(int))
                layer[attribute] = (
                    layer[attribute].at[index].set(row[petab.NOMINAL_VALUE])
                )
        # set values in model
        for net_id in model_pars:
            for layer_id in model_pars[net_id]:
                for attribute in model_pars[net_id][layer_id]:
                    model = eqx.tree_at(
                        lambda model: getattr(
                            model.nns[net_id].layers[layer_id], attribute
                        ),
                        model,
                        model_pars[net_id][layer_id][attribute],
                    )
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
        ), model

    def _get_inputs(self):
        inputs = {
            net: {} for net in self._petab_problem.mapping_df["netId"].unique()
        }
        for petab_id, row in self._petab_problem.mapping_df.iterrows():
            if (filepath := Path(petab_id)).is_file():
                data_flat = pd.read_csv(filepath, sep="\t").sort_values(
                    by="ix"
                )
                shape = tuple(
                    np.stack(
                        data_flat["ix"]
                        .astype(str)
                        .str.split(";")
                        .apply(np.array)
                    )
                    .astype(int)
                    .max(axis=0)
                    + 1
                )
                inputs[row["netId"]][row[petab.MODEL_ENTITY_ID]] = data_flat[
                    "value"
                ].values.reshape(shape)
        return inputs

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

    def _eval_nn(self, output_par: str):
        net_id = self._petab_problem.mapping_df.loc[output_par, "netId"]
        nn = self.model.nns[net_id]
        net_input = tuple(
            jax.lax.stop_gradient(self._inputs[net_id][input_id])
            for input_id in nn.inputs
        )
        return nn.forward(*net_input).squeeze()

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
                self._eval_nn(pname)
                if pname in self._petab_problem.mapping_df.index
                else pval
                if isinstance(pval := mapping.map_sim_var[pname], Number)
                else self.get_petab_parameter_by_id(pval)
                for pname in self.model.parameter_ids
            ]
        )
        pscale = tuple(
            [
                petab.LIN
                if pname in self._petab_problem.mapping_df.index
                else mapping.scale_map_sim_var[pname]
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
        ret: str = "llh",
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
        :param ret:
            which output to return. Valid values are
                - `llh`: log-likelihood (default)
                - `nllhs`: negative log-likelihood at each time point
                - `x0`: full initial state vector (after pre-equilibration)
                - `x0_solver`: reduced initial state vector (after pre-equilibration)
                - `x`: full state vector
                - `x_solver`: reduced state vector
                - `y`: observables
                - `sigmay`: standard deviations of the observables
                - `tcl`: total values for conservation laws (at final timepoint)
                - `res`: residuals (observed - simulated)
        :return:
            Tuple of output value and simulation statistics
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
            adjoint=diffrax.RecursiveCheckpointAdjoint()
            if ret == "llh"
            else diffrax.DirectAdjoint(),
            ret=ret,
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
    ret: str = "llh",
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
        which output to return. Valid values are
            - `llh`: log-likelihood (default)
            - `nllhs`: negative log-likelihood at each time point
            - `x0`: full initial state vector (after pre-equilibration)
            - `x0_solver`: reduced initial state vector (after pre-equilibration)
            - `x`: full state vector
            - `x_solver`: reduced state vector
            - `y`: observables
            - `sigmay`: standard deviations of the observables
            - `tcl`: total values for conservation laws (at final timepoint)
            - `res`: residuals (observed - simulated)
    :return:
        Overall negative log-likelihood and condition specific results and statistics.
    """
    if simulation_conditions is None:
        simulation_conditions = problem.get_all_simulation_conditions()

    results = {
        sc: problem.run_simulation(sc, solver, controller, max_steps, ret)
        for sc in simulation_conditions
    }
    if ret == "llh":
        output = sum(llh for llh, _ in results.values())
    else:
        output = {sc: res for sc, (res, _) in results.items()}
    return output, results
