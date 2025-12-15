"""PEtab wrappers for JAX models.""" ""

import copy
import logging
import re
import shutil
from collections.abc import Callable, Iterable, Sized
from numbers import Number
from pathlib import Path

import diffrax
import equinox as eqx
import h5py
import jax.lax
import jax.numpy as jnp
import jaxtyping as jt
import numpy as np
import optimistix
import pandas as pd
import petab.v1 as petabv1
import petab.v2 as petabv2
from optimistix import AbstractRootFinder

from amici import _module_from_path
from amici.importers.petab.v1.parameter_mapping import (
    ParameterMappingForCondition,
    create_parameter_mapping,
)
from amici.jax.model import JAXModel, ReturnValue
from amici.logging import get_logger
from amici.sim.jax import get_simulation_conditions_v2

DEFAULT_CONTROLLER_SETTINGS = {
    "atol": 1e-8,
    "rtol": 1e-8,
    "pcoeff": 0.4,
    "icoeff": 0.3,
    "dcoeff": 0.0,
}

DEFAULT_ROOT_FINDER_SETTINGS = {
    "atol": 1e-12,
    "rtol": 1e-12,
}

SCALE_TO_INT = {
    petabv1.LIN: 0,
    petabv1.LOG: 1,
    petabv1.LOG10: 2,
}

logger = get_logger(__name__, logging.WARNING)


def jax_unscale(
    parameter: jnp.float_,
    scale_str: str,
) -> jnp.float_:
    """Unscale parameter according to ``scale_str``.

    Arguments:
        parameter:
            Parameter to be unscaled.
        scale_str:
            One of ``petabv1.LIN``, ``petabv1.LOG``, ``petabv1.LOG10``.

    Returns:
        The unscaled parameter.
    """
    if scale_str == petabv1.LIN or not scale_str:
        return parameter
    if scale_str == petabv1.LOG:
        return jnp.exp(parameter)
    if scale_str == petabv1.LOG10:
        return jnp.power(10, parameter)
    raise ValueError(f"Invalid parameter scaling: {scale_str}")


# IDEA: Implement this class in petab-sciml instead?
class HybridProblem(petabv1.Problem):
    hybridization_df: pd.DataFrame

    def __init__(self, petab_problem: petabv1.Problem):
        self.__dict__.update(petab_problem.__dict__)
        self.hybridization_df = _get_hybridization_df(petab_problem)

# Implement v2 version of this class - way around the missing extensions config too?
class HybridV2Problem(petabv2.Problem):
    hybridization_df: pd.DataFrame
    extensions_config: dict

    def __init__(self, petab_problem: petabv2.Problem):
        if not hasattr(petab_problem, "extensions_config"):
            self.extensions_config = {}
        self.__dict__.update(petab_problem.__dict__)
        self.hybridization_df = _get_hybridization_df(petab_problem)


def _get_hybridization_df(petab_problem):
    if not hasattr(petab_problem, "extensions_config"):
        return None
    
    if "sciml" in petab_problem.extensions_config:
        hybridizations = [
            pd.read_csv(hf, sep="\t", index_col=0)
            for hf in petab_problem.extensions_config["sciml"][
                "hybridization_files"
            ]
        ]
        hybridization_df = pd.concat(hybridizations)
        return hybridization_df


def _get_hybrid_petab_problem(petab_problem: petabv1.Problem | petabv2.Problem):
    if isinstance(petab_problem, petabv2.Problem):
        return HybridV2Problem(petab_problem)
    return HybridProblem(petab_problem)


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
        Preprocessed arrays for each simulation condition.
    :ivar _petab_problem:
        PEtab problem to simulate.
    """

    parameters: jnp.ndarray
    model: JAXModel
    simulation_conditions: tuple[tuple[str, ...], ...]
    _parameter_mappings: dict[str, ParameterMappingForCondition]
    _ts_dyn: np.ndarray
    _ts_posteq: np.ndarray
    _my: np.ndarray
    _iys: np.ndarray
    _iy_trafos: np.ndarray
    _ts_masks: np.ndarray
    _op_numeric: np.ndarray
    _op_mask: np.ndarray
    _op_indices: np.ndarray
    _np_numeric: np.ndarray
    _np_mask: np.ndarray
    _np_indices: np.ndarray
    _petab_measurement_indices: np.ndarray
    _petab_problem: petabv1.Problem | HybridProblem | petabv2.Problem

    def __init__(self, model: JAXModel, petab_problem: petabv1.Problem | petabv2.Problem):
        """
        Initialize a JAXProblem instance with a model and a PEtab problem.

        :param model:
            JAXModel instance to use for simulation.
        :param petab_problem:
            PEtab problem to simulate.
        """
        if isinstance(petab_problem, petabv2.Problem):
            scs = get_simulation_conditions_v2(petab_problem)
            self.simulation_conditions = scs.simulationConditionId
        else:
            scs = petab_problem.get_simulation_conditions_from_measurement_df()
            self.simulation_conditions = tuple(tuple(sc) for sc in scs.values)
        self._petab_problem = _get_hybrid_petab_problem(petab_problem)
        self.parameters, self.model = (
            self._initialize_model_with_nominal_values(model)
        )
        self._parameter_mappings = self._get_parameter_mappings(scs)
        (
            self._ts_dyn,
            self._ts_posteq,
            self._my,
            self._iys,
            self._iy_trafos,
            self._ts_masks,
            self._petab_measurement_indices,
            self._op_numeric,
            self._op_mask,
            self._op_indices,
            self._np_numeric,
            self._np_mask,
            self._np_indices,
        ) = self._get_measurements(scs)

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
        petab_problem = petabv1.Problem.from_yaml(
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
            :meth:`petabv1.Problem.get_simulation_conditions_from_measurement_df`.
        :return:
            Dictionary mapping simulation conditions to parameter mappings.
        """
        scs = list(set(simulation_conditions.simulationConditionId))
        petab_problem = copy.deepcopy(self._petab_problem)
        # remove observable and noise parameters from measurement dataframe as we are mapping them elsewhere
        petab_problem.measurement_df.drop(
            columns=[petabv1.OBSERVABLE_PARAMETERS, petabv1.NOISE_PARAMETERS],
            inplace=True,
            errors="ignore",
        )
        mappings = create_parameter_mapping(
            petab_problem=petab_problem,
            simulation_conditions=[
                {petabv1.SIMULATION_CONDITION_ID: sc} for sc in scs
            ],
            scaled_parameters=False,
            allow_timepoint_specific_numeric_noise_parameters=True,
        )
        # fill in dummy variables
        for mapping in mappings:
            for sim_var, value in mapping.map_sim_var.items():
                if isinstance(value, Number) and not np.isfinite(value):
                    mapping.map_sim_var[sim_var] = 1.0
        return dict(zip(scs, mappings, strict=True))

    def _get_measurements(
        self, simulation_conditions: pd.DataFrame
    ) -> tuple[
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
    ]:
        """
        Get measurements for the model based on the provided simulation conditions.

        :param simulation_conditions:
            Simulation conditions to create parameter mappings for. Same format as returned by
            :meth:`petabv1.Problem.get_simulation_conditions_from_measurement_df`.
        :return:
            tuple of padded
             - dynamic time points
             - post-equilibrium time points
             - measurements
             - observable indices
             - observable transformations indices
             - measurement masks
             - data indices (index in petab measurement dataframe).
             - numeric values for observable parameter overrides
             - non-numeric mask for observable parameter overrides
             - parameter indices (problem parameters) for observable parameter overrides
             - numeric values for noise parameter overrides
             - non-numeric mask for noise parameter overrides
             - parameter indices (problem parameters) for noise parameter overrides
        """
        measurements = dict()
        petab_indices = dict()

        n_pars = dict()
        for col in [petabv1.OBSERVABLE_PARAMETERS, petabv1.NOISE_PARAMETERS]:
            n_pars[col] = 0
            if col in self._petab_problem.measurement_df:
                if np.issubdtype(
                    self._petab_problem.measurement_df[col].dtype, np.number
                ):
                    n_pars[col] = 1 - int(
                        self._petab_problem.measurement_df[col].isna().all()
                    )
                else:
                    n_pars[col] = (
                        self._petab_problem.measurement_df[col]
                        .str.split(petabv1.C.PARAMETER_SEPARATOR)
                        .apply(
                            lambda x: len(x)
                            if isinstance(x, Sized)
                            else 1 - int(pd.isna(x))
                        )
                        .max()
                    )

        for _, simulation_condition in simulation_conditions.iterrows():
            query = " & ".join(
                [f"{k} == '{v}'" for k, v in simulation_condition.items()]
            )
            m = self._petab_problem.measurement_df.query(query).sort_values(
                by=petabv1.TIME
            )

            ts = m[petabv1.TIME]
            ts_dyn = ts[np.isfinite(ts)]
            ts_posteq = ts[np.logical_not(np.isfinite(ts))]
            index = pd.concat([ts_dyn, ts_posteq]).index
            ts_dyn = ts_dyn.values
            ts_posteq = ts_posteq.values
            my = m[petabv1.MEASUREMENT].values
            iys = np.array(
                [
                    self.model.observable_ids.index(oid)
                    for oid in m[petabv1.OBSERVABLE_ID].values
                ]
            )
            if (
                petabv1.OBSERVABLE_TRANSFORMATION
                in self._petab_problem.observable_df
            ):
                iy_trafos = np.array(
                    [
                        SCALE_TO_INT[
                            self._petab_problem.observable_df.loc[
                                oid, petabv1.OBSERVABLE_TRANSFORMATION
                            ]
                        ]
                        for oid in m[petabv1.OBSERVABLE_ID].values
                    ]
                )
            else:
                iy_trafos = np.zeros_like(iys)

            parameter_overrides_par_indices = dict()
            parameter_overrides_numeric_vals = dict()
            parameter_overrides_mask = dict()

            def get_parameter_override(x):
                if (
                    x in self._petab_problem.parameter_df.index
                    and not self._petab_problem.parameter_df.loc[
                        x, petabv1.ESTIMATE
                    ]
                ):
                    return self._petab_problem.parameter_df.loc[
                        x, petabv1.NOMINAL_VALUE
                    ]
                return x

            for col in [petabv1.OBSERVABLE_PARAMETERS, petabv1.NOISE_PARAMETERS]:
                if col not in m or m[col].isna().all():
                    mat_numeric = jnp.ones((len(m), n_pars[col]))
                    par_mask = np.zeros_like(mat_numeric, dtype=bool)
                    par_index = np.zeros_like(mat_numeric, dtype=int)
                elif np.issubdtype(m[col].dtype, np.number):
                    mat_numeric = np.expand_dims(m[col].values, axis=1)
                    par_mask = np.zeros_like(mat_numeric, dtype=bool)
                    par_index = np.zeros_like(mat_numeric, dtype=int)
                else:
                    split_vals = m[col].str.split(petabv1.C.PARAMETER_SEPARATOR)
                    list_vals = split_vals.apply(
                        lambda x: [get_parameter_override(y) for y in x]
                        if isinstance(x, list)
                        else []
                        if pd.isna(x)
                        else [
                            x
                        ]  # every string gets transformed to lists, so this is already a float
                    )
                    vals = list_vals.apply(
                        lambda x: np.pad(
                            x,
                            (0, n_pars[col] - len(x)),
                            mode="constant",
                            constant_values=1.0,
                        )
                    )
                    mat = np.stack(vals)
                    # deconstruct such that we can reconstruct mapped parameter overrides via vectorized operations
                    # mat = np.where(par_mask, map(lambda ip: p.at[ip], par_index), mat_numeric)
                    par_index = np.vectorize(
                        lambda x: self.parameter_ids.index(x)
                        if x in self.parameter_ids
                        else -1
                    )(mat)
                    # map out numeric values
                    par_mask = par_index != -1
                    # remove non-numeric values
                    mat[par_mask] = 0.0
                    mat_numeric = mat.astype(float)
                    # replace dummy index with some valid index
                    par_index[~par_mask] = 0

                parameter_overrides_numeric_vals[col] = mat_numeric
                parameter_overrides_mask[col] = par_mask
                parameter_overrides_par_indices[col] = par_index

            measurements[tuple(simulation_condition)] = (
                ts_dyn,  # 0
                ts_posteq,  # 1
                my,  # 2
                iys,  # 3
                iy_trafos,  # 4
                parameter_overrides_numeric_vals[
                    petabv1.OBSERVABLE_PARAMETERS
                ],  # 5
                parameter_overrides_mask[petabv1.OBSERVABLE_PARAMETERS],  # 6
                parameter_overrides_par_indices[
                    petabv1.OBSERVABLE_PARAMETERS
                ],  # 7
                parameter_overrides_numeric_vals[petabv1.NOISE_PARAMETERS],  # 8
                parameter_overrides_mask[petabv1.NOISE_PARAMETERS],  # 9
                parameter_overrides_par_indices[petabv1.NOISE_PARAMETERS],  # 10
            )
            petab_indices[tuple(simulation_condition)] = tuple(index.tolist())

        # compute maximum lengths
        n_ts_dyn = max(len(mv[0]) for mv in measurements.values())
        n_ts_posteq = max(len(mv[1]) for mv in measurements.values())

        # pad with last value and stack
        ts_dyn = np.stack(
            [
                np.pad(mv[0], (0, n_ts_dyn - len(mv[0])), mode="edge")
                for mv in measurements.values()
            ]
        )
        ts_posteq = np.stack(
            [
                np.pad(mv[1], (0, n_ts_posteq - len(mv[1])), mode="edge")
                for mv in measurements.values()
            ]
        )

        def pad_measurement(x_dyn, x_peq):
            # only pad first axis
            pad_width_dyn = tuple(
                [(0, n_ts_dyn - len(x_dyn))] + [(0, 0)] * (x_dyn.ndim - 1)
            )
            pad_width_peq = tuple(
                [(0, n_ts_posteq - len(x_peq))] + [(0, 0)] * (x_peq.ndim - 1)
            )
            return np.concatenate(
                (
                    np.pad(x_dyn, pad_width_dyn, mode="edge"),
                    np.pad(x_peq, pad_width_peq, mode="edge"),
                )
            )

        def pad_and_stack(output_index: int):
            return np.stack(
                [
                    pad_measurement(
                        mv[output_index][: len(mv[0])],
                        mv[output_index][len(mv[0]) :],
                    )
                    for mv in measurements.values()
                ]
            )

        my = pad_and_stack(2)
        iys = pad_and_stack(3)
        iy_trafos = pad_and_stack(4)
        op_numeric = pad_and_stack(5)
        op_mask = pad_and_stack(6)
        op_indices = pad_and_stack(7)
        np_numeric = pad_and_stack(8)
        np_mask = pad_and_stack(9)
        np_indices = pad_and_stack(10)
        ts_masks = np.stack(
            [
                np.concatenate(
                    (
                        np.pad(
                            np.ones_like(mv[0]), (0, n_ts_dyn - len(mv[0]))
                        ),
                        np.pad(
                            np.ones_like(mv[1]), (0, n_ts_posteq - len(mv[1]))
                        ),
                    )
                )
                for mv in measurements.values()
            ]
        ).astype(bool)
        petab_indices = np.stack(
            [
                pad_measurement(
                    np.array(idx[: len(mv[0])]),
                    np.array(idx[len(mv[0]) :]),
                )
                for mv, idx in zip(
                    measurements.values(), petab_indices.values()
                )
            ]
        )

        return (
            ts_dyn,
            ts_posteq,
            my,
            iys,
            iy_trafos,
            ts_masks,
            petab_indices,
            op_numeric,
            op_mask,
            op_indices,
            np_numeric,
            np_mask,
            np_indices,
        )

    def get_all_simulation_conditions(self) -> tuple[tuple[str, ...], ...]:
        simulation_conditions = (
            self._petab_problem.get_simulation_conditions_from_measurement_df()
        )
        return tuple(tuple(row) for _, row in simulation_conditions.iterrows())

    def _initialize_model_parameters(self, model: JAXModel) -> dict:
        """
        Initialize model parameter structure with zeros.

        :param model:
            JAX model with neural networks

        :return:
            Nested dictionary structure for model parameters
        """
        return {
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

    def _load_parameter_arrays_from_files(self) -> dict:
        """
        Load neural network parameter arrays from HDF5 files.

        :return:
            Dictionary mapping network IDs to parameter arrays
        """
        if not self._petab_problem.extensions_config:
            return {}

        array_files = self._petab_problem.extensions_config["sciml"].get(
            "array_files", []
        )

        return {
            file_spec.split("_")[0]: h5py.File(file_spec, "r")["parameters"][
                file_spec.split("_")[0]
            ]
            for file_spec in array_files
            if "parameters" in h5py.File(file_spec, "r").keys()
        }

    def _load_input_arrays_from_files(self) -> dict:
        """
        Load neural network input arrays from HDF5 files.

        :return:
            Dictionary mapping network IDs to input arrays
        """
        if not self._petab_problem.extensions_config:
            return {}

        array_files = self._petab_problem.extensions_config["sciml"].get(
            "array_files", []
        )

        return {
            file_spec.split("_")[0]: h5py.File(file_spec, "r")["inputs"]
            for file_spec in array_files
            if "inputs" in h5py.File(file_spec, "r").keys()
        }

    def _parse_parameter_name(
        self, pname: str, model_pars: dict
    ) -> list[tuple[str, str]]:
        """
        Parse parameter name to determine which layers and attributes to set.

        :param pname:
            Parameter name from PEtab (format: net.layer.attribute)
        :param model_pars:
            Model parameters dictionary

        :return:
            List of (layer_name, attribute_name) tuples to set
        """
        net = pname.split("_")[0]
        nn = model_pars[net]
        to_set = []

        name_parts = pname.split(".")

        if len(name_parts) > 1:
            layer_name = name_parts[1]
            layer = nn[layer_name]
            if len(name_parts) > 2:
                # Specific attribute specified
                attribute_name = name_parts[2]
                to_set.append((layer_name, attribute_name))
            else:
                # All attributes of the layer
                to_set.extend(
                    [(layer_name, attribute) for attribute in layer.keys()]
                )
        else:
            # All layers and attributes
            to_set.extend(
                [
                    (layer_name, attribute)
                    for layer_name, layer in nn.items()
                    for attribute in layer.keys()
                ]
            )

        return to_set

    def _extract_nominal_values_from_petab(
        self, model: JAXModel, model_pars: dict, par_arrays: dict
    ) -> None:
        """
        Extract nominal parameter values from PEtab problem and populate model_pars.

        :param model:
            JAX model
        :param model_pars:
            Model parameters dictionary to populate (modified in place)
        :param par_arrays:
            Parameter arrays loaded from files
        """
        for pname, row in self._petab_problem.parameter_df.iterrows():
            net = pname.split("_")[0]
            if net not in model.nns:
                continue

            nn = model_pars[net]
            scalar = True

            # Determine value source (scalar from PEtab or array from file)
            if np.isnan(row[petabv1.NOMINAL_VALUE]):
                value = par_arrays[net]
                scalar = False
            else:
                value = float(row[petabv1.NOMINAL_VALUE])

            # Parse parameter name and set values
            to_set = self._parse_parameter_name(pname, model_pars)

            for layer, attribute in to_set:
                if scalar:
                    nn[layer][attribute] = value * jnp.ones_like(
                        getattr(model.nns[net].layers[layer], attribute)
                    )
                else:
                    nn[layer][attribute] = jnp.array(
                        value[layer][attribute][:]
                    )

    def _set_model_parameters(
        self, model: JAXModel, model_pars: dict
    ) -> JAXModel:
        """
        Set parameter values in the model using equinox tree_at.

        :param model:
            JAX model to update
        :param model_pars:
            Dictionary of parameter values to set

        :return:
            Updated JAX model
        """
        for net_id in model_pars:
            for layer_id in model_pars[net_id]:
                for attribute in model_pars[net_id][layer_id]:
                    logger.debug(
                        f"Setting {attribute} of layer {layer_id} in network "
                        f"{net_id} to {model_pars[net_id][layer_id][attribute]}"
                    )
                    model = eqx.tree_at(
                        lambda model: getattr(
                            model.nns[net_id].layers[layer_id], attribute
                        ),
                        model,
                        model_pars[net_id][layer_id][attribute],
                    )
        return model

    def _set_input_arrays(
        self, model: JAXModel, nn_input_arrays: dict, model_pars: dict
    ) -> JAXModel:
        """
        Set input arrays in the model if provided.

        :param model:
            JAX model to update
        :param nn_input_arrays:
            Input arrays loaded from files
        :param model_pars:
            Model parameters dictionary (for network IDs)

        :return:
            Updated JAX model
        """
        if len(nn_input_arrays) == 0:
            return model

        for net_id in model_pars:
            input_array = {
                input: {
                    k: jnp.array(
                        arr[:],
                        dtype=jnp.float64
                        if jax.config.jax_enable_x64
                        else jnp.float32,
                    )
                    for k, arr in nn_input_arrays[net_id][input].items()
                }
                for input in model.nns[net_id].inputs
            }
            model = eqx.tree_at(
                lambda model: model.nns[net_id].inputs, model, input_array
            )

        return model

    def _create_scaled_parameter_array(self) -> jt.Float[jt.Array, "np"]:
        """
        Create array of scaled nominal parameter values for estimation.

        :return:
            JAX array of scaled parameter values
        """
        return jnp.array(
            [
                petabv1.scale(
                    float(
                        self._petab_problem.parameter_df.loc[
                            pval, petabv1.NOMINAL_VALUE
                        ]
                    ),
                    self._petab_problem.parameter_df.loc[
                        pval, petabv1.PARAMETER_SCALE
                    ],
                )
                for pval in self.parameter_ids
            ]
        )

    def _initialize_model_with_nominal_values(
        self, model: JAXModel
    ) -> tuple[jt.Float[jt.Array, "np"], JAXModel]:
        """
        Initialize the model with nominal parameter values and inputs from the PEtab problem.

        This method:
        - Initializes model parameter structure
        - Loads parameter and input arrays from HDF5 files
        - Extracts nominal values from PEtab problem
        - Sets parameter values in the model
        - Sets input arrays in the model
        - Creates scaled parameter array to initialized to nominal values

        :param model:
            JAX model to initialize

        :return:
            Tuple of (scaled parameter array, initialized model)
        """
        # Initialize model parameters structure
        model_pars = self._initialize_model_parameters(model)

        # Load arrays from files (getters)
        par_arrays = self._load_parameter_arrays_from_files()
        nn_input_arrays = self._load_input_arrays_from_files()

        # Extract nominal values from PEtab problem
        self._extract_nominal_values_from_petab(model, model_pars, par_arrays)

        # Set values in model (setters)
        model = self._set_model_parameters(model, model_pars)
        model = self._set_input_arrays(model, nn_input_arrays, model_pars)

        # Create scaled parameter array
        parameter_array = self._create_scaled_parameter_array()

        return parameter_array, model

    def _get_inputs(self) -> dict:
        if self._petab_problem.mapping_df is None:
            return {}
        inputs = {net: {} for net in self.model.nns.keys()}
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
                inputs[row["netId"]][row[petabv1.MODEL_ENTITY_ID]] = data_flat[
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
            self._petab_problem.parameter_df[petabv1.ESTIMATE]
            == 1
            & pd.to_numeric(
                self._petab_problem.parameter_df[petabv1.NOMINAL_VALUE],
                errors="coerce",
            ).notna()
        ].index.tolist()

    @property
    def nn_output_ids(self) -> list[str]:
        """
        Parameter ids that are estimated in the PEtab problem. Same ordering as values in :attr:`parameters`.

        :return:
            PEtab parameter ids
        """
        if self._petab_problem.mapping_df is None:
            return []
        return self._petab_problem.mapping_df[
            self._petab_problem.mapping_df[petabv1.MODEL_ENTITY_ID]
            .str.split(".")
            .str[1]
            .str.startswith("output")
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

    def _eval_nn(self, output_par: str, condition_id: str):
        net_id = self._petab_problem.mapping_df.loc[
            output_par, petabv1.MODEL_ENTITY_ID
        ].split(".")[0]
        nn = self.model.nns[net_id]

        def _is_net_input(model_id):
            comps = model_id.split(".")
            return comps[0] == net_id and comps[1].startswith("inputs")

        model_id_map = (
            self._petab_problem.mapping_df[
                self._petab_problem.mapping_df[petabv1.MODEL_ENTITY_ID].apply(
                    _is_net_input
                )
            ]
            .reset_index()
            .set_index(petabv1.MODEL_ENTITY_ID)[petabv1.PETAB_ENTITY_ID]
            .to_dict()
        )

        condition_input_map = (
            dict(
                [
                    (
                        petab_id,
                        self._petab_problem.parameter_df.loc[
                            self._petab_problem.condition_df.loc[
                                condition_id, petab_id
                            ],
                            petabv1.NOMINAL_VALUE,
                        ],
                    )
                    if self._petab_problem.condition_df.loc[
                        condition_id, petab_id
                    ]
                    in self._petab_problem.parameter_df.index
                    else (
                        petab_id,
                        np.float64(
                            self._petab_problem.condition_df.loc[
                                condition_id, petab_id
                            ]
                        ),
                    )
                    for petab_id in model_id_map.values()
                ]
            )
            if not self._petab_problem.condition_df.empty
            else {}
        )

        hybridization_parameter_map = {
            petab_id: self._petab_problem.hybridization_df.loc[
                petab_id, "targetValue"
            ]
            for petab_id in model_id_map.values()
            if petab_id in set(self._petab_problem.hybridization_df.index)
        }

        # handle conditions
        if len(condition_input_map) > 0:
            net_input = jnp.array(
                [
                    condition_input_map[petab_id]
                    for _, petab_id in model_id_map.items()
                ]
            )
            return nn.forward(net_input).squeeze()

        # handle array inputs
        if isinstance(self.model.nns[net_id].inputs, dict):
            net_input = jnp.array(
                [
                    self.model.nns[net_id].inputs[petab_id][condition_id]
                    if condition_id in self.model.nns[net_id].inputs[petab_id]
                    else self.model.nns[net_id].inputs[petab_id]["0"]
                    for _, petab_id in model_id_map.items()
                ]
            )
            return nn.forward(net_input).squeeze()

        net_input = jnp.array(
            [
                jax.lax.stop_gradient(self.model.nns[net_id][model_id])
                if model_id in self.model.nns[net_id].inputs
                else self.get_petab_parameter_by_id(petab_id)
                if petab_id in self.parameter_ids
                else self._petab_problem.parameter_df.loc[
                    petab_id, petabv1.NOMINAL_VALUE
                ]
                if petab_id in set(self._petab_problem.parameter_df.index)
                else self._petab_problem.parameter_df.loc[
                    hybridization_parameter_map[petab_id], petabv1.NOMINAL_VALUE
                ]
                for model_id, petab_id in model_id_map.items()
            ]
        )
        return nn.forward(net_input).squeeze()

    def _map_model_parameter_value(
        self,
        mapping: ParameterMappingForCondition,
        pname: str,
        condition_id: str,
    ) -> jt.Float[jt.Scalar, ""] | float:  # noqa: F722
        pval = mapping.map_sim_var[pname]
        if hasattr(self, "nn_output_ids") and pval in self.nn_output_ids:
            nn_output = self._eval_nn(pval, condition_id)
            if nn_output.size > 1:
                entityId = self._petab_problem.mapping_df.loc[
                    pval, petabv1.MODEL_ENTITY_ID
                ]
                ind = int(re.search(r"\[\d+\]\[(\d+)\]", entityId).group(1))
                return nn_output[ind]
            else:
                return nn_output
        if isinstance(pval, Number):
            return pval
        return self.get_petab_parameter_by_id(pval)

    def load_model_parameters(
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
                self._map_model_parameter_value(
                    mapping, pname, simulation_condition
                )
                for pname in self.model.parameter_ids
            ]
        )
        pscale = tuple(
            [
                petabv1.LIN
                if self._petab_problem.mapping_df is not None
                and pname in self._petab_problem.mapping_df.index
                else mapping.scale_map_sim_var[pname]
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
        if state_id in self.nn_output_ids:
            return True

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
        if state_id in self.nn_output_ids:
            return self._eval_nn(state_id)

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
                    xval, petabv1.PARAMETER_SCALE
                ],
            )
        # only remaining option is nominal value for PEtab parameter
        # that is not estimated, return nominal value
        return self._petab_problem.parameter_df.loc[xval, petabv1.NOMINAL_VALUE]

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
            or hasattr(self, "nn_output_ids")
            and x_id in self.nn_output_ids
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

    def _prepare_conditions(
        self,
        conditions: list[str],
        op_numeric: np.ndarray | None = None,
        op_mask: np.ndarray | None = None,
        op_indices: np.ndarray | None = None,
        np_numeric: np.ndarray | None = None,
        np_mask: np.ndarray | None = None,
        np_indices: np.ndarray | None = None,
    ) -> tuple[
        jt.Float[jt.Array, "nc np"],  # noqa: F821, F722
        jt.Bool[jt.Array, "nx"],  # noqa: F821
        jt.Float[jt.Array, "nx"],  # noqa: F821
        jt.Float[jt.Array, "nc nt nop"],  # noqa: F821, F722
        jt.Float[jt.Array, "nc nt nnp"],  # noqa: F821, F722
    ]:
        """
        Prepare conditions for simulation.

        :param conditions:
            Simulation conditions to prepare.
        :param op_numeric:
            Numeric values for observable parameter overrides. If None, no overrides are used.
        :param op_mask:
            Mask for observable parameter overrides. True for free parameter overrides, False for numeric values.
        :param op_indices:
            Free parameter indices (wrt. `self.parameters`) for observable parameter overrides.
        :param np_numeric:
            Numeric values for noise parameter overrides. If None, no overrides are used.
        :param np_mask:
            Mask for noise parameter overrides. True for free parameter overrides, False for numeric values.
        :param np_indices:
            Free parameter indices (wrt. `self.parameters`) for noise parameter overrides.
        :return:
            Tuple of parameter arrays, reinitialisation masks and reinitialisation values, observable parameters and
            noise parameters.
        """
        p_array = jnp.stack(
            [self.load_model_parameters(sc) for sc in conditions]
        )

        if self.parameters.size:
            unscaled_parameters = jnp.stack(
                [
                    jax_unscale(
                        self.parameters[ip],
                        self._petab_problem.parameter_df.loc[
                            p_id, petabv1.PARAMETER_SCALE
                        ],
                    )
                    for ip, p_id in enumerate(self.parameter_ids)
                ]
            )
        else:
            unscaled_parameters = jnp.zeros((*self._ts_masks.shape[:2], 0))

        if op_numeric is not None and op_numeric.size:
            op_array = jnp.where(
                op_mask,
                jax.vmap(
                    jax.vmap(jax.vmap(lambda ip: unscaled_parameters[ip]))
                )(op_indices),
                op_numeric,
            )
        else:
            op_array = jnp.zeros((*self._ts_masks.shape[:2], 0))

        if np_numeric is not None and np_numeric.size:
            np_array = jnp.where(
                np_mask,
                jax.vmap(
                    jax.vmap(jax.vmap(lambda ip: unscaled_parameters[ip]))
                )(np_indices),
                np_numeric,
            )
        else:
            np_array = jnp.zeros((*self._ts_masks.shape[:2], 0))

        mask_reinit_array = jnp.stack(
            [
                self.load_reinitialisation(sc, p)[0]
                for sc, p in zip(conditions, p_array)
            ]
        )
        x_reinit_array = jnp.stack(
            [
                self.load_reinitialisation(sc, p)[1]
                for sc, p in zip(conditions, p_array)
            ]
        )
        return p_array, mask_reinit_array, x_reinit_array, op_array, np_array

    @eqx.filter_vmap(
        in_axes={
            "max_steps": None,
            "self": None,
        },  # only list arguments here where eqx.is_array(0) is not the right thing
    )
    def run_simulation(
        self,
        p: jt.Float[jt.Array, "np"],  # noqa: F821, F722
        ts_dyn: np.ndarray,
        ts_posteq: np.ndarray,
        my: np.ndarray,
        iys: np.ndarray,
        iy_trafos: np.ndarray,
        ops: jt.Float[jt.Array, "nt *nop"],  # noqa: F821, F722
        nps: jt.Float[jt.Array, "nt *nnp"],  # noqa: F821, F722
        mask_reinit: jt.Bool[jt.Array, "nx"],  # noqa: F821, F722
        x_reinit: jt.Float[jt.Array, "nx"],  # noqa: F821, F722
        init_override: jt.Float[jt.Array, "nx"],  # noqa: F821, F722
        init_override_mask: jt.Bool[jt.Array, "nx"],  # noqa: F821, F722
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        root_finder: AbstractRootFinder,
        steady_state_event: Callable[
            ..., diffrax._custom_types.BoolScalarLike
        ],
        max_steps: jnp.int_,
        x_preeq: jt.Float[jt.Array, "*nx"] = jnp.array([]),  # noqa: F821, F722
        ts_mask: np.ndarray = np.array([]),
        ret: ReturnValue = ReturnValue.llh,
    ) -> tuple[jnp.float_, dict]:
        """
        Run a simulation for a given simulation condition.

        :param p:
            Parameters for the simulation condition
        :param ts_dyn:
            (Padded) dynamic time points
        :param ts_posteq:
            (Padded) post-equilibrium time points
        :param my:
            (Padded) measurements
        :param iys:
            (Padded) observable indices
        :param iy_trafos:
            (Padded) observable transformations indices
        :param ops:
            (Padded) observable parameters
        :param nps:
            (Padded) noise parameters
        :param mask_reinit:
            Mask for states that need reinitialisation
        :param x_reinit:
            Reinitialisation values for states
        :param solver:
            ODE solver to use for simulation
        :param controller:
            Step size controller to use for simulation
        :param steady_state_event:
            Steady state event function to use for post-equilibration. Allows customisation of the steady state
            condition, see :func:`diffrax.steady_state_event` for details.
        :param max_steps:
            Maximum number of steps to take during simulation
        :param x_preeq:
            Pre-equilibration state. Can be empty if no pre-equilibration is available, in which case the states will
            be initialised to the model default values.
        :param ts_mask:
            padding mask, see :meth:`JAXModel.simulate_condition` for details.
        :param ret:
            which output to return. See :class:`ReturnValue` for available options.
        :return:
            Tuple of output value and simulation statistics
        """
        return self.model.simulate_condition(
            p=p,
            ts_dyn=jax.lax.stop_gradient(jnp.array(ts_dyn)),
            ts_posteq=jax.lax.stop_gradient(jnp.array(ts_posteq)),
            my=jax.lax.stop_gradient(jnp.array(my)),
            iys=jax.lax.stop_gradient(jnp.array(iys)),
            iy_trafos=jax.lax.stop_gradient(jnp.array(iy_trafos)),
            nps=nps,
            ops=ops,
            x_preeq=x_preeq,
            mask_reinit=jax.lax.stop_gradient(mask_reinit),
            x_reinit=x_reinit,
            init_override=init_override,
            init_override_mask=jax.lax.stop_gradient(init_override_mask),
            ts_mask=jax.lax.stop_gradient(jnp.array(ts_mask)),
            solver=solver,
            controller=controller,
            root_finder=root_finder,
            max_steps=max_steps,
            steady_state_event=steady_state_event,
            adjoint=diffrax.RecursiveCheckpointAdjoint()
            if ret in (ReturnValue.llh, ReturnValue.chi2)
            else diffrax.DirectAdjoint(),
            ret=ret,
        )

    def run_simulations(
        self,
        simulation_conditions: list[str],
        preeq_array: jt.Float[jt.Array, "ncond *nx"],  # noqa: F821, F722
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        root_finder: AbstractRootFinder,
        steady_state_event: Callable[
            ..., diffrax._custom_types.BoolScalarLike
        ],
        max_steps: jnp.int_,
        ret: ReturnValue = ReturnValue.llh,
    ):
        """
        Run simulations for a list of simulation conditions.

        :param simulation_conditions:
            List of simulation conditions to run simulations for.
        :param preeq_array:
            Matrix of pre-equilibrated states for the simulation conditions. Ordering must match the simulation
            conditions. If no pre-equilibration is available for a condition, the corresponding row must be empty.
        :param solver:
            ODE solver to use for simulation.
        :param controller:
            Step size controller to use for simulation.
        :param steady_state_event:
            Steady state event function to use for post-equilibration. Allows customisation of the steady state
            condition, see :func:`diffrax.steady_state_event` for details.
        :param max_steps:
            Maximum number of steps to take during simulation.
        :param ret:
            which output to return. See :class:`ReturnValue` for available options.
        :return:
            Output value and condition specific results and statistics. Results and statistics are returned as a dict
            with arrays with the leading dimension corresponding to the simulation conditions.
        """
        p_array, mask_reinit_array, x_reinit_array, op_array, np_array = (
            self._prepare_conditions(
                simulation_conditions,
                self._op_numeric,
                self._op_mask,
                self._op_indices,
                self._np_numeric,
                self._np_mask,
                self._np_indices,
            )
        )

        init_override_mask = jnp.stack(
            [
                jnp.array(
                    [
                        p
                        in set(self._parameter_mappings[sc].map_sim_var.keys())
                        for p in self.model.state_ids
                    ]
                )
                for sc in simulation_conditions
            ]
        )
        init_override = jnp.stack(
            [
                jnp.array(
                    [
                        self._eval_nn(
                            self._parameter_mappings[sc].map_sim_var[p], sc
                        )
                        if p
                        in set(self._parameter_mappings[sc].map_sim_var.keys())
                        else 1.0
                        for p in self.model.state_ids
                    ]
                )
                for sc in simulation_conditions
            ]
        )

        return self.run_simulation(
            p_array,
            self._ts_dyn,
            self._ts_posteq,
            self._my,
            self._iys,
            self._iy_trafos,
            op_array,
            np_array,
            mask_reinit_array,
            x_reinit_array,
            init_override,
            init_override_mask,
            solver,
            controller,
            root_finder,
            steady_state_event,
            max_steps,
            preeq_array,
            self._ts_masks,
            ret,
        )

    @eqx.filter_vmap(
        in_axes={
            "max_steps": None,
            "self": None,
        },  # only list arguments here where eqx.is_array(0) is not the right thing
    )
    def run_preequilibration(
        self,
        p: jt.Float[jt.Array, "np"],  # noqa: F821, F722
        mask_reinit: jt.Bool[jt.Array, "nx"],  # noqa: F821, F722
        x_reinit: jt.Float[jt.Array, "nx"],  # noqa: F821, F722
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        root_finder: AbstractRootFinder,
        steady_state_event: Callable[
            ..., diffrax._custom_types.BoolScalarLike
        ],
        max_steps: jnp.int_,
    ) -> tuple[jt.Float[jt.Array, "nx"], dict]:  # noqa: F821
        """
        Run a pre-equilibration simulation for a given simulation condition.

        :param p:
            Parameters for the simulation condition
        :param mask_reinit:
            Mask for states that need reinitialisation
        :param x_reinit:
            Reinitialisation values for states
        :param solver:
            ODE solver to use for simulation
        :param controller:
            Step size controller to use for simulation
        :param steady_state_event:
            Steady state event function to use for pre-equilibration. Allows customisation of the steady state
            condition, see :func:`diffrax.steady_state_event` for details.
        :param max_steps:
            Maximum number of steps to take during simulation
        :return:
            Pre-equilibration state
        """
        return self.model.preequilibrate_condition(
            p=p,
            mask_reinit=mask_reinit,
            x_reinit=x_reinit,
            solver=solver,
            controller=controller,
            root_finder=root_finder,
            max_steps=max_steps,
            steady_state_event=steady_state_event,
        )

    def run_preequilibrations(
        self,
        simulation_conditions: list[str],
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        root_finder: AbstractRootFinder,
        steady_state_event: Callable[
            ..., diffrax._custom_types.BoolScalarLike
        ],
        max_steps: jnp.int_,
    ):
        p_array, mask_reinit_array, x_reinit_array, _, _ = (
            self._prepare_conditions(simulation_conditions, None, None)
        )
        return self.run_preequilibration(
            p_array,
            mask_reinit_array,
            x_reinit_array,
            solver,
            controller,
            root_finder,
            steady_state_event,
            max_steps,
        )


def run_simulations(
    problem: JAXProblem,
    simulation_conditions: Iterable[tuple[str, ...]] | None = None,
    solver: diffrax.AbstractSolver = diffrax.Kvaerno5(),
    controller: diffrax.AbstractStepSizeController = diffrax.PIDController(
        **DEFAULT_CONTROLLER_SETTINGS
    ),
    root_finder: AbstractRootFinder = optimistix.Newton(
        **DEFAULT_ROOT_FINDER_SETTINGS
    ),
    steady_state_event: Callable[
        ..., diffrax._custom_types.BoolScalarLike
    ] = diffrax.steady_state_event(),
    max_steps: int = 2**10,
    ret: ReturnValue | str = ReturnValue.llh,
):
    """
    Run simulations for a problem.

    :param problem:
        Problem to run simulations for.
    :param simulation_conditions:
        Simulation conditions to run simulations for. This is a series of tuples, where each tuple contains the
        simulation condition or the pre-equilibration condition followed by the simulation condition. Default is to run
        simulations for all conditions.
    :param solver:
        ODE solver to use for simulation.
    :param controller:
        Step size controller to use for simulation.
    :param root_finder:
        Root finder to use for event detection.
    :param steady_state_event:
        Steady state event function to use for pre-/post-equilibration. Allows customisation of the steady state
        condition, see :func:`diffrax.steady_state_event` for details.
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

    dynamic_conditions = [sc[0] for sc in simulation_conditions]
    preequilibration_conditions = list(
        {sc[1] for sc in simulation_conditions if len(sc) > 1}
    )

    conditions = {
        "dynamic_conditions": dynamic_conditions,
        "preequilibration_conditions": preequilibration_conditions,
        "simulation_conditions": simulation_conditions,
    }

    if preequilibration_conditions:
        preeqs, preresults = problem.run_preequilibrations(
            preequilibration_conditions,
            solver,
            controller,
            root_finder,
            steady_state_event,
            max_steps,
        )
    else:
        preresults = {
            "stats_preeq": None,
        }

    if dynamic_conditions:
        preeq_array = jnp.stack(
            [
                preeqs[preequilibration_conditions.index(sc[1]), :]
                if len(sc) > 1
                else jnp.array([])
                for sc in simulation_conditions
            ]
        )
        output, results = problem.run_simulations(
            dynamic_conditions,
            preeq_array,
            solver,
            controller,
            root_finder,
            steady_state_event,
            max_steps,
            ret,
        )
    else:
        output = jnp.array(0.0)
        results = {
            "llh": jnp.array([]),
            "stats_dyn": None,
            "stats_posteq": None,
            "ts": jnp.array([]),
            "x": jnp.array([]),
        }

    if ret in (ReturnValue.llh, ReturnValue.chi2):
        output = jnp.sum(output)

    return output, results | preresults | conditions


def petab_simulate(
    problem: JAXProblem,
    solver: diffrax.AbstractSolver = diffrax.Kvaerno5(),
    controller: diffrax.AbstractStepSizeController = diffrax.PIDController(
        **DEFAULT_CONTROLLER_SETTINGS
    ),
    steady_state_event: Callable[
        ..., diffrax._custom_types.BoolScalarLike
    ] = diffrax.steady_state_event(),
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
    :param steady_state_event:
        Steady state event function to use for pre-/post-equilibration. Allows customisation of the steady state
        condition, see :func:`diffrax.steady_state_event` for details.
    :return:
        petab simulation dataframe.
    """
    y, r = run_simulations(
        problem,
        solver=solver,
        controller=controller,
        steady_state_event=steady_state_event,
        max_steps=max_steps,
        ret=ReturnValue.y,
    )
    dfs = []
    for ic, sc in enumerate(r["dynamic_conditions"]):
        obs = [
            problem.model.observable_ids[io]
            for io in problem._iys[ic, problem._ts_masks[ic, :]]
        ]
        t = jnp.concat(
            (
                problem._ts_dyn[ic, :],
                problem._ts_posteq[ic, :],
            )
        )
        df_sc = pd.DataFrame(
            {
                petabv1.SIMULATION: y[ic, problem._ts_masks[ic, :]],
                petabv1.TIME: t[problem._ts_masks[ic, :]],
                petabv1.OBSERVABLE_ID: obs,
                petabv1.SIMULATION_CONDITION_ID: [sc] * len(t),
            },
            index=problem._petab_measurement_indices[ic, :],
        )
        if (
            petabv1.OBSERVABLE_PARAMETERS
            in problem._petab_problem.measurement_df
        ):
            df_sc[petabv1.OBSERVABLE_PARAMETERS] = (
                problem._petab_problem.measurement_df.query(
                    f"{petabv1.SIMULATION_CONDITION_ID} == '{sc}'"
                )[petabv1.OBSERVABLE_PARAMETERS]
            )
        if petabv1.NOISE_PARAMETERS in problem._petab_problem.measurement_df:
            df_sc[petabv1.NOISE_PARAMETERS] = (
                problem._petab_problem.measurement_df.query(
                    f"{petabv1.SIMULATION_CONDITION_ID} == '{sc}'"
                )[petabv1.NOISE_PARAMETERS]
            )
        if (
            petabv1.PREEQUILIBRATION_CONDITION_ID
            in problem._petab_problem.measurement_df
        ):
            df_sc[petabv1.PREEQUILIBRATION_CONDITION_ID] = (
                problem._petab_problem.measurement_df.query(
                    f"{petabv1.SIMULATION_CONDITION_ID} == '{sc}'"
                )[petabv1.PREEQUILIBRATION_CONDITION_ID]
            )
        dfs.append(df_sc)
    return pd.concat(dfs).sort_index()
