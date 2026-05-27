import os
from contextlib import contextmanager
from pathlib import Path

import equinox as eqx
import jax
import jax.numpy as jnp
import numpy as np
import pandas as pd
import pytest
from amici.importers.petab import *
from amici.sim.jax import run_simulations
from amici.sim.jax.petab import _try_float
from petab import v1, v2
from yaml import safe_load


@contextmanager
def change_directory(destination):
    # Save the current working directory
    original_directory = os.getcwd()
    try:
        # Change to the new directory
        os.chdir(destination)
        yield
    finally:
        # Change back to the original directory
        os.chdir(original_directory)


jax.config.update("jax_enable_x64", True)

cases_dir = Path(__file__).parent / "testsuite" / "Benchmark-Models"


@pytest.mark.parametrize(
    "test", sorted([d.name for d in cases_dir.iterdir() if d.is_dir()])
)
def test_benchmark(test):
    test_dir = cases_dir / test

    for split in ("train", "validation"):
        with open(test_dir / f"{split}_{test}.yaml") as f:
            petab_yaml = safe_load(f)
        with open(test_dir / f"expected_{test}.yaml") as f:
            solutions = safe_load(f)

        with change_directory(test_dir):
            # HACK!! Again!! Around "array" in parameters table
            petab_problem = _v2_sciml_problem_helper(petab_yaml, test_dir)

            pi = PetabImporter(
                petab_problem=petab_problem,
                module_name="hybrid" + test,
                compile_=True,
                jax=jax,
                validate=False,  # And again...around "array" in parameters table
            )

            jax_problem = pi.create_simulator(
                force_import=True,
            )

        # llh
        llh, _ = run_simulations(jax_problem)
        solution_split = "training" if split == "train" else "validation"
        np.testing.assert_allclose(
            llh,
            solutions["llh"][solution_split],
            atol=1e-3,
            rtol=1e-3,
        )


def _v2_sciml_problem_helper(yaml_config, base_path):
    config = v2.ProblemConfig(**yaml_config)

    parameter_tables = []
    for f in config.parameter_files:
        df = pd.read_csv(f, sep="\t")
        df.nominalValue = df.nominalValue.apply(_try_float)
        if "priorParameters" in df.columns:
            df.priorParameters = df.priorParameters.apply(
                _process_prior_params
            )
        parameters = [
            v2.Parameter.model_construct(**row.to_dict())
            for _, row in df.reset_index().iterrows()
        ]
        parameter_tables.append(v2.ParameterTable(elements=parameters))

    models = [
        v1.models.model.model_factory(
            model_info.location,
            base_path=base_path,
            model_language=model_info.language,
            model_id=model_id,
        )
        for model_id, model_info in (config.model_files or {}).items()
    ]

    measurement_tables = (
        [
            v2.MeasurementTable.from_tsv(f, base_path)
            for f in config.measurement_files
        ]
        if config.measurement_files
        else None
    )

    experiment_tables = (
        [
            v2.ExperimentTable.from_tsv(f, base_path)
            for f in config.experiment_files
        ]
        if config.experiment_files
        else None
    )

    condition_tables = (
        [
            v2.ConditionTable.from_tsv(f, base_path)
            for f in config.condition_files
        ]
        if config.condition_files
        else None
    )

    if condition_tables is None and experiment_tables is not None:
        cond_ids = [
            cid
            for exp_table in experiment_tables
            for exp in exp_table.elements
            for p in exp.periods
            for cid in p.condition_ids
        ]
        condition_tables = [
            v2.ConditionTable(elements=[v2.Condition(id=cid, changes=[])])
            for cid in set(cond_ids)
        ]

    observable_tables = (
        [
            v2.ObservableTable.from_tsv(f, base_path)
            for f in config.observable_files
        ]
        if config.observable_files
        else None
    )

    mapping_tables = (
        [v2.MappingTable.from_tsv(f, base_path) for f in config.mapping_files]
        if config.mapping_files
        else None
    )

    return v2.Problem(
        config=config,
        models=models,
        condition_tables=condition_tables,
        experiment_tables=experiment_tables,
        observable_tables=observable_tables,
        measurement_tables=measurement_tables,
        parameter_tables=parameter_tables,
        mapping_tables=mapping_tables,
    )


def _process_prior_params(prior_params):
    if isinstance(prior_params, float):
        return prior_params
    else:
        return [float(param) for param in prior_params.split(";")]


def _normal_logpdf(x: jnp.ndarray, mean: float, std: float) -> jnp.ndarray:
    var = std**2
    return jnp.sum(
        -0.5 * jnp.log(2.0 * jnp.pi * var) - 0.5 * ((x - mean) ** 2) / var
    )


def _uniform_logpdf(x: jnp.ndarray, low: float, high: float) -> jnp.ndarray:
    return jnp.sum(
        jnp.where(
            (x >= low) & (x <= high),
            -jnp.log(high - low),
            -jnp.inf,
        )
    )


def _tree_array_lognormprior(tree, mean: float, std: float) -> jnp.ndarray:
    arrays, _ = eqx.partition(tree, eqx.is_inexact_array)
    leaves = jax.tree_util.tree_leaves(arrays)

    total = jnp.array(0.0)
    for leaf in leaves:
        if leaf is not None:
            total = total + _normal_logpdf(leaf, mean, std)
    return total


def _tree_array_loguniformprior(tree, low: float, high: float) -> jnp.ndarray:
    arrays, _ = eqx.partition(tree, eqx.is_inexact_array)
    leaves = jax.tree_util.tree_leaves(arrays)

    total = jnp.array(0.0)
    for leaf in leaves:
        if leaf is not None:
            total = total + _uniform_logpdf(leaf, low, high)
    return total


def _model_logprior(
    model, layer1_bias_std=1.0, layer1_weight_std=1.0
) -> jnp.ndarray:
    mech = model.parameters
    layer1_bias = model.model.nns["net1"].layers["layer1"].bias
    layer1_weight = model.model.nns["net1"].layers["layer1"].weight
    rest = eqx.tree_at(
        lambda m: m["net1"].layers["layer1"], model.model.nns, replace=None
    )

    return (
        _tree_array_loguniformprior(mech, low=0.0, high=15.0)
        + _tree_array_lognormprior(layer1_bias, mean=0.0, std=layer1_bias_std)
        + _tree_array_lognormprior(
            layer1_weight, mean=0.0, std=layer1_weight_std
        )
        + _tree_array_lognormprior(rest, mean=0.0, std=1.0)
    )
