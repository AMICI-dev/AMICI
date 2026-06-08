import os
from contextlib import contextmanager
from pathlib import Path

import amici
import diffrax
import equinox as eqx
import h5py
import jax
import jax.numpy as jnp
import jax.random as jr
import numpy as np
import pandas as pd
import pytest
from amici.exporters.jax import generate_equinox
from amici.importers.petab import *
from amici.sim.jax import petab_simulate, run_simulations
from petab import v2
from petab_sciml import NNModelStandard
from sciml_helpers import _v2_sciml_problem_helper
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


# pip install git+https://github.com/sebapersson/petab_sciml@add_standard#egg=petab_sciml\&subdirectory=src/python

cases_dir = Path(__file__).parent / "testsuite" / "test_cases"
net_cases_dir = cases_dir / "ml_model_import"
ude_cases_dir = cases_dir / "sciml_problem_import"
initialization_cases_dir = cases_dir / "initialization"


def _reshape_flat_array(array_flat):
    array_flat["ix"] = array_flat["ix"].astype(str)
    ix_cols = [
        f"ix_{i}" for i in range(len(array_flat["ix"].values[0].split(";")))
    ]
    if len(ix_cols) == 1:
        array_flat[ix_cols[0]] = array_flat["ix"].apply(int)
    else:
        array_flat[ix_cols] = pd.DataFrame(
            array_flat["ix"].str.split(";").apply(np.array).to_list(),
            index=array_flat.index,
        ).astype(int)
    array_flat.sort_values(by=ix_cols, inplace=True)
    array_shape = tuple(array_flat[ix_cols].max().astype(int) + 1)
    array = np.array(array_flat["value"].values).reshape(array_shape)
    return array


@pytest.mark.parametrize(
    "test", sorted([d.stem for d in initialization_cases_dir.glob("[0-9]*")])
)
def test_initialization(test):
    test_dir = initialization_cases_dir / test

    with open(test_dir / "solutions.yaml") as f:
        solutions = safe_load(f)

    with open(test_dir / "petab" / "problem.yaml") as f:
        petab_yaml = safe_load(f)

    with change_directory(test_dir / "petab"):
        petab_problem = _v2_sciml_problem_helper(
            petab_yaml, test_dir / "petab"
        )

        pi = PetabImporter(
            petab_problem=petab_problem,
            module_name="init" + test,
            compile_=True,
            jax=jax,
            validate=False,
        )

        jax_problem = pi.create_simulator(force_import=True)

    tol = solutions["tol"]

    for net_id, ref_file in solutions["parameter_files"].items():
        nn = jax_problem.model.nns[net_id]
        with h5py.File(test_dir / ref_file, "r") as expected:
            for layer_name in nn.layers.keys():
                for attribute in ("weight", "bias"):
                    actual = getattr(nn.layers[layer_name], attribute)
                    np.testing.assert_allclose(
                        np.squeeze(actual),
                        np.squeeze(
                            expected["parameters"][net_id][layer_name][
                                attribute
                            ][:]
                        ),
                        atol=tol,
                        rtol=tol,
                    )


@pytest.mark.parametrize(
    "test", sorted(d.stem for d in net_cases_dir.glob("[0-9]*"))
)
def test_ml_model_import(test):
    test_dir = net_cases_dir / test
    with open(test_dir / "solutions.yaml") as f:
        solutions = safe_load(f)

    if test.endswith("_alt"):
        net_file = cases_dir / test.replace("_alt", "") / solutions["net_file"]
    else:
        net_file = test_dir / solutions["net_file"]
    ml_model = NNModelStandard.load_data(net_file)

    nets = {}
    outdir = Path(__file__).parent / "models" / test
    module_dir = outdir / f"{ml_model.nn_model_id}.py"
    if test in (
        "009",
        "018",
    ):
        with pytest.raises(NotImplementedError):
            generate_equinox(ml_model, module_dir)
        return
    generate_equinox(ml_model, module_dir)
    nets[ml_model.nn_model_id] = amici._module_from_path(
        ml_model.nn_model_id, module_dir
    ).net

    if test in ("052", "053"):
        input_files = [
            (i1, i2)
            for i1, i2 in zip(
                solutions["net_input_arg0"], solutions["net_input_arg1"]
            )
        ]
    else:
        input_files = solutions["net_input"]

    for input_file, par_file, output_file in zip(
        input_files,
        solutions.get("net_ps", input_files),
        solutions["net_output"],
    ):
        if test in ("052", "053"):
            input = tuple(
                [
                    h5py.File(test_dir / in_file, "r")["inputs"]["input0"][
                        "data"
                    ][:]
                    for in_file in input_file
                ]
            )
        else:
            input = h5py.File(test_dir / input_file, "r")["inputs"]["input0"][
                "data"
            ][:]
        output = h5py.File(test_dir / output_file, "r")["outputs"]["output0"][
            "data"
        ][:]

        if "net_ps" in solutions:
            par = h5py.File(test_dir / par_file, "r")
            net = nets[ml_model.nn_model_id](jr.PRNGKey(0))
            for layer in net.layers.keys():
                if (
                    isinstance(net.layers[layer], eqx.Module)
                    and hasattr(net.layers[layer], "weight")
                    and net.layers[layer].weight is not None
                ):
                    w = par["parameters"][ml_model.nn_model_id][layer][
                        "weight"
                    ][:]
                    if isinstance(net.layers[layer], eqx.nn.ConvTranspose):
                        # see FAQ in https://docs.kidger.site/equinox/api/nn/conv/#equinox.nn.ConvTranspose
                        w = np.flip(w, axis=tuple(range(2, w.ndim))).swapaxes(
                            0, 1
                        )
                    assert w.shape == net.layers[layer].weight.shape
                    net = eqx.tree_at(
                        lambda x: x.layers[layer].weight,
                        net,
                        jnp.array(w),
                    )
                if (
                    isinstance(net.layers[layer], eqx.Module)
                    and hasattr(net.layers[layer], "bias")
                    and net.layers[layer].bias is not None
                ):
                    b = par["parameters"][ml_model.nn_model_id][layer]["bias"][
                        :
                    ]
                    if isinstance(
                        net.layers[layer],
                        eqx.nn.Conv | eqx.nn.ConvTranspose,
                    ):
                        b = np.expand_dims(
                            b,
                            tuple(
                                range(
                                    1,
                                    net.layers[layer].num_spatial_dims + 1,
                                )
                            ),
                        )
                    assert b.shape == net.layers[layer].bias.shape
                    net = eqx.tree_at(
                        lambda x: x.layers[layer].bias,
                        net,
                        jnp.array(b),
                    )

            if test in ("019", "020", "021", "022"):
                iters = solutions.get("dropout", 40000)
                res = net.forward(
                    input, key=jr.PRNGKey(np.random.randint(0, 2**31))
                )
                for _ in range(iters):
                    res += net.forward(
                        input, key=jr.PRNGKey(np.random.randint(0, 2**31))
                    )
                np.testing.assert_allclose(
                    res / iters,
                    output,
                    atol=2e-2,
                    rtol=2e-2,
                )
            else:
                net = eqx.nn.inference_mode(net)

                if test == "net_004_alt":
                    return  # skipping, no support for non-cross-correlation in equinox

                np.testing.assert_allclose(
                    net.forward(input),
                    output,
                    atol=1e-3,
                    rtol=1e-3,
                )


@pytest.mark.parametrize(
    "test", sorted([d.stem for d in ude_cases_dir.glob("[0-9]*")])
)
def test_sciml_problem_import(test):
    test_dir = ude_cases_dir / test

    with open(test_dir / "petab" / "problem.yaml") as f:
        petab_yaml = safe_load(f)
    with open(test_dir / "solutions.yaml") as f:
        solutions = safe_load(f)

    with change_directory(test_dir / "petab"):
        # HACK!! Again!! Around "array" in parameters table
        petab_problem = _v2_sciml_problem_helper(
            petab_yaml, test_dir / "petab"
        )

        if test in ("003",):
            with pytest.raises(NotImplementedError):
                pi = PetabImporter(
                    petab_problem=petab_problem,
                    module_name="hybrid" + test,
                    compile_=True,
                    jax=jax,
                    validate=False,  # And again...around "array" in parameters table
                )
            return

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
    if test in (
        "032",
        "033",
        "034",
    ):
        configs = {
            "032": {},
            "033": {"layer1_weight_std": 2.0, "layer1_bias_std": 2.0},
            "034": {"layer1_weight_std": 2.0},
        }
        logposterior = llh + _model_logprior(jax_problem, **configs[test])
        np.testing.assert_allclose(
            logposterior,
            solutions["log_posterior"],
            atol=solutions["tol_log_posterior"],
            rtol=solutions["tol_log_posterior"],
        )
    else:
        np.testing.assert_allclose(
            llh,
            solutions["llh"],
            atol=solutions["tol_llh"],
            rtol=solutions["tol_llh"],
        )
    simulations = pd.concat(
        [
            pd.read_csv(test_dir / simulation, sep="\t")
            for simulation in solutions["simulation_files"]
        ]
    )

    # simulations
    sort_by = [v2.C.OBSERVABLE_ID, v2.C.TIME, v2.C.EXPERIMENT_ID]
    actual = petab_simulate(jax_problem).sort_values(by=sort_by)
    expected = simulations.sort_values(by=sort_by)
    np.testing.assert_allclose(
        actual[v2.C.SIMULATION].values,
        expected[v2.C.SIMULATION].values,
        atol=solutions["tol_simulations"],
        rtol=solutions["tol_simulations"],
    )

    # gradient
    sllh, aux = eqx.filter_grad(run_simulations, has_aux=True)(
        jax_problem,
        solver=diffrax.Kvaerno5(),
        controller=diffrax.PIDController(atol=1e-14, rtol=1e-14),
        max_steps=2**16,
    )
    for component, file in solutions["grad_files"].items():
        actual_dict = {}
        if component == "mech":
            expected = pd.read_csv(test_dir / file, sep="\t").set_index(
                v2.C.PARAMETER_ID
            )

            for ip in expected.index:
                if ip in jax_problem.parameter_ids:
                    actual_dict[ip] = sllh.parameters[
                        jax_problem.parameter_ids.index(ip)
                    ].item()
            actual = pd.Series(actual_dict).loc[expected.index].values
            np.testing.assert_allclose(
                actual,
                expected["value"].values,
                atol=solutions["tol_grad"],
                rtol=solutions["tol_grad"],
            )
        else:
            expected = h5py.File(test_dir / file, "r")
            for layer_name, layer in jax_problem.model.nns[
                component
            ].layers.items():
                for attribute in dir(layer):
                    if not isinstance(
                        getattr(layer, attribute), jax.numpy.ndarray
                    ):
                        continue
                    actual = getattr(
                        sllh.model.nns[component].layers[layer_name], attribute
                    )
                    if (
                        isinstance(layer, eqx.nn.ConvTranspose)
                        and attribute == "weight"
                    ):
                        actual = np.flip(
                            actual.swapaxes(0, 1),
                            axis=tuple(range(2, actual.ndim)),
                        )
                    if (
                        np.squeeze(
                            expected["parameters"][component][layer_name][
                                attribute
                            ][:]
                        ).size
                        == 0
                    ):
                        assert np.all(actual == 0.0)
                    else:
                        np.testing.assert_allclose(
                            np.squeeze(actual),
                            np.squeeze(
                                expected["parameters"][component][layer_name][
                                    attribute
                                ][:]
                            ),
                            atol=solutions["tol_grad"],
                            rtol=solutions["tol_grad"],
                        )


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
