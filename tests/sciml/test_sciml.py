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
import petab.v1 as petab
import pytest
from amici.importers.petab.v1 import import_petab_problem
from amici.jax import (
    generate_equinox,
    petab_simulate,
    run_simulations,
)
from petab_sciml import NNModelStandard
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
net_cases_dir = cases_dir / "net_import"
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
    "test", sorted(d.stem for d in net_cases_dir.glob("[0-9]*"))
)
def test_net(test):
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
        "002",
        "009",
        "018",
        "019",
        "020",
        "021",
        "022",
        "042",
        "043",
        "044",
        "045",
        "046",
        "047",
        "048",
    ):
        with pytest.raises(NotImplementedError):
            generate_equinox(ml_model, module_dir)
        return
    generate_equinox(ml_model, module_dir)
    nets[ml_model.nn_model_id] = amici._module_from_path(
        ml_model.nn_model_id, module_dir
    ).net

    for input_file, par_file, output_file in zip(
        solutions["net_input"],
        solutions.get("net_ps", solutions["net_input"]),
        solutions["net_output"],
    ):
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
def test_ude(test):
    test_dir = ude_cases_dir / test
    with open(test_dir / "petab" / "problem.yaml") as f:
        petab_yaml = safe_load(f)
    with open(test_dir / "solutions.yaml") as f:
        solutions = safe_load(f)

    with change_directory(test_dir / "petab"):
        from petab.v2 import Problem

        petab_yaml["format_version"] = "2.0.0"  # TODO: fixme
        petab_problem = Problem.from_yaml(petab_yaml)
        jax_problem = import_petab_problem(
            petab_problem,
            output_dir=Path(__file__).parent / "models" / test,
            compile_=True,
            jax=True,
        )

    # llh
    llh, r = run_simulations(jax_problem)
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
    sort_by = [petab.OBSERVABLE_ID, petab.TIME, petab.SIMULATION_CONDITION_ID]
    actual = petab_simulate(jax_problem).sort_values(by=sort_by)
    expected = simulations.sort_values(by=sort_by)
    np.testing.assert_allclose(
        actual[petab.SIMULATION].values,
        expected[petab.SIMULATION].values,
        atol=solutions["tol_simulations"],
        rtol=solutions["tol_simulations"],
    )

    # gradient
    sllh, _ = eqx.filter_grad(run_simulations, has_aux=True)(
        jax_problem,
        solver=diffrax.Kvaerno5(),
        controller=diffrax.PIDController(atol=1e-14, rtol=1e-14),
        max_steps=2**16,
    )
    for component, file in solutions["grad_files"].items():
        actual_dict = {}
        if component == "mech":
            expected = pd.read_csv(test_dir / file, sep="\t").set_index(
                petab.PARAMETER_ID
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
                atol=solutions["tol_grad_llh"],
                rtol=solutions["tol_grad_llh"],
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
                            atol=solutions["tol_grad_llh"],
                            rtol=solutions["tol_grad_llh"],
                        )
