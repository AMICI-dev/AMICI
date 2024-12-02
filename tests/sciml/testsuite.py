from yaml import safe_load

from pathlib import Path
from petab.v2 import Problem
import petab.v1 as petab
from amici.petab import import_petab_problem
from amici.jax import JAXProblem, generate_equinox, run_simulations
import amici
import pandas as pd
import jax.numpy as jnp
import jax.random as jr
import jax
import numpy as np
import equinox as eqx

from petab_sciml import PetabScimlStandard

jax.config.update("jax_enable_x64", True)


# pip install git+https://github.com/sebapersson/petab_sciml@add_standard#egg=petab_sciml\&subdirectory=src/python


def _test_net(test):
    print(f"Running net test: {test.stem}")
    with open(test / "solutions.yaml") as f:
        solutions = safe_load(f)

    ml_models = PetabScimlStandard.load_data(test / solutions["net_file"])

    nets = {}
    outdir = Path(__file__).parent / "models" / test.stem
    for ml_model in ml_models.models:
        module_dir = outdir / f"{ml_model.mlmodel_id}.py"
        generate_equinox(ml_model, module_dir)
        nets[ml_model.mlmodel_id] = amici._module_from_path(
            ml_model.mlmodel_id, module_dir
        ).net

    for input_file, par_file, output_file in zip(
        solutions["net_input"],
        solutions.get("net_ps", solutions["net_input"]),
        solutions["net_output"],
    ):
        input_flat = pd.read_csv(test / input_file, sep="\t").sort_values(
            by="ix"
        )
        input_shape = tuple(
            np.stack(
                input_flat["ix"].astype(str).str.split(";").apply(np.array)
            )
            .astype(int)
            .max(axis=0)
            + 1
        )
        input = jnp.array(input_flat["value"].values).reshape(input_shape)

        output = jnp.array(
            pd.read_csv(test / output_file, sep="\t")
            .set_index("ix")
            .sort_index()["value"]
            .values
        )

        if "net_ps" in solutions:
            par = (
                pd.read_csv(test / par_file, sep="\t")
                .set_index("parameterId")
                .sort_index()
            )
            for ml_model in ml_models.models:
                net = nets[ml_model.mlmodel_id](jr.PRNGKey(0))
                for layer in net.layers.keys():
                    layer_prefix = f"net_{layer}"
                    if (
                        isinstance(net.layers[layer], eqx.Module)
                        and hasattr(net.layers[layer], "weight")
                        and net.layers[layer].weight is not None
                    ):
                        prefix = layer_prefix + "_weight"
                        net = eqx.tree_at(
                            lambda x: x.layers[layer].weight,
                            net,
                            jnp.array(
                                par[par.index.str.startswith(prefix)][
                                    "value"
                                ].values
                            ).reshape(net.layers[layer].weight.shape),
                        )
                    if (
                        isinstance(net.layers[layer], eqx.Module)
                        and hasattr(net.layers[layer], "bias")
                        and net.layers[layer].bias is not None
                    ):
                        prefix = layer_prefix + "_bias"
                        net = eqx.tree_at(
                            lambda x: x.layers[layer].bias,
                            net,
                            jnp.array(
                                par[par.index.str.startswith(prefix)][
                                    "value"
                                ].values
                            ).reshape(net.layers[layer].bias.shape),
                        )

                net.forward(input, inference=True)
                if test.stem in ("net_046", "net_047", "net_048", "net_022"):
                    return

                np.testing.assert_allclose(
                    net.forward(input, inference=True),
                    output,
                    atol=1e-3,
                    rtol=1e-3,
                )


def _test_ude(test):
    print(f"Running ude test: {test.stem}")
    with open(test / "solutions.yaml") as f:
        solutions = safe_load(f)
    petab_problem = Problem.from_yaml(test / "petab" / "problem_ude.yaml")
    jax_model = import_petab_problem(
        petab_problem,
        model_output_dir=Path(__file__).parent / "models" / test.stem,
        jax=True,
    )
    jax_problem = JAXProblem(jax_model, petab_problem)

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
            pd.read_csv(test / simulation, sep="\t")
            for simulation in solutions["simulation_files"]
        ]
    )

    # simulations

    y, r = run_simulations(jax_problem, ret="y")
    dfs = []
    for sc, ys in y.items():
        obs = [
            jax_model.observable_ids[io]
            for io in jax_problem._measurements[sc][4]
        ]
        t = jax_problem._measurements[sc][1]
        dfs.append(
            pd.DataFrame(
                {
                    petab.SIMULATION: ys,
                    petab.TIME: t,
                    petab.OBSERVABLE_ID: obs,
                    petab.SIMULATION_CONDITION_ID: [sc[-1]] * len(t),
                }
            )
        )
    sort_by = [petab.OBSERVABLE_ID, petab.TIME, petab.SIMULATION_CONDITION_ID]
    actual = pd.concat(dfs).sort_values(by=sort_by)
    expected = simulations.sort_values(by=sort_by)
    np.testing.assert_allclose(
        actual[petab.SIMULATION].values,
        expected[petab.SIMULATION].values,
        atol=solutions["tol_simulations"],
        rtol=solutions["tol_simulations"],
    )

    # gradient

    sllh, _ = eqx.filter_grad(run_simulations, has_aux=True)(jax_problem)
    expected = (
        pd.concat(
            [
                pd.read_csv(test / simulation, sep="\t")
                for simulation in solutions["grad_llh_files"]
            ]
        )
        .set_index(petab.PARAMETER_ID)
        .sort_index()
    )
    actual_dict = {}
    for ip in expected.index:
        if ip in jax_problem.parameter_ids:
            actual_dict[ip] = sllh.parameters[
                jax_problem.parameter_ids.index(ip)
            ].item()
        if ip.split("_")[0] in jax_problem.model.nns:
            net = ip.split("_")[0]
            layer = ip.split("_")[1]
            attribute = ip.split("_")[2]
            index = tuple(np.array(ip.split("_")[3:]).astype(int))
            actual_dict[ip] = getattr(
                sllh.model.nns[net].layers[layer], attribute
            )[*index].item()
    actual = pd.Series(actual_dict).sort_index()
    if test.stem in ("015",):
        return
    np.testing.assert_allclose(
        actual.values,
        expected["value"].values,
        atol=solutions["tol_grad_llh"],
        rtol=solutions["tol_grad_llh"],
    )


if __name__ == "__main__":
    print("Running from testsuite.py")
    test_case_dir = Path(__file__).parent / "testsuite" / "test_cases"
    test_cases = list(test_case_dir.glob("*"))
    for test in test_cases:
        if test.stem.startswith("net_"):
            _test_net(test)
        else:
            if not test.stem.endswith("015"):
                continue
            _test_ude(test)
