import logging
from functools import partial

import numpy as np
import pytest
import jax
import jax.numpy as jnp
import equinox as eqx
from beartype import beartype

import amici
from amici.jax.petab import run_simulations, JAXProblem
from amici.petab.petab_import import import_petab_problem
from amici.petab.simulations import simulate_petab, LLH, SLLH
from test_petab_benchmark import (
    benchmark_outdir,
    problems_for_gradient_check,
    settings,
)

jax.config.update("jax_enable_x64", True)


@pytest.mark.filterwarnings(
    "ignore:The following problem parameters were not used *",
    "ignore: The environment variable *",
    "ignore:Adjoint sensitivity analysis for models with discontinuous ",
)
def test_jax_llh(benchmark_problem):
    problem_id, flat_petab_problem, petab_problem, amici_model = (
        benchmark_problem
    )
    if problem_id == "Smith_BMCSystBiol2013":
        pytest.skip(
            "Skipping Smith_BMCSystBiol2013 due to non-supported events in JAX."
        )

    amici_solver = amici_model.getSolver()
    cur_settings = settings[problem_id]
    amici_solver.setAbsoluteTolerance(1e-8)
    amici_solver.setRelativeTolerance(1e-8)
    amici_solver.setMaxSteps(10_000)

    simulate_amici = partial(
        simulate_petab,
        petab_problem=flat_petab_problem,
        amici_model=amici_model,
        solver=amici_solver,
        scaled_parameters=True,
        scaled_gradients=True,
        log_level=logging.DEBUG,
    )

    np.random.seed(cur_settings.rng_seed)

    problem_parameters = None
    if problem_id in problems_for_gradient_check:
        point = flat_petab_problem.x_nominal_free_scaled
        for _ in range(20):
            amici_solver.setSensitivityMethod(amici.SensitivityMethod.adjoint)
            amici_solver.setSensitivityOrder(amici.SensitivityOrder.first)
            amici_model.setSteadyStateSensitivityMode(
                cur_settings.ss_sensitivity_mode
            )
            point_noise = (
                np.random.randn(len(point)) * cur_settings.noise_level
            )
            point += point_noise  # avoid small gradients at nominal value

            problem_parameters = dict(
                zip(flat_petab_problem.x_free_ids, point)
            )

            r_amici = simulate_amici(
                problem_parameters=problem_parameters,
            )
            if np.isfinite(r_amici[LLH]):
                break
        else:
            raise RuntimeError("Could not compute expected derivative.")
    else:
        r_amici = simulate_amici()
    llh_amici = r_amici[LLH]

    jax_model = import_petab_problem(
        petab_problem,
        model_output_dir=benchmark_outdir / (problem_id + "_jax"),
        jax=True,
    )
    jax_problem = JAXProblem(jax_model, petab_problem)
    if problem_parameters:
        jax_problem = eqx.tree_at(
            lambda x: x.parameters,
            jax_problem,
            jnp.array(
                [problem_parameters[pid] for pid in jax_problem.parameter_ids]
            ),
        )

    if problem_id in problems_for_gradient_check:
        beartype(run_simulations)(jax_problem)
        (llh_jax, _), sllh_jax = eqx.filter_value_and_grad(
            run_simulations, has_aux=True
        )(jax_problem)
    else:
        llh_jax, _ = beartype(run_simulations)(jax_problem)

    np.testing.assert_allclose(
        llh_jax,
        llh_amici,
        rtol=1e-3,
        atol=1e-3,
        err_msg=f"LLH mismatch for {problem_id}",
    )

    if problem_id in problems_for_gradient_check:
        sllh_amici = r_amici[SLLH]
        np.testing.assert_allclose(
            sllh_jax.parameters,
            np.array([sllh_amici[pid] for pid in jax_problem.parameter_ids]),
            rtol=1e-2,
            atol=1e-2,
            err_msg=f"SLLH mismatch for {problem_id}, {dict(zip(jax_problem.parameter_ids, sllh_jax.parameters))}",
        )
