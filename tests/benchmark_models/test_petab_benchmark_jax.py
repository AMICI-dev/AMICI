import logging
from functools import partial

import diffrax
import equinox as eqx
import jax
import jax.numpy as jnp
import numpy as np
import petab.v1 as petab
import pytest
from amici.importers.petab.v1 import (
    import_petab_problem,
)
from amici.sim.jax.petab import run_simulations
from amici.sim.sundials import SensitivityMethod, SensitivityOrder
from amici.sim.sundials.petab.v1 import (
    LLH,
    SLLH,
    simulate_petab,
)
from beartype import beartype
from petab.v1.parameters import unscale
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
    # v1->v2 auto-upgrade rewrites unsupported log10-normal priors/noise to
    # log-normal and drops initialisation priors; both are benign for the data
    # likelihood checked here (priors/init are not part of the objective).
    "ignore:.*Using `log-normal` instead.*:UserWarning",
    "ignore:.*Initialisation priors in parameter table are not supported.*:UserWarning",
)
def test_jax_llh(benchmark_problem):
    problem_id, flat_petab_problem, petab_problem, amici_model = (
        benchmark_problem
    )

    # Models with events not yet supported by the JAX backend.
    events_not_supported = [
        "Liu_IFACPapersOnLine2025",
        "Oliveira_NatCommun2021",
        "SalazarCavazos_MBoC2020",
        "Smith_BMCSystBiol2013",
    ]
    if problem_id in events_not_supported:
        pytest.skip(
            f"Skipping {problem_id} due to non-supported events in JAX."
        )

    # Skipped only to stay within the CI time budget, not for correctness:
    # Weber uses max_steps=4e7 under reverse-mode AD, which takes hours (it
    # dominated a full-suite run). The preequilibration + parameter-dependent-
    # event handling it exercises is also covered by Brannmark_JBC2010 and
    # Isensee_JCB2018, which stay enabled.
    slow_models = [
        "Weber_BMC2015",
    ]
    if problem_id in slow_models:
        pytest.skip(f"Skipping {problem_id}: excessive runtime for CI.")

    # PEtab v2 has no log10-normal noise: the v1->v2 upgrade rewrites
    # log10-normal to log-normal (a genuinely different likelihood, since the
    # noise parameter is not rescaled by ln(10)), so the JAX (v2) result cannot
    # match the log10-normal SUNDIALS reference. Mirrors the SUNDIALS v2 test
    # (test_petab_benchmark.py::test_nominal_parameters_llh_v2).
    obs_df = flat_petab_problem.observable_df
    if (
        petab.C.OBSERVABLE_TRANSFORMATION in obs_df
        and petab.C.LOG10 in obs_df[petab.C.OBSERVABLE_TRANSFORMATION].unique()
    ):
        pytest.skip(
            f"Skipping {problem_id}: log10-normal noise is not supported in "
            "PEtab v2."
        )

    amici_solver = amici_model.create_solver()
    cur_settings = settings[problem_id]
    amici_solver.set_absolute_tolerance(1e-8)
    amici_solver.set_relative_tolerance(1e-8)
    amici_solver.set_max_steps(10_000)

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
            amici_solver.set_sensitivity_method(SensitivityMethod.forward)
            amici_solver.set_sensitivity_order(SensitivityOrder.first)
            amici_model.set_steady_state_sensitivity_mode(
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

    try:
        jax_problem = import_petab_problem(
            petab_problem,
            output_dir=benchmark_outdir / (problem_id + "_jax"),
            jax=True,
        )
        if problem_parameters:
            # ``problem_parameters`` holds scaled (e.g. log10) values, as fed to
            # the SUNDIALS reference with ``scaled_parameters=True``. The PEtab
            # v2 JAX backend operates on *linear* parameters, so unscale before
            # injecting to evaluate both backends at the same physical point.
            jax_problem = eqx.tree_at(
                lambda x: x.parameters,
                jax_problem,
                jnp.array(
                    [
                        unscale(
                            problem_parameters[pid],
                            flat_petab_problem.parameter_df.loc[
                                pid, petab.C.PARAMETER_SCALE
                            ],
                        )
                        for pid in jax_problem.parameter_ids
                    ]
                ),
            )

        if problem_id in problems_for_gradient_check:
            if problem_id == "Weber_BMC2015":
                atol = cur_settings.atol_sim
                rtol = cur_settings.rtol_sim
                max_steps = 4 * 10**7
            else:
                atol = 1e-8
                rtol = 1e-8
                max_steps = 2 * 10**5
            beartype(run_simulations)(jax_problem)
            (llh_jax, _), sllh_jax = eqx.filter_value_and_grad(
                run_simulations, has_aux=True
            )(
                jax_problem,
                max_steps=max_steps,
                controller=diffrax.PIDController(
                    atol=atol,
                    rtol=rtol,
                ),
            )
        else:
            llh_jax, _ = beartype(run_simulations)(
                jax_problem,
                max_steps=2 * 10**5,
            )

        np.testing.assert_allclose(
            llh_jax,
            llh_amici,
            rtol=1e-3,
            atol=1e-3,
            err_msg=f"LLH mismatch for {problem_id}",
        )

        if problem_id in problems_for_gradient_check:
            sllh_amici = r_amici[SLLH]
            # SUNDIALS returns the gradient w.r.t. the *scaled* parameters
            # (scaled_gradients=True), while JAX differentiates w.r.t. the
            # *linear* parameters. Convert the JAX gradient to the estimation
            # scale via the chain rule
            #   d(llh)/d(scaled) = d(llh)/d(linear) * d(linear)/d(scaled),
            # with d(linear)/d(scaled) = p*ln(10) for log10, p for log, 1 for lin.
            scales = [
                flat_petab_problem.parameter_df.loc[
                    pid, petab.C.PARAMETER_SCALE
                ]
                for pid in jax_problem.parameter_ids
            ]
            p_lin = np.asarray(jax_problem.parameters)
            chain = np.array(
                [
                    p * np.log(10.0)
                    if s == petab.C.LOG10
                    else p
                    if s == petab.C.LOG
                    else 1.0
                    for s, p in zip(scales, p_lin)
                ]
            )
            sllh_jax_scaled = np.asarray(sllh_jax.parameters) * chain
            np.testing.assert_allclose(
                sllh_jax_scaled,
                np.array(
                    [sllh_amici[pid] for pid in jax_problem.parameter_ids]
                ),
                rtol=1e-2,
                atol=1e-2,
                err_msg=f"SLLH mismatch for {problem_id}, {dict(zip(jax_problem.parameter_ids, sllh_jax_scaled))}",
            )
    except (NotImplementedError, TypeError) as err:
        if "JAXProblem does not support PEtab v1 problems" in str(err):
            pytest.skip(str(err))
        elif "The JAX backend does not support" in str(err):
            pytest.skip(str(err))
        raise err
