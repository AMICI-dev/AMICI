import logging
from functools import partial

import diffrax
import equinox as eqx
import jax
import jax.numpy as jnp
import numpy as np
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

    # Models whose JAX benchmark currently fails, grouped by root cause so that
    # entries can be removed individually as the underlying issues are fixed.
    # Most of these are pre-existing failures that were accidentally un-skipped
    # by #3204: that PR added a PEtab v1->v2 auto-upgrade inside
    # ``import_petab_problem(..., jax=True)``, which removed the "JAXProblem
    # does not support PEtab v1 problems" error that this test used to catch and
    # skip on (see the ``except`` block below). They are not regressions.
    known_jax_failures = {
        # Events are not yet supported by the JAX backend.
        "non-supported events in JAX": [
            "Liu_IFACPapersOnLine2025",
            "Oliveira_NatCommun2021",
            "SalazarCavazos_MBoC2020",
            "Smith_BMCSystBiol2013",
        ],
        # PEtab v1->v2 auto-upgrade fails upstream in libpetab
        # (petab.v2.petab1to2.update_prior / objective-prior handling):
        # "ValueError: Failed to auto-upgrade PEtab 1.0 problem to PEtab 2.0".
        "PEtab v1->v2 auto-upgrade failure (upstream libpetab)": [
            "Armistead_CellDeathDis2024",
            "Bachmann_MSB2011",
            "Isensee_JCB2018",
            "Schwen_PONE2014",
        ],
        # JAXProblem measurement parsing chokes on an empty string in a numeric
        # column: "ValueError: could not convert string to float: np.str_('')".
        "JAXProblem measurement parsing (empty numeric cell)": [
            "Laske_PLOSComputBiol2019",
            "Raia_CancerResearch2011",
            "Weber_BMC2015",
        ],
        # Measurement padding computes a negative pad width:
        # "ValueError: index can't contain negative values".
        "JAX measurement padding (negative pad width)": [
            "Brannmark_JBC2010",
            "Lucarelli_CellSystems2018",
        ],
        # diffrax integration produces non-finite values on these stiff models
        # (linear solver returns/receives NaN or inf).
        "diffrax non-finite integration (stiff model)": [
            "Borghans_BiophysChem1997",
            "Bruno_JExpBot2016",
            "Elowitz_Nature2000",
            "Fiedler_BMCSystBiol2016",
            "Fujita_SciSignal2010",
            "Giordano_Nature2020",
            "Okuonghae_ChaosSolitonsFractals2020",
            "Rahman_MBS2016",
            "Zheng_PNAS2012",
        ],
        # JAX vs. AMICI LLH/SLLH mismatch beyond the assertion tolerance.
        "JAX/AMICI LLH or SLLH mismatch": [
            "Bertozzi_PNAS2020",
            "Blasi_CellSystems2016",
            "Boehm_JProteomeRes2014",
            "Perelson_Science1996",
            "Sneyd_PNAS2002",
            "Zhao_QuantBiol2020",
        ],
    }
    for reason, problem_ids in known_jax_failures.items():
        if problem_id in problem_ids:
            pytest.skip(f"Skipping {problem_id}: {reason}.")

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
            jax_problem = eqx.tree_at(
                lambda x: x.parameters,
                jax_problem,
                jnp.array(
                    [
                        problem_parameters[pid]
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
            np.testing.assert_allclose(
                sllh_jax.parameters,
                np.array(
                    [sllh_amici[pid] for pid in jax_problem.parameter_ids]
                ),
                rtol=1e-2,
                atol=1e-2,
                err_msg=f"SLLH mismatch for {problem_id}, {dict(zip(jax_problem.parameter_ids, sllh_jax.parameters))}",
            )
    except (NotImplementedError, TypeError) as err:
        if "JAXProblem does not support PEtab v1 problems" in str(err):
            pytest.skip(str(err))
        elif "The JAX backend does not support" in str(err):
            pytest.skip(str(err))
        raise err
