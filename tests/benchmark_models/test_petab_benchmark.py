"""Tests based on the PEtab benchmark problems.

Tests simulate_petab, correctness of the log-likelihood computation at nominal
parameters, correctness of the gradient computation, and simulation times
for a subset of the benchmark problems.
"""

import contextlib
import logging
import os
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import benchmark_models_petab
import fiddy
import numpy as np
import pandas as pd
import petab.v1 as petab
import pytest
import yaml
from amici import get_model_root_dir
from amici.adapters.fiddy import (
    simulate_petab_to_cached_functions,
    simulate_petab_v2_to_cached_functions,
)
from amici.importers.petab.v1 import (
    import_petab_problem,
)
from amici.logging import get_logger
from amici.sim.sundials import (
    AMICI_SUCCESS,
    SensitivityMethod,
    SensitivityOrder,
    SteadyStateComputationMode,
    SteadyStateSensitivityMode,
)
from amici.sim.sundials.petab.v1 import (
    LLH,
    RDATAS,
    rdatas_to_measurement_df,
    simulate_petab,
)
from fiddy import MethodId, get_derivative
from fiddy.derivative_check import NumpyIsCloseDerivativeCheck
from fiddy.success import Consistency
from petab.v1.lint import measurement_table_has_timepoint_specific_mappings
from petab.v1.visualize import plot_problem

# Enable various debug output
debug = False

logger = get_logger(
    f"amici.{__name__}", logging.DEBUG if debug else logging.INFO
)

script_dir = Path(__file__).parent.absolute()
benchmark_outdir = get_model_root_dir() / "test_bmc"
debug_path = script_dir / "debug"
if debug:
    debug_path.mkdir(exist_ok=True, parents=True)

# reference values for simulation times and log-likelihoods
references_yaml = script_dir / "benchmark_models.yaml"
with open(references_yaml) as f:
    reference_values = yaml.full_load(f)

# problem IDs for which to check the gradient
# TODO: extend
problems_for_gradient_check = set(benchmark_models_petab.MODELS) - {
    # excluded due to excessive runtime
    "Bachmann_MSB2011",
    "Chen_MSB2009",
    "Froehlich_CellSystems2018",
    "Raimundez_PCB2020",
    "Lucarelli_CellSystems2018",
    "Isensee_JCB2018",
    "Beer_MolBioSystems2014",
    "Alkan_SciSignal2018",
    "Lang_PLOSComputBiol2024",
    "Smith_BMCSystBiol2013",
    # excluded due to excessive numerical failures
    "Crauste_CellSystems2017",
}
problems_for_gradient_check = list(sorted(problems_for_gradient_check))


# Problems for checking the log-likelihood computation at nominal parameters
#
# not merged:
# Becker_Science2010 (multiple models)             https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/tree/model_Becker_Science2010
# Hass_PONE2017 (???)                              https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/tree/model_Hass_PONE2017
# Korkut_eLIFE2015 (???)                           https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/tree/model_Korkut_eLIFE2015
# Casaletto_PNAS2019 (yaml missing)                https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/tree/model_Casaletto_PNAS2019
# Merkle_PCB2016 (model missing)                   https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/tree/model_Merkle_PCB2016
# Parmar_PCB2019 (SBML extensions)                 https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/tree/model_Parmar_PCB2019
# Swameye_PNAS2003 (splines)                       https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/tree/model_Swameye_PNAS2003
# Sobotta_Frontiers2017 (???)                      https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/tree/model_Sobotta_Frontiers2017
# Raia_CancerResearch2011 (state dependent sigmas) https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/tree/model_Raia_CancerResearch2011
#
# no reference value:
# Alkan_SciSignal2018 (d2d: Alkan_DDRP_SciSignal2018)
# Bertozzi_PNAS2020 (gh: vanako, doi: https://doi.org/10.1073/pnas.2006520117, code/data: https://github.com/gomohler/pnas2020 (R))
# Blasi_CellSystems2016 (gh: Leonard Schmiester, doi: https://doi.org/10.1016/j.cels.2016.01.002, code/data: not available)
# Giordano_Nature2020 (gh: Paul Jonas Jost, doi: https://doi.org/10.1038/s41591-020-0883-7, code/data: http://users.dimi.uniud.it/~giulia.giordano/docs/papers/SIDARTHEcode.zip (MATLAB))
# Laske_PLOSComputBiol2019 (gh: Clemens Peiter, doi: https://doi.org/10.1128/JVI.00080-12 (?), code/data: ???)
# Okuonghae_ChaosSolitonsFractals2020 (gh: Paul Jonas Jost, doi: https://doi.org/10.1016/j.chaos.2020.110032, code/data: ???)
# Oliveira_NatCommun2021 (gh: lorenapohl, doi: https://doi.org/10.1038/s41467-020-19798-3, code: https://github.com/cidacslab/Mathematical-and-Statistical-Modeling-of-COVID19-in-Brazil (python) data: https://infovis.sei.ba.gov.br/covid19/ )
# Perelson_Science1996 (gh: Philipp Staedter, doi: https://doi.org/10.1126/science.271.5255.1582, code/data: ???)
# Rahman_MBS2016 (gh: Yannik Schaelte, doi: https://doi.org/10.1016/j.mbs.2016.07.009, code: not available, data: table in paper ...)
# Raimundez_PCB2020 (gh: Elba Raimundez, doi: https://doi.org/10.1371/journal.pcbi.1007147, code/data: https://zenodo.org/record/2908234#.Y5hUUS8w3yw (d2d))
# SalazarCavazos_MBoC2020 (gh: Dilan Pathirana, doi: https://doi.org/10.1091/mbc.E19-09-0548, code/data: supplement (BNGL))
# Zhao_QuantBiol2020 (gh: Iva Ewert, doi: https://doi.org/10.1007/s40484-020-0199-0, code: not available, data: table in supp)
#
# covered by performance test:
# Froehlich_CellSystems2018
#
# Unknown reasons:
# Chen_MSB2009
#
# Confirmed to be working:
# TODO: extend
problems_for_llh_check = [
    "Bachmann_MSB2011",
    "Beer_MolBioSystems2014",
    "Boehm_JProteomeRes2014",
    "Borghans_BiophysChem1997",
    "Brannmark_JBC2010",
    "Bruno_JExpBot2016",
    "Crauste_CellSystems2017",
    "Elowitz_Nature2000",
    "Fiedler_BMCSystBiol2016",
    "Fujita_SciSignal2010",
    # Excluded until https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/pull/253
    #  is sorted out
    # "Isensee_JCB2018",
    "Lucarelli_CellSystems2018",
    "Schwen_PONE2014",
    "Smith_BMCSystBiol2013",
    "Sneyd_PNAS2002",
    "Weber_BMC2015",
    "Zheng_PNAS2012",
]

# all PEtab problems we need to import
problems = list(
    sorted(set(problems_for_gradient_check + problems_for_llh_check))
)


@dataclass
class GradientCheckSettings:
    """Problem-specific settings for gradient checks."""

    # Absolute and relative tolerances for simulation
    atol_sim: float = 1e-16
    rtol_sim: float = 1e-12
    # Absolute and relative tolerances for finite difference gradient checks.
    atol_check: float = 1e-3
    rtol_check: float = 1e-2
    # Absolute and relative tolerances for fiddy consistency check between
    # forward/backward/central differences.
    atol_consistency: float = 1e-5
    rtol_consistency: float = 1e-1
    # Step sizes for finite difference gradient checks.
    step_sizes: list[float] = field(
        default_factory=lambda: [
            2e-1,
            1e-1,
            5e-2,
            1e-2,
            5e-1,
            1e-3,
            1e-4,
            1e-5,
        ]
    )
    rng_seed: int = 0
    ss_computation_mode: SteadyStateComputationMode = (
        SteadyStateComputationMode.integrationOnly
    )
    ss_sensitivity_mode: SteadyStateSensitivityMode = (
        SteadyStateSensitivityMode.integrateIfNewtonFails
    )
    noise_level: float = 0.05


settings = defaultdict(GradientCheckSettings)
# NOTE: Newton method fails badly with ASA for Blasi_CellSystems2016
settings["Blasi_CellSystems2016"] = GradientCheckSettings(
    atol_check=1e-12,
    rtol_check=1e-4,
    ss_sensitivity_mode=SteadyStateSensitivityMode.integrationOnly,
)
settings["Borghans_BiophysChem1997"] = GradientCheckSettings(
    rng_seed=2,
    atol_check=1e-5,
    rtol_check=1e-3,
)
settings["Brannmark_JBC2010"] = GradientCheckSettings(
    ss_sensitivity_mode=SteadyStateSensitivityMode.integrationOnly,
)
settings["Fujita_SciSignal2010"] = GradientCheckSettings(
    atol_check=1e-7,
    rtol_check=5e-4,
)
settings["Giordano_Nature2020"] = GradientCheckSettings(
    atol_check=1e-6, rtol_check=1e-3, rng_seed=1
)
settings["Okuonghae_ChaosSolitonsFractals2020"] = GradientCheckSettings(
    atol_sim=1e-14,
    rtol_sim=1e-14,
    noise_level=0.01,
    atol_consistency=1e-3,
)
settings["Oliveira_NatCommun2021"] = GradientCheckSettings(
    # Avoid "root after reinitialization"
    atol_sim=1e-12,
    rtol_sim=1e-12,
)
settings["Raia_CancerResearch2011"] = GradientCheckSettings(
    atol_check=1e-10,
    rtol_check=1e-3,
)
settings["Smith_BMCSystBiol2013"] = GradientCheckSettings(
    atol_sim=1e-10,
    rtol_sim=1e-10,
)
settings["Sneyd_PNAS2002"] = GradientCheckSettings(
    atol_sim=1e-15,
    rtol_sim=1e-12,
    atol_check=1e-5,
    rtol_check=1e-4,
    rng_seed=7,
)
settings["Weber_BMC2015"] = GradientCheckSettings(
    atol_sim=1e-12,
    rtol_sim=1e-12,
    atol_check=1e-6,
    rtol_check=1e-2,
    rng_seed=1,
)
settings["Zheng_PNAS2012"] = GradientCheckSettings(
    rng_seed=1,
    rtol_sim=1e-15,
    atol_check=5e-4,
    rtol_check=4e-3,
    noise_level=0.01,
    ss_sensitivity_mode=SteadyStateSensitivityMode.integrationOnly,
    step_sizes=[
        3e-1,
        2e-1,
        1e-1,
        5e-2,
        1e-2,
        5e-1,
        1e-3,
        1e-4,
        1e-5,
    ],
)


@pytest.mark.filterwarnings(
    "ignore:divide by zero encountered in log",
    # https://github.com/AMICI-dev/AMICI/issues/18
    "ignore:Adjoint sensitivity analysis for models with discontinuous "
    "right hand sides .*:UserWarning",
)
def test_nominal_parameters_llh(benchmark_problem):
    """Test the log-likelihood computation at nominal parameters.

    Also check that the simulation time is within the reference range.
    """
    problem_id, petab_problem, _, amici_model = benchmark_problem
    if problem_id not in problems_for_llh_check:
        pytest.skip("Excluded from log-likelihood check.")

    amici_solver = amici_model.create_solver()
    amici_solver.set_absolute_tolerance(1e-8)
    amici_solver.set_relative_tolerance(1e-8)
    amici_solver.set_max_steps(10_000)
    if problem_id in ("Brannmark_JBC2010", "Isensee_JCB2018"):
        amici_model.set_steady_state_sensitivity_mode(
            SteadyStateSensitivityMode.integrationOnly
        )

    times = dict()

    for label, sensi_mode in {
        "t_sim": SensitivityMethod.none,
        "t_fwd": SensitivityMethod.forward,
        "t_adj": SensitivityMethod.adjoint,
    }.items():
        amici_solver.set_sensitivity_method(sensi_mode)
        if sensi_mode == SensitivityMethod.none:
            amici_solver.set_sensitivity_order(SensitivityOrder.none)
        else:
            amici_solver.set_sensitivity_order(SensitivityOrder.first)

        res_repeats = [
            simulate_petab(
                petab_problem=petab_problem,
                amici_model=amici_model,
                solver=amici_solver,
                log_level=logging.DEBUG,
            )
            for _ in range(3)  # repeat to get more stable timings
        ]
        res = res_repeats[0]

        times[label] = np.min(
            [
                sum(r.cpu_time + r.cpu_time_b for r in res[RDATAS]) / 1000
                # only forwards/backwards simulation
                for res in res_repeats
            ]
        )

        if sensi_mode == SensitivityMethod.none:
            rdatas = res[RDATAS]
            llh = res[LLH]
    # TODO: check that all llh match, check that all sllh match
    times["np"] = sum(petab_problem.parameter_df[petab.ESTIMATE])

    pd.Series(times).to_csv(script_dir / f"{problem_id}_benchmark.csv")
    for rdata in rdatas:
        assert rdata.status == AMICI_SUCCESS, (
            f"Simulation failed for {rdata.id}"
        )

    if debug:
        # create simulation PEtab table
        sim_df = rdatas_to_measurement_df(
            rdatas=rdatas,
            model=amici_model,
            measurement_df=petab_problem.measurement_df,
        )
        sim_df.rename(
            columns={petab.MEASUREMENT: petab.SIMULATION}, inplace=True
        )
        sim_df.to_csv(
            benchmark_outdir / problem_id / f"{problem_id}_simulation.tsv",
            index=False,
            sep="\t",
        )

        # visualize fit, save to file
        with contextlib.suppress(NotImplementedError):
            axs = plot_problem(
                petab_problem=petab_problem, simulations_df=sim_df
            )
            for plot_id, ax in axs.items():
                fig_path = (
                    benchmark_outdir
                    / problem_id
                    / f"{problem_id}_{plot_id}_vis.png"
                )
                logger.info(f"Saving figure to {fig_path}")
                ax.get_figure().savefig(fig_path, dpi=150)

    compare_to_reference(problem_id, llh)

    for label, key in {
        "simulation": "t_sim",
        "adjoint sensitivity": "t_adj",
        "forward sensitivity": "t_fwd",
    }.items():
        try:
            ref = reference_values[problem_id][key]
            if times[key] > ref:
                pytest.fail(
                    f"Computation time for {label} ({times[key]:.2e}) "
                    f"exceeds reference ({ref:.2e})."
                )
            else:
                logger.info(
                    f"Computation time for {label} ({times[key]:.2e}) "
                    f"within reference ({ref:.2e})."
                )
        except KeyError:
            logger.error(
                f"No reference time for {label} found for "
                f"{problem_id} in {references_yaml}"
            )


@pytest.mark.filterwarnings(
    "ignore:divide by zero encountered in log",
    # https://github.com/AMICI-dev/AMICI/issues/18
    "ignore:Adjoint sensitivity analysis for models with discontinuous "
    "right hand sides .*:UserWarning",
)
@pytest.mark.parametrize("scale", (True, False), ids=["scaled", "unscaled"])
@pytest.mark.parametrize(
    "sensitivity_method",
    (SensitivityMethod.forward, SensitivityMethod.adjoint),
    ids=["forward", "adjoint"],
)
def test_benchmark_gradient(
    benchmark_problem, scale, sensitivity_method, request
):
    problem_id, petab_problem, _, amici_model = benchmark_problem
    if problem_id not in problems_for_gradient_check:
        pytest.skip("Excluded from gradient check.")

    if not scale and problem_id in (
        "Smith_BMCSystBiol2013",
        "Brannmark_JBC2010",
        "Elowitz_Nature2000",
        "Borghans_BiophysChem1997",
        "Sneyd_PNAS2002",
        "Bertozzi_PNAS2020",
        "Zheng_PNAS2012",
    ):
        # not really worth the effort trying to fix these cases if they
        # only fail on linear scale
        pytest.skip("scale=False disabled for this problem")

    petab_problem = benchmark_models_petab.get_problem(problem_id)
    if measurement_table_has_timepoint_specific_mappings(
        petab_problem.measurement_df,
    ):
        petab.flatten_timepoint_specific_output_overrides(petab_problem)

    # Only compute gradient for estimated parameters.
    parameter_ids = petab_problem.x_free_ids
    cur_settings = settings[problem_id]

    # Setup AMICI objects.
    amici_model = import_petab_problem(
        petab_problem,
        model_output_dir=benchmark_outdir / problem_id,
    )
    amici_solver = amici_model.create_solver()
    amici_solver.set_absolute_tolerance(cur_settings.atol_sim)
    amici_solver.set_relative_tolerance(cur_settings.rtol_sim)
    amici_solver.set_max_steps(2 * 10**5)
    amici_solver.set_sensitivity_method(sensitivity_method)
    # TODO: we should probably test all sensitivity modes
    amici_model.set_steady_state_sensitivity_mode(
        cur_settings.ss_sensitivity_mode
    )

    amici_function, amici_derivative = simulate_petab_to_cached_functions(
        petab_problem=petab_problem,
        free_parameter_ids=parameter_ids,
        amici_model=amici_model,
        solver=amici_solver,
        scaled_parameters=scale,
        scaled_gradients=scale,
        # FIXME: there is some issue with caching in fiddy
        #  e.g. Elowitz_Nature2000-True fails with cache=True,
        #  but not with cache=False
        # cache=not debug,
        cache=False,
        num_threads=os.cpu_count(),
    )
    np.random.seed(cur_settings.rng_seed)

    # find a point where the derivative can be computed
    for _ in range(5):
        if scale:
            point = petab_problem.x_nominal_free_scaled
            point_noise = (
                np.random.randn(len(point)) * cur_settings.noise_level
            )
        else:
            point = petab_problem.x_nominal_free
            point_noise = (
                np.random.randn(len(point)) * point * cur_settings.noise_level
            )
        point += point_noise  # avoid small gradients at nominal value

        try:
            expected_derivative = amici_derivative(point)
            break
        except RuntimeError as e:
            print(e)
            continue
    else:
        raise RuntimeError("Could not compute expected derivative.")

    derivative = get_derivative(
        function=amici_function,
        point=point,
        sizes=cur_settings.step_sizes,
        direction_ids=parameter_ids,
        method_ids=[MethodId.CENTRAL, MethodId.FORWARD, MethodId.BACKWARD],
        success_checker=Consistency(
            rtol=cur_settings.rtol_consistency,
            atol=cur_settings.atol_consistency,
        ),
        expected_result=expected_derivative,
        relative_sizes=not scale,
    )

    print()
    print("Testing at:", point)
    print("Expected derivative (amici):", expected_derivative)
    print("Print actual derivative (fiddy):", derivative.series.values)

    if debug:
        write_debug_output(
            debug_path / f"{request.node.callspec.id}.tsv",
            derivative,
            expected_derivative,
            parameter_ids,
        )

    assert_gradient_check_success(
        derivative,
        expected_derivative,
        point,
        rtol=cur_settings.rtol_check,
        atol=cur_settings.atol_check,
        always_print=True,
    )


def assert_gradient_check_success(
    derivative: fiddy.Derivative,
    expected_derivative: np.ndarray,
    point: np.ndarray,
    atol: float,
    rtol: float,
    always_print: bool = False,
) -> None:
    if not derivative.df.success.all():
        raise AssertionError(
            f"Failed to compute finite differences:\n{derivative.df}"
        )
    check = NumpyIsCloseDerivativeCheck(
        derivative=derivative,
        expectation=expected_derivative,
        point=point,
    )
    check_result = check(rtol=rtol, atol=atol)

    if check_result.success is True and not always_print:
        return

    df = check_result.df
    df["abs_diff"] = np.abs(df["expectation"] - df["test"])
    df["rel_diff"] = df["abs_diff"] / np.abs(df["expectation"])
    df["atol_success"] = df["abs_diff"] <= atol
    df["rtol_success"] = df["rel_diff"] <= rtol
    max_adiff = df["abs_diff"].max()
    max_rdiff = df["rel_diff"].max()

    success_fail = "succeeded" if check_result.success else "failed"
    with pd.option_context(
        "display.max_columns",
        None,
        "display.width",
        None,
        "display.max_rows",
        None,
    ):
        message = (
            f"Gradient check {success_fail}:\n{df}\n\n"
            f"Maximum absolute difference: {max_adiff} (tolerance: {atol})\n"
            f"Maximum relative difference: {max_rdiff} (tolerance: {rtol})"
        )

    if check_result.success is False:
        raise AssertionError(message)

    if always_print:
        print(message)


def write_debug_output(
    file_name, derivative, expected_derivative, parameter_ids
):
    df = pd.DataFrame(
        [
            {
                (
                    "fd",
                    r.metadata["size_absolute"],
                    str(r.method_id),
                ): r.value
                for c in d.computers
                for r in c.results
            }
            for d in derivative.directional_derivatives
        ],
        index=parameter_ids,
    )
    df[("fd", "full", "")] = derivative.series.values
    df[("amici", "", "")] = expected_derivative
    df["abs_diff"] = np.abs(df[("fd", "full", "")] - df[("amici", "", "")])
    df["rel_diff"] = df["abs_diff"] / np.abs(df[("amici", "", "")])

    df.to_csv(file_name, sep="\t")


@pytest.mark.filterwarnings(
    "ignore:divide by zero encountered in log",
    # https://github.com/AMICI-dev/AMICI/issues/18
    "ignore:Adjoint sensitivity analysis for models with discontinuous "
    "right hand sides .*:UserWarning",
    "ignore:.*has `useValuesFromTriggerTime=true'.*:UserWarning",
    "ignore:.*Using `log-normal` instead.*:UserWarning",
)
@pytest.mark.parametrize("problem_id", problems_for_llh_check)
def test_nominal_parameters_llh_v2(problem_id):
    """Test the log-likelihood computation at nominal parameters
    after auto-conversion of PEtab v1 benchmark problems to PEtab v2.

    Also check that the simulation time is within the reference range.
    """
    from amici.importers.petab import (
        PetabImporter,
        flatten_timepoint_specific_output_overrides,
        has_timepoint_specific_overrides,
        rdatas_to_simulation_df,
        unflatten_simulation_df,
    )
    from petab.v2 import Problem

    # TODO: differences in llh, due to log10-normal -> log-normal noise:
    #  max abs and rel diff in simulations <1e-4
    #   Bachmann_MSB2011
    #   Borghans_BiophysChem1997
    #   Elowitz_Nature2000
    #   Lucarelli_CellSystems2018
    #   Schwen_PONE2014
    #  store new reference values or recompute llh with log10-normal noise?
    v1_problem = benchmark_models_petab.get_problem(problem_id)
    if (
        petab.C.OBSERVABLE_TRANSFORMATION in v1_problem.observable_df
        and petab.C.LOG10
        in v1_problem.observable_df[petab.C.OBSERVABLE_TRANSFORMATION].unique()
    ):
        pytest.skip(
            "Problem uses log10-normal noise, not supported in PEtab v2."
        )

    if problem_id not in problems_for_llh_check:
        pytest.skip("Excluded from log-likelihood check.")

    benchmark_outdir = get_model_root_dir() / "test_bmc_v2"
    model_output_dir = benchmark_outdir / problem_id

    try:
        # Load PEtab v1 problem. This will be upgraded to v2 automatically.
        problem = Problem.from_yaml(
            benchmark_models_petab.get_problem_yaml_path(problem_id)
        )
    except ValueError as e:
        cause = f": {e.__cause__}" if e.__cause__ else ""
        pytest.skip(f"Could not load problem {problem_id}: {e}{cause}")

    was_flattened = False
    if has_timepoint_specific_overrides(problem):
        flatten_timepoint_specific_output_overrides(problem)
        was_flattened = True

    jax = False

    pi = PetabImporter(
        petab_problem=problem,
        module_name=f"{problem_id}_v2",
        outdir=model_output_dir,
        compile_=True,
        jax=jax,
    )

    ps = pi.create_simulator(force_import=True)
    ps.solver.set_absolute_tolerance(1e-8)
    ps.solver.set_relative_tolerance(1e-8)
    ps.solver.set_max_steps(10_000)
    ps.num_threads = os.cpu_count()

    if problem_id in ("Brannmark_JBC2010", "Isensee_JCB2018"):
        ps.model.set_steady_state_sensitivity_mode(
            SteadyStateSensitivityMode.integrationOnly
        )
    problem_parameters = problem.get_x_nominal_dict(free=True, fixed=True)
    res = ps.simulate(problem_parameters=problem_parameters)

    simulation_df = rdatas_to_simulation_df(
        res.rdatas, ps.model, pi.petab_problem
    )
    if was_flattened:
        simulation_df = unflatten_simulation_df(simulation_df, problem)
    print("v2 simulations:")
    print(simulation_df)
    simulation_df.to_csv("benchmark_sim.tsv", sep="\t")

    # compare to benchmark_models_petab simulations, if available
    if (
        (
            simulation_df_bm := benchmark_models_petab.get_simulation_df(
                problem_id
            )
        )
        is not None
        # https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/issues/275
        and problem_id != "Lucarelli_CellSystems2018"
    ):
        print("benchmark collection simulations:")

        assert len(simulation_df_bm) == len(simulation_df)
        # caveat: not necessarily in the same order
        simulation_df_bm["actual_simulation"] = simulation_df["simulation"]
        simulation_df_bm["abs_diff"] = np.abs(
            simulation_df_bm["actual_simulation"]
            - simulation_df_bm["simulation"]
        )
        simulation_df_bm["rel_diff"] = simulation_df_bm["abs_diff"] / np.abs(
            simulation_df_bm["simulation"]
        )
        with pd.option_context(
            "display.max_columns",
            None,
            "display.width",
            None,
            "display.max_rows",
            None,
        ):
            print(simulation_df_bm)
        print("max abs diff:", simulation_df_bm["abs_diff"].max())
        print("max rel diff:", simulation_df_bm["rel_diff"].max())

    # check llh
    compare_to_reference(problem_id, res.llh)

    # check gradient
    if problem_id not in problems_for_gradient_check:
        return None
        # pytest.skip("Excluded from gradient check.")

    # sensitivities computed w.r.t. the expected parameters? (`plist` correct?)
    ps.solver.set_sensitivity_order(SensitivityOrder.first)
    ps.solver.set_sensitivity_method(SensitivityMethod.forward)
    ps.model.set_always_check_finite(True)
    result = ps.simulate(problem_parameters=problem_parameters)
    assert result.sllh is not None
    actual_sens_pars = set(result.sllh.keys())
    expected_sens_pars = set(problem.x_free_ids)
    assert actual_sens_pars == expected_sens_pars

    # TODO
    scale = False

    # also excluded from v1 test
    if not scale and problem_id in (
        "Smith_BMCSystBiol2013",
        "Brannmark_JBC2010",
        "Elowitz_Nature2000",
        "Borghans_BiophysChem1997",
        "Sneyd_PNAS2002",
        "Bertozzi_PNAS2020",
        # "Zheng_PNAS2012",
    ):
        # not really worth the effort trying to fix these cases if they
        # only fail on linear scale
        pytest.skip("scale=False disabled for this problem")

    cur_settings = settings[problem_id]
    ps.solver.set_absolute_tolerance(cur_settings.atol_sim)
    ps.solver.set_relative_tolerance(cur_settings.rtol_sim)
    ps.solver.set_max_steps(200_000)

    # TODO: ASA + FSA
    sensitivity_method = SensitivityMethod.forward
    ps.solver.set_sensitivity_method(sensitivity_method)
    ps.model.set_steady_state_computation_mode(
        cur_settings.ss_computation_mode
    )
    # TODO: we should probably test all sensitivity modes
    ps.model.set_steady_state_sensitivity_mode(
        cur_settings.ss_sensitivity_mode
    )

    parameter_ids = ps._petab_problem.x_free_ids
    amici_function, amici_derivative = simulate_petab_v2_to_cached_functions(
        ps,
        free_parameter_ids=parameter_ids,
        cache=False,
    )
    np.random.seed(cur_settings.rng_seed)

    # find a point where the derivative can be computed
    for _ in range(5):
        if scale:
            point = ps._petab_problem.x_nominal_free_scaled
            point_noise = (
                np.random.randn(len(point)) * cur_settings.noise_level
            )
        else:
            point = ps._petab_problem.x_nominal_free
            point_noise = (
                np.random.randn(len(point)) * point * cur_settings.noise_level
            )
        point += point_noise  # avoid small gradients at nominal value

        try:
            expected_derivative = amici_derivative(point)
            break
        except RuntimeError as e:
            print(e)
            continue
    else:
        raise RuntimeError("Could not compute expected derivative.")

    derivative = get_derivative(
        function=amici_function,
        point=point,
        sizes=cur_settings.step_sizes,
        direction_ids=parameter_ids,
        method_ids=[MethodId.CENTRAL, MethodId.FORWARD, MethodId.BACKWARD],
        success_checker=Consistency(
            rtol=cur_settings.rtol_consistency,
            atol=cur_settings.atol_consistency,
        ),
        expected_result=expected_derivative,
        relative_sizes=not scale,
    )

    print()
    print("Testing at:", point)
    print("Expected derivative (amici):", expected_derivative)
    print("Print actual derivative (fiddy):", derivative.series.values)

    # if debug:
    #     write_debug_output(
    #         debug_path / f"{request.node.callspec.id}.tsv",
    #         derivative,
    #         expected_derivative,
    #         parameter_ids,
    #     )

    assert_gradient_check_success(
        derivative,
        expected_derivative,
        point,
        rtol=cur_settings.rtol_check,
        atol=cur_settings.atol_check,
        always_print=True,
    )


def compare_to_reference(problem_id: str, llh: float):
    """Compare simulation to reference value.

    For now, only the log-likelihood is checked.
    """
    try:
        ref_llh = reference_values[problem_id]["llh"]
    except KeyError:
        logger.error(
            "No reference likelihood found for "
            f"{problem_id} in {references_yaml}"
        )
        return

    if (
        simulation_df_bm := benchmark_models_petab.get_simulation_df(
            problem_id
        )
    ) is not None:
        from petab.v1.calculate import calculate_llh

        petab_problem = benchmark_models_petab.get_problem(problem_id)
        if problem_id in (
            "Bachmann_MSB2011",
            "Borghans_BiophysChem1997",
            "Elowitz_Nature2000",
            "Lucarelli_CellSystems2018",
            "Schwen_PONE2014",
        ):
            print(
                petab_problem.observable_df[
                    petab.C.OBSERVABLE_TRANSFORMATION
                ].unique()
            )
            # differences due to log10-normal -> log-normal noise
            petab_problem.observable_df.loc[
                petab_problem.observable_df[petab.C.OBSERVABLE_TRANSFORMATION]
                == petab.C.LOG10,
                petab.C.OBSERVABLE_TRANSFORMATION,
            ] = petab.C.LOG
            print(
                petab_problem.observable_df[
                    petab.C.OBSERVABLE_TRANSFORMATION
                ].unique()
            )

        try:
            ref_llh_bm = calculate_llh(
                measurement_dfs=petab_problem.measurement_df,
                simulation_dfs=simulation_df_bm,
                observable_dfs=petab_problem.observable_df,
                parameter_dfs=petab_problem.parameter_df,
            )
            print(
                "Reference llh from benchmark collection simulation table:",
                ref_llh_bm,
            )
            # TODO https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/issues/278
            # assert np.isclose(
            #     ref_llh, ref_llh_bm
            # ), f"Stored Reference llh {ref_llh} differs from the value computed "\
            #    f"from the benchmark collection simulation table {ref_llh_bm}"
        except Exception as e:
            print(
                "Could not compute reference llh from benchmark collection"
                " simulation table:",
                e,
            )
    rdiff = np.abs((llh - ref_llh) / ref_llh)
    rtol = 1e-3
    adiff = np.abs(llh - ref_llh)
    atol = 1e-3
    tolstr = (
        f" Absolute difference is {adiff:.2e} "
        f"(tol {atol:.2e}) and relative difference is "
        f"{rdiff:.2e} (tol {rtol:.2e})."
    )

    if np.isclose(llh, ref_llh, rtol=rtol, atol=atol):
        logger.info(
            f"Computed llh {llh:.4e} matches reference {ref_llh:.4e}." + tolstr
        )
    else:
        pytest.fail(
            f"Computed llh {llh:.4e} does not match reference "
            f"{ref_llh:.4e}." + tolstr
        )
