import copy
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import pytest
import petab.v1 as petab
from petab.v1.lint import measurement_table_has_timepoint_specific_mappings

import benchmark_models_petab
import amici
from amici.petab.petab_import import import_petab_problem

script_dir = Path(__file__).parent
repo_root = script_dir.parent.parent
benchmark_outdir = repo_root / "test_bmc"

# problem IDs for which to check the gradient
# TODO: extend
problems_for_gradient_check = list(
    sorted(
        set(benchmark_models_petab.MODELS)
        - {
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
    )
)

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
    "Lucarelli_CellSystems2018",
    "Schwen_PONE2014",
    "Smith_BMCSystBiol2013",
    "Sneyd_PNAS2002",
    "Weber_BMC2015",
    "Zheng_PNAS2012",
]

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
    ss_sensitivity_mode: amici.SteadyStateSensitivityMode = (
        amici.SteadyStateSensitivityMode.integrateIfNewtonFails
    )
    noise_level: float = 0.05


settings = defaultdict(GradientCheckSettings)
# NOTE: Newton method fails badly with ASA for Blasi_CellSystems2016
settings["Blasi_CellSystems2016"] = GradientCheckSettings(
    atol_check=1e-12,
    rtol_check=1e-4,
    ss_sensitivity_mode=amici.SteadyStateSensitivityMode.integrationOnly,
)
settings["Borghans_BiophysChem1997"] = GradientCheckSettings(
    rng_seed=2,
    atol_check=1e-5,
    rtol_check=1e-3,
)
settings["Brannmark_JBC2010"] = GradientCheckSettings(
    ss_sensitivity_mode=amici.SteadyStateSensitivityMode.integrationOnly,
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
)


@pytest.fixture(scope="session", params=problems, ids=problems)
def benchmark_problem(request):
    """Fixture providing model and PEtab problem for a benchmark model."""
    problem_id = request.param
    petab_problem = benchmark_models_petab.get_problem(problem_id)
    flat_petab_problem = copy.deepcopy(petab_problem)
    if measurement_table_has_timepoint_specific_mappings(
        petab_problem.measurement_df,
    ):
        petab.flatten_timepoint_specific_output_overrides(flat_petab_problem)

    amici_model = import_petab_problem(
        flat_petab_problem,
        model_output_dir=benchmark_outdir / problem_id,
    )
    return problem_id, flat_petab_problem, petab_problem, amici_model
