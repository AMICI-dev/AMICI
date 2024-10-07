"""Tests for simulate_petab on PEtab benchmark problems."""

from pathlib import Path

import fiddy

import amici
import numpy as np
import pandas as pd
import petab.v1 as petab
import pytest
from amici.petab.petab_import import import_petab_problem
import benchmark_models_petab
from collections import defaultdict
from dataclasses import dataclass
from amici import SensitivityMethod
from fiddy import MethodId, get_derivative
from fiddy.derivative_check import NumpyIsCloseDerivativeCheck
from fiddy.extensions.amici import simulate_petab_to_cached_functions
from fiddy.success import Consistency

repo_root = Path(__file__).parent.parent.parent

# reuse compiled models from test_benchmark_collection.sh
benchmark_outdir = repo_root / "test_bmc"
models = set(benchmark_models_petab.MODELS) - {
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
    "Fujita_SciSignal2010",
    # FIXME: re-enable
    "Raia_CancerResearch2011",
    "Zheng_PNAS2012",
}
models = list(sorted(models))

debug = False
if debug:
    debug_path = Path(__file__).parent / "debug"
    debug_path.mkdir(exist_ok=True, parents=True)


@dataclass
class GradientCheckSettings:
    # Absolute and relative tolerances for simulation
    atol_sim: float = 1e-16
    rtol_sim: float = 1e-12
    # Absolute and relative tolerances for finite difference gradient checks.
    atol_check: float = 1e-3
    rtol_check: float = 1e-2
    # Step sizes for finite difference gradient checks.
    step_sizes = [
        1e-1,
        5e-2,
        1e-2,
        1e-3,
        1e-4,
        1e-5,
    ]
    rng_seed: int = 0


settings = defaultdict(GradientCheckSettings)
settings["Smith_BMCSystBiol2013"] = GradientCheckSettings(
    atol_sim=1e-10,
    rtol_sim=1e-10,
)
settings["Oliveira_NatCommun2021"] = GradientCheckSettings(
    # Avoid "root after reinitialization"
    atol_sim=1e-12,
    rtol_sim=1e-12,
)
settings["Okuonghae_ChaosSolitonsFractals2020"] = GradientCheckSettings(
    atol_sim=1e-14,
    rtol_sim=1e-14,
)
settings["Okuonghae_ChaosSolitonsFractals2020"].step_sizes.insert(0, 0.2)
settings["Giordano_Nature2020"] = GradientCheckSettings(
    atol_check=1e-6, rtol_check=1e-3, rng_seed=1
)


def assert_gradient_check_success(
    derivative: fiddy.Derivative,
    expected_derivative: np.array,
    point: np.array,
    atol: float,
    rtol: float,
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

    if check_result.success is True:
        return

    raise AssertionError(f"Gradient check failed:\n{check_result.df}")


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
@pytest.mark.parametrize("model", models)
def test_benchmark_gradient(model, scale, sensitivity_method, request):
    if not scale and model in (
        "Smith_BMCSystBiol2013",
        "Brannmark_JBC2010",
        "Elowitz_Nature2000",
        "Borghans_BiophysChem1997",
        "Sneyd_PNAS2002",
        "Bertozzi_PNAS2020",
        "Okuonghae_ChaosSolitonsFractals2020",
    ):
        # not really worth the effort trying to fix these cases if they
        # only fail on linear scale
        pytest.skip()

    if (
        model
        in (
            "Blasi_CellSystems2016",
            # events with parameter-dependent triggers
            #  https://github.com/AMICI-dev/AMICI/issues/18
            "Oliveira_NatCommun2021",
        )
        and sensitivity_method == SensitivityMethod.adjoint
    ):
        # FIXME
        pytest.skip()

    if (
        model
        in (
            "Weber_BMC2015",
            "Sneyd_PNAS2002",
        )
        and sensitivity_method == SensitivityMethod.forward
    ):
        # FIXME
        pytest.skip()

    petab_problem = benchmark_models_petab.get_problem(model)
    petab.flatten_timepoint_specific_output_overrides(petab_problem)

    # Only compute gradient for estimated parameters.
    parameter_df_free = petab_problem.parameter_df.loc[
        petab_problem.x_free_ids
    ]
    parameter_ids = list(parameter_df_free.index)
    cur_settings = settings[model]

    # Setup AMICI objects.
    amici_model = import_petab_problem(
        petab_problem,
        model_output_dir=benchmark_outdir / model,
    )
    amici_solver = amici_model.getSolver()
    amici_solver.setAbsoluteTolerance(cur_settings.atol_sim)
    amici_solver.setRelativeTolerance(cur_settings.rtol_sim)
    amici_solver.setMaxSteps(int(1e5))
    amici_solver.setSensitivityMethod(sensitivity_method)

    if model in ("Brannmark_JBC2010",):
        amici_model.setSteadyStateSensitivityMode(
            amici.SteadyStateSensitivityMode.integrationOnly
        )

    amici_function, amici_derivative = simulate_petab_to_cached_functions(
        petab_problem=petab_problem,
        parameter_ids=parameter_ids,
        amici_model=amici_model,
        solver=amici_solver,
        scaled_parameters=scale,
        scaled_gradients=scale,
        # FIXME: there is some issue with caching in fiddy
        #  e.g. Elowitz_Nature2000-True fails with cache=True,
        #  but not with cache=False
        # cache=not debug,
        cache=False,
    )

    noise_level = 0.1
    np.random.seed(cur_settings.rng_seed)

    # find a point where the derivative can be computed
    for _ in range(5):
        if scale:
            point = np.asarray(
                list(
                    petab_problem.scale_parameters(
                        dict(parameter_df_free.nominalValue)
                    ).values()
                )
            )
            point_noise = np.random.randn(len(point)) * noise_level
        else:
            point = parameter_df_free.nominalValue.values
            point_noise = np.random.randn(len(point)) * point * noise_level
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
        success_checker=Consistency(atol=1e-5, rtol=1e-1),
        expected_result=expected_derivative,
        relative_sizes=not scale,
    )

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
    )


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
