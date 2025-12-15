#!/usr/bin/env python3
import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import amici
import amici.exporters.sundials.de_export
import amici.importers.sbml
import petab.v1 as petab
from amici.importers.petab.v1.petab_import import import_model_sbml
from amici.sim.sundials import (
    AMICI_SUCCESS,
    ExpData,
    SensitivityMethod,
    SensitivityOrder,
    SteadyStateSensitivityMode,
    run_simulation,
)


def parse_args():
    arg = sys.argv[1]
    if "_" in arg and re.match(r"O[0-2]", arg.split("_")[-1]):
        optim = arg.split("_")[-1]
        os.environ["AMICI_CXXFLAGS"] = f"-{optim}"
        print(f"{os.environ['AMICI_CXXFLAGS']=}")
        suffix = f"_{optim}"
        arg = "_".join(arg.split("_")[:-1])
    else:
        suffix = ""

    return arg, suffix


def check_results(rdata):
    diagnostics = [
        "numsteps",
        "numsteps_b",
        "num_rhs_evals",
        "num_rhs_evals_b",
        "num_err_test_fails",
        "num_err_test_fails_b",
        "num_non_lin_solv_conv_fails",
        "num_non_lin_solv_conv_fails_b",
        "preeq_status",
        "preeq_numsteps",
        "preeq_numsteps_b",
        "preeq_cpu_time",
        "preeq_cpu_time_b",
        "cpu_time_total",
        "cpu_time",
        "cpu_time_b",
        "posteq_status",
        "posteq_cpu_time",
        "posteq_cpu_time_b",
        "posteq_numsteps",
        "posteq_numsteps_b",
    ]
    for d in diagnostics:
        print(d, rdata[d])
    assert rdata["status"] == AMICI_SUCCESS


def run_import(model_name, model_dir: Path):
    git_dir = Path("CS_Signalling_ERBB_RAS_AKT")
    if not git_dir.exists():
        subprocess.run(
            [
                "git",
                "clone",
                "--depth",
                "1",
                "https://github.com/ICB-DCM/CS_Signalling_ERBB_RAS_AKT",
            ]
        )

    pp = petab.Problem.from_yaml(
        git_dir / "FroehlichKes2018" / "PEtab" / "FroehlichKes2018.yaml"
    )
    petab.lint_problem(pp)
    amici.exporters.sundials.de_export.logger.setLevel(logging.DEBUG)
    amici.importers.sbml.logger.setLevel(logging.DEBUG)
    import_model_sbml(
        model_name=model_name,
        output_dir=model_dir,
        petab_problem=pp,
        compile=False,
        verbose=True,
    )


def compile_model(model_dir_source: Path, model_dir_compiled: Path):
    if model_dir_source != model_dir_compiled:
        shutil.copytree(model_dir_source, model_dir_compiled)

    cmd = ["python", "setup.py", "build_ext", "--build-lib=.", "--force"]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, cwd=model_dir_compiled, check=True)


def prepare_simulation(arg, model, solver, edata):
    if arg == "forward_simulation":
        solver.set_sensitivity_method(SensitivityMethod.none)
        solver.set_sensitivity_order(SensitivityOrder.none)
    elif arg == "forward_sensitivities":
        model.set_parameter_list(list(range(100)))
        solver.set_sensitivity_method(SensitivityMethod.forward)
        solver.set_sensitivity_order(SensitivityOrder.first)
    elif arg == "adjoint_sensitivities":
        solver.set_sensitivity_method(SensitivityMethod.adjoint)
        solver.set_sensitivity_order(SensitivityOrder.first)
    elif arg == "forward_simulation_non_optimal_parameters":
        tmp_par = model.get_free_parameters()
        model.set_free_parameters([0.1 for _ in tmp_par])
        solver.set_sensitivity_method(SensitivityMethod.none)
        solver.set_sensitivity_order(SensitivityOrder.none)
    elif arg == "adjoint_sensitivities_non_optimal_parameters":
        tmp_par = model.get_free_parameters()
        model.set_free_parameters([0.1 for _ in tmp_par])
        solver.set_sensitivity_method(SensitivityMethod.adjoint)
        solver.set_sensitivity_order(SensitivityOrder.first)
    elif arg == "forward_steadystate_sensitivities_non_optimal_parameters":
        tmp_par = model.get_free_parameters()
        model.set_free_parameters([0.1 for _ in tmp_par])
        solver.set_sensitivity_method(SensitivityMethod.forward)
        solver.set_sensitivity_order(SensitivityOrder.first)
        model.set_steady_state_sensitivity_mode(
            SteadyStateSensitivityMode.newtonOnly
        )
        edata.set_timepoints([float("inf")])
    elif arg == "adjoint_steadystate_sensitivities_non_optimal_parameters":
        tmp_par = model.get_free_parameters()
        model.set_free_parameters([0.1 for _ in tmp_par])
        solver.set_sensitivity_method(SensitivityMethod.adjoint)
        solver.set_sensitivity_order(SensitivityOrder.first)
        edata.set_timepoints([float("inf")])
    else:
        print("Unknown argument:", arg)
        sys.exit(1)


def main():
    arg, suffix = parse_args()

    # Model is imported once to this directory
    model_dir_source = Path("model_performance_test")
    # and copied to and compiled in this directory with different compiler
    #  options
    model_dir_compiled = Path(f"model_performance_test_{suffix}")
    model_name = "model_performance_test"

    if arg == "import":
        run_import(model_name, model_dir_source)
        return
    elif arg == "compile":
        compile_model(model_dir_source, model_dir_compiled)
        return
    else:
        model_module = amici.import_model_module(
            model_name, model_dir_compiled
        )
        model = model_module.get_model()
        solver = model.create_solver()
        # TODO
        edata = ExpData(model)
        edata.set_timepoints([1e8])
        edata.set_measurements([1.0])
        edata.set_noise_scales([1.0])

    prepare_simulation(arg, model, solver, edata)

    rdata = run_simulation(model, solver, edata)

    check_results(rdata)


if __name__ == "__main__":
    main()
