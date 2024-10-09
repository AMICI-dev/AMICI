#!/usr/bin/env python3
import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import amici
import petab.v1 as petab
from amici.petab.petab_import import import_model


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
        "numstepsB",
        "numrhsevals",
        "numrhsevalsB",
        "numerrtestfails",
        "numerrtestfailsB",
        "numnonlinsolvconvfails",
        "numnonlinsolvconvfailsB",
        "preeq_status",
        "preeq_numsteps",
        "preeq_numstepsB",
        "preeq_cpu_time",
        "preeq_cpu_timeB",
        "cpu_time_total",
        "cpu_time",
        "cpu_timeB",
        "posteq_status",
        "posteq_cpu_time",
        "posteq_cpu_timeB",
        "posteq_numsteps",
        "posteq_numstepsB",
    ]
    for d in diagnostics:
        print(d, rdata[d])
    assert rdata["status"] == amici.AMICI_SUCCESS


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
    amici.de_export.logger.setLevel(logging.DEBUG)
    amici.sbml_import.logger.setLevel(logging.DEBUG)
    import_model(
        model_name=model_name,
        model_output_dir=model_dir,
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
        solver.setSensitivityMethod(amici.SensitivityMethod.none)
        solver.setSensitivityOrder(amici.SensitivityOrder.none)
    elif arg == "forward_sensitivities":
        model.setParameterList(list(range(100)))
        solver.setSensitivityMethod(amici.SensitivityMethod.forward)
        solver.setSensitivityOrder(amici.SensitivityOrder.first)
    elif arg == "adjoint_sensitivities":
        solver.setSensitivityMethod(amici.SensitivityMethod.adjoint)
        solver.setSensitivityOrder(amici.SensitivityOrder.first)
    elif arg == "forward_simulation_non_optimal_parameters":
        tmp_par = model.getParameters()
        model.setParameters([0.1 for _ in tmp_par])
        solver.setSensitivityMethod(amici.SensitivityMethod.none)
        solver.setSensitivityOrder(amici.SensitivityOrder.none)
    elif arg == "adjoint_sensitivities_non_optimal_parameters":
        tmp_par = model.getParameters()
        model.setParameters([0.1 for _ in tmp_par])
        solver.setSensitivityMethod(amici.SensitivityMethod.adjoint)
        solver.setSensitivityOrder(amici.SensitivityOrder.first)
    elif arg == "forward_steadystate_sensitivities_non_optimal_parameters":
        tmp_par = model.getParameters()
        model.setParameters([0.1 for _ in tmp_par])
        solver.setSensitivityMethod(amici.SensitivityMethod.forward)
        solver.setSensitivityOrder(amici.SensitivityOrder.first)
        model.setSteadyStateSensitivityMode(
            amici.SteadyStateSensitivityMode.newtonOnly
        )
        edata.setTimepoints([float("inf")])
    elif arg == "adjoint_steadystate_sensitivities_non_optimal_parameters":
        tmp_par = model.getParameters()
        model.setParameters([0.1 for _ in tmp_par])
        solver.setSensitivityMethod(amici.SensitivityMethod.adjoint)
        solver.setSensitivityOrder(amici.SensitivityOrder.first)
        edata.setTimepoints([float("inf")])
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
        model = model_module.getModel()
        solver = model.getSolver()
        # TODO
        edata = amici.ExpData(model)
        edata.setTimepoints([1e8])
        edata.setObservedData([1.0])
        edata.setObservedDataStdDev([1.0])

    prepare_simulation(arg, model, solver, edata)

    rdata = amici.runAmiciSimulation(model, solver, edata)

    check_results(rdata)


if __name__ == "__main__":
    main()
