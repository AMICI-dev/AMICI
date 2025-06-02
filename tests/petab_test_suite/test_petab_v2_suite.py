#!/usr/bin/env python3
"""Test PEtab v2 import"""

import logging
import shutil
import sys

import numpy as np
import pandas as pd
import petabtests
import pytest
from _pytest.outcomes import Skipped

# from amici.gradient_check import check_derivatives as amici_check_derivatives
from amici.logging import get_logger, set_log_level
from amici.petab.petab_importer import *
from petab import v2

logger = get_logger(__name__, logging.DEBUG)
set_log_level(get_logger("amici.petab_import"), logging.DEBUG)
stream_handler = logging.StreamHandler()
logger.addHandler(stream_handler)


@pytest.mark.filterwarnings(
    "ignore:Event `_E0` has `useValuesFromTriggerTime=true'"
)
def test_case(case, model_type, version, jax):
    """Wrapper for _test_case for handling test outcomes"""
    try:
        _test_case(case, model_type, version, jax)
    except Exception as e:
        if isinstance(
            e, NotImplementedError
        ) or "Timepoint-specific parameter overrides" in str(e):
            logger.info(
                f"Case {case} expectedly failed. "
                "Required functionality is not yet "
                f"implemented: {e}"
            )
            pytest.skip(str(e))
        else:
            raise e


def _test_case(case, model_type, version, jax):
    """Run a single PEtab test suite case"""

    case = petabtests.test_id_str(case)
    logger.debug(f"Case {case} [{model_type}] [{version}] [{jax}]")

    # load
    case_dir = petabtests.get_case_dir(case, model_type, version)
    yaml_file = case_dir / petabtests.problem_yaml_name(case)
    problem = v2.Problem.from_yaml(yaml_file)

    # compile amici model
    if case.startswith("0006") and not jax:
        # TODO: petab.flatten_timepoint_specific_output_overrides(problem)
        # petab.flatten_timepoint_specific_output_overrides(problem)
        pytest.skip("Timepoint-specific output overrides not yet supported")

    model_name = (
        f"petab_{model_type}_test_case_{case}_{version.replace('.', '_')}"
    )
    model_output_dir = f"amici_models/{model_name}" + ("_jax" if jax else "")

    pi = PetabImporter(
        petab_problem=problem,
        model_id=model_name,
        outdir=model_output_dir,
        compile_=True,
        jax=jax,
        force_import=True,
    )
    # TODO force re-import
    shutil.rmtree(pi.outdir, ignore_errors=True)

    ps = pi.create_simulator()
    ps._solver.set_steady_state_tolerance_factor(1.0)

    problem_parameters = dict(
        zip(problem.x_free_ids, problem.x_nominal_free, strict=True)
    )

    ret = ps.simulate(problem_parameters=problem_parameters)

    rdatas = ret["rdatas"]
    chi2 = sum(rdata["chi2"] for rdata in rdatas)
    llh = ret["llh"]
    simulation_df = rdatas_to_measurement_df(
        rdatas, ps._model, pi.petab_problem
    )
    # TODO petab.check_measurement_df(simulation_df, problem.observable_df)
    simulation_df = simulation_df.rename(
        columns={v2.C.MEASUREMENT: v2.C.SIMULATION}
    )
    # revert setting default experiment Id
    simulation_df.loc[
        simulation_df[v2.C.EXPERIMENT_ID] == "__default__", v2.C.EXPERIMENT_ID
    ] = np.nan
    # FIXME: why int?? can be inf
    # simulation_df[v2.C.TIME] = simulation_df[v2.C.TIME].astype(int)
    solution = petabtests.load_solution(case, model_type, version=version)
    gt_chi2 = solution[petabtests.CHI2]
    gt_llh = solution[petabtests.LLH]
    gt_simulation_dfs = solution[petabtests.SIMULATION_DFS]
    if case.startswith("0006") and not jax:
        # account for flattening
        gt_simulation_dfs[0].loc[:, v2.C.OBSERVABLE_ID] = (
            "obs_a__10__c0",
            "obs_a__15__c0",
        )
    tol_chi2 = solution[petabtests.TOL_CHI2]
    tol_llh = solution[petabtests.TOL_LLH]
    tol_simulations = solution[petabtests.TOL_SIMULATIONS]

    chi2s_match = petabtests.evaluate_chi2(chi2, gt_chi2, tol_chi2)
    llhs_match = petabtests.evaluate_llh(llh, gt_llh, tol_llh)
    simulations_match = petabtests.evaluate_simulations(
        [simulation_df], gt_simulation_dfs, tol_simulations
    )

    logger.log(
        logging.DEBUG if simulations_match else logging.ERROR,
        f"Simulations: match = {simulations_match}",
    )
    if not simulations_match:
        with pd.option_context(
            "display.max_rows",
            None,
            "display.max_columns",
            None,
            "display.width",
            200,
        ):
            logger.log(
                logging.DEBUG,
                f"x_ss: {ps._model.get_state_ids()} "
                f"{[rdata.x_ss for rdata in rdatas]}"
                f"preeq_t: {[rdata.preeq_t for rdata in rdatas]}"
                f"posteq_t: {[rdata.posteq_t for rdata in rdatas]}",
            )
            logger.log(
                logging.ERROR, f"Expected simulations:\n{gt_simulation_dfs}"
            )
            logger.log(logging.ERROR, f"Actual simulations:\n{simulation_df}")
    logger.log(
        logging.DEBUG if chi2s_match else logging.ERROR,
        f"CHI2: simulated: {chi2}, expected: {gt_chi2}, match = {chi2s_match}",
    )
    logger.log(
        logging.DEBUG if simulations_match else logging.ERROR,
        f"LLH: simulated: {llh}, expected: {gt_llh}, match = {llhs_match}",
    )

    if jax:
        pass  # skip derivative checks for now
    else:
        pass
        # FIXME: later
        # check_derivatives(ps._petab_problem, ps._model, ps._solver, problem_parameters)

    if not all([llhs_match, simulations_match]) or not chi2s_match:
        logger.error(f"Case {case} failed.")
        raise AssertionError(
            f"Case {case}: Test results do not match expectations"
        )

    logger.info(f"Case {case} passed.")


#
# def check_derivatives(
#     problem: petab.Problem,
#     model: amici.Model,
#     solver: amici.Solver,
#     problem_parameters: dict[str, float],
# ) -> None:
#     """Check derivatives using finite differences for all experimental
#     conditions
#
#     Arguments:
#         problem: PEtab problem
#         model: AMICI model matching ``problem``
#         solver: AMICI solver
#         problem_parameters: Dictionary of problem parameters
#     """
#     solver.setSensitivityMethod(amici.SensitivityMethod.forward)
#     solver.setSensitivityOrder(amici.SensitivityOrder.first)
#     # Required for case 9 to not fail in
#     #  amici::NewtonSolver::computeNewtonSensis
#     model.setSteadyStateSensitivityMode(
#         SteadyStateSensitivityMode.integrateIfNewtonFails
#     )
#
#     for edata in create_parameterized_edatas(
#         amici_model=model,
#         petab_problem=problem,
#         problem_parameters=problem_parameters,
#     ):
#         # check_derivatives does currently not support parameters in ExpData
#         # set parameter scales before setting parameter values!
#         model.setParameterScale(edata.pscale)
#         model.setParameters(edata.parameters)
#         edata.parameters = []
#         edata.pscale = amici.parameterScalingFromIntVector([])
#         amici_check_derivatives(model, solver, edata)


def run():
    """Run the full PEtab test suite"""
    n_success = 0
    n_skipped = 0
    n_total = 0
    version = "v2.0.0"
    jax = False

    cases = petabtests.get_cases("sbml", version=version)
    # FIXME
    #  0019: "PEtab v2.0.0 mapping tables are not yet supported."
    #   -> not compliant with v2.0.0?! unnecessary aliasing
    # cases = [
    #     "0019",
    # ]

    n_total += len(cases)
    for case in cases:
        try:
            test_case(case, "sbml", version=version, jax=jax)
            n_success += 1
        except Skipped:
            n_skipped += 1
        except Exception as e:
            # run all despite failures
            logger.error(f"Case {case} failed.")
            logger.exception(e)

    logger.info(f"{n_success} / {n_total} successful, {n_skipped} skipped")
    if n_success != len(cases):
        sys.exit(1)


if __name__ == "__main__":
    run()
