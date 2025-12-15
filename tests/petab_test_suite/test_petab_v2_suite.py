#!/usr/bin/env python3
"""Test PEtab v2 import"""

import logging
import sys

import pandas as pd
import petabtests
import pytest
from _pytest.outcomes import Skipped
from amici.importers.petab import *
from amici.logging import get_logger, set_log_level
from amici.sim.sundials import (
    AMICI_SUCCESS,
    SensitivityMethod,
    SensitivityOrder,
)
from amici.sim.sundials.gradient_check import (
    check_derivatives as amici_check_derivatives,
)
from amici.sim.sundials.petab import PetabSimulator
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
                f"Required functionality is not yet implemented: {e}"
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
        flatten_timepoint_specific_output_overrides(problem)

    model_name = (
        f"petab_{model_type}_test_case_{case}_{version.replace('.', '_')}"
    )

    pi = PetabImporter(
        petab_problem=problem,
        module_name=model_name,
        compile_=True,
        jax=jax,
    )

    ps = pi.create_simulator(
        force_import=True,
    )
    ps.solver.set_steady_state_tolerance_factor(1.0)

    problem_parameters = problem.get_x_nominal_dict(free=True, fixed=True)
    res = ps.simulate(problem_parameters=problem_parameters)
    rdatas = res.rdatas
    for rdata in rdatas:
        assert rdata.status == AMICI_SUCCESS, (
            f"Simulation failed for {rdata.id}"
        )
    chi2 = sum(rdata.chi2 for rdata in rdatas)
    llh = res.llh
    simulation_df = rdatas_to_simulation_df(rdatas, ps.model, pi.petab_problem)

    solution = petabtests.load_solution(case, model_type, version=version)
    gt_chi2 = solution[petabtests.CHI2]
    gt_llh = solution[petabtests.LLH]
    gt_simulation_dfs = solution[petabtests.SIMULATION_DFS]
    if case.startswith("0006") and not jax:
        unflattened_problem = v2.Problem.from_yaml(yaml_file)
        simulation_df = unflatten_simulation_df(
            simulation_df, unflattened_problem
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
            if not jax:
                logger.log(
                    logging.DEBUG,
                    f"x_ss: {ps.model.get_state_ids()} "
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
        if (case, model_type, version) in (
            ("0016", "sbml", "v2.0.0"),
            ("0013", "pysb", "v2.0.0"),
        ):
            # FIXME: issue with events and sensitivities
            ...
        else:
            check_derivatives(ps, problem_parameters)

    if not all([llhs_match, simulations_match]) or not chi2s_match:
        logger.error(f"Case {case} failed.")
        raise AssertionError(
            f"Case {case}: Test results do not match expectations"
        )

    logger.info(f"Case {case} passed.")


def check_derivatives(
    petab_simulator: PetabSimulator,
    problem_parameters: dict[str, float],
) -> None:
    """Check derivatives using finite differences for each experimental
    conditions.

    :param petab_simulator: PetabSimulator object
    :param problem_parameters: Dictionary of problem parameters
    """
    solver = petab_simulator.solver
    model = petab_simulator.model
    edatas = petab_simulator.exp_man.create_edatas()

    solver.set_sensitivity_method(SensitivityMethod.forward)
    solver.set_sensitivity_order(SensitivityOrder.first)

    for edata in edatas:
        petab_simulator.exp_man.apply_parameters(edata, problem_parameters)
        amici_check_derivatives(model, solver=solver, edata=edata)

    # TODO check aggregated sensitivities over all conditions


def run():
    """Run the full PEtab test suite"""
    n_success = 0
    n_skipped = 0
    n_total = 0
    version = "v2.0.0"
    jax = False

    cases = list(petabtests.get_cases("sbml", version=version))
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
