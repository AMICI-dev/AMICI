#!/usr/bin/env python3
"""Run PEtab v1 test suite (https://github.com/PEtab-dev/petab_test_suite)"""

import logging
import sys

import diffrax
import pandas as pd
import petab.v1 as petab
import petabtests
import pytest
from _pytest.outcomes import Skipped
from amici.importers.petab.v1 import (
    import_petab_problem,
)
from amici.logging import get_logger, set_log_level
from amici.sim.sundials import (
    Model,
    SensitivityMethod,
    SensitivityOrder,
    Solver,
    SteadyStateSensitivityMode,
)
from amici.sim.sundials.gradient_check import (
    check_derivatives as amici_check_derivatives,
)
from amici.sim.sundials.petab.v1 import (
    create_parameterized_edatas,
    rdatas_to_measurement_df,
    simulate_petab,
)

logger = get_logger(__name__, logging.DEBUG)
set_log_level(get_logger("amici.petab_import"), logging.DEBUG)
stream_handler = logging.StreamHandler()
logger.addHandler(stream_handler)


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
    problem = petab.Problem.from_yaml(yaml_file)

    if problem.mapping_df is not None:
        pytest.skip(
            "PEtab test suite cases with mapping_df are not supported yet."
        )

    # compile amici model
    if case.startswith("0006") and not jax:
        petab.flatten_timepoint_specific_output_overrides(problem)
    model_name = (
        f"petab_{model_type}_test_case_{case}_{version.replace('.', '_')}"
    )
    imported = import_petab_problem(
        petab_problem=problem,
        model_name=model_name,
        compile_=True,
        jax=jax,
    )
    if jax:
        from amici.jax import petab_simulate, run_simulations

        steady_state_event = diffrax.steady_state_event(rtol=1e-6, atol=1e-6)
        jax_problem = (
            imported  # import_petab_problem returns JAXProblem when jax=True
        )
        llh, ret = run_simulations(
            jax_problem, steady_state_event=steady_state_event
        )
        chi2, _ = run_simulations(
            jax_problem, ret="chi2", steady_state_event=steady_state_event
        )
        simulation_df = petab_simulate(
            jax_problem, steady_state_event=steady_state_event
        )
        simulation_df.rename(
            columns={petab.SIMULATION: petab.MEASUREMENT}, inplace=True
        )
    else:
        model = imported  # import_petab_problem returns Model when jax=False
        solver = model.create_solver()
        solver.set_steady_state_tolerance_factor(1.0)
        problem_parameters = dict(
            zip(problem.x_free_ids, problem.x_nominal_free, strict=True)
        )

        # simulate
        ret = simulate_petab(
            problem,
            model,
            problem_parameters=problem_parameters,
            solver=solver,
            log_level=logging.DEBUG,
        )

        rdatas = ret["rdatas"]
        chi2 = sum(rdata["chi2"] for rdata in rdatas)
        llh = ret["llh"]
        simulation_df = rdatas_to_measurement_df(
            rdatas, model, problem.measurement_df
        )
    petab.check_measurement_df(simulation_df, problem.observable_df)
    simulation_df = simulation_df.rename(
        columns={petab.MEASUREMENT: petab.SIMULATION}
    )
    simulation_df[petab.TIME] = simulation_df[petab.TIME].astype(int)
    solution = petabtests.load_solution(case, model_type, version=version)
    gt_chi2 = solution[petabtests.CHI2]
    gt_llh = solution[petabtests.LLH]
    gt_simulation_dfs = solution[petabtests.SIMULATION_DFS]
    if case.startswith("0006") and not jax:
        # account for flattening
        gt_simulation_dfs[0].loc[:, petab.OBSERVABLE_ID] = (
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
            if not jax:
                logger.log(
                    logging.DEBUG,
                    f"x_ss: {model.get_state_ids()} "
                    f"{[rdata.x_ss for rdata in rdatas]}",
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
        check_derivatives(problem, model, solver, problem_parameters)

    if not all([llhs_match, simulations_match]) or not chi2s_match:
        logger.error(f"Case {case} failed.")
        raise AssertionError(
            f"Case {case}: Test results do not match expectations"
        )

    logger.info(f"Case {case} passed.")


def check_derivatives(
    problem: petab.Problem,
    model: Model,
    solver: Solver,
    problem_parameters: dict[str, float],
) -> None:
    """Check derivatives using finite differences for all experimental
    conditions

    Arguments:
        problem: PEtab problem
        model: AMICI model matching ``problem``
        solver: AMICI solver
        problem_parameters: Dictionary of problem parameters
    """
    solver.set_sensitivity_method(SensitivityMethod.forward)
    solver.set_sensitivity_order(SensitivityOrder.first)
    # Required for case 9 to not fail in
    #  amici::NewtonSolver::computeNewtonSensis
    model.set_steady_state_sensitivity_mode(
        SteadyStateSensitivityMode.integrateIfNewtonFails
    )

    for edata in create_parameterized_edatas(
        amici_model=model,
        petab_problem=problem,
        problem_parameters=problem_parameters,
    ):
        amici_check_derivatives(model, solver, edata)


def run():
    """Run the full PEtab test suite"""
    n_success = 0
    n_skipped = 0
    n_total = 0
    for version in ("v1.0.0", "v2.0.0"):
        for jax in (False, True):
            cases = petabtests.get_cases("sbml", version=version)
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
                    logger.error(e)

    logger.info(f"{n_success} / {n_total} successful, {n_skipped} skipped")
    if n_success != len(cases):
        sys.exit(1)


if __name__ == "__main__":
    run()
