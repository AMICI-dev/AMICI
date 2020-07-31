#!/usr/bin/env python3

"""Run PEtab test suite (https://github.com/PEtab-dev/petab_test_suite)"""

import logging
import os
import sys

import amici
import petab
import petabtests
import pytest
from _pytest.outcomes import Skipped
from amici.gradient_check import check_derivatives as amici_check_derivatives
from amici.logging import get_logger, set_log_level
from amici.petab_import import import_petab_problem, PysbPetabProblem
from amici.petab_objective import (
    simulate_petab, rdatas_to_measurement_df, create_parameterized_edatas)
from amici import SteadyStateSensitivityMode_simulationFSA

logger = get_logger(__name__, logging.DEBUG)
set_log_level(get_logger("amici.petab_import"), logging.DEBUG)
stream_handler = logging.StreamHandler()
logger.addHandler(stream_handler)


def test_case(case, model_type):
    """Wrapper for _test_case for handling test outcomes"""
    try:
        _test_case(case, model_type)
    except Exception as e:
        if isinstance(e, NotImplementedError) \
                or "Timepoint-specific parameter overrides" in str(e):
            logger.info(f"Case {case} expectedly failed. "
                        "Required functionality is not yet "
                        f"implemented: {e}")
            pytest.skip(str(e))
        else:
            raise e


def _test_case(case, model_type):
    """Run a single PEtab test suite case"""
    case = petabtests.test_id_str(case)
    logger.debug(f"Case {case} [{model_type}]")

    # load
    if model_type == "sbml":
        case_dir = os.path.join(petabtests.SBML_DIR, case)
        # import petab problem
        yaml_file = os.path.join(case_dir, petabtests.problem_yaml_name(case))
        problem = petab.Problem.from_yaml(yaml_file)
    elif model_type == "pysb":
        import pysb
        pysb.SelfExporter.cleanup()
        pysb.SelfExporter.do_export = True
        case_dir = os.path.join(petabtests.PYSB_DIR, case)
        # import petab problem
        yaml_file = os.path.join(case_dir, petabtests.problem_yaml_name(case))
        problem = PysbPetabProblem.from_yaml(yaml_file)
    else:
        raise ValueError()

    # compile amici model
    model_output_dir = f'amici_models/model_{case}'
    model = import_petab_problem(
        problem, model_output_dir=model_output_dir,
        force_compile=True)

    # simulate
    ret = simulate_petab(problem, model, log_level=logging.DEBUG)

    rdatas = ret['rdatas']
    chi2 = sum(rdata['chi2'] for rdata in rdatas)
    llh = ret['llh']
    simulation_df = rdatas_to_measurement_df(rdatas, model,
                                             problem.measurement_df)
    petab.check_measurement_df(simulation_df, problem.observable_df)
    simulation_df = simulation_df.rename(
        columns={petab.MEASUREMENT: petab.SIMULATION})
    simulation_df[petab.TIME] = simulation_df[petab.TIME].astype(int)
    solution = petabtests.load_solution(case, model_type)
    gt_chi2 = solution[petabtests.CHI2]
    gt_llh = solution[petabtests.LLH]
    gt_simulation_dfs = solution[petabtests.SIMULATION_DFS]
    tol_chi2 = solution[petabtests.TOL_CHI2]
    tol_llh = solution[petabtests.TOL_LLH]
    tol_simulations = solution[petabtests.TOL_SIMULATIONS]

    chi2s_match = petabtests.evaluate_chi2(chi2, gt_chi2, tol_chi2)
    llhs_match = petabtests.evaluate_llh(llh, gt_llh, tol_llh)
    simulations_match = petabtests.evaluate_simulations(
        [simulation_df], gt_simulation_dfs, tol_simulations)

    logger.log(logging.DEBUG if chi2s_match else logging.ERROR,
               f"CHI2: simulated: {chi2}, expected: {gt_chi2},"
               f" match = {chi2s_match}")
    logger.log(logging.DEBUG if simulations_match else logging.ERROR,
               f"LLH: simulated: {llh}, expected: {gt_llh}, "
               f"match = {llhs_match}")
    logger.log(logging.DEBUG if simulations_match else logging.ERROR,
               f"Simulations: match = {simulations_match}")

    # FIXME case 7 fails due to #963
    if case not in ['0007', '0016']:
        check_derivatives(problem, model)

    # FIXME case 7 fails due to #963
    if not all([llhs_match, simulations_match]) \
            or (not chi2s_match and case not in ['0007', '0016']):
        # chi2s_match ignored until fixed in amici
        logger.error(f"Case {case} failed.")
        raise AssertionError(f"Case {case}: Test results do not match "
                             "expectations")

    logger.info(f"Case {case} passed.")


def check_derivatives(problem: petab.Problem, model: amici.Model) -> None:
    """Check derivatives using finite differences for all experimental
    conditions

    Arguments:
        problem: PEtab problem
        model: AMICI model matching ``problem``
    """
    problem_parameters = {t.Index: getattr(t, petab.NOMINAL_VALUE) for t in
                          problem.parameter_df.itertuples()}
    solver = model.getSolver()
    solver.setSensitivityMethod(amici.SensitivityMethod_forward)
    solver.setSensitivityOrder(amici.SensitivityOrder_first)
    # Required for case 9 to not fail in
    #  amici::NewtonSolver::computeNewtonSensis
    model.setSteadyStateSensitivityMode(
        SteadyStateSensitivityMode_simulationFSA)

    def assert_true(x):
        assert x

    for edata in create_parameterized_edatas(
            amici_model=model, petab_problem=problem,
            problem_parameters=problem_parameters):
        # check_derivatives does currently not support parameters in ExpData
        model.setParameters(edata.parameters)
        model.setParameterScale(edata.pscale)
        edata.parameters = []
        edata.pscale = amici.parameterScalingFromIntVector([])
        amici_check_derivatives(model, solver, edata, assert_true)


def run():
    """Run the full PEtab test suite"""

    n_success = 0
    n_skipped = 0
    for case in petabtests.CASES_LIST:
        try:
            test_case(case)
            n_success += 1
        except Skipped:
            n_skipped += 1
        except Exception as e:
            # run all despite failures
            logger.error(f"Case {case} failed.")
            logger.error(e)

    logger.info(f"{n_success} / {len(petabtests.CASES_LIST)} successful, "
                f"{n_skipped} skipped")
    if n_success != len(petabtests.CASES_LIST):
        sys.exit(1)


if __name__ == '__main__':
    run()
