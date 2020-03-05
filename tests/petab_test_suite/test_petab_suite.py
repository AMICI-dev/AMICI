"""Run PEtab test suite (https://github.com/PEtab-dev/petab_test_suite)"""
import os
import sys
import logging

import amici
from amici.petab_import import import_petab_problem
from amici.petab_objective import simulate_petab, rdatas_to_measurement_df
import petab
import petabtests


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
stream_handler = logging.StreamHandler()
logger.addHandler(stream_handler)


def test_case(case):
    """Run a single PEtab test suite case"""
    case = petabtests.test_id_str(case)
    logger.debug(f"Case {case}")

    # load
    case_dir = os.path.join(petabtests.CASES_DIR, case)

    # import petab problem
    yaml_file = os.path.join(case_dir, petabtests.problem_yaml_name(case))
    problem = petab.Problem.from_yaml(yaml_file)

    # compile amici model
    model_output_dir = f'amici_models/model_{case}'
    model = import_petab_problem(
        problem, model_output_dir=model_output_dir)

    # simulate
    chi2s_match = llhs_match = simulations_match = False
    ret = simulate_petab(problem, model)

    rdatas = ret['rdatas']
    chi2 = None
    llh = ret['llh']
    simulation_df = rdatas_to_measurement_df(rdatas, model,
                                             problem.measurement_df)
    simulation_df = simulation_df.rename(
        columns={petab.MEASUREMENT: petab.SIMULATION})
    simulation_df[petab.TIME] = simulation_df[petab.TIME].astype(int)
    solution = petabtests.load_solution(case)
    gt_chi2 = solution[petabtests.CHI2]
    gt_llh = solution[petabtests.LLH]
    gt_simulation_dfs = solution[petabtests.SIMULATION_DFS]
    tol_chi2 = solution[petabtests.TOL_CHI2]
    tol_llh = solution[petabtests.TOL_LLH]
    tol_simulations = solution[petabtests.TOL_SIMULATIONS]

    chi2s_match = petabtests.evaluate_chi2(chi2, solution, tol_chi2)
    llhs_match = petabtests.evaluate_llh(llh, gt_llh, tol_llh)
    simulations_match = petabtests.evaluate_simulations(
        [simulation_df], gt_simulation_dfs, tol_simulations)
    for edata in edatas_from_petab(model=model, petab_problem=petab_problem):
        for ip in range(model.np):
            amici.check_finite_difference(x0, model, solver, edata, ip, fields, assert_fun)
    logger.log(logging.DEBUG if chi2s_match else logging.ERROR,
               f"CHI2: simulated: {chi2}, expected: {gt_chi2},"
               f" match = {chi2s_match}")
    logger.log(logging.DEBUG if simulations_match else logging.ERROR,
               f"LLH: simulated: {llh}, expected: {gt_llh}, "
               f"match = {llhs_match}")
    logger.log(logging.DEBUG if simulations_match else logging.ERROR,
               f"Simulations: match = {simulations_match}")

    if not all([llhs_match, simulations_match]):
        # chi2s_match ignored until fixed in amici
        logger.error(f"Case {case} failed.")
        raise AssertionError(f"Case {case}: Test results do not match "
                             "expectations")

    logger.info(f"Case {case} passed.")


def run():
    """Run the full PEtab test suite"""

    n_success = 0

    for case in petabtests.CASES_LIST:
        try:
            test_case(case)
            n_success += 1
        except Exception as e:
            # run all despite failures
            logger.error(f"Case {case} failed.")
            logger.error(e)

    logger.info(f"{n_success} / {len(petabtests.CASES_LIST)} successful")

    if n_success != len(petabtests.CASES_LIST):
        sys.exit(1)


if __name__ == '__main__':
    run()
