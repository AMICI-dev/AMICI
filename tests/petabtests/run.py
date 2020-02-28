import os
import sys
import logging

import amici
from amici.petab_import import import_petab_problem
from amici.petab_objective import simulate_petab
import petab
import petabtests

logger = logging.getLogger(__name__)


def run():
    for case in petabtests.CASES_LIST:
        logger.info(f"Case {case}")

        # load
        case_dir = os.path.join(petabtests.CASES_DIR, case)

        # import petab problem
        yaml_file = os.path.join(case_dir, f'_{case}.yaml')
        problem = petab.Problem.from_yaml(yaml_file)

        # compile amici model
        model_output_dir = f'amici_models/model_{case}'
        model = import_petab_problem(
            problem, model_output_dir=model_output_dir)

        # simulate
        ret = simulate_petab(problem, model)
        print(ret['llh'])

if __name__ == '__main__':
    run()
