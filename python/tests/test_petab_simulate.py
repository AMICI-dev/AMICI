"""Test for petab.yaml"""

import os
import petab

from amici.petab_simulate import PetabSimulator

def test_petab_simulate():
    model_name = 'conversion_reaction'
    petab_yaml_filepath = os.path.join(
        'petab_test_models',
        model_name,
        model_name + '.yaml',
    )
    petab_problem = petab.Problem.from_yaml(petab_yaml_filepath)
    simulator = PetabSimulator(petab_problem)
    simulator.simulate(noise=True)
    simulator.clean_working_dir()

