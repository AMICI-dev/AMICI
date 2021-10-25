"""Tests for petab_objective.py."""

import copy

import amici.petab_objective
from amici.petab_objective import (
    LLH,
    SLLH,
)
import petab


"""
models with preequilibration
Weber_BMC2015       1
Isensee_JCB2018     1
Zheng_PNAS2012      1
Brannmark_JBC2010   1
Raimundez_PCB2020   4
"""


def test_simulate_petab_with_sensis():
    """Test simulation of PEtab problems with sensitivities."""
    model_id = 'Weber_BMC2015'
    petab_problem = \
        petab.Problem.from_yaml(f'benchmark_petab/{model_id}/{model_id}.yaml')
    amici_model = amici.petab_import.import_petab_problem(petab_problem)

    problem_parameters0 = dict(zip(
        petab_problem.x_ids,
        petab_problem.x_nominal,
    ))

    problem_parameters = {
        id: value
        for id, value in problem_parameters0.items()
        if id in amici_model.getParameterIds()
    }
    print(set(problem_parameters0).difference(problem_parameters))

    def simulate(
        vector,
        petab_problem=petab_problem,
        amici_model=amici_model,
    ):
        return amici.petab_objective.simulate_petab(
            petab_problem=petab_problem,
            amici_model=amici_model,
            problem_parameters=vector,
        )

    result = simulate(problem_parameters)

    epsilon = 1e-5
    finite_differences = {}
    for x_id in problem_parameters:
        backward_parameters = copy.deepcopy(problem_parameters)
        forward_parameters = copy.deepcopy(problem_parameters)

        backward_parameters[x_id] -= epsilon
        forward_parameters[x_id] += epsilon

        backward_llh = simulate(backward_parameters)[LLH]
        forward_llh = simulate(forward_parameters)[LLH]

        central_differences[x_id] = (
            (forward_llh - backward_llh)
            /
            (2 * epsilon)
        )

    breakpoint()
