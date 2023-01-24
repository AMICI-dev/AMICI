#!/usr/bin/python
import os
from contextlib import suppress

import benchmark_models_petab
import numpy as np

import amici
from amici.petab_import import import_petab_problem
from amici.petab_objective import RDATAS, simulate_petab

with suppress(KeyError):
    del os.environ["AMICI_EXPERIMENTAL_SBML_NONCONST_CLS"]
    del os.environ["AMICI_EXTRACT_CSE"]
#petab_problem = benchmark_models_petab.get_problem("Blasi_CellSystems2016")
petab_problem = benchmark_models_petab.get_problem("Brannmark_JBC2010")
amici_model = import_petab_problem(petab_problem, verbose=False, force_compile=True, compute_conservation_laws=False)
amici_solver = amici_model.getSolver()
#amici_solver.setSensitivityOrder(amici.SensitivityOrder.first)
#amici_solver.setSensitivityMethod(amici.SensitivityMethod.forward)

for i in range(1850, 4000):
    print(f"{'-' * 10} {i} {'-' * 10}")

    np.random.seed(i)
    problem_parameters = dict(
        zip(
            petab_problem.x_free_ids,
            petab_problem.sample_parameter_startpoints(n_starts=1)[0],
        )
    )
    res = simulate_petab(
        petab_problem=petab_problem,
        amici_model=amici_model,
        problem_parameters=problem_parameters,
        scaled_parameters=True,
        solver=amici_solver,
    )
    if not all(rdata.status == amici.AMICI_SUCCESS for rdata in res[RDATAS]):
        print([amici.simulation_status_to_str(rdata.status) for rdata in res[RDATAS]])
#
# print(f"{res[RDATAS][0].status=}")
# print(f"{res[RDATAS][0].posteq_status=}")
# print(f"{res[RDATAS][0].preeq_status=}")
# print(f"{res[RDATAS][0].numsteps=}")
# print(f"{res[RDATAS][0].numstepsB=}")
# print(f"{res[RDATAS][0].x0=}")
# print(f"{res[RDATAS][0].x[-1]=}")
# assert res[RDATAS][0].status != amici.AMICI_SUCCESS
