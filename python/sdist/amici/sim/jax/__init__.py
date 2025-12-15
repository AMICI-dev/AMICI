"""Functionality for simulating JAX-based AMICI models."""

import petab.v1 as petabv1
import petab.v2 as petabv2

def get_simulation_conditions_v2(petab_problem) -> dict[str, tuple[str, ...]]:
    """Get simulation conditions from PEtab v2 measurement DataFrame.

    Returns:
        A dictionary mapping simulation condition IDs to tuples of
        preequilibration and simulation condition IDs.
    """
    experiment_df = petab_problem.experiment_df
    exps = {}
    for exp_id in experiment_df[petabv2.C.EXPERIMENT_ID].unique():
        exps[exp_id] = experiment_df[
            experiment_df[petabv2.C.EXPERIMENT_ID] == exp_id
        ][petabv2.C.CONDITION_ID].unique()

    experiment_df = experiment_df.rename(columns={"conditionId": "simulationConditionId"})
    return experiment_df

def reformat_for_v2(petab_problem: petabv2.Problem):
    """Reformat PEtab v2 problem for parameter mapping.

    Args:
        petab_problem: PEtab v2 problem to reformat.
    """
    measurement_df = petab_problem.measurement_df.rename(
        columns={petabv2.C.EXPERIMENT_ID: petabv1.C.SIMULATION_CONDITION_ID}
    ) if isinstance(petab_problem, petabv2.Problem) else petab_problem.measurement_df

    return measurement_df