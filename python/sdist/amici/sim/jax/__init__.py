"""Functionality for simulating JAX-based AMICI models."""

import petab.v1 as petabv1
import petab.v2 as petabv2

import pandas as pd
import jax.numpy as jnp

def add_default_experiment_names_to_v2_problem(petab_problem: petabv2.Problem):
    """Add default experiment names to PEtab v2 problem.

    Args:
        petab_problem: PEtab v2 problem to modify.
    """
    if not hasattr(petab_problem, "extensions_config"):
        petab_problem.extensions_config = {}

    petab_problem.visualization_df = None

    if petab_problem.condition_df is None:
        default_condition = petabv2.core.Condition(id="__default__", changes=[], conditionId="__default__")
        petab_problem.condition_tables[0].elements = [default_condition]

    if petab_problem.experiment_df is None or petab_problem.experiment_df.empty:
        condition_ids = petab_problem.condition_df[petabv2.C.CONDITION_ID].values
        condition_ids = [c for c in condition_ids if "preequilibration" not in c]
        default_experiment = petabv2.core.Experiment(
            id="__default__",
            periods=[
                petabv2.core.ExperimentPeriod(
                    time=0.0, 
                    condition_ids=condition_ids
                )
            ],
        )
        petab_problem.experiment_tables[0].elements = [default_experiment]

        measurement_tables = petab_problem.measurement_tables.copy()
        for mt in measurement_tables:
            for m in mt.elements:
                m.experiment_id = "__default__"

        petab_problem.measurement_tables = measurement_tables

    return petab_problem

def get_simulation_conditions_v2(petab_problem) -> pd.DataFrame:
    """Get simulation conditions from PEtab v2 measurement DataFrame.

    Returns:
        A pandas DataFrame mapping experiment_ids to condition ids.
    """
    experiment_df = petab_problem.experiment_df
    exps = {}
    for exp_id in experiment_df[petabv2.C.EXPERIMENT_ID].unique():
        exps[exp_id] = experiment_df[
            experiment_df[petabv2.C.EXPERIMENT_ID] == exp_id
        ][petabv2.C.CONDITION_ID].unique()

    experiment_df = experiment_df.rename(columns={"conditionId": "simulationConditionId"})
    experiment_df = experiment_df.drop(columns=[petabv2.C.TIME])
    return experiment_df

def _build_simulation_df_v2(problem, y, dyn_conditions):
    """Build petab simulation DataFrame of similation results from a PEtab v2 problem."""
    dfs = []
    for ic, sc in enumerate(dyn_conditions):
        experiment_id = _conditions_to_experiment_map(
            problem._petab_problem.experiment_df
        )[sc]

        if experiment_id == "__default__":
            experiment_id = jnp.nan

        obs = [
            problem.model.observable_ids[io]
            for io in problem._iys[ic, problem._ts_masks[ic, :]]
        ]
        t = jnp.concat(
            (
                problem._ts_dyn[ic, :],
                problem._ts_posteq[ic, :],
            )
        )
        df_sc = pd.DataFrame(
            {
                petabv2.C.MODEL_ID: [float("nan")] * len(t),
                petabv1.OBSERVABLE_ID: obs,
                petabv2.C.EXPERIMENT_ID: [experiment_id] * len(t),
                petabv1.TIME: t[problem._ts_masks[ic, :]],
                petabv1.SIMULATION: y[ic, problem._ts_masks[ic, :]],
            },
            index=problem._petab_measurement_indices[ic, :],
        )
        if (
            petabv1.OBSERVABLE_PARAMETERS
            in problem._petab_problem.measurement_df
        ):                
            df_sc[petabv2.C.OBSERVABLE_PARAMETERS] = (
                problem._petab_problem.measurement_df.query(
                    f"{petabv2.C.EXPERIMENT_ID} == '{experiment_id}'"
                )[petabv2.C.OBSERVABLE_PARAMETERS]
            )
        if petabv1.NOISE_PARAMETERS in problem._petab_problem.measurement_df:
            df_sc[petabv2.C.NOISE_PARAMETERS] = (
                problem._petab_problem.measurement_df.query(
                    f"{petabv2.C.EXPERIMENT_ID} == '{experiment_id}'"
                )[petabv2.C.NOISE_PARAMETERS]
            )
        dfs.append(df_sc)
    return pd.concat(dfs).sort_index()

def _conditions_to_experiment_map(experiment_df: pd.DataFrame) -> dict[str, str]:
    condition_to_experiment = {
        row.conditionId: row.experimentId
        for row in experiment_df.itertuples()
    }
    return condition_to_experiment

def _try_float(value):
    try:
        return float(value)
    except Exception as e:
        msg = str(e).lower()
        if isinstance(e, ValueError) and "could not convert" in msg:
            return value
        raise