"""Functionality for simulating JAX-based AMICI models."""

import petab.v1 as petabv1
import petab.v2 as petabv2

import pandas as pd
import jax.numpy as jnp

import libsbml

def add_events_to_sbml(petab_problem: petabv2.Problem):
    """Add events to SBML models based on PEtab v2 experiments."""
    ok = petab_problem.model.sbml_document.setLevelAndVersion(3, 1)
    if not ok:
        raise RuntimeError("Could not set SBML level to 3.")

    for experiment in petab_problem.experiments:
        for i_period, period in enumerate(experiment.sorted_periods):
            # TODO: use function from libpetab.ExperimentsToSbmlConverter if possible
            ev = petab_problem.model.sbml_model.createEvent()
            ev.setId(f"event_exp{experiment.id}_period{i_period}")
            ev.setUseValuesFromTriggerTime(True)
            trigger = ev.createTrigger()
            trigger.setInitialValue(False)
            trigger.setPersistent(True)
            trig_expr = libsbml.parseL3Formula(f"time >= {period.time}")
            trigger.setMath(trig_expr)

            for condition_id in period.condition_ids:
                cond = petab_problem.condition_df.set_index(petabv2.C.CONDITION_ID)
                target_id = cond.at[condition_id, petabv2.C.TARGET_ID]
                target_value = cond.at[condition_id, petabv2.C.TARGET_VALUE]

                if not jnp.issubdtype(type(target_value), jnp.number):
                    continue

                # Don't handle compartment size as event?
                if target_id in [
                    c.id for c in petab_problem.model.sbml_model.getListOfCompartments()
                ]:
                    for comp in petab_problem.model.sbml_model.getListOfCompartments():
                        if comp.id == target_id:
                            comp.setSize(float(target_value))
                    continue
                
                assignment = ev.createEventAssignment()
                assignment.setVariable(target_id)
                assignment_expr = libsbml.parseL3Formula(str(target_value))
                assignment.setMath(assignment_expr)

                _set_assignment_entities_false(
                    petab_problem.model.sbml_model, target_id
                )

    return petab_problem

def _set_assignment_entities_false(sbml_model, target_id: str):
    """Set a parameter used in event assignments to be non-constant."""
    for param in sbml_model.getListOfParameters():
        if param.id == target_id:
            param.setConstant(False)

    for species in sbml_model.getListOfSpecies():
        if species.id == target_id:
            species.setConstant(False)

def reformat_petab_v2_to_v1(petab_problem: petabv2.Problem) -> petabv1.Problem:
    """Reformat a PEtab v2 problem to PEtab v1 format.

    Args:
        petab_problem: PEtab v2 problem to reformat.

    Returns:
        A PEtab v1 problem.
    """
    if not hasattr(petab_problem, "extensions_config"):
        petab_problem.extensions_config = {}

    petab_problem.visualization_df = None

    if petab_problem.condition_df is None:
        default_condition = petabv2.core.Condition(id="__default__", changes=[], conditionId="__default__")
        petab_problem.condition_tables[0].elements = [default_condition]

    if petab_problem.experiment_df is None:
        default_experiment = petabv2.core.Experiment(
            id="__default__",
            periods=[petabv2.core.ExperimentPeriod(time=0.0, condition_ids=["__default__"])],
        )
        petab_problem.experiment_tables[0].elements = [default_experiment]

        measurement_tables = petab_problem.measurement_tables.copy()
        for mt in measurement_tables:
            for m in mt.elements:
                m.experiment_id = "__default__"

        petab_problem.measurement_tables = measurement_tables

    return petab_problem

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
    experiment_df = experiment_df.drop(columns=[petabv2.C.TIME])
    return experiment_df

def reformat_for_v2(petab_problem: petabv2.Problem):
    """Reformat PEtab v2 problem for parameter mapping.

    Args:
        petab_problem: PEtab v2 problem to reformat.
    """
    measurement_df = petab_problem.measurement_df.rename(
        columns={petabv2.C.EXPERIMENT_ID: petabv1.C.SIMULATION_CONDITION_ID}
    ) if isinstance(petab_problem, petabv2.Problem) else petab_problem.measurement_df

    condition_df = _long_to_wide(petab_problem.condition_df)

    return measurement_df, condition_df

def _long_to_wide(condition_df):
    wide_df = condition_df.pivot_table(
        index="conditionId",
        columns="targetId",
        values="targetValue",
        aggfunc="first"
    )
    wide_df.columns.name = None
    return wide_df

def _build_simulation_df_v2(problem, y, dyn_conditions):
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

def fixup_v2_parameter_mapping(prelim_mapping_for_condition, petab_problem):
    # No preequilibration in PEtab v2
    # Replace preeq mappings with sim mappings to satisfy existing code
    (
        _,
        condition_map_sim,
        _,
        condition_scale_map_sim,
    ) = prelim_mapping_for_condition

    for k, _ in condition_map_sim.items():
        if k not in condition_scale_map_sim:
            condition_scale_map_sim[k] = "lin"

    for placeholder_col, param_col in (
        (petabv2.C.OBSERVABLE_PLACEHOLDERS, petabv2.C.OBSERVABLE_PARAMETERS),
        (petabv2.C.NOISE_PLACEHOLDERS, petabv2.C.NOISE_PARAMETERS),
    ):
        placeholders = petab_problem.observable_df[
            placeholder_col
        ].unique()

        for placeholders in placeholders:
            placeholder_list = placeholders.split(";")
            params_list = petab_problem.measurement_df[param_col][0].split(";")
            for i, p in enumerate(placeholder_list):
                if p == '':
                    continue
                condition_map_sim[p] = _try_float(params_list[i])

    return (
        condition_map_sim,
        condition_map_sim,
        condition_scale_map_sim,
        condition_scale_map_sim,
    )

def _try_float(value):
    try:
        return float(value)
    except Exception as e:
        msg = str(e).lower()
        if isinstance(e, ValueError) and "could not convert" in msg:
            return value
        raise