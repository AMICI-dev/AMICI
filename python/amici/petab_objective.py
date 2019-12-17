"""Functionality related to evaluating the objective function for a PEtab
problem"""

import numbers
import numpy as np
import pandas as pd

from typing import List, Sequence
import amici
import petab


def edatas_from_petab(
        model: amici.Model, measurement_df: pd.DataFrame,
        condition_df: pd.DataFrame,
        simulation_conditions=None) -> List[amici.ExpData()]:
    """
    Create list of amici.ExpData objects for PEtab problem.

    Arguments:
        model:
            AMICI model
        measurement_df:
            PEtab measurement table
        condition_df:
            PEtab condition table
        simulation_conditions:
            Result of petab.get_simulation_conditions. Can be provided to save
            time if this has be obtained before.

    Returns:
        List with one ExpData per simulation condition.
    """

    condition_df = condition_df.reset_index()

    # number of amici simulations will be number of unique
    # (preequilibrationConditionId, simulationConditionId) pairs.
    # Can be improved by checking for identical condition vectors.
    if simulation_conditions is None:
        simulation_conditions = petab.get_simulation_conditions(
            measurement_df)

    observable_ids = model.getObservableIds()

    fixed_parameter_ids = model.getFixedParameterIds()

    edatas = []
    for _, condition in simulation_conditions.iterrows():
        # amici.ExpData for each simulation

        # extract rows for condition
        df_for_condition = petab.get_rows_for_condition(
            measurement_df, condition)

        # make list of all timepoints for which measurements exist
        timepoints = sorted(
            df_for_condition.time.unique().astype(float))

        # init edata object
        edata = amici.ExpData(model)

        # find replicate numbers of time points
        timepoints_w_reps = []
        for time in timepoints:
            # subselect for time
            df_for_time = df_for_condition[df_for_condition.time == time]
            # rep number is maximum over rep numbers for observables
            n_reps = max(df_for_time.groupby(
                ['observableId', 'time']).size())
            # append time point n_rep times
            timepoints_w_reps.extend([time] * n_reps)

        # set time points in edata
        edata.setTimepoints(timepoints_w_reps)

        # handle fixed parameters
        _fixed_parameters_to_edata(edata, condition_df,
                                   fixed_parameter_ids, condition)

        # prepare measurement matrix
        y = np.full(shape=(edata.nt(), edata.nytrue()), fill_value=np.nan)
        # prepare sigma matrix
        sigma_y = y.copy()

        # add measurements and sigmas
        # iterate over time points
        for time in timepoints:
            # subselect for time
            df_for_time = df_for_condition[df_for_condition.time == time]
            time_ix_0 = timepoints_w_reps.index(time)

            # remember used time indices for each observable
            time_ix_for_obs_ix = {}

            # iterate over measurements
            for _, measurement in df_for_time.iterrows():
                # extract observable index
                observable_ix = observable_ids.index(
                    f'observable_{measurement.observableId}')

                # update time index for observable
                if observable_ix in time_ix_for_obs_ix:
                    time_ix_for_obs_ix[observable_ix] += 1
                else:
                    time_ix_for_obs_ix[observable_ix] = time_ix_0

                # fill observable and possibly noise parameter
                y[time_ix_for_obs_ix[observable_ix],
                  observable_ix] = measurement.measurement
                if isinstance(measurement.noiseParameters, numbers.Number):
                    sigma_y[time_ix_for_obs_ix[observable_ix],
                            observable_ix] = measurement.noiseParameters

        # fill measurements and sigmas into edata
        edata.setObservedData(y.flatten())
        edata.setObservedDataStdDev(sigma_y.flatten())

        # append edata to edatas list
        edatas.append(edata)

    return edatas


def _fixed_parameters_to_edata(
        edata: amici.ExpData, condition_df: pd.DataFrame,
        fixed_parameter_ids: Sequence[str], condition) -> None:
    """
    Set fixed parameters for a given simulation condition to the corresponding
    ExpData.

    Parameters:
        edata:
            Current edata.

        condition_df:
            The conditions table.

        fixed_parameter_ids:
            Ids of parameters that are to be considered constant (in correct
            AMICI order).

        condition:
            The current condition, as created by
            petab.get_simulation_conditions.
    """

    if len(fixed_parameter_ids) == 0:
        # nothing to be done
        return

    # find fixed parameter values
    fixed_parameter_vals = condition_df.loc[
        condition_df.conditionId == condition.simulationConditionId,
        fixed_parameter_ids].values
    # fill into edata
    edata.fixedParameters = fixed_parameter_vals.astype(
        float).flatten()

    # same for preequilibration if necessary
    if all(('preequilibrationConditionId' in condition,
            condition.preequilibrationConditionId,
            'preequilibrationConditionId' in condition_df)):
        fixed_preequilibration_parameter_vals = condition_df.loc[
            condition_df.conditionId == condition.preequilibrationConditionId,
            fixed_parameter_ids].values
        edata.fixedParametersPreequilibration = \
            fixed_preequilibration_parameter_vals.astype(float) \
                                                 .flatten()
