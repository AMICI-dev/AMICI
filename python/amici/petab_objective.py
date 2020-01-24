"""Functionality related to running simulations or evaluating the objective
function as defined by a PEtab problem"""

import copy
import logging
import numbers
from typing import List, Sequence, Optional, Dict, Tuple

import amici
import numpy as np
import pandas as pd
import petab
from amici.logging import get_logger
from petab.C import *


logger = get_logger(__name__)


def simulate_petab(petab_problem: petab.Problem,
                   amici_model: amici.Model,
                   solver: Optional[amici.Solver] = None,
                   problem_parameters: Optional[Dict[str, float]] = None,
                   log_level: int = logging.WARNING
                   ) -> Tuple[float, List[amici.ReturnData]]:
    """Simulate PEtab model

    Arguments:
        petab_problem: PEtab problem to work on
        amici_model: AMICI Model assumed to be compatible with
            ``petab_problem``
        solver: An AMICI solver. Will use default options if None.
        problem_parameters: Run simulation with these parameters. If None,
            PEtab `nominalValues` will be used). To be provided as dict,
            mapping PEtab problem parameters to SBML IDs.

    Returns:
        Tuple of cost function value and a list of `ReturnData`s, corresponding
        to the different simulation conditions. For ordering of simulation
        conditions, see
        `petab.Problem.get_simulation_conditions_from_measurement_df`.
    """
    logger.setLevel(log_level)

    if solver is None:
        solver = amici_model.getSolver()

    if problem_parameters is None:
        # Use PEtab nominal values as default
        problem_parameters = {t.Index: getattr(t, NOMINAL_VALUE) for t in
                              petab_problem.parameter_df.itertuples()}


    edatas = edatas_from_petab(model=amici_model,
                               petab_problem=petab_problem,
                               problem_parameters=problem_parameters)


    # Simulate
    rdatas = amici.runAmiciSimulations(amici_model, solver, edata_list=edatas)

    # Compute total llh
    llh = sum(rdata['llh'] for rdata in rdatas)

    # log results
    sim_cond = petab_problem.get_simulation_conditions_from_measurement_df()

    for i, rdata in enumerate(rdatas):
        logger.debug(f"Condition: {sim_cond.iloc[i, :].values}, status: "
                     f"{rdata['status']}, llh: {rdata['llh']}")

    return llh, rdatas


def edatas_from_petab(
        model: amici.Model, petab_problem: petab.Problem,
        problem_parameters,
        simulation_conditions=None,
        parameter_mapping: List[petab.ParMappingDictTuple] = None,
        parameter_scale_mapping: List[petab.ScaleMappingDictTuple] = None
) -> List[amici.ExpData]:
    """
    Create list of ``amici.ExpData`` objects for PEtab problem.

    Sets timepoints, fixed parameters (including preequilibration),
    non-fixed parameters, and observed data and sigmas.

    Arguments:
        model:
            AMICI model.
        petab_problem:
            PEtab problem
        simulation_conditions:
            Result of petab.get_simulation_conditions. Can be provided to save
            time if this has be obtained before.
        parameter_mapping:
            Optional precomputed mapping for efficiency.

    Returns:
        List with one ``ExpData`` per simulation condition.
    """

    # TODO Rem condition_df = petab_problem.condition_df.reset_index()

    # number of amici simulations will be number of unique
    # (preequilibrationConditionId, simulationConditionId) pairs.
    # Can be optimized by checking for identical condition vectors.
    if simulation_conditions is None:
        simulation_conditions = \
            petab_problem.get_simulation_conditions_from_measurement_df()

    if parameter_mapping is None:
        parameter_mapping = \
            petab_problem.get_optimization_to_simulation_parameter_mapping(
                warn_unmapped=False)
    if parameter_scale_mapping is None:
        parameter_scale_mapping = \
            petab.get_optimization_to_simulation_scale_mapping(
                mapping_par_opt_to_par_sim=parameter_mapping,
                parameter_df=petab_problem.parameter_df,
                measurement_df=petab_problem.measurement_df
            )

    observable_ids = model.getObservableIds()

    fixed_parameter_ids = model.getFixedParameterIds()

    edatas = []
    for (_, condition), cur_parameter_mapping, cur_parameter_scale_mapping \
            in zip(simulation_conditions.iterrows(),
                   parameter_mapping, parameter_scale_mapping):
        # Create amici.ExpData for each simulation
        edata = get_edata_for_condition(
            condition=condition, amici_model=model, petab_problem=petab_problem,
            problem_parameters=problem_parameters,
            fixed_parameter_ids=fixed_parameter_ids,
            observable_ids=observable_ids,
            parameter_mapping=cur_parameter_mapping,
            parameter_scale_mapping=cur_parameter_scale_mapping
        )
        edatas.append(edata)

    return edatas


def get_edata_for_condition(
        condition,
        problem_parameters,
        amici_model: amici.Model,
        petab_problem: petab.Problem,
        fixed_parameter_ids,
        observable_ids,
        parameter_mapping: Optional[petab.ParMappingDictTuple] = None,
        parameter_scale_mapping: Optional[petab.ScaleMappingDictTuple] = None
) -> amici.ExpData:
    """Get ``amici.ExpData`` for the given PEtab condition

    Sets timepoints, fixed parameters (including preequilibration),
    non-fixed parameters, and observed data and sigmas.

    Arguments:
        condition:
            pandas.DataFrame row with preequilibrationConditionId and
            simulationConditionId.
        problem_parameters:
            PEtab problem parameters as parameterId=>value dict. Only
            parameters included here will be set. Remaining parameters will
            be used currently set in `amici_model`.
        amici_model:
            AMICI model
        measurement_df:
            PEtab measurement table
        condition_df:
            PEtab condition table

    Returns:
        ExpData instance
    """

    # extract measurement table rows for condition
    measurement_df = petab.get_rows_for_condition(
        measurement_df=petab_problem.measurement_df, condition=condition)

    # create sorted list of all timepoints for which measurements exist
    timepoints = sorted(measurement_df[TIME].unique().astype(float))

    # find replicate numbers of time points
    timepoints_w_reps = _get_timepoints_with_replicates(
        timepoints=timepoints, df_for_condition=measurement_df)

    edata = amici.ExpData(amici_model)

    ##########################################################################
    # variable parameters

    # get mapping if required
    if parameter_mapping is None:
        # TODO petab.get_parameter_mapping_for_condition
        raise NotImplementedError()

    if parameter_scale_mapping is None:
        # TODO petab.get_parameter_scale_mapping_for_condition
        raise NotImplementedError()

    condition_map_preeq, condition_map_sim = parameter_mapping
    condition_scale_map_preeq, condition_scale_map_sim = \
        parameter_scale_mapping

    # petab parameter mapping may contain parameter_ids as values, these *must*
    # be replaced
    def _get_par(model_par, value):
        """Replace parameter IDs in mapping dicts by values from
        problem_parameters where necessary"""
        if isinstance(value, str):
            return problem_parameters[value]
        if  model_par in problem_parameters:
            return problem_parameters[model_par]
        return model_par

    condition_map_sim = {key: _get_par(key, val)
                         for key, val in condition_map_sim.items()}
    condition_map_preeq = {key: _get_par(key, val)
                           for key, val in condition_map_preeq.items()}

    # separate fixed and variable AMICI parameters, because we may have
    # different fixed parameters for preeq and sim condition, but we cannot
    # have different variable parameters. without splitting,
    # merge_preeq_and_sim_pars_condition below may fail.
    variable_par_ids = amici_model.getParameterIds()
    fixed_par_ids = amici_model.getFixedParameterIds()

    condition_map_preeq_var = {
        key: val for (key, val) in condition_map_preeq.items()
        if key in variable_par_ids}
    condition_map_preeq_fix = {
        key: val for (key, val) in condition_map_preeq.items()
        if key in fixed_par_ids}

    condition_scale_map_preeq_var = {
        key: val for (key, val) in condition_scale_map_preeq.items()
        if key in variable_par_ids}
    condition_scale_map_preeq_fix = {
        key: val for (key, val) in condition_scale_map_preeq.items()
        if key in fixed_par_ids}

    condition_map_sim_var = {
        key: val for (key, val) in condition_map_sim.items()
        if key in variable_par_ids}
    condition_map_sim_fix = {
        key: val for (key, val) in condition_map_sim.items()
        if key in fixed_par_ids}

    condition_scale_map_sim_var = {
        key: val for (key, val) in condition_scale_map_sim.items()
        if key in variable_par_ids}
    condition_scale_map_sim_fix = {
        key: val for (key, val) in condition_scale_map_sim.items()
        if key in fixed_par_ids}

    petab.merge_preeq_and_sim_pars_condition(
        condition_map_preeq_var, condition_map_sim_var,
        condition_scale_map_preeq_var, condition_scale_map_sim_var,
        condition)

    # variable parameters default to those set on the model, and are
    # overwritten by those specified in problem_parameters
    parameters = [condition_map_sim_var[par_id]
                  for par_id in amici_model.getParameterIds()]

    edata.parameters = parameters

    ##########################################################################
    # parameter scale
    # TODO: argument for problem_parameter scale, for now assume linear scale
    edata.pscale = amici.parameterScalingFromIntVector(
        [amici.ParameterScaling_none] * len(parameters))
    # pscale = [petab_scale_to_amici_scale(condition_scale_map_sim[model_par])
    #          for model_par in amici_model.getParameterIds()]

    ##########################################################################
    # timepoints
    edata.setTimepoints(timepoints_w_reps)

    ##########################################################################
    # initial states
    # initial states are set in the model. if they were overwritten in the
    # PEtab condition table, they are now fixed model parameters and handled
    # below

    ##########################################################################
    # fixed parameters preequilibration
    fixed_pars_preeq = [condition_map_preeq_fix[par_id]
                      for par_id in amici_model.getFixedParameterIds()]
    edata.fixedParameters = fixed_pars_preeq

    ##########################################################################
    # fixed parameters simulation
    fixed_pars_sim = [condition_map_sim_fix[par_id]
                      for par_id in amici_model.getFixedParameterIds()]
    edata.fixedParameters = fixed_pars_sim

    ##########################################################################
    # measurements and sigmas
    y, sigma_y = _get_measurements_and_sigmas(
        measurement_df, timepoints, timepoints_w_reps, observable_ids, edata)
    edata.setObservedData(y.flatten())
    edata.setObservedDataStdDev(sigma_y.flatten())

    return edata


def _get_timepoints_with_replicates(timepoints, df_for_condition):
    """
    # TODO
    :param timepoints:
    :param df_for_condition:
    :return:
    """
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
    return timepoints_w_reps


def _get_measurements_and_sigmas(df_for_condition,
                                 timepoints, timepoints_w_reps,
                                    observable_ids, edata
                                 ) -> Tuple[np.array, np.array]:
    """Get measurements and sigmas

    Generate arrays with measurements and sigmas in AMICI format from a
    PEtab measurement table subset for a single condition.

    # TODO
    """
    # prepare measurement matrix
    y = np.full(shape=(edata.nt(), edata.nytrue()), fill_value=np.nan)
    # prepare sigma matrix
    sigma_y = y.copy()

    for time in timepoints:
        # subselect for time
        df_for_time = df_for_condition[df_for_condition[TIME] == time]
        time_ix_0 = timepoints_w_reps.index(time)

        # remember used time indices for each observable
        time_ix_for_obs_ix = {}

        # iterate over measurements
        for _, measurement in df_for_time.iterrows():
            # extract observable index
            observable_ix = observable_ids.index(
                f'observable_{measurement[OBSERVABLE_ID]}')

            # update time index for observable
            if observable_ix in time_ix_for_obs_ix:
                time_ix_for_obs_ix[observable_ix] += 1
            else:
                time_ix_for_obs_ix[observable_ix] = time_ix_0

            # fill observable and possibly noise parameter
            y[time_ix_for_obs_ix[observable_ix],
              observable_ix] = measurement[MEASUREMENT]
            if isinstance(measurement[NOISE_PARAMETERS], numbers.Number):
                sigma_y[time_ix_for_obs_ix[observable_ix],
                        observable_ix] = measurement[NOISE_PARAMETERS]
    return y, sigma_y

'''
def _fixed_parameters_to_edata(
        edata: amici.ExpData, condition_df: pd.DataFrame,
        fixed_parameter_ids: Sequence[str], condition) -> None:
    """
    Apply fixed parameters for a given simulation condition to the
    corresponding ExpData.

    Parameters:
        edata:
            ExpData to set fixed parameters on.
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
    if (PREEQUILIBRATION_CONDITION_ID in condition
           and condition[PREEQUILIBRATION_CONDITION_ID]):
        fixed_preequilibration_parameter_vals = \
            condition_df.loc[ \
                condition_df[CONDITION_ID] \
                == condition[PREEQUILIBRATION_CONDITION_ID],
                fixed_parameter_ids].values
        edata.fixedParametersPreequilibration = \
            fixed_preequilibration_parameter_vals.astype(float) \
                                                 .flatten()
'''

def rdatas_to_measurement_df(
        rdatas: List[amici.ReturnData], model: amici.Model,
        measurement_df: pd.DataFrame) -> pd.DataFrame:
    """
    Create a measurement dataframe in the PEtab format from the passed
    `rdatas` and own information.

    Parameters:
        rdatas:
            A list of rdatas with the ordering of
            `petab.get_simulation_conditions`.
        model:
            AMICI model used to generate `rdatas`.
        measurement_df:
            PEtab measurement table used to generate `rdatas`.

    Returns:
        A dataframe built from the rdatas in the format of `measurement_df`.
    """

    df = pd.DataFrame(columns=list(measurement_df.columns))

    simulation_conditions = petab.get_simulation_conditions(
        measurement_df)

    observable_ids = model.getObservableIds()

    # iterate over conditions
    for data_idx, condition in simulation_conditions.iterrows():
        # current rdata
        rdata = rdatas[data_idx]
        # current simulation matrix
        y = rdata['y']
        # time array used in rdata
        t = list(rdata['t'])

        # extract rows for condition
        cur_measurement_df = petab.get_rows_for_condition(
            measurement_df, condition)

        # iterate over entries for the given condition
        # note: this way we only generate a dataframe entry for every
        # row that existed in the original dataframe. if we want to
        # e.g. have also timepoints non-existent in the original file,
        # we need to instead iterate over the rdata['y'] entries
        for _, row in cur_measurement_df.iterrows():
            # copy row
            row_sim = copy.deepcopy(row)

            # extract simulated measurement value
            timepoint_idx = t.index(row[TIME])
            observable_idx = observable_ids.index(
                "observable_" + row[OBSERVABLE_ID])
            measurement_sim = y[timepoint_idx, observable_idx]

            # change measurement entry
            row_sim.measurement = measurement_sim

            # append to dataframe
            df = df.append(row_sim, ignore_index=True)

    return df
