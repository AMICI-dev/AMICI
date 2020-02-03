"""Functionality related to running simulations or evaluating the objective
function as defined by a PEtab problem"""

import copy
import logging
import numbers
from typing import (List, Sequence, Optional, Dict, Tuple, Union, Any,
                    Collection, Iterator)

import amici
import numpy as np
import pandas as pd
import petab
from amici.logging import get_logger, log_execution_time
from petab.C import *


logger = get_logger(__name__)


@log_execution_time('Simulating PEtab model', logger)
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
        log_level: Log level, see `logging` module
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

    # Generate ExpData with all condition-specific information
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
        model: amici.Model,
        petab_problem: petab.Problem,
        problem_parameters: Dict[str, numbers.Number],
        simulation_conditions: Union[pd.DataFrame, Dict] = None,
        parameter_mapping: List[petab.ParMappingDictTuple] = None,
        parameter_scale_mapping: List[petab.ScaleMappingDictTuple] = None,
        scaled_parameters: Optional[bool] = False
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
        problem_parameters:
            Dictionary mapping parameter names of the PEtab problem to
            parameter values
        simulation_conditions:
            Result of petab.get_simulation_conditions. Can be provided to save
            time if this has be obtained before.
        parameter_mapping:
            Optional precomputed PEtab parameter mapping for efficiency.
        parameter_scale_mapping:
            Optional precomputed PEtab parameter scale mapping for efficiency.
        scaled_parameters:
            If True, problem_parameters are assumed to be on the scale provided
            in the PEtab parameter table and will be unscaled. If False, they
            are assumed to be in linear scale.

    Returns:
        List with one ``ExpData`` per simulation condition.
    """

    # number of amici simulations will be number of unique
    # (preequilibrationConditionId, simulationConditionId) pairs.
    # Can be optimized by checking for identical condition vectors.
    if simulation_conditions is None:
        simulation_conditions = \
            petab_problem.get_simulation_conditions_from_measurement_df()

    # Get parameter mapping if not user-provided
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

    logger.debug(f"Problem parameters: {problem_parameters}")

    edatas = []
    for (_, condition), cur_parameter_mapping, cur_parameter_scale_mapping \
            in zip(simulation_conditions.iterrows(),
                   parameter_mapping, parameter_scale_mapping):
        # Create amici.ExpData for each simulation
        edata = get_edata_for_condition(
            condition=condition, amici_model=model, petab_problem=petab_problem,
            problem_parameters=problem_parameters,
            observable_ids=observable_ids,
            parameter_mapping=cur_parameter_mapping,
            parameter_scale_mapping=cur_parameter_scale_mapping,
            scaled_parameters=scaled_parameters
        )
        edatas.append(edata)

    return edatas


def subset_dict(full: Dict[Any, Any],
                *args: Collection[Any]) -> Iterator[Dict[Any, Any]]:
    """Get subset of dictionary based on provides keys

    Arguments:
        full: Dictionary to subset
        *args: Collections of keys to be contained in the different subsets
    """

    for keys in args:
        yield {key: val for (key, val) in full.items() if key in keys}


def get_edata_for_condition(
        condition: Dict,
        problem_parameters: Dict[str, numbers.Number],
        amici_model: amici.Model,
        petab_problem: petab.Problem,
        observable_ids: List[str],
        parameter_mapping: Optional[petab.ParMappingDictTuple] = None,
        parameter_scale_mapping: Optional[petab.ScaleMappingDictTuple] = None,
        scaled_parameters: Optional[bool] = False
) -> amici.ExpData:
    """Get ``amici.ExpData`` for the given PEtab condition

    Sets timepoints, fixed parameters (including preequilibration),
    variable parameters, and observed data and sigmas.

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
        petab_problem:
            Underlying PEtab problem
        observable_ids:
            List of observable IDs
        parameter_mapping:
            PEtab parameter mapping for current condition
        parameter_scale_mapping:
            PEtab parameter scale mapping for current condition
        scaled_parameters:
            If True, problem_parameters are assumed to be on the scale provided
            in the PEtab parameter table and will be unscaled. If False, they
            are assumed to be in linear scale.

    Returns:
        ExpData instance
    """

    # extract measurement table rows for condition
    measurement_df = petab.get_rows_for_condition(
        measurement_df=petab_problem.measurement_df, condition=condition)

    if amici_model.nytrue != len(observable_ids):
        raise AssertionError("Number of AMICI model observables does not"
                             "match number of PEtab observables.")

    edata = amici.ExpData(amici_model)

    ##########################################################################
    # parameter mapping

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

    logger.debug(f"PEtab mapping: {parameter_mapping}")

    if len(condition_map_preeq) != len(condition_scale_map_preeq) \
            or len(condition_map_sim) != len(condition_scale_map_sim):
        raise AssertionError("Number of parameters and number of parameter "
                             "scales do not match.")
    if len(condition_map_preeq) \
            and len(condition_map_preeq) != len(condition_map_sim):
        logger.debug(f"Preequilibration parameter map: {condition_map_preeq}")
        logger.debug(f"Simulation parameter map: {condition_map_sim}")
        raise AssertionError("Number of parameters for preequilbration "
                             "and simulation do not match.")

    # PEtab parameter mapping may contain parameter_ids as values, these *must*
    # be replaced

    def _get_par(model_par, value):
        """Replace parameter IDs in mapping dicts by values from
        problem_parameters where necessary"""
        if isinstance(value, str):
            # estimated parameter
            # (condition table overrides have been handled by PEtab
            # parameter mapping)
            return problem_parameters[value]
        if model_par in problem_parameters:
            # user-provided
            return problem_parameters[model_par]
        # constant value
        return value

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

    condition_map_preeq_var, condition_map_preeq_fix = \
        subset_dict(condition_map_preeq, variable_par_ids, fixed_par_ids)

    condition_scale_map_preeq_var, condition_scale_map_preeq_fix = \
        subset_dict(condition_scale_map_preeq, variable_par_ids, fixed_par_ids)

    condition_map_sim_var, condition_map_sim_fix = \
        subset_dict(condition_map_sim, variable_par_ids, fixed_par_ids)

    condition_scale_map_sim_var, condition_scale_map_sim_fix = \
        subset_dict(condition_scale_map_sim, variable_par_ids, fixed_par_ids)

    logger.debug("Fixed parameters preequilibration: "
                 f"{condition_map_preeq_fix}")
    logger.debug("Fixed parameters simulation: "
                 f"{condition_map_sim_fix}")
    logger.debug("Variable parameters preequilibration: "
                 f"{condition_map_preeq_var}")
    logger.debug("Variable parameters simulation: "
                 f"{condition_map_sim_var}")

    petab.merge_preeq_and_sim_pars_condition(
        condition_map_preeq_var, condition_map_sim_var,
        condition_scale_map_preeq_var, condition_scale_map_sim_var,
        condition)

    logger.debug(f"Merged: {condition_map_sim_var}")

    # If necessary, bring parameters to linear scale
    if scaled_parameters:
        # For dynamic parameters we could also change ExpData.pscale, but since
        # we need to do it for fixed parameters anyways, we just do it for all
        # and set pscale to linear. we can skip preequilibration parameters,
        # because they are identical with simulation parameters, and only the
        # latter are used from here on
        unscale_parameters_dict(condition_map_preeq_fix,
                                condition_scale_map_preeq_fix)
        unscale_parameters_dict(condition_map_sim_fix,
                                condition_scale_map_sim_fix)
        unscale_parameters_dict(condition_map_sim_var,
                                condition_scale_map_sim_var)

    ##########################################################################
    # variable parameters and parameter scale

    # parameter list from mapping dict
    parameters = [condition_map_sim_var[par_id]
                  for par_id in amici_model.getParameterIds()]

    edata.parameters = parameters

    edata.pscale = amici.parameterScalingFromIntVector(
        [amici.ParameterScaling_none] * len(parameters))

    ##########################################################################
    # timepoints

    # find replicate numbers of time points
    timepoints_w_reps = _get_timepoints_with_replicates(
        df_for_condition=measurement_df)

    edata.setTimepoints(timepoints_w_reps)

    ##########################################################################
    # initial states
    # initial states have been set during model import. if they were
    # overwritten in the PEtab condition table, they would be handled as fixed
    # model parameters below

    ##########################################################################
    # fixed parameters preequilibration
    if condition_map_preeq:
        fixed_pars_preeq = [condition_map_preeq_fix[par_id]
                            for par_id in amici_model.getFixedParameterIds()]
        edata.fixedParametersPreequilibration = fixed_pars_preeq

    ##########################################################################
    # fixed parameters simulation
    fixed_pars_sim = [condition_map_sim_fix[par_id]
                      for par_id in amici_model.getFixedParameterIds()]
    edata.fixedParameters = fixed_pars_sim

    ##########################################################################
    # measurements and sigmas
    y, sigma_y = _get_measurements_and_sigmas(
        df_for_condition=measurement_df, timepoints_w_reps=timepoints_w_reps,
        observable_ids=observable_ids)
    edata.setObservedData(y.flatten())
    edata.setObservedDataStdDev(sigma_y.flatten())

    return edata


def unscale_parameter(value: numbers.Number,
                      petab_scale: str) -> numbers.Number:
    """Parameter to linear scale

    Arguments:
        value:
            Value to unscale
        petab_scale:
            Current scale of ``value``

    Returns:
        ``value`` on linear scale
    """
    if petab_scale == LIN:
        return value
    if petab_scale == LOG10:
        return np.power(10, value)
    if petab_scale == LOG:
        return np.exp(value)
    raise ValueError(f"Unknown parameter scale {petab_scale}. "
                     f"Must be from {(LIN, LOG, LOG10)}")


def unscale_parameters(values: Sequence[numbers.Number],
                       petab_scales: Sequence[str]) -> List[numbers.Number]:
    """Parameters to linear scale

    Arguments:
        values:
            Values to unscale
        petab_scales:
            Current scales of ``values``

    Returns:
        List of ``values`` on linear scale
    """
    return [unscale_parameter(value, scale)
            for value, scale in zip(values, petab_scales)]


def unscale_parameters_dict(
        value_dict: Dict[Any, numbers.Number],
        petab_scale_dict: Dict[Any, str]) -> None:
    """Parameters to linear scale

    Bring values in ``value_dict`` from current scale provided in
    ``petab_scale_dict`` to linear scale (in-place).
    Both arguments are expected to have the same length and matching keys.

    Arguments:
        value_dict:
            Values to unscale
        petab_scale_dict:
            Current scales of ``values``
    """
    if not value_dict.keys() == petab_scale_dict.keys():
        raise AssertionError("Keys don't match.")

    print(value_dict)
    for key, value in value_dict.items():
        value_dict[key] = unscale_parameter(value, petab_scale_dict[key])
    print(value_dict)


def _get_timepoints_with_replicates(
        df_for_condition: pd.DataFrame) -> List[numbers.Number]:
    """
    Get list of timepoints including replicate measurements

    Arguments:
        df_for_condition:
            PEtab measurement table subset for a single condition.

    Returns:
        Sorted list of timepoints, including multiple timepoints accounting
        for replicate measurements.
    """
    # create sorted list of all timepoints for which measurements exist
    timepoints = sorted(df_for_condition[TIME].unique().astype(float))

    # find replicate numbers of time points
    timepoints_w_reps = []
    for time in timepoints:
        # subselect for time
        df_for_time = df_for_condition[df_for_condition.time == time]
        # rep number is maximum over rep numbers for observables
        n_reps = max(df_for_time.groupby(
            [OBSERVABLE_ID, TIME]).size())
        # append time point n_rep times
        timepoints_w_reps.extend([time] * n_reps)

    return timepoints_w_reps


def _get_measurements_and_sigmas(
        df_for_condition: pd.DataFrame,
        timepoints_w_reps: Sequence[numbers.Number],
        observable_ids: Sequence[str]) -> Tuple[np.array, np.array]:
    """Get measurements and sigmas

    Generate arrays with measurements and sigmas in AMICI format from a
    PEtab measurement table subset for a single condition.

    Arguments:
        df_for_condition:
            Subset of PEtab measurement table for one condition
        timepoints_w_reps:
            Timepoints for which there exist measurements, including replicates
        observable_ids:
            List of observable IDs for mapping IDs to indices.
    """
    # prepare measurement matrix
    y = np.full(shape=(len(timepoints_w_reps), len(observable_ids)),
                fill_value=np.nan)
    # prepare sigma matrix
    sigma_y = y.copy()

    timepoints = sorted(df_for_condition[TIME].unique().astype(float))

    for time in timepoints:
        # subselect for time
        df_for_time = df_for_condition[df_for_condition[TIME] == time]
        time_ix_0 = timepoints_w_reps.index(time)

        # remember used time indices for each observable
        time_ix_for_obs_ix = {}

        # iterate over measurements
        for _, measurement in df_for_time.iterrows():
            # extract observable index
            observable_ix = observable_ids.index(measurement[OBSERVABLE_ID])

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


def rdatas_to_measurement_df(
        rdatas: Sequence[amici.ReturnData],
        model: amici.Model,
        measurement_df: pd.DataFrame) -> pd.DataFrame:
    """
    Create a measurement dataframe in the PEtab format from the passed
    `rdatas` and own information.

    Parameters:
        rdatas:
            A sequence of rdatas with the ordering of
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
    for (_, condition), rdata in zip(simulation_conditions.iterrows(), rdatas):
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
            observable_idx = observable_ids.index(row[OBSERVABLE_ID])
            measurement_sim = y[timepoint_idx, observable_idx]

            # change measurement entry
            row_sim.measurement = measurement_sim

            # append to dataframe
            df = df.append(row_sim, ignore_index=True)

    return df
