"""
PEtab Objective
---------------
Functionality related to running simulations or evaluating the objective
function as defined by a PEtab problem
"""

import copy
import logging
import numbers
from typing import (List, Sequence, Optional, Dict, Tuple, Union, Any,
                    Collection, Iterator)

import amici
import libsbml
import numpy as np
import pandas as pd
import petab
from .logging import get_logger, log_execution_time
from petab.C import *
from .petab_import import PREEQ_INDICATOR_ID

LLH = 'llh'
SLLH = 'sllh'
FIM = 'fim'
S2LLH = 's2llh'
RES = 'res'
SRES = 'sres'
RDATAS = 'rdatas'

logger = get_logger(__name__)


@log_execution_time('Simulating PEtab model', logger)
def simulate_petab(
        petab_problem: petab.Problem,
        amici_model: amici.Model,
        solver: Optional[amici.Solver] = None,
        problem_parameters: Optional[Dict[str, float]] = None,
        simulation_conditions: Union[pd.DataFrame, Dict] = None,
        parameter_mapping: List[petab.ParMappingDictQuadruple] = None,
        scaled_parameters: Optional[bool] = False,
        log_level: int = logging.WARNING
) -> Dict[str, Any]:
    """
    Simulate PEtab model

    :param petab_problem:
        PEtab problem to work on.

    :param amici_model:
        AMICI Model assumed to be compatible with ``petab_problem``.

    :param solver:
        An AMICI solver. Will use default options if None.

    :param problem_parameters:
        Run simulation with these parameters. If None, PEtab `nominalValues`
        will be used). To be provided as dict,  mapping PEtab problem
        parameters to SBML IDs.

    :param simulation_conditions:
        Result of `petab.get_simulation_conditions`. Can be provided to save
        time if this has be obtained before.

    :param parameter_mapping:
        Optional precomputed PEtab parameter mapping for efficiency.

    :param scaled_parameters:
        If True, problem_parameters are assumed to be on the scale provided
        in the PEtab parameter table and will be unscaled. If False, they
        are assumed to be in linear scale.

    :param log_level:
        Log level, see :mod:`amici.logging` module.

    :return:
        Dictionary of

        * cost function value (LLH),
        * const function sensitivity w.r.t. parameters (SLLH),
          (**NOTE**: Sensitivities are computed for the non-scaled parameters)
        * list of `ReturnData` (RDATAS),

        corresponding to the different simulation conditions.
        For ordering of simulation conditions, see
        :meth:`petab.Problem.get_simulation_conditions_from_measurement_df`.
    """
    logger.setLevel(log_level)

    if solver is None:
        solver = amici_model.getSolver()

    # Get parameters
    if problem_parameters is None:
        # Use PEtab nominal values as default
        problem_parameters = {t.Index: getattr(t, NOMINAL_VALUE) for t in
                              petab_problem.parameter_df.itertuples()}

    # Because AMICI globalizes all local parameters during model import,
    # we need to do that here as well to prevent parameter mapping errors
    # (PEtab does currently not care about SBML LocalParameters)
    if petab_problem.sbml_document:
        converter_config = libsbml.SBMLLocalParameterConverter()\
            .getDefaultProperties()
        petab_problem.sbml_document.convert(converter_config)
    else:
        logger.debug("No petab_problem.sbml_document is set. Cannot convert "
                     "SBML LocalParameters. If the model contains "
                     "LocalParameters, parameter mapping will fail.")

    # Get parameter mapping
    if parameter_mapping is None:
        parameter_mapping = \
            petab_problem.get_optimization_to_simulation_parameter_mapping(
                warn_unmapped=False, scaled_parameters=scaled_parameters)

    # Generate ExpData with all condition-specific information
    edatas = edatas_from_petab(
        model=amici_model,
        petab_problem=petab_problem,
        problem_parameters=problem_parameters,
        simulation_conditions=simulation_conditions,
        parameter_mapping=parameter_mapping,
        scaled_parameters=scaled_parameters)

    # only required if we have preequilibration AND species or
    #  compartments in condition table, but shouldn't hurt to enable anyways.
    amici_model.setReinitializeFixedParameterInitialStates(True)

    # Simulate
    rdatas = amici.runAmiciSimulations(amici_model, solver, edata_list=edatas)

    # Compute total llh
    llh = sum(rdata['llh'] for rdata in rdatas)
    # Compute total sllh
    sllh = aggregate_sllh(amici_model=amici_model, rdatas=rdatas,
                          parameter_mapping=parameter_mapping)

    # TODO: implement me
    # # Compute total fim
    # fim = None
    # # Compute total s2llh
    # s2llh = None
    # # Compute total res
    # res = None
    # # Compute total sres
    # sres = None

    # log results
    sim_cond = petab_problem.get_simulation_conditions_from_measurement_df()

    for i, rdata in enumerate(rdatas):
        logger.debug(f"Condition: {sim_cond.iloc[i, :].values}, status: "
                     f"{rdata['status']}, llh: {rdata['llh']}")

    return {
        LLH: llh,
        SLLH: sllh,
        # FIM: fim,
        # S2LLH: s2llh,
        # RES: res,
        # SRES: sres,
        RDATAS: rdatas
    }


def edatas_from_petab(
        model: amici.Model,
        petab_problem: petab.Problem,
        problem_parameters: Dict[str, numbers.Number],
        simulation_conditions: Union[pd.DataFrame, Dict] = None,
        parameter_mapping: List[petab.ParMappingDictTuple] = None,
        scaled_parameters: Optional[bool] = False
) -> List[amici.ExpData]:
    """
    Create list of :class:`amici.amici.ExpData` objects for PEtab problem.

    Sets timepoints, fixed parameters (including preequilibration),
    non-fixed parameters, and observed data and sigmas.

    :param model:
        AMICI model.

    :param petab_problem:
        PEtab problem

    :param problem_parameters:
        Dictionary mapping parameter names of the PEtab problem to
        parameter values

    :param simulation_conditions:
        Result of petab.get_simulation_conditions. Can be provided to save
        time if this has be obtained before.

    :param parameter_mapping:
        Optional precomputed PEtab parameter mapping for efficiency.

    :param scaled_parameters:
        If True, problem_parameters are assumed to be on the scale provided
        in the PEtab parameter table and will be unscaled. If False,
        they are assumed to be in linear scale.

    :return:
        List with one :class:`amici.amici.ExpData` per simulation condition.
    """
    # number of amici simulations will be number of unique
    # (preequilibrationConditionId, simulationConditionId) pairs.
    # Can be optimized by checking for identical condition vectors.
    if simulation_conditions is None:
        simulation_conditions = \
            petab_problem.get_simulation_conditions_from_measurement_df()

    # Get parameter mapping
    if parameter_mapping is None:
        parameter_mapping = \
            petab_problem.get_optimization_to_simulation_parameter_mapping(
                warn_unmapped=False, scaled_parameters=scaled_parameters)

    observable_ids = model.getObservableIds()

    logger.debug(f"Problem parameters: {problem_parameters}")

    edatas = []
    for (_, condition), cur_parameter_mapping \
            in zip(simulation_conditions.iterrows(), parameter_mapping):
        # Create amici.ExpData for each simulation
        edata = get_edata_for_condition(
            condition=condition, amici_model=model, petab_problem=petab_problem,
            problem_parameters=problem_parameters,
            observable_ids=observable_ids,
            parameter_mapping=cur_parameter_mapping,
            scaled_parameters=scaled_parameters
        )
        edatas.append(edata)

    return edatas


def subset_dict(full: Dict[Any, Any],
                *args: Collection[Any]) -> Iterator[Dict[Any, Any]]:
    """
    Get subset of dictionary based on provided keys

    :param full:
        Dictionary to subset
    :param args:
        Collections of keys to be contained in the different subsets

    :return:
        subsetted dictionary
    """

    for keys in args:
        yield {key: val for (key, val) in full.items() if key in keys}


def get_edata_for_condition(
        condition: Union[Dict, pd.Series],
        problem_parameters: Dict[str, numbers.Number],
        amici_model: amici.Model,
        petab_problem: petab.Problem,
        observable_ids: List[str],
        parameter_mapping: Optional[petab.ParMappingDictTuple] = None,
        scaled_parameters: Optional[bool] = False
) -> amici.ExpData:
    """Get :class:`amici.amici.ExpData` for the given PEtab condition

    Sets timepoints, fixed parameters (including preequilibration),
    variable parameters, and observed data and sigmas.

    :param condition:
        pandas.DataFrame row with preequilibrationConditionId and
        simulationConditionId.

    :param problem_parameters:
        PEtab problem parameters as parameterId=>value dict. Only
        parameters included here will be set. Remaining parameters will
        be used currently set in `amici_model`.

    :param amici_model:
        AMICI model

    :param petab_problem:
        Underlying PEtab problem

    :param observable_ids:
        List of observable IDs

    :param parameter_mapping:
        PEtab parameter mapping for current condition

    :param scaled_parameters:
        If True, problem_parameters are assumed to be on the scale provided
        in the PEtab parameter table and will be unscaled. If False, they
        are assumed to be in linear scale.

    :return:
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

    (condition_map_preeq, condition_map_sim, condition_scale_map_preeq,
        condition_scale_map_sim) = parameter_mapping

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

    ##########################################################################
    # initial states
    # Initial states have been set during model import based on the SBML model.
    # If initial states were overwritten in the PEtab condition table, they are
    # applied here.
    # During model generation, parameter for initial concentrations and
    # respective initial assignments have been created for the
    # relevant species, here we add these parameters to the parameter mapping.
    # In absence of preequilibration this could also be handled via
    # ExpData.x0, but in the case of preequilibration this would not allow for
    # resetting initial states.

    species = [col for col in petab_problem.condition_df
               if petab_problem.sbml_model.getSpecies(col) is not None]
    if species:
        # set indicator fixed parameter for preeq
        # (we expect here, that this parameter was added during import and
        # that it was not added by the user with a different meaning...)
        if condition_map_preeq:
            condition_map_preeq[PREEQ_INDICATOR_ID] = 1.0
            condition_scale_map_preeq[PREEQ_INDICATOR_ID] = LIN

        condition_map_sim[PREEQ_INDICATOR_ID] = 0.0
        condition_scale_map_sim[PREEQ_INDICATOR_ID] = LIN

        def _set_initial_concentration(condition_id, species_id, init_par_id,
                                      par_map, scale_map):
            value = petab.to_float_if_float(
                petab_problem.condition_df.loc[condition_id, species_id])
            if isinstance(value, float):
                # numeric initial state
                par_map[init_par_id] = value
                scale_map[init_par_id] = petab.LIN
            else:
                # parametric initial state
                try:
                    # try find in mapping
                    par_map[init_par_id] = par_map[value]
                    scale_map[init_par_id] = scale_map[value]
                except KeyError:
                    # otherwise look up in parameter table
                    par_map[init_par_id] = problem_parameters[value]
                    if (scaled_parameters == False
                            or PARAMETER_SCALE
                            not in petab_problem.parameter_df
                            or not petab_problem.parameter_df.loc[
                                value, PARAMETER_SCALE]):
                        scale_map[init_par_id] = LIN
                    else:
                        scale_map[init_par_id] = \
                            petab_problem.parameter_df.loc[
                                value, PARAMETER_SCALE]

        for species_id in species:
            # for preequilibration
            init_par_id = f'initial_{species_id}_preeq'
            if PREEQUILIBRATION_CONDITION_ID in condition \
                    and condition[PREEQUILIBRATION_CONDITION_ID]:
                condition_id = condition[PREEQUILIBRATION_CONDITION_ID]
                _set_initial_concentration(
                    condition_id, species_id, init_par_id, condition_map_preeq,
                    condition_scale_map_preeq)
            else:
                # need to set dummy value for preeq parameter anyways, as it
                #  expected below (set to 0, not nan, because will be
                #  multiplied with indicator variable in initial assignment)
                condition_map_sim[init_par_id] = 0.0
                condition_scale_map_sim[init_par_id] = LIN

            # for simulation
            condition_id = condition[SIMULATION_CONDITION_ID]
            init_par_id = f'initial_{species_id}_sim'
            _set_initial_concentration(
                condition_id, species_id, init_par_id, condition_map_sim,
                condition_scale_map_sim)

    ##########################################################################

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

    # If necessary, scale parameters
    if not scaled_parameters:
        # We scale all parameters to the scale they are estimated on, and pass
        # that information to amici via edata.{parameters,pscale}.
        # The scaling is necessary to obtain correct derivatives.
        scale_parameters_dict(condition_map_preeq_fix,
                              condition_scale_map_preeq_fix)
        scale_parameters_dict(condition_map_sim_fix,
                              condition_scale_map_sim_fix)
        scale_parameters_dict(condition_map_sim_var,
                              condition_scale_map_sim_var)
        # We can skip preequilibration parameters, because they are identical
        # with simulation parameters, and only the latter are used from here
        # on.

    ##########################################################################
    # variable parameters and parameter scale

    # parameter list from mapping dict
    parameters = [condition_map_sim_var[par_id]
                  for par_id in amici_model.getParameterIds()]

    # scales list from mapping dict
    scales = [_to_amici_scale(condition_scale_map_sim_var[par_id])
              for par_id in amici_model.getParameterIds()]

    edata.parameters = parameters

    edata.pscale = amici.parameterScalingFromIntVector(scales)

    ##########################################################################
    # timepoints

    # find replicate numbers of time points
    timepoints_w_reps = _get_timepoints_with_replicates(
        df_for_condition=measurement_df)

    edata.setTimepoints(timepoints_w_reps)

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


def _to_amici_scale(petab_scale: str) -> int:
    """Convert petab scale id to amici scale id."""
    if petab_scale == LIN:
        return amici.ParameterScaling_none
    if petab_scale == LOG10:
        return amici.ParameterScaling_log10
    if petab_scale == LOG:
        return amici.ParameterScaling_ln


def scale_parameter(value: numbers.Number,
                    petab_scale: str) -> numbers.Number:
    """Bring parameter from linear scale to target scale.

    :param value:
        Value to scale
    :param petab_scale:
        Target scale of ``value``

    :return:
        ``value`` on linear scale
    """
    if petab_scale == LIN:
        return value
    if petab_scale == LOG10:
        return np.log10(value)
    if petab_scale == LOG:
        return np.log(value)
    raise ValueError(f"Unknown parameter scale {petab_scale}. "
                     f"Must be from {(LIN, LOG, LOG10)}")


def scale_parameters_dict(
        value_dict: Dict[Any, numbers.Number],
        petab_scale_dict: Dict[Any, str]) -> None:
    """
    Bring parameters from linear scale to target scale.

    Bring values in ``value_dict`` from linear scale to the scale
    provided in ``petab_scale_dict`` (in-place).
    Both arguments are expected to have the same length and matching keys.

    :param value_dict:
        Values to scale

    :param petab_scale_dict:
        Target scales of ``values``

    """
    if not value_dict.keys() == petab_scale_dict.keys():
        raise AssertionError("Keys don't match.")

    for key, value in value_dict.items():
        value_dict[key] = scale_parameter(value, petab_scale_dict[key])


def _get_timepoints_with_replicates(
        df_for_condition: pd.DataFrame) -> List[numbers.Number]:
    """
    Get list of timepoints including replicate measurements

    :param df_for_condition:
        PEtab measurement table subset for a single condition.

    :return:
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
    """
    Get measurements and sigmas

    Generate arrays with measurements and sigmas in AMICI format from a
    PEtab measurement table subset for a single condition.

    :param df_for_condition:
        Subset of PEtab measurement table for one condition

    :param timepoints_w_reps:
        Timepoints for which there exist measurements, including replicates

    :param observable_ids:
        List of observable IDs for mapping IDs to indices.

    :return:
        arrays for measurement and sigmas
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
            if isinstance(measurement.get(NOISE_PARAMETERS, None),
                    numbers.Number):
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

    :param rdatas:
        A sequence of rdatas with the ordering of
        `petab.get_simulation_conditions`.

    :param model:
        AMICI model used to generate `rdatas`.

    :param measurement_df:
        PEtab measurement table used to generate `rdatas`.

    :return:
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


def aggregate_sllh(
        amici_model: amici.Model,
        rdatas: Sequence[amici.ReturnDataView],
        parameter_mapping: Optional[List[petab.ParMappingDictTuple]],
) -> Union[None, Dict[str, float]]:
    """
    Aggregate likelihood gradient for all conditions, according to PEtab
    parameter mapping.

    :param amici_model:
        AMICI model from which ``rdatas`` were obtained.
    :param rdatas:
        Simulation results.
    :param parameter_mapping:
        PEtab parameter mapping to condition-specific
            simulation parameters.
    """
    sllh = {}
    model_par_ids = amici_model.getParameterIds()
    for (_, par_map_sim, _, _), rdata in zip(parameter_mapping, rdatas):
        if rdata['status'] != amici.AMICI_SUCCESS \
                or 'sllh' not in rdata\
                or rdata['sllh'] is None:
            return None

        for model_par_id, problem_par_id in par_map_sim.items():
            if isinstance(problem_par_id, str):
                model_par_idx  = model_par_ids.index(model_par_id)
                cur_par_sllh = rdata['sllh'][model_par_idx]
                try:
                    sllh[problem_par_id] += cur_par_sllh
                except KeyError:
                    sllh[problem_par_id] = cur_par_sllh
    return sllh
