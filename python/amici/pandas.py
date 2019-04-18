import pandas as pd
import numpy as np
import math
import copy

from .numpy import ExpDataView
import amici
from amici import ExpData


def getDataObservablesAsDataFrame(model, edata_list, by_id=False):
    """ Write Observables from experimental data as DataFrame.

    Arguments:
        model: Model instance.
        edata_list: list of ExpData instances with experimental data.
            May also be a single ExpData instance.
        by_id: bool (optional, default = False)
            If True, uses observable ids as identifiers in dataframe,
            otherwise the possibly more descriptive observable names
            are used.

    Returns:
        pandas DataFrame with conditions and observables.

    Raises:

    """
    if isinstance(edata_list, (amici.amici.ExpData, amici.amici.ExpDataPtr)):
        edata_list = [edata_list]

    # list of all column names using either ids or names
    cols = _get_extended_observable_cols(model, by_id=by_id)

    # initialize dataframe with columns
    df_edata = pd.DataFrame(columns=cols)

    # append all converted edatas
    for edata in edata_list:
        npdata = ExpDataView(edata)
        for i_time, timepoint in enumerate(edata.getTimepoints()):
            datadict = {
                'time': timepoint,
                'datatype': 'data'
            }
            # add observables and noises
            for i_obs, obs in enumerate(_get_names_or_ids(
                    model, 'Observable', by_id=by_id)):
                datadict[obs] = npdata['observedData'][i_time, i_obs]
                datadict[obs + '_std'] = \
                    npdata['observedDataStdDev'][i_time, i_obs]

            # add conditions
            _fill_conditions_dict(datadict, model, edata, by_id=by_id)

            df_edata.loc[len(df_edata)] = datadict

    return df_edata


def getSimulationObservablesAsDataFrame(
        model, edata_list, rdata_list, by_id=False):
    """ Write Observables from simulation results as DataFrame.

    Arguments:
        model: Model instance.
        edata_list: list of ExpData instances with experimental data.
            May also be a single ExpData instance.
        rdata_list: list of ReturnData instances corresponding to ExpData.
            May also be a single ReturnData instance.
        by_id: bool, optional (default = False)
            If True, ids are used as identifiers, otherwise the possibly more
            descriptive names.

    Returns:
        pandas DataFrame with conditions and observables.

    Raises:

    """
    if isinstance(edata_list, (amici.amici.ExpData, amici.amici.ExpDataPtr)):
        edata_list = [edata_list]
    if isinstance(rdata_list, (amici.amici.ReturnData, amici.amici.ReturnDataPtr)):
        rdata_list = [rdata_list]

    # list of all column names using either names or ids
    cols = _get_extended_observable_cols(model, by_id=by_id)

    # initialize dataframe with columns
    df_rdata = pd.DataFrame(columns=cols)

    # append all converted rdatas
    for edata, rdata in zip(edata_list, rdata_list):
        for i_time, timepoint in enumerate(rdata['t']):
            datadict = {
                'time': timepoint,
                'datatype': 'simulation',
            }
            # append simulations
            for i_obs, obs in enumerate(_get_names_or_ids(
                    model, 'Observable', by_id=by_id)):
                datadict[obs] = rdata['y'][i_time, i_obs]
                datadict[obs + '_std'] = rdata['sigmay'][i_time, i_obs]

            # use edata to fill conditions columns
            _fill_conditions_dict(datadict, model, edata, by_id=by_id)

            # append to dataframe
            df_rdata.loc[len(df_rdata)] = datadict

    return df_rdata


def getSimulationStatesAsDataFrame(
        model, edata_list, rdata_list, by_id=False):
    """ Compute model residuals according to lists of ReturnData and ExpData.

    Arguments:
        model: Model instance.
        edata_list: list of ExpData instances with experimental data.
            May also be a single ExpData instance.
        rdata_list: list of ReturnData instances corresponding to ExpData.
            May also be a single ReturnData instance.
        by_id: bool, optional (default = False)
            If True, ids are used as identifiers, otherwise the possibly more
            descriptive names.

    Returns:
        pandas DataFrame with conditions and observables.

    Raises:

    """
    if isinstance(edata_list, (amici.amici.ExpData, amici.amici.ExpDataPtr)):
        edata_list = [edata_list]
    if isinstance(rdata_list, (amici.amici.ReturnData, amici.amici.ReturnDataPtr)):
        rdata_list = [rdata_list]

    # get conditions and state column names by name or id
    cols = _get_state_cols(model, by_id=by_id)

    # initialize dataframe with columns
    df_rdata = pd.DataFrame(columns=cols)

    # append states
    for edata, rdata in zip(edata_list, rdata_list):
        for i_time, timepoint in enumerate(rdata['t']):
            datadict = {
                'time': timepoint,
            }

            # append states
            for i_state, state in enumerate(
                    _get_names_or_ids(model, 'State', by_id=by_id)):
                datadict[state] = rdata['x'][i_time, i_state]

            # use data to fill condition columns
            _fill_conditions_dict(datadict, model, edata, by_id=by_id)

            # append to dataframe
            df_rdata.loc[len(df_rdata)] = datadict

    return df_rdata


def getResidualsAsDataFrame(model, edata_list, rdata_list, by_id=False):
    """ Convert a list of ExpData to pandas DataFrame.

    Arguments:
        model: Model instance.
        edata_list: list of ExpData instances with experimental data.
            May also be a single ExpData instance.
        rdata_list: list of ReturnData instances corresponding to ExpData.
            May also be a single ReturnData instance.
        by_id: bool, optional (default = False)
            If True, ids are used as identifiers, otherwise the possibly more
            descriptive names.

    Returns:
        pandas DataFrame with conditions and observables.

    Raises:

    """
    if isinstance(edata_list, (amici.amici.ExpData, amici.amici.ExpDataPtr)):
        edata_list = [edata_list]
    if isinstance(rdata_list, (amici.amici.ReturnData, amici.amici.ReturnDataPtr)):
        rdata_list = [rdata_list]

    # create observable and simulation dataframes
    df_edata = getDataObservablesAsDataFrame(
        model, edata_list, by_id=by_id)
    df_rdata = getSimulationObservablesAsDataFrame(
        model, edata_list, rdata_list, by_id=by_id)

    # get all column names using names or ids
    cols = _get_observable_cols(model, by_id=by_id)

    #  initialize dataframe with columns
    df_res = pd.DataFrame(columns=cols)

    # iterate over rdata rows
    for row in df_rdata.index:
        datadict = {
            'time': df_rdata.loc[row]['time'],
            't_presim': df_rdata.loc[row]['t_presim']
        }

        # iterate over observables
        for obs in _get_names_or_ids(model, 'Observable', by_id=by_id):
            # compute residual and append to dict
            datadict[obs] = abs(
                (df_edata.loc[row][obs] - df_rdata.loc[row][obs]) /
                df_rdata.loc[row][obs + '_std'])

        # iterate over fixed parameters
        for par in _get_names_or_ids(model, 'FixedParameter', by_id=by_id):
            # fill in conditions
            datadict[par] = df_rdata.loc[row][par]
            datadict[par + '_preeq'] = df_rdata.loc[row][par + '_preeq']
            datadict[par + '_presim'] = df_rdata.loc[row][par + '_presim']

        # append to dataframe
        df_res.loc[len(df_res)] = datadict

    return df_res


def _fill_conditions_dict(datadict, model, edata, by_id) -> dict:
    """
    Helper function that fills in condition parameters from model and
    edata.

    Arguments:
        datadict: dictionary in which condition parameters will be inserted
            as key value pairs.
        model: Model instance.
        edata: ExpData instance.
        by_id: bool
            If True, ids are used as identifiers, otherwise the possibly more
            descriptive names.

    Returns:
        dictionary with filled condition parameters.

    Raises:

    """
    datadict['t_presim'] = edata.t_presim

    for i_par, par in enumerate(
            _get_names_or_ids(model, 'FixedParameter', by_id=by_id)):
        if len(edata.fixedParameters):
            datadict[par] = edata.fixedParameters[i_par]
        else:
            datadict[par] = model.getFixedParameters()[i_par]

        if len(edata.fixedParametersPreequilibration):
            datadict[par + '_preeq'] = \
                edata.fixedParametersPreequilibration[i_par]
        else:
            datadict[par + '_preeq'] = math.nan

        if len(edata.fixedParametersPresimulation):
            datadict[par + '_presim'] = \
                edata.fixedParametersPresimulation[i_par]
        else:
            datadict[par + '_presim'] = math.nan
    return datadict


def _get_extended_observable_cols(model, by_id) -> list:
    """ Construction helper for extended observable dataframe headers.

    Arguments:
        model: Model instance.
        by_id: bool
            If True, ids are used as identifiers, otherwise the possibly more
            descriptive names.

    Returns:
        column names as list.

    Raises:

    """
    return \
        ['time', 'datatype', 't_presim'] + \
        _get_names_or_ids(model, 'FixedParameter', by_id=by_id) + \
        [name + '_preeq' for name in
            _get_names_or_ids(model, 'FixedParameter', by_id=by_id)] + \
        [name + '_presim' for name in
            _get_names_or_ids(model, 'FixedParameter', by_id=by_id)] + \
        _get_names_or_ids(model, 'Observable', by_id=by_id) + \
        [name + '_std' for name in
            _get_names_or_ids(model, 'Observable', by_id=by_id)]


def _get_observable_cols(model, by_id):
    """ Construction helper for observable dataframe headers.

    Arguments:
        model: Model instance.
        by_id: bool
            If True, ids are used as identifiers, otherwise the possibly more
            descriptive names.

    Returns:
        column names as list.

    Raises:

    """
    return \
        ['time', 't_presim'] + \
        _get_names_or_ids(model, 'FixedParameter', by_id=by_id) + \
        [name + '_preeq' for name in
         _get_names_or_ids(model, 'FixedParameter', by_id=by_id)] + \
        [name + '_presim' for name in
         _get_names_or_ids(model, 'FixedParameter', by_id=by_id)] + \
        _get_names_or_ids(model, 'Observable', by_id=by_id)


def _get_state_cols(model, by_id):
    """ Construction helper for state dataframe headers.

    Arguments:
        model: Model instance.
        by_id: bool
            If True, ids are used as identifiers, otherwise the possibly more
            descriptive names.

    Returns:
        column names as list.

    Raises:

    """
    return \
        ['time', 't_presim'] + \
        _get_names_or_ids(model, 'FixedParameter', by_id=by_id) + \
        [name + '_preeq' for name in
            _get_names_or_ids(model, 'FixedParameter', by_id=by_id)] + \
        [name + '_presim' for name in
            _get_names_or_ids(model, 'FixedParameter', by_id=by_id)] + \
        _get_names_or_ids(model, 'State', by_id=by_id)


def _get_names_or_ids(model, variable, by_id):
    """
    Obtains a unique list of identifiers for the specified variable.
    First tries model.getVariableNames and then uses model.getVariableIds.

    Arguments:
        model: Model instance.
        variable: variable name.
        by_id: bool
            If True, ids are used as identifiers, otherwise first the possibly
            more descriptive names are used.

    Returns:
        column names as list.

    Raises:

    """
    # check whether variable type permitted
    variable_options = ['Parameter', 'FixedParameter', 'Observable', 'State']
    if variable not in variable_options:
        raise ValueError('Variable must be in ' + str(variable_options))

    # extract attributes
    names = list(getattr(model, f'get{variable}Names')())
    ids = list(getattr(model, f'get{variable}Ids')())

    # find out if model has names and ids
    has_names = getattr(model, f'has{variable}Names')()
    has_ids = getattr(model, f'has{variable}Ids')()

    # extract labels
    if not by_id and has_names and len(set(names)) == len(names):
        # use variable names
        return names
    elif has_ids:
        # use variable ids
        return ids
    else:
        # unable to create unique labels
        if by_id:
            msg = f"Model {variable} ids are not set."
        else:
            msg = f"Model {variable} names are not unique and " \
                  f"{variable} ids are not set."
        raise ValueError(msg)


def _get_specialized_fixed_parameters(
        model, condition, overwrite, by_id) -> list:
    """
    Copies values in condition and overwrites them according to key
    value pairs specified in overwrite.

    Arguments:
        model: Model instance.
        condition: dict/pd.Series containing FixedParameter values.
        overwrite: dict specifying which values in condition are to be replaced.
        by_id: bool
            If True, ids are used as identifiers, otherwise the possibly more
            descriptive names.

    Returns:
        overwritten FixedParameter as list.

    Raises:

    """
    cond = copy.deepcopy(condition)
    for field in overwrite:
        cond[field] = overwrite[field]
    return [float(cond[name]) for name in _get_names_or_ids(
        model, 'FixedParameter', by_id=by_id)]


def constructEdataFromDataFrame(df, model, condition, by_id=False):
    """ Constructs an ExpData instance according to the provided Model and DataFrame.

    Arguments:
        df: pd.DataFrame with Observable Names/Ids as columns.
            Standard deviations may be specified by appending '_std' as suffix.
        model: Model instance.
        condition: pd.Series with FixedParameter Names/Ids as columns.
            Preequilibration conditions may be specified by appending '_preeq' as suffix.
            Presimulation conditions may be specified by appending '_presim' as suffix.
        by_id: bool, optional (default = False)
            Indicate whether in the arguments, column headers are based on ids or names.
            This should correspond to the way `df` and `condition` was created in the
            first place.

    Returns:
        ExpData instance.

    Raises:

    """
    # initialize edata
    edata = ExpData(model.get())

    # timepoints
    df = df.sort_values(by='time', ascending=True)
    edata.setTimepoints(df['time'].values.astype(float))

    # get fixed parameters from condition
    overwrite_preeq = {}
    overwrite_presim = {}
    for par in list(_get_names_or_ids(model, 'FixedParameter', by_id=by_id)):
        if par + '_preeq' in condition.keys() \
                and not math.isnan(condition[par + '_preeq']):
            overwrite_preeq[par] = condition[par + '_preeq']
        if par + '_presim' in condition.keys() \
                and not math.isnan(condition[par + '_presim']):
            overwrite_presim[par] = condition[par + '_presim']

    # fill in fixed parameters
    edata.fixedParameters = \
        condition[_get_names_or_ids(model, 'FixedParameter', by_id=by_id)].values

    # fill in preequilibration parameters
    if any([overwrite_preeq[key] != condition[key] for key in
            overwrite_preeq.keys()]):
        edata.fixedParametersPreequilibration = \
            _get_specialized_fixed_parameters(
                model, condition, overwrite_preeq, by_id=by_id)
    elif len(overwrite_preeq.keys()):
        edata.fixedParametersPreequilibration = copy.deepcopy(
            edata.fixedParameters
        )

    # fill in presimulation parameters
    if any([overwrite_presim[key] != condition[key] for key in
            overwrite_presim.keys()]):
        edata.fixedParametersPresimulation = _get_specialized_fixed_parameters(
            model, condition, overwrite_presim, by_id=by_id
        )
    elif len(overwrite_presim.keys()):
        edata.fixedParametersPresimulation = copy.deepcopy(
            edata.fixedParameters
        )

    # fill in presimulation time
    if 't_presim' in condition.keys():
        edata.t_presim = float(condition['t_presim'])

    # fill in data and stds
    for obs_index, obs in enumerate(
            _get_names_or_ids(model, 'Observable', by_id=by_id)):
        if obs in df.keys():
            edata.setObservedData(df[obs].values.astype(float), obs_index)
        if obs + '_std' in df.keys():
            edata.setObservedDataStdDev(
                df[obs + '_std'].values.astype(float), obs_index
            )

    return edata


def getEdataFromDataFrame(model, df, by_id=False):
    """ Constructs a ExpData instance according to the provided Model and DataFrame.

    Arguments:
        df: pd.DataFrame with Observable Names/Ids, FixedParameter Names/Ids and time as columns.
            Standard deviations may be specified by appending '_std' as suffix.
            Preequilibration fixedParameters may be specified by appending '_preeq' as suffix.
            Presimulation fixedParameters may be specified by appending '_presim' as suffix.
            Presimulation time may be specified as 't_presim' column.
        model: Model instance.
        by_id: bool, optional (default = False)
            Whether the column names in `df` are based on ids or names,
            corresponding to how the dataframe was created in the first place.

    Returns:
        ExpData instance.

    Raises:

    """
    edata_list = []

    # aggregate features that define a condition

    # fixed parameters
    condition_parameters = _get_names_or_ids(model, 'FixedParameter', by_id=by_id)
    # preeq and presim parameters
    for par in _get_names_or_ids(model, 'FixedParameter', by_id=by_id):
        if par + '_preeq' in df.columns:
            condition_parameters.append(par + '_preeq')
        if par + '_presim' in df.columns:
            condition_parameters.append(par + '_presim')
    # presimulation time
    if 't_presim' in df.columns:
        condition_parameters.append('t_presim')
    # drop duplicates to create final conditions
    conditions = df[condition_parameters].drop_duplicates()

    # iterate over conditions
    for row in conditions.iterrows():
        # subselect rows that match condition
        selected = np.ones((len(df),), dtype=bool)
        for par_label, par in row[1].iteritems():
            if math.isnan(par):
                selected = selected & np.isnan(df[par_label].values)
            else:
                selected = selected & (df[par_label] == par)
        edata_df = df[selected]

        edata_list.append(constructEdataFromDataFrame(edata_df, model, row[1], by_id=by_id))

    return edata_list
