"""
Pandas Wrappers
---------------
This module contains convenience wrappers that allow for easy interconversion
between C++ objects from :mod:`amici.amici` and pandas DataFrames
"""

from __future__ import annotations

import copy
import math
from typing import SupportsFloat

import numpy as np
import pandas as pd

from . import (
    ExpData,
    ExpDataPtr,
    ExpDataView,
    Model,
    ModelPtr,
    ParameterScaling,
    ReturnDataView,
    parameter_scaling_from_int_vector,
)

__all__ = [
    "get_expressions_as_dataframe",
    "get_edata_from_data_frame",
    "get_data_observables_as_data_frame",
    "get_simulation_observables_as_data_frame",
    "get_simulation_states_as_data_frame",
    "get_residuals_as_data_frame",
]

ExpDatas = list[ExpData] | list[ExpDataPtr] | ExpData | ExpDataPtr

ReturnDatas = list[ReturnDataView] | ReturnDataView

AmiciModel = ModelPtr | Model


def _process_edata_list(edata_list: ExpDatas) -> list[ExpData]:
    """
    Maps single instances of :class:`ExpData` to lists of
    :class:`ExpData`

    :param edata_list:
        list of instances or single instance

    :return:
        list of instance(s)
    """
    if isinstance(edata_list, (ExpData | ExpDataPtr)):
        return [edata_list]
    else:
        return edata_list


def _process_rdata_list(rdata_list: ReturnDatas) -> list[ReturnDataView]:
    """
    Maps single instances of :class:`ReturnData` to lists of
    :class:`ReturnData`

    :param rdata_list:
        list of instances or single instance

    :return:
        list of instance(s)
    """
    if isinstance(rdata_list, ReturnDataView):
        return [rdata_list]
    else:
        return rdata_list


def get_data_observables_as_data_frame(
    model: AmiciModel, edata_list: ExpDatas, by_id: bool | None = False
) -> pd.DataFrame:
    """
    Write Observables from experimental data as DataFrame.

    :param model:
        Model instance.

    :param edata_list:
        list of ExpData instances with experimental data.
        May also be a single ExpData instance.

    :param by_id:
        If True, uses observable ids as column names in the generated
        DataFrame, otherwise the possibly more descriptive observable names
        are used.

    :return:
        pandas DataFrame with conditions/timepoints as rows and observables as
        columns.
    """
    edata_list = _process_edata_list(edata_list)

    # list of all column names using either ids or names
    cols = _get_extended_observable_cols(model, by_id=by_id)

    # aggregate records
    dicts = []
    for edata in edata_list:
        npdata = ExpDataView(edata)
        for i_time, timepoint in enumerate(edata.get_timepoints()):
            datadict = {"time": timepoint, "datatype": "data"}
            # add observables and noises
            for i_obs, obs in enumerate(
                _get_names_or_ids(model, "observable", by_id=by_id)
            ):
                datadict[obs] = npdata["measurements"][i_time, i_obs]
                datadict[obs + "_std"] = npdata["noise_scales"][i_time, i_obs]

            # add conditions
            _fill_conditions_dict(datadict, model, edata, by_id=by_id)

            dicts.append(datadict)

    return pd.DataFrame.from_records(dicts, columns=cols)


def get_simulation_observables_as_data_frame(
    model: Model,
    edata_list: ExpDatas,
    rdata_list: ReturnDatas,
    by_id: bool | None = False,
) -> pd.DataFrame:
    """
    Write Observables from simulation results as DataFrame.

    :param model:
        Model instance.

    :param edata_list:
        list of ExpData instances with experimental data.
        May also be a single ExpData instance.

    :param rdata_list:
        list of ReturnData instances corresponding to ExpData.
        May also be a single ReturnData instance.

    :param by_id:
        If True, ids are used as identifiers, otherwise the possibly more
        descriptive names.

    :return:
        pandas DataFrame with conditions/timepoints as rows and observables as
        columns.
    """
    edata_list = _process_edata_list(edata_list)
    rdata_list = _process_rdata_list(rdata_list)

    # list of all column names using either names or ids
    cols = _get_extended_observable_cols(model, by_id=by_id)

    # aggregate records
    dicts = []
    for edata, rdata in zip(edata_list, rdata_list, strict=True):
        for i_time, timepoint in enumerate(rdata["t"]):
            datadict = {
                "time": timepoint,
                "datatype": "simulation",
            }
            # append simulations
            for i_obs, obs in enumerate(
                _get_names_or_ids(model, "observable", by_id=by_id)
            ):
                datadict[obs] = rdata["y"][i_time, i_obs]
                datadict[obs + "_std"] = rdata["sigmay"][i_time, i_obs]

            # use edata to fill conditions columns
            _fill_conditions_dict(datadict, model, edata, by_id=by_id)

            # append to dataframe
            dicts.append(datadict)

    return pd.DataFrame.from_records(dicts, columns=cols)


def get_simulation_states_as_data_frame(
    model: Model,
    edata_list: ExpDatas,
    rdata_list: ReturnDatas,
    by_id: bool | None = False,
) -> pd.DataFrame:
    """
    Get model state according to lists of ReturnData and ExpData.

    :param model:
        Model instance.

    :param edata_list:
        list of ExpData instances with experimental data.
        May also be a single ExpData instance.

    :param rdata_list:
        list of ReturnData instances corresponding to ExpData.
        May also be a single ReturnData instance.

    :param by_id:
        If True, ids are used as identifiers, otherwise the possibly more
        descriptive names.

    :return: pandas DataFrame with conditions/timepoints as rows and
        state variables as columns.
    """
    edata_list = _process_edata_list(edata_list)
    rdata_list = _process_rdata_list(rdata_list)

    # get conditions and state column names by name or id
    cols = _get_state_cols(model, by_id=by_id)

    # aggregate records
    dicts = []
    for edata, rdata in zip(edata_list, rdata_list, strict=True):
        for i_time, timepoint in enumerate(rdata["t"]):
            datadict = {
                "time": timepoint,
            }

            # append states
            for i_state, state in enumerate(
                _get_names_or_ids(model, "state", by_id=by_id)
            ):
                datadict[state] = rdata["x"][i_time, i_state]

            # use data to fill condition columns
            _fill_conditions_dict(datadict, model, edata, by_id=by_id)

            # append to dataframe
            dicts.append(datadict)

    return pd.DataFrame.from_records(dicts, columns=cols)


def get_expressions_as_dataframe(
    model: Model,
    edata_list: ExpDatas,
    rdata_list: ReturnDatas,
    by_id: bool | None = False,
) -> pd.DataFrame:
    """
    Get values of model expressions from lists of ReturnData as DataFrame.

    :param model:
        Model instance.

    :param edata_list:
        list of ExpData instances with experimental data.
        May also be a single ExpData instance.

    :param rdata_list:
        list of ReturnData instances corresponding to ExpData.
        May also be a single ReturnData instance.

    :param by_id:
        If True, ids are used as identifiers, otherwise the possibly more
        descriptive names.

    :return: pandas DataFrame with conditions/timepoints as rows and
        model expressions as columns.
    """
    edata_list = _process_edata_list(edata_list)
    rdata_list = _process_rdata_list(rdata_list)

    # get conditions and state column names by name or id
    cols = _get_expression_cols(model, by_id=by_id)

    # aggregate records
    dicts = []
    for edata, rdata in zip(edata_list, rdata_list, strict=True):
        for i_time, timepoint in enumerate(rdata["t"]):
            datadict = {
                "time": timepoint,
            }

            # append expressions
            for i_expr, expr in enumerate(
                _get_names_or_ids(model, "expression", by_id=by_id)
            ):
                datadict[expr] = rdata["w"][i_time, i_expr]

            # use data to fill condition columns
            _fill_conditions_dict(datadict, model, edata, by_id=by_id)

            # append to dataframe
            dicts.append(datadict)

    return pd.DataFrame.from_records(dicts, columns=cols)


def get_residuals_as_data_frame(
    model: Model,
    edata_list: ExpDatas,
    rdata_list: ReturnDatas,
    by_id: bool | None = False,
) -> pd.DataFrame:
    """
    Convert a list of ReturnData and ExpData to pandas DataFrame with
    residuals.

    :param model:
        Model instance.

    :param edata_list:
        list of ExpData instances with experimental data. May also be a
        single ExpData instance.

    :param rdata_list:
        list of ReturnData instances corresponding to ExpData. May also be a
        single ReturnData instance.

    :param by_id: bool, optional (default = False)
        If True, ids are used as identifiers, otherwise the possibly more
        descriptive names.

    :return:
        pandas DataFrame with conditions and residuals.
    """
    edata_list = _process_edata_list(edata_list)
    rdata_list = _process_rdata_list(rdata_list)

    # create observable and simulation dataframes
    df_edata = get_data_observables_as_data_frame(
        model, edata_list, by_id=by_id
    )
    df_rdata = get_simulation_observables_as_data_frame(
        model, edata_list, rdata_list, by_id=by_id
    )

    # get all column names using names or ids
    cols = _get_observable_cols(model, by_id=by_id)

    # aggregate records
    dicts = []
    for row in df_rdata.index:
        datadict = {
            "time": df_rdata.loc[row]["time"],
            "t_presim": df_rdata.loc[row]["t_presim"],
        }

        # iterate over observables
        for obs in _get_names_or_ids(model, "observable", by_id=by_id):
            # compute residual and append to dict
            datadict[obs] = abs(
                (df_edata.loc[row][obs] - df_rdata.loc[row][obs])
                / df_rdata.loc[row][obs + "_std"]
            )

        # iterate over fixed parameters
        for par in _get_names_or_ids(model, "fixed_parameter", by_id=by_id):
            # fill in conditions
            datadict[par] = df_rdata.loc[row][par]
            datadict[par + "_preeq"] = df_rdata.loc[row][par + "_preeq"]
            datadict[par + "_presim"] = df_rdata.loc[row][par + "_presim"]

        # append to dataframe
        dicts.append(datadict)

    return pd.DataFrame.from_records(dicts, columns=cols)


def _fill_conditions_dict(
    datadict: dict[str, float],
    model: AmiciModel,
    edata: ExpData,
    by_id: bool,
) -> dict[str, float]:
    """
    Helper function that fills in condition parameters from model and
    edata.

    :param datadict:
        dictionary in which condition parameters will be inserted
        as key value pairs.

    :param model:
        Model instance.

    :param edata:
        ExpData instance.

    :param by_id:
        If True, ids are used as identifiers, otherwise the possibly more
        descriptive names.

    :return:
        dictionary with filled condition parameters.

    """
    datadict["condition_id"] = edata.id
    datadict["t_presim"] = edata.t_presim

    for i_par, par in enumerate(
        _get_names_or_ids(model, "fixed_parameter", by_id=by_id)
    ):
        if len(edata.fixed_parameters):
            datadict[par] = edata.fixed_parameters[i_par]
        else:
            datadict[par] = model.get_fixed_parameters()[i_par]

        if len(edata.fixed_parameters_pre_equilibration):
            datadict[par + "_preeq"] = (
                edata.fixed_parameters_pre_equilibration[i_par]
            )
        else:
            datadict[par + "_preeq"] = np.nan

        if len(edata.fixed_parameters_presimulation):
            datadict[par + "_presim"] = edata.fixed_parameters_presimulation[
                i_par
            ]
        else:
            datadict[par + "_presim"] = np.nan

    for i_par, par in enumerate(
        _get_names_or_ids(model, "free_parameter", by_id=by_id)
    ):
        if len(edata.free_parameters):
            datadict[par] = edata.free_parameters[i_par]
        else:
            datadict[par] = model.get_free_parameters()[i_par]

        if len(edata.pscale):
            datadict[par + "_scale"] = edata.pscale[i_par]
        else:
            datadict[par + "_scale"] = model.get_parameter_scale()[i_par]

    return datadict


def _get_extended_observable_cols(model: AmiciModel, by_id: bool) -> list[str]:
    """
    Construction helper for extended observable dataframe headers.

    :param model:
        Model instance.

    :param by_id:
        If True, ids are used as identifiers, otherwise the possibly more
        descriptive names.

    :return:
        column names as list.
    """
    return (
        ["condition_id", "time", "datatype", "t_presim"]
        + _get_names_or_ids(model, "fixed_parameter", by_id=by_id)
        + [
            name + "_preeq"
            for name in _get_names_or_ids(
                model, "fixed_parameter", by_id=by_id
            )
        ]
        + [
            name + "_presim"
            for name in _get_names_or_ids(
                model, "fixed_parameter", by_id=by_id
            )
        ]
        + _get_names_or_ids(model, "free_parameter", by_id=by_id)
        + [
            name + "_scale"
            for name in _get_names_or_ids(model, "free_parameter", by_id=by_id)
        ]
        + _get_names_or_ids(model, "observable", by_id=by_id)
        + [
            name + "_std"
            for name in _get_names_or_ids(model, "observable", by_id=by_id)
        ]
    )


def _get_observable_cols(model: AmiciModel, by_id: bool) -> list[str]:
    """
    Construction helper for observable dataframe headers.

    :param model:
        Model instance.

    :param by_id:
        If True, ids are used as identifiers, otherwise the possibly more
        descriptive names.

    :return:
        column names as list.
    """
    return (
        ["condition_id", "time", "t_presim"]
        + _get_names_or_ids(model, "fixed_parameter", by_id=by_id)
        + [
            name + "_preeq"
            for name in _get_names_or_ids(
                model, "fixed_parameter", by_id=by_id
            )
        ]
        + [
            name + "_presim"
            for name in _get_names_or_ids(
                model, "fixed_parameter", by_id=by_id
            )
        ]
        + _get_names_or_ids(model, "free_parameter", by_id=by_id)
        + [
            name + "_scale"
            for name in _get_names_or_ids(model, "free_parameter", by_id=by_id)
        ]
        + _get_names_or_ids(model, "observable", by_id=by_id)
    )


def _get_state_cols(model: AmiciModel, by_id: bool) -> list[str]:
    """
    Construction helper for state dataframe headers.

    :param model:
        Model instance.

    :param by_id:
        If True, ids are used as identifiers, otherwise the possibly more
        descriptive names.

    :return:
        column names as list.
    """
    return (
        ["condition_id", "time", "t_presim"]
        + _get_names_or_ids(model, "fixed_parameter", by_id=by_id)
        + [
            name + "_preeq"
            for name in _get_names_or_ids(
                model, "fixed_parameter", by_id=by_id
            )
        ]
        + [
            name + "_presim"
            for name in _get_names_or_ids(
                model, "fixed_parameter", by_id=by_id
            )
        ]
        + _get_names_or_ids(model, "free_parameter", by_id=by_id)
        + [
            name + "_scale"
            for name in _get_names_or_ids(model, "free_parameter", by_id=by_id)
        ]
        + _get_names_or_ids(model, "state", by_id=by_id)
    )


def _get_expression_cols(model: AmiciModel, by_id: bool) -> list[str]:
    """Construction helper for expression dataframe headers.

    :param model:
        Model instance.

    :param by_id:
        If True, ids are used as identifiers, otherwise the possibly more
        descriptive names.

    :return:
        column names as list.
    """
    return (
        ["condition_id", "time", "t_presim"]
        + _get_names_or_ids(model, "fixed_parameter", by_id=by_id)
        + [
            name + "_preeq"
            for name in _get_names_or_ids(
                model, "fixed_parameter", by_id=by_id
            )
        ]
        + [
            name + "_presim"
            for name in _get_names_or_ids(
                model, "fixed_parameter", by_id=by_id
            )
        ]
        + _get_names_or_ids(model, "free_parameter", by_id=by_id)
        + [
            name + "_scale"
            for name in _get_names_or_ids(model, "free_parameter", by_id=by_id)
        ]
        + _get_names_or_ids(model, "expression", by_id=by_id)
    )


def _get_names_or_ids(
    model: AmiciModel, variable: str, by_id: bool
) -> list[str]:
    """
    Obtains a unique list of identifiers for the specified variable.
    First tries model.getVariableNames and then uses model.getVariableIds.

    :param model:
        Model instance.

    :param variable:
        variable name.

    :param by_id:
        If True, ids are used as identifiers, otherwise first the possibly
        more descriptive names are used.

    :return:
        column names as list.
    """
    # check whether variable type permitted
    variable_options = [
        "free_parameter",
        "fixed_parameter",
        "observable",
        "state",
        "expression",
    ]
    if variable not in variable_options:
        raise ValueError("Variable must be in " + str(variable_options))

    # extract attributes
    names = list(getattr(model, f"get_{variable}_names")())
    ids = list(getattr(model, f"get_{variable}_ids")())

    # find out if model has names and ids
    has_names = getattr(model, f"has_{variable}_names")()
    has_ids = getattr(model, f"has_{variable}_ids")()

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
            msg = (
                f"Model {variable} names are not unique and "
                f"{variable} ids are not set."
            )
        raise ValueError(msg)


def _get_specialized_fixed_parameters(
    model: AmiciModel,
    condition: dict[str, SupportsFloat] | pd.Series,
    overwrite: dict[str, SupportsFloat] | pd.Series,
    by_id: bool,
) -> list[float]:
    """
    Copies values in condition and overwrites them according to key
    value pairs specified in overwrite.

    :param model:
        Model instance.
    :param condition:
        fixedParameter values.
    :param overwrite:
        dict specifying which values in condition are to be replaced.
    :param by_id:
        If True, ids are used as identifiers, otherwise the possibly more
        descriptive names.

    :return:
        overwritten FixedParameter as list.
    """
    cond = copy.deepcopy(condition)
    for field in overwrite:
        cond[field] = overwrite[field]
    return [
        float(cond[name])
        for name in _get_names_or_ids(model, "fixed_parameter", by_id=by_id)
    ]


def construct_edata_from_data_frame(
    df: pd.DataFrame,
    model: AmiciModel,
    condition: pd.Series,
    by_id: bool | None = False,
) -> ExpData:
    """
    Constructs an ExpData instance according to the provided Model
    and DataFrame.

    :param df:
        pd.DataFrame with Observable Names/Ids as columns.
        Standard deviations may be specified by appending '_std' as suffix.

    :param model:
        Model instance.

    :param condition:
        pd.Series with (Fixed)Parameter Names/Ids as columns.
        Preequilibration conditions may be specified by appending
        '_preeq' as suffix. Presimulation conditions may be specified by
        appending '_presim' as suffix. Parameter scales may be specified by
        appending '_scale' as suffix.

    :param by_id:
        Indicate whether in the arguments, column headers are based on ids or
        names. This should correspond to the way `df` and `condition` was
        created in the first place.

    :return:
        ExpData instance.
    """
    # initialize edata
    edata = ExpData(model.get())

    # timepoints
    df = df.sort_values(by="time", ascending=True)
    edata.set_timepoints(df["time"].values.astype(float))

    # get fixed parameters from condition
    overwrite_preeq = {}
    overwrite_presim = {}
    for par in list(_get_names_or_ids(model, "fixed_parameter", by_id=by_id)):
        if par + "_preeq" in condition.keys() and not math.isnan(
            condition[par + "_preeq"].astype(float)
        ):
            overwrite_preeq[par] = condition[par + "_preeq"].astype(float)
        if par + "_presim" in condition.keys() and not math.isnan(
            condition[par + "_presim"].astype(float)
        ):
            overwrite_presim[par] = condition[par + "_presim"].astype(float)

    # fill in fixed parameters
    edata.fixed_parameters = (
        condition[_get_names_or_ids(model, "fixed_parameter", by_id=by_id)]
        .astype(float)
        .values
    )

    # fill in parameters
    edata.free_parameters = (
        condition[_get_names_or_ids(model, "free_parameter", by_id=by_id)]
        .astype(float)
        .values
    )

    edata.pscale = parameter_scaling_from_int_vector(
        [
            ParameterScaling(condition[par + "_scale"].astype(int))
            for par in list(
                _get_names_or_ids(model, "free_parameter", by_id=by_id)
            )
        ]
    )

    # fill in preequilibration parameters
    if any(
        [overwrite_preeq[key] != condition[key] for key in overwrite_preeq]
    ):
        edata.fixed_parameters_pre_equilibration = (
            _get_specialized_fixed_parameters(
                model, condition, overwrite_preeq, by_id=by_id
            )
        )
    elif len(overwrite_preeq):
        edata.fixed_parameters_pre_equilibration = copy.deepcopy(
            edata.fixed_parameters
        )

    # fill in presimulation parameters
    if any(
        [
            overwrite_presim[key] != condition[key]
            for key in overwrite_presim.keys()
        ]
    ):
        edata.fixed_parameters_presimulation = (
            _get_specialized_fixed_parameters(
                model, condition, overwrite_presim, by_id=by_id
            )
        )
    elif len(overwrite_presim.keys()):
        edata.fixed_parameters_presimulation = copy.deepcopy(
            edata.fixed_parameters
        )

    # fill in presimulation time
    if "t_presim" in condition.keys():
        edata.t_presim = float(condition["t_presim"])

    # fill in data and stds
    for obs_index, obs in enumerate(
        _get_names_or_ids(model, "observable", by_id=by_id)
    ):
        if obs in df.keys():
            edata.set_measurements(df[obs].values.astype(float), obs_index)
        if obs + "_std" in df.keys():
            edata.set_noise_scales(
                df[obs + "_std"].values.astype(float), obs_index
            )

    return edata


def get_edata_from_data_frame(
    model: AmiciModel, df: pd.DataFrame, by_id: bool | None = False
) -> list[ExpData]:
    """
    Constructs a ExpData instances according to the provided Model and
    DataFrame.

    :param df:
        dataframe with Observable Names/Ids, FixedParameter Names/Ids
        and time as columns. Standard deviations may be specified by
        appending '_std' as suffix. Preequilibration fixedParameters may be
        specified by appending '_preeq' as suffix. Presimulation
        fixedParameters may be specified by appending '_presim' as suffix.
        Presimulation time may be specified as 't_presim' column.

    :param model:
        Model instance.

    :param by_id:
        Whether the column names in `df` are based on ids or names,
        corresponding to how the dataframe was created in the first place.

    :return:
        list of ExpData instances.
    """
    edata_list = []

    # aggregate features that define a condition

    # fixed parameters
    condition_parameters = _get_names_or_ids(
        model, "fixed_parameter", by_id=by_id
    )
    # preeq and presim parameters
    for par in _get_names_or_ids(model, "fixed_parameter", by_id=by_id):
        if par + "_preeq" in df.columns:
            condition_parameters.append(par + "_preeq")
        if par + "_presim" in df.columns:
            condition_parameters.append(par + "_presim")
    # parameters & scales
    for par in _get_names_or_ids(model, "free_parameter", by_id=by_id):
        condition_parameters.append(par)
        condition_parameters.append(par + "_scale")
    # presimulation time
    if "t_presim" in df.columns:
        condition_parameters.append("t_presim")
    # drop duplicates to create final conditions
    conditions = df[condition_parameters].drop_duplicates()

    # iterate over conditions
    for ir, row in conditions.iterrows():
        # subselect rows that match condition
        selected = np.ones((len(df),), dtype=bool)
        for par_label, par in row.items():
            if math.isnan(par):
                selected = selected & np.isnan(
                    df[par_label].astype(float).values
                )
            else:
                selected = selected & (df[par_label] == par)
        edata_df = df[selected]

        edata_list.append(
            construct_edata_from_data_frame(edata_df, model, row, by_id=by_id)
        )

    return edata_list
