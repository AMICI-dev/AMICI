"""Importer-independent processing of observation models.

Translates :class:`amici.importers.utils.MeasurementChannel` definitions into
observables, sigmas, and log-likelihood contributions, stored in the
``SymbolId``-keyed symbol dictionaries used by the importers.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from typing import Any

import sympy as sp

from amici._symbolic.de_model import DEModel
from amici._symbolic.de_model_components import symbol_to_type
from amici.constants import SymbolId
from amici.importers.utils import (
    MeasurementChannel,
    noise_distribution_to_cost_function,
    symbol_with_assumptions,
)

__all__ = [
    "add_components",
    "add_log_likelihoods",
    "add_observables",
    "check_event_id",
    "get_placeholders",
]

#: ``SymbolId`` -> {symbol -> component attributes}
SymbolDict = dict[SymbolId, dict[sp.Symbol, dict[str, Any]]]


def add_components(de_model: DEModel, symbols: SymbolDict) -> None:
    """Add the entries of a symbol dictionary to a model as components.

    :param de_model:
        The model to add the components to.
    :param symbols:
        The symbol dictionary. Each entry is added as a component of the
        type corresponding to its ``SymbolId``.
    """
    for symbol_name in symbols:
        # transform dict of lists into a list of dicts
        args = ["name", "symbol"]

        if symbol_name == SymbolId.SPECIES:
            args += ["dt", "init"]
        elif symbol_name == SymbolId.ALGEBRAIC_STATE:
            args += ["init"]
        else:
            args += ["value"]

        if symbol_name == SymbolId.EVENT:
            args += [
                "assignments",
                "initial_value",
                "priority",
                "use_values_from_trigger_time",
            ]
        elif symbol_name == SymbolId.OBSERVABLE:
            args += ["noise_distribution", "measurement_symbol"]
        elif symbol_name == SymbolId.EVENT_OBSERVABLE:
            args += [
                "event",
                "measurement_symbol",
                "regularization_symbol",
                "noise_distribution",
            ]

        comp_kwargs = [
            {
                "symbol": var_id,
                **{k: v for k, v in var.items() if k in args},
            }
            for var_id, var in symbols[symbol_name].items()
        ]

        for comp_kwarg in comp_kwargs:
            de_model.add_component(symbol_to_type[symbol_name](**comp_kwarg))


def check_event_id(symbols: SymbolDict, channel: MeasurementChannel) -> None:
    """Check that the event referenced by a measurement channel exists.

    :param symbols:
        The importer's symbol dictionary, including all events.
    :param channel:
        The measurement channel to check.
    :raises ValueError: If the referenced event does not exist.
    """
    if (
        channel.event_id is not None
        and sp.Symbol(channel.event_id) not in symbols[SymbolId.EVENT]
    ):
        raise ValueError(
            "Could not find an event with the event identifier "
            f"{channel.event_id} for the event observable with "
            f"name {channel.name or channel.id}."
        )


def get_placeholders(
    observation_model: Iterable[MeasurementChannel],
) -> tuple[list[str], list[str]]:
    """Get the placeholder parameters declared by the measurement channels.

    :param observation_model:
        The measurement channels.
    :return:
        The identifiers of the noise and of the observable placeholders,
        each without duplicates, in order of first declaration.
    :raises NotImplementedError:
        If an event-resolved channel declares placeholders.
    """
    noise_placeholders = {}
    observable_placeholders = {}
    for channel in observation_model:
        if channel.is_event_resolved and (
            channel.noise_placeholders or channel.observable_placeholders
        ):
            raise NotImplementedError(
                "Placeholder parameters are only supported for "
                f"time-resolved observables, but {channel.id} is "
                "event-resolved."
            )
        noise_placeholders.update(dict.fromkeys(channel.noise_placeholders))
        observable_placeholders.update(
            dict.fromkeys(channel.observable_placeholders)
        )
    return list(noise_placeholders), list(observable_placeholders)


def add_observables(
    symbols: SymbolDict,
    observation_model: list[MeasurementChannel],
    sympify: Callable[[Any], sp.Expr],
) -> None:
    """Add observables, event observables, and placeholder parameters.

    :param symbols:
        The importer's symbol dictionary, including all events.
        Observable-related entries are written to it.
    :param observation_model:
        The measurement channels.
    :param sympify:
        Function to convert formulas to sympy expressions, resolving
        identifiers in the importer's namespace.
    """
    for channel in observation_model:
        check_event_id(symbols, channel)

    noise_parameters, observable_parameters = get_placeholders(
        observation_model
    )
    symbols[SymbolId.NOISE_PARAMETER] = {
        symbol_with_assumptions(noise_par): {"name": noise_par}
        for noise_par in noise_parameters
    }
    symbols[SymbolId.OBSERVABLE_PARAMETER] = {
        symbol_with_assumptions(obs_par): {"name": obs_par}
        for obs_par in observable_parameters
    }

    # Add time-resolved observables
    symbols[SymbolId.OBSERVABLE] = {
        symbol_with_assumptions(channel.id): {
            "name": channel.name or f"y{iobs}",
            "value": sympify(channel.formula),
            "noise_distribution": channel.noise_distribution,
        }
        for iobs, channel in enumerate(
            filter(lambda c: c.is_time_resolved, observation_model)
        )
    }
    # check for nesting of observables (unsupported)
    _check_symbol_nesting(symbols[SymbolId.OBSERVABLE], "observable")

    # Add event observables
    symbols[SymbolId.EVENT_OBSERVABLE] = {
        symbol_with_assumptions(channel.id): {
            "name": channel.name or f"z{iobs}",
            "value": sympify(channel.formula),
            "event": sp.Symbol(channel.event_id),
            "noise_distribution": channel.noise_distribution,
        }
        for iobs, channel in enumerate(
            filter(lambda c: c.is_event_resolved, observation_model)
        )
    }
    # check for nesting of observables (unsupported)
    _check_symbol_nesting(
        symbols[SymbolId.EVENT_OBSERVABLE], "eventObservable"
    )


def add_log_likelihoods(
    symbols: SymbolDict,
    observation_model: list[MeasurementChannel] | None,
    sympify: Callable[[Any], sp.Expr],
    make_derived_name_unique: Callable[[str], str],
    str_symbols: dict[str, str] | None = None,
) -> None:
    """Add sigmas and log-likelihood contributions for all observables.

    Observables and event observables must already be present in
    ``symbols``.

    :param symbols:
        The importer's symbol dictionary. Sigma and log-likelihood entries
        are written to it.
    :param observation_model:
        The measurement channels, or ``None`` if the observables were not
        generated from measurement channels.
    :param sympify:
        See :func:`add_observables`.
    :param make_derived_name_unique:
        Function to disambiguate the name of a derived symbol against all
        names already used in the model.
    :param str_symbols:
        The ``str_symbol`` to pass to custom noise distribution callables,
        by observable ID. Defaults to the observable ID.
    """
    channels = {channel.id: channel for channel in observation_model or []}
    str_symbols = str_symbols or {}
    _process_log_likelihood(
        symbols,
        channels,
        sympify,
        make_derived_name_unique,
        str_symbols,
        obs_type=SymbolId.OBSERVABLE,
        sigma_type=SymbolId.SIGMAY,
        llh_type=SymbolId.LLHY,
    )
    _process_log_likelihood(
        symbols,
        channels,
        sympify,
        make_derived_name_unique,
        str_symbols,
        obs_type=SymbolId.EVENT_OBSERVABLE,
        sigma_type=SymbolId.SIGMAZ,
        llh_type=SymbolId.LLHZ,
        llh_reg_type=SymbolId.LLHRZ,
    )


def _process_log_likelihood(
    symbols: SymbolDict,
    channels: dict[str, MeasurementChannel],
    sympify: Callable[[Any], sp.Expr],
    make_derived_name_unique: Callable[[str], str],
    str_symbols: dict[str, str],
    obs_type: SymbolId,
    sigma_type: SymbolId,
    llh_type: SymbolId,
    llh_reg_type: SymbolId | None = None,
):
    """Add sigmas and log-likelihood contributions for one observable type.

    :param symbols:
        See :func:`add_log_likelihoods`.
    :param channels:
        Measurement channels by observable ID. Observables without a channel
        get a normal noise distribution with unit sigma.
    :param sympify:
        See :func:`add_observables`.
    :param make_derived_name_unique:
        See :func:`add_log_likelihoods`.
    :param str_symbols:
        See :func:`add_log_likelihoods`.
    :param obs_type:
        The type of the observables to process.
    :param sigma_type:
        The type of the sigmas to add.
    :param llh_type:
        The type of the log-likelihood contributions to add.
    :param llh_reg_type:
        If not ``None``, additionally add log-likelihood contributions of
        this type for event observable regularization: the measurement is
        replaced by 0, and the observable by its regularization symbol.
    """
    observables = symbols[obs_type]
    for obs_id, obs in observables.items():
        obs["measurement_symbol"] = symbol_with_assumptions(
            make_derived_name_unique(f"m{obs_id}")
        )

    sigmas = {
        str(obs_id): sympify(channel.sigma)
        for obs_id in observables
        if (channel := channels.get(str(obs_id))) and channel.sigma is not None
    }
    # check that only the corresponding observable ID is used in the
    #  sigma formula, but not any other observable ID
    #  https://github.com/AMICI-dev/AMICI/issues/2561
    obs_syms = set(observables.keys())
    for obs_id in observables:
        if (sigma := sigmas.get(str(obs_id))) and sigma.has(
            *(obs_syms - {obs_id})
        ):
            raise ValueError(
                f"Sigma expression for {obs_id} ({sigma}) must not "
                f"contain any observable symbols other than {obs_id}."
            )
    symbols[sigma_type] = {
        symbol_with_assumptions(make_derived_name_unique(f"sigma_{obs_id}")): {
            "name": f"sigma_{obs['name']}",
            "value": sigmas.get(str(obs_id), sp.Float(1.0)),
        }
        for obs_id, obs in observables.items()
    }

    symbols[llh_type] = _log_likelihoods(
        symbols,
        channels,
        str_symbols,
        obs_type,
        sigma_type,
        make_derived_name_unique,
    )
    if llh_reg_type is None:
        return

    for obs_id, obs in observables.items():
        obs["regularization_symbol"] = symbol_with_assumptions(
            make_derived_name_unique(f"r{obs_id}")
        )
    symbols[llh_reg_type] = _log_likelihoods(
        symbols,
        channels,
        str_symbols,
        obs_type,
        sigma_type,
        make_derived_name_unique,
        regularization=True,
    )


def _log_likelihoods(
    symbols: SymbolDict,
    channels: dict[str, MeasurementChannel],
    str_symbols: dict[str, str],
    obs_type: SymbolId,
    sigma_type: SymbolId,
    make_derived_name_unique: Callable[[str], str],
    regularization: bool = False,
) -> dict[sp.Symbol, dict[str, Any]]:
    """Create the log-likelihood contributions for one observable type.

    See :func:`_process_log_likelihood`.
    """
    llhs = {}
    for (obs_id, obs), sigma_id in zip(
        symbols[obs_type].items(),
        symbols[sigma_type].keys(),
        strict=True,
    ):
        symbol = symbol_with_assumptions(
            make_derived_name_unique(f"J{obs_id}")
        )
        channel = channels.get(str(obs_id))
        dist = channel.noise_distribution if channel else "normal"
        value = noise_distribution_to_cost_function(dist)(
            obs_id,
            obs["measurement_symbol"],
            sigma_id,
            str_symbol=str_symbols.get(str(obs_id)),
        )
        if regularization:
            value = value.subs(obs["measurement_symbol"], 0.0)
            value = value.subs(obs_id, obs["regularization_symbol"])
        llhs[symbol] = {
            "name": f"J{obs['name']}",
            "value": value,
            "dist": dist,
        }
    return llhs


def _check_symbol_nesting(
    symbols: dict[sp.Symbol, dict[str, sp.Expr]], symbol_type: str
):
    observable_syms = set(symbols.keys())
    for obs in symbols.values():
        if any(sym in observable_syms for sym in obs["value"].free_symbols):
            raise ValueError(
                "Nested observables are not supported, "
                f"but {symbol_type} `{obs['name']} = {obs['value']}` "
                "references another observable."
            )
