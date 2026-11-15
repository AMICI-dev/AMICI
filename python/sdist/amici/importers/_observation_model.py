"""Importer-independent processing of observation models.

Translates :class:`amici.importers.utils.MeasurementChannel` definitions into
observables, sigmas, and log-likelihood contributions, stored in the
``SymbolId``-keyed symbol dictionaries used by the importers.
"""

from __future__ import annotations

import warnings
from collections.abc import Callable, Iterable
from typing import Any

import sympy as sp

from amici.constants import SymbolId
from amici.importers.utils import (
    MeasurementChannel,
    _get_str_symbol_identifiers,
    noise_distribution_to_cost_function,
    noise_distribution_to_observable_transformation,
    symbol_with_assumptions,
)

__all__ = ["add_log_likelihoods", "add_observables", "check_event_id"]

#: ``SymbolId`` -> {symbol -> component attributes}
SymbolDict = dict[SymbolId, dict[sp.Symbol, dict[str, Any]]]


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


def add_observables(
    symbols: SymbolDict,
    observation_model: list[MeasurementChannel],
    sympify: Callable[[Any], sp.Expr],
    noise_parameters: Iterable[str] = (),
    observable_parameters: Iterable[str] = (),
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
    :param noise_parameters:
        Names of the noise parameter placeholders.
    :param observable_parameters:
        Names of the observable parameter placeholders.
    """
    for channel in observation_model:
        check_event_id(symbols, channel)

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
            "transformation": noise_distribution_to_observable_transformation(
                channel.noise_distribution
            ),
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
            "transformation": noise_distribution_to_observable_transformation(
                channel.noise_distribution
            ),
        }
        for iobs, channel in enumerate(
            filter(lambda c: c.is_event_resolved, observation_model)
        )
    }
    wrong_t = sp.Symbol("t")
    for eo in symbols[SymbolId.EVENT_OBSERVABLE].values():
        if eo["value"].has(wrong_t):
            warnings.warn(
                f"Event observable {eo['name']} uses `t` in "
                "it's formula which is not the time variable. "
                "For the time variable, please use `time` "
                "instead!",
                stacklevel=1,
            )
    # check for nesting of observables (unsupported)
    _check_symbol_nesting(
        symbols[SymbolId.EVENT_OBSERVABLE], "eventObservable"
    )


def add_log_likelihoods(
    symbols: SymbolDict,
    observation_model: list[MeasurementChannel] | None,
    sympify: Callable[[Any], sp.Expr],
    make_derived_name_unique: Callable[[str], str],
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
    """
    sigmas = {
        channel.id: channel.sigma
        for channel in observation_model or []
        if channel.is_time_resolved and channel.sigma is not None
    }
    noise_distributions = {
        channel.id: channel.noise_distribution
        for channel in observation_model or []
        if channel.is_time_resolved
    }
    _process_log_likelihood(
        symbols, sigmas, noise_distributions, sympify, make_derived_name_unique
    )

    sigmas = {
        channel.id: channel.sigma
        for channel in observation_model or []
        if channel.is_event_resolved and channel.sigma is not None
    }
    noise_distributions = {
        channel.id: channel.noise_distribution
        for channel in observation_model or []
        if channel.is_event_resolved
    }
    _process_log_likelihood(
        symbols,
        sigmas,
        noise_distributions,
        sympify,
        make_derived_name_unique,
        events=True,
    )
    _process_log_likelihood(
        symbols,
        sigmas,
        noise_distributions,
        sympify,
        make_derived_name_unique,
        events=True,
        event_reg=True,
    )


def _process_log_likelihood(
    symbols: SymbolDict,
    sigmas: dict[str, str | float],
    noise_distributions: dict[str, str],
    sympify: Callable[[Any], sp.Expr],
    make_derived_name_unique: Callable[[str], str],
    events: bool = False,
    event_reg: bool = False,
):
    """
    Perform symbolic computations required for objective function
    evaluation.

    :param symbols:
        See :func:`add_log_likelihoods`.

    :param sigmas:
        Sigma expressions by observable ID.

    :param noise_distributions:
        Noise distributions by observable ID.

    :param sympify:
        See :func:`add_observables`.

    :param make_derived_name_unique:
        See :func:`add_log_likelihoods`.

    :param events:
        indicates whether the passed definitions are for observables
        (False) or for event observables (True).

    :param event_reg:
        indicates whether log-likelihood definitions should be processed
        for event observable regularization (Jrz). If this is activated,
        measurements are substituted by 0 and the observable by the
        respective regularization symbol.
    """

    if events:
        if event_reg:
            obs_symbol = SymbolId.EVENT_OBSERVABLE
            sigma_symbol = SymbolId.SIGMAZ
            llh_symbol = SymbolId.LLHRZ
        else:
            obs_symbol = SymbolId.EVENT_OBSERVABLE
            sigma_symbol = SymbolId.SIGMAZ
            llh_symbol = SymbolId.LLHZ
    else:
        assert not event_reg
        obs_symbol = SymbolId.OBSERVABLE
        sigma_symbol = SymbolId.SIGMAY
        llh_symbol = SymbolId.LLHY

    for obs_id, obs in symbols[obs_symbol].items():
        # event observables go through this method twice (once with
        # `event_reg=False`, once with `event_reg=True`, both against
        # the same `symbols[SymbolId.EVENT_OBSERVABLE]` dict) --
        # only mint the measurement symbol once, so a second pass
        # doesn't see its own first choice as an already-taken name and
        # needlessly disambiguate again
        if "measurement_symbol" not in obs:
            obs["measurement_symbol"] = symbol_with_assumptions(
                make_derived_name_unique(f"m{obs_id}")
            )
        if event_reg:
            obs["regularization_symbol"] = symbol_with_assumptions(
                make_derived_name_unique(f"r{obs_id}")
            )

    if not event_reg:
        sigmas = {obs_id: sympify(sigma) for obs_id, sigma in sigmas.items()}
        obs_syms = set(symbols[obs_symbol].keys())
        for obs_id in symbols[obs_symbol]:
            if (sigma := sigmas.get(str(obs_id))) and sigma.has(
                *(obs_syms - {obs_id})
            ):
                raise ValueError(
                    f"Sigma expression for {obs_id} ({sigma}) must not "
                    f"contain any observable symbols other than {obs_id}."
                )
        # check that only the corresponding observable ID is used in the
        #  sigma formula, but not any other observable ID
        #  https://github.com/AMICI-dev/AMICI/issues/2561
        symbols[sigma_symbol] = {
            symbol_with_assumptions(
                make_derived_name_unique(f"sigma_{obs_id}")
            ): {
                "name": f"sigma_{obs['name']}",
                "value": sigmas.get(str(obs_id), sp.Float(1.0)),
            }
            for obs_id, obs in symbols[obs_symbol].items()
        }

    symbols[llh_symbol] = {}
    for (obs_id, obs), (sigma_id, sigma) in zip(
        symbols[obs_symbol].items(),
        symbols[sigma_symbol].items(),
        strict=True,
    ):
        symbol = symbol_with_assumptions(
            make_derived_name_unique(f"J{obs_id}")
        )
        dist = noise_distributions.get(str(obs_id), "normal")
        cost_fun = noise_distribution_to_cost_function(dist)(obs_id)
        # TODO: clarify expected grammar for cost_fun
        #   python vs SBML L3 vs PEtab
        value = sp.sympify(
            cost_fun,
            locals=dict(
                zip(
                    _get_str_symbol_identifiers(obs_id),
                    (obs_id, obs["measurement_symbol"], sigma_id),
                    strict=True,
                )
            ),
        )
        if event_reg:
            value = value.subs(obs["measurement_symbol"], 0.0)
            value = value.subs(obs_id, obs["regularization_symbol"])
        symbols[llh_symbol][symbol] = {
            "name": f"J{obs['name']}",
            "value": value,
            "dist": dist,
        }


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
