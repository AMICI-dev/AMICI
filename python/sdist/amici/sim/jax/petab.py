"""PEtab wrappers for JAX models.""" ""

import logging
import os
import re
import shutil
from collections.abc import Callable, Iterable, Sized
from numbers import Number
from pathlib import Path
from typing import NamedTuple

import diffrax
import equinox as eqx
import jax.lax
import jax.numpy as jnp
import jaxtyping as jt
import numpy as np
import optimistix
import pandas as pd
import petab.v1 as petabv1
import petab.v2 as petabv2
from optimistix import AbstractRootFinder

from amici import _module_from_path
from amici.logging import get_logger
from amici.sim.jax.model import JAXModel, ReturnValue

DEFAULT_CONTROLLER_SETTINGS = {
    "atol": 1e-8,
    "rtol": 1e-8,
    "pcoeff": 0.4,
    "icoeff": 0.3,
    "dcoeff": 0.0,
}

DEFAULT_ROOT_FINDER_SETTINGS = {
    "atol": 1e-12,
    "rtol": 1e-12,
}

SCALE_TO_INT = {
    petabv2.C.LIN: 0,
    petabv2.C.LOG: 1,
    petabv2.C.LOG10: 2,
}

logger = get_logger(__name__, logging.WARNING)


class SteadyStateEvent(eqx.Module):
    """Steady-state termination event with value-based equality.

    :func:`diffrax.steady_state_event` returns a fresh closure on every
    call, which Python compares by identity. Since ``steady_state_event`` is
    passed as a static argument into :func:`equinox.filter_jit`-compiled
    functions (:meth:`JAXModel.simulate_condition`,
    :meth:`JAXModel.preequilibrate_condition`), constructing a new closure
    with identical settings on every call (e.g. once per iteration of an
    optimization loop) silently defeats the JIT cache and forces a full
    recompilation, even though only numeric inputs (parameters) changed.
    Wrapping the settings in an :class:`equinox.Module` gives value-based
    ``__eq__``/``__hash__``, so equivalent instances share the compiled
    executable regardless of how many times they are (re-)constructed.
    """

    rtol: float | None = None
    atol: float | None = None

    def __call__(self, *args, **kwargs):
        return diffrax.steady_state_event(rtol=self.rtol, atol=self.atol)(
            *args, **kwargs
        )


def jax_unscale(
    parameter: jnp.float_,
    scale_str: str,
) -> jnp.float_:
    """Unscale parameter according to ``scale_str``.

    Arguments:
        parameter:
            Parameter to be unscaled.
        scale_str:
            One of ``petabv2.C.LIN``, ``petabv2.C.LOG``, ``petabv2.C.LOG10``.

    Returns:
        The unscaled parameter.
    """
    if scale_str == petabv2.C.LIN or not scale_str:
        return parameter
    if scale_str == petabv2.C.LOG:
        return jnp.exp(parameter)
    if scale_str == petabv2.C.LOG10:
        return jnp.power(10, parameter)
    raise ValueError(f"Invalid parameter scaling: {scale_str}")


class OverrideColumn(NamedTuple):
    """Numeric values, free-parameter mask, and free-parameter indices for
    one observable/noise parameter override column (each of shape
    ``(n_rows, n_pars)``). Recombine via
    ``jnp.where(mask, p[index], numeric)``.
    """

    numeric: np.ndarray
    mask: np.ndarray
    index: np.ndarray

    @classmethod
    def placeholder(cls, n_rows: int, n_pars: int) -> "OverrideColumn":
        """All-numeric-one column, used for an absent/empty override
        column, or for a synthetic/padding measurement row."""
        numeric = np.ones((n_rows, n_pars))
        return cls(
            numeric,
            np.zeros_like(numeric, dtype=bool),
            np.zeros_like(numeric, dtype=int),
        )

    @classmethod
    def concatenate(cls, *columns: "OverrideColumn") -> "OverrideColumn":
        return cls(
            np.concatenate([c.numeric for c in columns]),
            np.concatenate([c.mask for c in columns]),
            np.concatenate([c.index for c in columns]),
        )


def _get_fixed_parameter_values(
    petab_problem: petabv2.Problem,
) -> dict[str, float]:
    """Nominal (linear) values of fixed (non-estimated) parameters, used to
    resolve observable/noise parameter overrides that reference them.

    Array-valued (PEtab-SciML) fixed parameters are excluded: their nominal
    value is the literal string ``"array"``, not a substitutable number.

    Built once (from :attr:`petabv2.Problem.parameter_tables` directly)
    rather than via the :attr:`petabv2.Problem.parameter_df` property on
    every override column, since building the latter can be expensive and,
    for some problems (e.g. array-valued PEtab-SciML parameters), emits
    pydantic serialization warnings.
    """
    return {
        p.id: p.nominal_value
        for pt in petab_problem.parameter_tables
        for p in pt.elements
        if not p.estimate and p.nominal_value != "array"
    }


def _resolve_override_symbol(value, fixed_parameter_values: dict[str, float]):
    """Replace a reference to a non-estimated PEtab parameter by its
    nominal value; pass through anything else (numeric literals, estimated
    parameter ids, model entity ids, array-valued fixed parameters)
    unchanged."""
    return fixed_parameter_values.get(value, value)


def _split_override_column(
    col_values: pd.Series,
    n_pars: int,
    fixed_parameter_values: dict[str, float],
) -> np.ndarray:
    """Split ``;``-separated observable/noise parameter override strings
    into a ``(n_rows, n_pars)`` matrix of numeric values / free parameter
    ids, resolving non-estimated parameter references to their nominal
    value and right-padding with ``1.0``."""

    def resolve_row(entry) -> list:
        # `col_values` may have `object` dtype with a mix of `;`-separated
        # override strings and already-numeric entries (e.g. a column
        # where only some rows use string-encoded overrides); splitting
        # via the `.str` accessor on the whole column would silently turn
        # every non-string entry into NaN; handle each entry's type here
        # instead.
        if pd.isna(entry):
            return []
        values = (
            # an empty cell (as opposed to a missing/NaN one) splits to a
            # single empty-string token; drop it rather than trying to
            # resolve/convert "" as if it were a real override
            [v for v in entry.split(petabv2.C.PARAMETER_SEPARATOR) if v]
            if isinstance(entry, str)
            else [entry]
        )
        return [
            _resolve_override_symbol(v, fixed_parameter_values)
            for v in values
        ]

    rows = col_values.apply(resolve_row)
    padded = rows.apply(
        lambda row: np.pad(
            row, (0, n_pars - len(row)), mode="constant", constant_values=1.0
        )
    )
    return np.stack(padded)


def _override_triple_from_matrix(
    mat: np.ndarray, parameter_ids: tuple[str, ...]
) -> OverrideColumn:
    """Split a raw (numeric-or-parameter-id) override matrix into an
    :class:`OverrideColumn`."""
    par_index = np.vectorize(
        lambda x: parameter_ids.index(x) if x in parameter_ids else -1
    )(mat)
    par_mask = par_index != -1
    # in-place assignment (rather than e.g. `np.where`) is required here:
    # `mat` may be a fixed-width numpy string array (not `object` dtype) if
    # every entry is a parameter reference, and `np.where(mask, 0.0, mat)`
    # then fails to find a common dtype for the float/string mix.
    mat = mat.copy()
    mat[par_mask] = 0.0
    mat = mat.astype(float)
    par_index[~par_mask] = 0
    return OverrideColumn(mat, par_mask, par_index)


def _column_overrides(
    m: pd.DataFrame,
    col: str,
    n_pars: int,
    fixed_parameter_values: dict[str, float],
    parameter_ids: tuple[str, ...],
) -> OverrideColumn:
    """Numeric values, non-numeric mask and parameter indices for one
    observable/noise parameter override column of the rows in ``m``."""
    if col not in m or m[col].isna().all() or (m[col] == "").all():
        return OverrideColumn.placeholder(len(m), n_pars)
    if pd.api.types.is_numeric_dtype(m[col].dtype):
        mat_numeric = np.expand_dims(m[col].values, axis=1)
        return OverrideColumn(
            mat_numeric,
            np.zeros_like(mat_numeric, dtype=bool),
            np.zeros_like(mat_numeric, dtype=int),
        )
    mat = _split_override_column(m[col], n_pars, fixed_parameter_values)
    return _override_triple_from_matrix(mat, parameter_ids)


def _get_overrides(
    m: pd.DataFrame,
    n_pars: dict[str, int],
    fixed_parameter_values: dict[str, float],
    parameter_ids: tuple[str, ...],
) -> dict[str, OverrideColumn]:
    """Numeric values, non-numeric mask and parameter indices for
    observable/noise parameter overrides of the rows in ``m``."""
    return {
        col: _column_overrides(
            m, col, n_pars[col], fixed_parameter_values, parameter_ids
        )
        for col in (petabv2.C.OBSERVABLE_PARAMETERS, petabv2.C.NOISE_PARAMETERS)
    }


def _get_iy_trafos(
    iys: np.ndarray, petab_problem: petabv2.Problem
) -> np.ndarray:
    """Observable transformation index (see ``SCALE_TO_INT``) for each
    observable index in ``iys``."""
    if petabv2.C.NOISE_DISTRIBUTION in petab_problem.observable_df:
        return np.array(
            [
                SCALE_TO_INT[petabv2.C.LOG]
                if obs.noise_distribution == petabv2.C.LOG_NORMAL
                else SCALE_TO_INT[petabv2.C.LIN]
                for obs in petab_problem.observables
            ]
        )
    return np.zeros_like(iys)


class _PeriodMeasurements(NamedTuple):
    """One experiment period's bucketed measurement data, as built by
    :meth:`JAXProblem._get_measurements` and consumed by its
    padding/stacking step."""

    ts_dyn: np.ndarray
    ts_posteq: np.ndarray
    my: np.ndarray
    iys: np.ndarray
    iy_trafos: np.ndarray
    op_overrides: OverrideColumn
    noise_overrides: OverrideColumn
    valid: np.ndarray


def _masked_placeholder_period(
    t: float, n_pars: dict[str, int]
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    dict[str, OverrideColumn],
]:
    """A single masked, zero-information timepoint at time ``t``.

    Used both for a real period that has no actual measurements in its
    time window (e.g. a pure post-equilibration period), and for a padding
    period that doesn't exist for a given experiment -- in both cases the
    only requirement is a valid, non-backward-in-time, fully-masked-out
    entry so downstream padding/chaining has something to work with.

    :return:
        Tuple of ``(ts_dyn, dyn_valid, my, iys, iy_trafos, overrides)``.
    """
    return (
        np.array([t]),
        np.array([False]),
        np.array([0.0]),
        np.array([0]),
        np.array([0]),
        {
            col: OverrideColumn.placeholder(1, n_pars[col])
            for col in (petabv2.C.OBSERVABLE_PARAMETERS, petabv2.C.NOISE_PARAMETERS)
        },
    )


def _pad_measurement(
    x_dyn: np.ndarray, x_peq: np.ndarray, n_ts_dyn: int, n_ts_posteq: int
) -> np.ndarray:
    """Right-pad ``x_dyn``/``x_peq`` (edge mode: repeat the last value) to
    ``n_ts_dyn``/``n_ts_posteq`` along the first axis, then concatenate."""
    pad_width_dyn = tuple(
        [(0, n_ts_dyn - len(x_dyn))] + [(0, 0)] * (x_dyn.ndim - 1)
    )
    pad_width_peq = tuple(
        [(0, n_ts_posteq - len(x_peq))] + [(0, 0)] * (x_peq.ndim - 1)
    )
    return np.concatenate(
        (
            np.pad(x_dyn, pad_width_dyn, mode="edge")
            if len(x_dyn)
            else np.zeros((n_ts_dyn, *x_dyn.shape[1:]), dtype=x_dyn.dtype),
            np.pad(x_peq, pad_width_peq, mode="edge")
            if len(x_peq)
            else np.zeros(
                (n_ts_posteq, *x_peq.shape[1:]), dtype=x_peq.dtype
            ),
        )
    )


def _pad_and_stack(
    measurements: dict[tuple[str, int], _PeriodMeasurements],
    extractor: Callable[[_PeriodMeasurements], np.ndarray],
    n_ts_dyn: int,
    n_ts_posteq: int,
) -> np.ndarray:
    """Apply ``extractor`` to every bucketed period, split each result at
    its own dynamic/post-equilibrium boundary, pad, and stack."""
    return np.stack(
        [
            _pad_measurement(
                extractor(mv)[: len(mv.ts_dyn)],
                extractor(mv)[len(mv.ts_dyn) :],
                n_ts_dyn,
                n_ts_posteq,
            )
            for mv in measurements.values()
        ]
    )


class JAXProblem(eqx.Module):
    """
    PEtab problem wrapper for JAX models.

    :ivar parameters:
        Values for the model parameters. Do not change dimensions, values may be changed during, e.g. model training.
    :ivar model:
        JAXModel instance to use for simulation.
    :ivar _parameter_mappings:
        :class:`ParameterMappingForCondition` instances for each simulation condition.
    :ivar _measurements:
        Preprocessed arrays for each simulation condition.
    :ivar _petab_problem:
        PEtab problem to simulate.
    """

    parameters: jnp.ndarray
    model: JAXModel
    simulation_conditions: tuple[tuple[str, ...], ...]
    _parameter_mappings: dict[str, ...]
    _max_periods: int
    _ts_dyn: np.ndarray
    _ts_posteq: np.ndarray
    _my: np.ndarray
    _iys: np.ndarray
    _iy_trafos: np.ndarray
    _ts_masks: np.ndarray
    _op_numeric: np.ndarray
    _op_mask: np.ndarray
    _op_indices: np.ndarray
    _np_numeric: np.ndarray
    _np_mask: np.ndarray
    _np_indices: np.ndarray
    _petab_measurement_indices: np.ndarray
    _petab_problem: petabv2.Problem
    _unconverted_problem: petabv2.Problem | None
    _all_condition_targets: frozenset[str]

    def __init__(
        self,
        model: JAXModel,
        petab_problem: petabv1.Problem | petabv2.Problem,
        unconverted_problem: petabv2.Problem | None = None,
    ):
        """
        Initialize a JAXProblem instance with a model and a PEtab problem.

        :param model:
            JAXModel instance to use for simulation.
        :param petab_problem:
            PEtab problem to simulate.
        """
        if isinstance(petab_problem, petabv1.Problem):
            raise TypeError(
                "JAXProblem does not support PEtab v1 problems. Upgrade the problem to PEtab v2."
            )
        petab_problem = add_default_experiment_names_to_v2_problem(
            petab_problem
        )
        scs = get_simulation_conditions_v2(petab_problem)
        self.simulation_conditions = scs.conditionId.to_list()
        self._petab_problem = petab_problem
        self._unconverted_problem = unconverted_problem
        self._all_condition_targets = frozenset(
            change.target_id
            for condition in self._petab_problem.conditions
            for change in condition.changes
        )
        self.parameters, self.model = (
            self._initialize_model_with_nominal_values(model)
        )
        self._parameter_mappings = self._get_parameter_mappings()
        (
            self._max_periods,
            self._ts_dyn,
            self._ts_posteq,
            self._my,
            self._iys,
            self._iy_trafos,
            self._ts_masks,
            self._petab_measurement_indices,
            self._op_numeric,
            self._op_mask,
            self._op_indices,
            self._np_numeric,
            self._np_mask,
            self._np_indices,
        ) = self._get_measurements(self._petab_problem.experiments)

    def save(self, directory: Path):
        """
        Save the problem to a directory.

        :param directory:
            Directory to save the problem to.
        """
        if self._petab_problem.config is None:
            self._petab_problem.config = petabv2.ProblemConfig(
                format_version="2.0.0"
            )
        self._petab_problem.config.filepath = "problem.yaml"
        self._petab_problem.to_files(base_path=directory)
        shutil.copy(self.model.jax_py_file, directory / "jax_py_file.py")
        with open(directory / "parameters.pkl", "wb") as f:
            eqx.tree_serialise_leaves(f, self)

    @classmethod
    def load(cls, directory: Path):
        """
        Load a problem from a directory.

        :param directory:
            Directory to load the problem from.

        :return:
            Loaded problem instance.
        """
        petab_problem = petabv2.Problem.from_yaml(
            directory / "problem.yaml",
        )
        model = _module_from_path("jax", directory / "jax_py_file.py").Model()
        problem = cls(model, petab_problem)
        with open(directory / "parameters.pkl", "rb") as f:
            return eqx.tree_deserialise_leaves(f, problem)

    def _get_measurements(
        self, experiments: list[petabv2.Experiment]
    ) -> tuple[
        int,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
    ]:
        """
        Get measurements for the model, bucketed per experiment and per
        (non-pre-equilibration) period, and padded to a common shape
        ``(n_experiments, max_periods, n_timepoints, ...)`` suitable for
        :func:`eqx.filter_vmap` over experiments with an inner per-period
        chain (see :meth:`run_simulation`).

        Periods beyond an experiment's own real period count, and (for
        non-terminal real periods) one synthetic time point at the start of
        the following period, are added so that state can always be handed
        off at exactly the period boundary; those entries are marked
        ``False`` in the returned mask and never contribute to the
        log-likelihood.

        :param experiments:
            Experiments to build measurement arrays for.
        :return:
            tuple of padded
             - dynamic time points
             - post-equilibrium time points
             - measurements
             - observable indices
             - observable transformations indices
             - measurement masks
             - data indices (index in petab measurement dataframe, -1 for
               synthetic/padding entries).
             - numeric values for observable parameter overrides
             - non-numeric mask for observable parameter overrides
             - parameter indices (problem parameters) for observable parameter overrides
             - numeric values for noise parameter overrides
             - non-numeric mask for noise parameter overrides
             - parameter indices (problem parameters) for noise parameter overrides
        """
        measurements: dict[tuple[str, int], _PeriodMeasurements] = {}
        petab_indices = dict()

        # Nominal (linear) values of fixed (non-estimated) parameters, used to
        # resolve observable/noise parameter overrides that reference them.
        # Built once up front (rather than lazily per override column) since
        # it is otherwise rebuilt once per period.
        fixed_parameter_values = _get_fixed_parameter_values(
            self._petab_problem
        )

        n_pars = dict()
        for col in [
            petabv2.C.OBSERVABLE_PARAMETERS,
            petabv2.C.NOISE_PARAMETERS,
        ]:
            n_pars[col] = 0
            if col in self._petab_problem.measurement_df:
                if pd.api.types.is_numeric_dtype(
                    self._petab_problem.measurement_df[col].dtype
                ):
                    n_pars[col] = 1 - int(
                        self._petab_problem.measurement_df[col].isna().all()
                    )
                else:
                    n_pars[col] = (
                        self._petab_problem.measurement_df[col]
                        .str.split(petabv2.C.PARAMETER_SEPARATOR)
                        .apply(
                            lambda x: (
                                len(x)
                                if isinstance(x, Sized)
                                else 1 - int(pd.isna(x))
                            )
                        )
                        .max()
                    )

        dyn_periods_by_exp = {
            exp.id: [
                period
                for period in exp.sorted_periods
                if not period.is_preequilibration
            ]
            for exp in experiments
        }
        max_periods = max(
            (len(periods) for periods in dyn_periods_by_exp.values()),
            default=0,
        )
        max_periods = max(max_periods, 1)

        for exp in experiments:
            dyn_periods = dyn_periods_by_exp[exp.id]

            # All measurements for this experiment. Periods are
            # distinguished by time window below since PEtab v2
            # measurements reference an experiment, not an individual
            # condition/period.
            query = f"{petabv2.C.EXPERIMENT_ID} == '{exp.id}'"
            m_full = self._petab_problem.measurement_df.query(
                query
            ).sort_values(by=petabv2.C.TIME)
            m_full_times = m_full[petabv2.C.TIME].to_numpy()

            last_period_time = dyn_periods[-1].time if dyn_periods else 0.0

            for i_period in range(max_periods):
                if i_period >= len(dyn_periods):
                    # Padding slot: this experiment has no period here at
                    # all. A single masked, zero-duration integration step
                    # that leaves the carried-over state unchanged.
                    ts_dyn, dyn_valid, my, iys, iy_trafos, overrides = (
                        _masked_placeholder_period(last_period_time, n_pars)
                    )
                    ts_posteq = np.array([])
                    posteq_valid = np.array([])
                    index = [-1]
                else:
                    is_own_last = i_period == len(dyn_periods) - 1
                    t_lo = dyn_periods[i_period].time
                    t_hi = (
                        np.inf
                        if is_own_last
                        else dyn_periods[i_period + 1].time
                    )
                    in_window = (
                        (m_full_times >= t_lo)
                        & (m_full_times < t_hi)
                        & np.isfinite(m_full_times)
                    )
                    m = m_full[in_window]

                    ts_dyn_real = m[petabv2.C.TIME].values
                    my_real = m[petabv2.C.MEASUREMENT].values
                    iys_real = np.array(
                        [
                            self.model.observable_ids.index(oid)
                            for oid in m[petabv2.C.OBSERVABLE_ID].values
                        ]
                    )
                    iy_trafos_real = _get_iy_trafos(
                        iys_real, self._petab_problem
                    )
                    overrides_real = _get_overrides(
                        m, n_pars, fixed_parameter_values, self.parameter_ids
                    )
                    index_dyn = list(m.index)

                    if is_own_last:
                        posteq_mask = np.isfinite(m_full_times) == False  # noqa: E712
                        m_posteq = m_full[posteq_mask]
                        ts_posteq = m_posteq[petabv2.C.TIME].values
                        index_posteq = list(m_posteq.index)

                        # No further period to hand off to within this
                        # experiment's own chain; padding slots (if any)
                        # after this one are pure no-ops.
                        if len(ts_dyn_real):
                            ts_dyn = ts_dyn_real
                            dyn_valid = np.ones(len(ts_dyn_real), dtype=bool)
                            my = my_real
                            iys = iys_real
                            iy_trafos = iy_trafos_real
                            overrides = overrides_real
                            index_dyn_full = index_dyn
                        else:
                            # No real dyn measurements in this period (e.g.
                            # a pure post-equilibration period): still need
                            # a valid (masked), non-backward-in-time entry
                            # so that padding never produces a time point
                            # before this period's start.
                            (
                                ts_dyn,
                                dyn_valid,
                                my,
                                iys,
                                iy_trafos,
                                overrides,
                            ) = _masked_placeholder_period(t_lo, n_pars)
                            index_dyn_full = [-1]
                    else:
                        ts_posteq = np.array([])
                        index_posteq = []

                        # Append a synthetic, masked boundary time point at
                        # the start of the next period so that the ODE
                        # state is always available at exactly the period
                        # boundary, regardless of whether there is a real
                        # measurement there.
                        ts_dyn = np.append(ts_dyn_real, t_hi)
                        dyn_valid = np.append(
                            np.ones(len(ts_dyn_real), dtype=bool), False
                        )
                        my = np.append(my_real, 0.0)
                        iys = np.append(iys_real, 0)
                        iy_trafos = np.append(iy_trafos_real, 0)
                        overrides = {
                            col: OverrideColumn.concatenate(
                                overrides_real[col],
                                OverrideColumn.placeholder(1, n_pars[col]),
                            )
                            for col in overrides_real
                        }
                        index_dyn_full = [*index_dyn, -1]

                    index = [*index_dyn_full, *index_posteq]
                    posteq_valid = np.ones(len(ts_posteq), dtype=bool)

                valid = np.concatenate([dyn_valid, posteq_valid]).astype(
                    bool
                )

                measurements[(exp.id, i_period)] = _PeriodMeasurements(
                    ts_dyn=ts_dyn,
                    ts_posteq=ts_posteq,
                    my=my,
                    iys=iys,
                    iy_trafos=iy_trafos,
                    op_overrides=overrides[petabv2.C.OBSERVABLE_PARAMETERS],
                    noise_overrides=overrides[petabv2.C.NOISE_PARAMETERS],
                    valid=valid,
                )
                petab_indices[(exp.id, i_period)] = tuple(index)

        # compute maximum lengths
        n_ts_dyn = max(len(mv.ts_dyn) for mv in measurements.values())
        n_ts_posteq = max(len(mv.ts_posteq) for mv in measurements.values())

        # pad with last value and stack
        ts_dyn = np.stack(
            [
                np.pad(mv.ts_dyn, (0, n_ts_dyn - len(mv.ts_dyn)), mode="edge")
                if len(mv.ts_dyn)
                else np.zeros(n_ts_dyn, dtype=mv.ts_dyn.dtype)
                for mv in measurements.values()
            ]
        )
        ts_posteq = np.stack(
            [
                np.pad(
                    mv.ts_posteq,
                    (0, n_ts_posteq - len(mv.ts_posteq)),
                    mode="edge",
                )
                if len(mv.ts_posteq)
                else np.zeros(n_ts_posteq, dtype=mv.ts_posteq.dtype)
                for mv in measurements.values()
            ]
        )

        my = _pad_and_stack(
            measurements, lambda mv: mv.my, n_ts_dyn, n_ts_posteq
        )
        iys = _pad_and_stack(
            measurements, lambda mv: mv.iys, n_ts_dyn, n_ts_posteq
        )
        iy_trafos = _pad_and_stack(
            measurements, lambda mv: mv.iy_trafos, n_ts_dyn, n_ts_posteq
        )
        op_numeric = _pad_and_stack(
            measurements,
            lambda mv: mv.op_overrides.numeric,
            n_ts_dyn,
            n_ts_posteq,
        )
        op_mask = _pad_and_stack(
            measurements,
            lambda mv: mv.op_overrides.mask,
            n_ts_dyn,
            n_ts_posteq,
        )
        op_indices = _pad_and_stack(
            measurements,
            lambda mv: mv.op_overrides.index,
            n_ts_dyn,
            n_ts_posteq,
        )
        np_numeric = _pad_and_stack(
            measurements,
            lambda mv: mv.noise_overrides.numeric,
            n_ts_dyn,
            n_ts_posteq,
        )
        np_mask = _pad_and_stack(
            measurements,
            lambda mv: mv.noise_overrides.mask,
            n_ts_dyn,
            n_ts_posteq,
        )
        np_indices = _pad_and_stack(
            measurements,
            lambda mv: mv.noise_overrides.index,
            n_ts_dyn,
            n_ts_posteq,
        )
        # mask padding must stay `False` (not repeat the last real mask
        # value), so this is stacked directly rather than via
        # `_pad_and_stack` (which pads in "edge" mode).
        ts_masks = np.stack(
            [
                np.concatenate(
                    (
                        np.pad(
                            mv.valid[: len(mv.ts_dyn)],
                            (0, n_ts_dyn - len(mv.ts_dyn)),
                        ),
                        np.pad(
                            mv.valid[len(mv.ts_dyn) :],
                            (0, n_ts_posteq - len(mv.ts_posteq)),
                        ),
                    )
                )
                for mv in measurements.values()
            ]
        ).astype(bool)
        petab_indices = np.stack(
            [
                _pad_measurement(
                    np.array(idx[: len(mv.ts_dyn)]),
                    np.array(idx[len(mv.ts_dyn) :]),
                    n_ts_dyn,
                    n_ts_posteq,
                )
                for mv, idx in zip(
                    measurements.values(), petab_indices.values()
                )
            ]
        )

        n_exp = len(experiments)
        outputs = (
            ts_dyn,
            ts_posteq,
            my,
            iys,
            iy_trafos,
            ts_masks,
            petab_indices,
            op_numeric,
            op_mask,
            op_indices,
            np_numeric,
            np_mask,
            np_indices,
        )
        return (
            max_periods,
            *(
                arr.reshape(n_exp, max_periods, *arr.shape[1:])
                for arr in outputs
            ),
        )

    def _get_parameter_mappings(self) -> dict[str, ...]:
        # `targets_map` intentionally stores each value only lightly parsed
        # (a numeric literal, or a parameter id `str`, via
        # `_resolve_petab_change_value`) rather than eagerly resolved
        # against `self.parameters`: this dict is cached once at
        # construction time (`self._parameter_mappings`, see `__init__`),
        # and `update_parameters` (`eqx.tree_at(lambda p: p.parameters,
        # self, p)`) does not recompute it, so an eagerly-resolved value
        # would go stale after a parameter update.
        # `_resolve_parameter_reference` re-resolves the actual value live
        # against `self.parameters` at use time instead.
        targets_map = {
            c.id: {
                ch.target_id: _resolve_petab_change_value(ch.target_value)
                for ch in c.changes
            }
            for c in self._petab_problem.conditions
        }

        hybrid_map = {}
        if self._petab_problem.config.extensions.get("sciml"):
            hybridization_df = (
                self._petab_problem.extensions.sciml.hybridization_df
            )
            hybrid_map = (
                hybridization_df.set_index("targetId")["targetValue"]
                .astype(str)
                .to_dict()
            )

        return {"targets_map": targets_map, "hybrid_map": hybrid_map}

    def _resolve_parameter_reference(
        self, value: float | str
    ) -> jt.Float[jt.Scalar, ""]:  # noqa: F722
        """
        Resolve a value from ``targets_map`` (see :meth:`_get_parameter_mappings`)
        to a JAX scalar. Numeric values pass through; a parameter id
        resolves to that parameter's current (estimated) or nominal
        (fixed) value. Symbolic references to estimated parameters keep
        their dependence on :attr:`parameters` so gradients flow through
        correctly.

        :param value:
            A numeric literal or PEtab parameter id, as returned by
            :func:`_resolve_petab_change_value`.
        """
        if isinstance(value, str):
            if value in self.parameter_ids:
                return self.parameters[self.parameter_ids.index(value)]
            return jnp.asarray(
                self._petab_problem.parameter_df.loc[
                    value, petabv2.C.NOMINAL_VALUE
                ],
                dtype=self.model.parameters.dtype,
            )
        return jnp.asarray(value, dtype=self.model.parameters.dtype)

    def get_all_simulation_conditions(self) -> tuple[tuple[str, ...], ...]:
        simulation_conditions = get_simulation_conditions_v2(
            self._petab_problem
        )
        return tuple(
            tuple([row.conditionId])
            for _, row in simulation_conditions.iterrows()
        )

    def _initialize_model_parameters(self, model: JAXModel) -> dict:
        """
        Initialize model parameter structure with zeros.

        :param model:
            JAX model with neural networks

        :return:
            Nested dictionary structure for model parameters
        """
        return {
            net_id: {
                layer_id: {
                    attribute: jnp.zeros_like(getattr(layer, attribute))
                    for attribute in ["weight", "bias"]
                    if hasattr(layer, attribute)
                }
                for layer_id, layer in nn.layers.items()
            }
            for net_id, nn in model.nns.items()
        }

    def _load_parameter_arrays_from_files(self) -> dict:
        """
        Load neural network parameter arrays from HDF5 files.

        :return:
            Dictionary mapping network IDs to parameter arrays
        """
        if not self._petab_problem.config.extensions:
            return {}

        array_files = self._petab_problem.config.extensions.get(
            "sciml", {}
        ).array_files

        import h5py

        net_ids = list(
            self._petab_problem.config.extensions[
                "sciml"
            ].neural_networks.keys()
        )

        # TODO(performance): Avoid opening each file multiple times
        return {
            net_id: h5py.File(file_spec, "r")["parameters"][net_id]
            for file_spec in array_files
            for net_id in net_ids
            if "parameters" in h5py.File(file_spec, "r").keys()
            and net_id in h5py.File(file_spec, "r")["parameters"].keys()
        }

    def _load_input_arrays_from_files(self) -> dict:
        """
        Load neural network input arrays from HDF5 files.

        :return:
            Dictionary mapping network IDs to input arrays
        """
        if not self._petab_problem.config.extensions:
            return {}

        array_files = self._petab_problem.config.extensions.get(
            "sciml", {}
        ).array_files

        import h5py

        net_ids = list(
            self._petab_problem.config.extensions[
                "sciml"
            ].neural_networks.keys()
        )

        # TODO(performance): Avoid opening each file multiple times
        return {
            net_id: h5py.File(file_spec, "r")["inputs"]
            for file_spec in array_files
            for net_id in net_ids
            if "inputs" in h5py.File(file_spec, "r").keys()
        }

    def _parse_parameter_name(
        self, pname: str, model_pars: dict
    ) -> list[tuple[str, str]]:
        """
        Parse parameter name to determine which layers and attributes to set.

        :param pname:
            Parameter name from PEtab (format: net.layer.attribute)
        :param model_pars:
            Model parameters dictionary

        :return:
            List of (layer_name, attribute_name) tuples to set
        """
        net = pname.split("_")[0]
        nn = model_pars[net]
        to_set = []

        name_parts = pname.split(".")

        if len(name_parts) > 1:
            layer_name = name_parts[1]
            layer = nn[layer_name]
            if len(name_parts) > 2:
                # Specific attribute specified
                attribute_name = name_parts[2]
                to_set.append((layer_name, attribute_name))
            else:
                # All attributes of the layer
                to_set.extend(
                    [(layer_name, attribute) for attribute in layer.keys()]
                )
        else:
            # All layers and attributes
            to_set.extend(
                [
                    (layer_name, attribute)
                    for layer_name, layer in nn.items()
                    for attribute in layer.keys()
                ]
            )

        return to_set

    def _extract_nominal_values_from_petab(
        self, model: JAXModel, model_pars: dict, par_arrays: dict
    ) -> None:
        """
        Extract nominal parameter values from PEtab problem and populate model_pars.

        :param model:
            JAX model
        :param model_pars:
            Model parameters dictionary to populate (modified in place)
        :param par_arrays:
            Parameter arrays loaded from files
        """
        mapping_df = self._petab_problem.mapping_df

        def _lookup_mid(pname: str) -> str:
            if mapping_df is not None and pname in mapping_df.index:
                return mapping_df.loc[pname, petabv2.C.MODEL_ENTITY_ID]
            return ""

        for table in self._petab_problem.parameter_tables:
            # Array-indexed params (e.g. net.parameters[layer]) must be processed
            # after scalar params to avoid overwriting a full-layer assignment.
            def _sort_key(p):
                model_id = str(_lookup_mid(p.id))
                if "parameters[" not in model_id:
                    return 0
                return model_id.count(".") + model_id.count("[")

            sorted_params = sorted(table.elements, key=_sort_key)

            for parameter in sorted_params:
                pname = parameter.id
                net = pname.split("_")[0]
                if net not in model.nns:
                    continue

                nn = model_pars[net]
                scalar = True

                if parameter.nominal_value == "array":
                    value = par_arrays[net]
                    scalar = False
                else:
                    value = float(parameter.nominal_value)

                model_entity_id = _lookup_mid(pname)
                if model_entity_id != "" and "parameters[" in str(
                    model_entity_id
                ):
                    to_set = _parse_model_entity_id(model_entity_id, nn)
                else:
                    to_set = self._parse_parameter_name(pname, model_pars)

                for layer, attribute in to_set:
                    if scalar:
                        nn[layer][attribute] = value * jnp.ones_like(
                            getattr(model.nns[net].layers[layer], attribute)
                        )
                    else:
                        nn[layer][attribute] = jnp.array(
                            value[layer][attribute][:]
                        )

    def _set_model_parameters(
        self, model: JAXModel, model_pars: dict
    ) -> JAXModel:
        """
        Set parameter values in the model using equinox tree_at.

        :param model:
            JAX model to update
        :param model_pars:
            Dictionary of parameter values to set

        :return:
            Updated JAX model
        """
        for net_id in model_pars:
            for layer_id in model_pars[net_id]:
                for attribute in model_pars[net_id][layer_id]:
                    logger.debug(
                        f"Setting {attribute} of layer {layer_id} in network "
                        f"{net_id} to {model_pars[net_id][layer_id][attribute]}"
                    )
                    model = eqx.tree_at(
                        lambda model: getattr(
                            model.nns[net_id].layers[layer_id], attribute
                        ),
                        model,
                        model_pars[net_id][layer_id][attribute],
                    )
        return model

    def _set_input_arrays(
        self, model: JAXModel, nn_input_arrays: dict, model_pars: dict
    ) -> JAXModel:
        """
        Set input arrays in the model if provided.

        :param model:
            JAX model to update
        :param nn_input_arrays:
            Input arrays loaded from files
        :param model_pars:
            Model parameters dictionary (for network IDs)

        :return:
            Updated JAX model
        """
        if len(nn_input_arrays) == 0:
            return model

        array_inputs = self._petab_problem.extensions.sciml.hybridization_df[
            self._petab_problem.extensions.sciml.hybridization_df[
                "targetValue"
            ]
            == "array"
        ].index.tolist()

        for net_id in model_pars:
            net_id_inputs = array_inputs + model.nns[net_id].inputs
            input_array = {
                input: {
                    k: jnp.array(
                        arr[:],
                        dtype=jnp.float64
                        if jax.config.jax_enable_x64
                        else jnp.float32,
                    )
                    for k, arr in nn_input_arrays[net_id][input].items()
                }
                for input in net_id_inputs
                if input in nn_input_arrays[net_id]
            }
            model = eqx.tree_at(
                lambda model: model.nns[net_id].inputs, model, input_array
            )

        return model

    def _initialize_model_with_nominal_values(
        self, model: JAXModel
    ) -> tuple[jt.Float[jt.Array, "np"], JAXModel]:
        """
        Initialize the model with nominal parameter values and inputs from the PEtab problem.

        This method:
        - Initializes model parameter structure
        - Loads parameter and input arrays from HDF5 files
        - Extracts nominal values from PEtab problem
        - Sets parameter values in the model
        - Sets input arrays in the model
        - Creates scaled parameter array to initialized to nominal values

        :param model:
            JAX model to initialize

        :return:
            Tuple of (scaled parameter array, initialized model)
        """
        # Initialize model parameters structure
        model_pars = self._initialize_model_parameters(model)

        # Load arrays from files (getters)
        par_arrays = self._load_parameter_arrays_from_files()
        nn_input_arrays = self._load_input_arrays_from_files()

        # Extract nominal values from PEtab problem
        self._extract_nominal_values_from_petab(model, model_pars, par_arrays)

        # Set values in model (setters)
        model = self._set_model_parameters(model, model_pars)
        model = self._set_input_arrays(model, nn_input_arrays, model_pars)

        # Create scaled parameter array
        param_map = self._petab_problem.get_x_nominal_dict()
        parameter_array = jnp.array(
            [float(param_map[pval]) for pval in self.parameter_ids]
        )

        return parameter_array, model

    @property
    def parameter_ids(self) -> list[str]:
        """
        Parameter ids that are estimated in the PEtab problem. Same ordering as values in :attr:`parameters`.

        :return:
            PEtab parameter ids
        """
        return [
            p.id
            for pt in self._petab_problem.parameter_tables
            for p in pt.elements
            if p.estimate and p.nominal_value != "array"
        ]

    @property
    def nn_output_ids(self) -> list[str]:
        """
        Parameter ids that are estimated in the PEtab problem. Same ordering as values in :attr:`parameters`.

        :return:
            PEtab parameter ids
        """
        if self._petab_problem.mapping_df is None:
            return []
        # A mapping table is not exclusive to neural nets: e.g. PySB problems
        # map PEtab entities to model expressions like ``A_() ** compartment``
        # (no ``.``). Match only ``<net>.output...``-style entities and skip
        # everything else, rather than relying on vectorized ``.str`` accessors
        # that break when an entity has no ``.`` (``str[1]`` yields a float NaN
        # series).
        return [
            petab_entity_id
            for petab_entity_id, model_entity_id in self._petab_problem.mapping_df[
                petabv2.C.MODEL_ENTITY_ID
            ].items()
            if isinstance(model_entity_id, str)
            and len(parts := model_entity_id.split(".")) > 1
            and parts[1].startswith("output")
        ]

    def get_petab_parameter_by_id(self, name: str) -> jnp.float_:
        """
        Get the value of a PEtab parameter by name.

        :param name:
            PEtab parameter id, as returned by :attr:`parameter_ids`.
        :return:
            Value of the parameter
        """
        return self.parameters[self.parameter_ids.index(name)]

    def _unscale(
        self, p: jt.Float[jt.Array, "np"], scales: tuple[str, ...]
    ) -> jt.Float[jt.Array, "np"]:
        """
        Unscaling of parameters.

        :param p:
            Parameter values
        :param scales:
            Parameter scalings
        :return:
            Unscaled parameter values
        """
        return jnp.array(
            [jax_unscale(pval, scale) for pval, scale in zip(p, scales)]
        )

    def _resolve_original_condition_id(self, condition_id: str) -> str:
        """Map a converted condition ID back to its original unconverted ID."""
        if self._unconverted_problem is None:
            return condition_id
        for orig_exp, conv_exp in zip(
            self._unconverted_problem.experiments,
            self._petab_problem.experiments,
        ):
            for orig_period, conv_period in zip(
                orig_exp.sorted_periods, conv_exp.sorted_periods
            ):
                if (
                    condition_id in conv_period.condition_ids
                    and orig_period.condition_ids
                ):
                    return orig_period.condition_ids[0]
        return condition_id

    def _eval_nn(self, output_par: str, condition_id: str):
        entity_id = self._petab_problem.mapping_df.loc[
            output_par, petabv2.C.MODEL_ENTITY_ID
        ]
        net_id = entity_id.split(".")[0]
        ind = int(re.search(r"\[\d+\]\[(\d+)\]", entity_id).group(1))
        nn = self.model.nns[net_id]
        original_condition_id = self._resolve_original_condition_id(
            condition_id
        )

        def _is_net_input(model_id):
            comps = model_id.split(".")
            return comps[0] == net_id and comps[1].startswith("inputs")

        model_id_map = (
            self._petab_problem.mapping_df[
                self._petab_problem.mapping_df[
                    petabv2.C.MODEL_ENTITY_ID
                ].apply(_is_net_input)
            ]
            .reset_index()
            .set_index(petabv2.C.MODEL_ENTITY_ID)[petabv2.C.PETAB_ENTITY_ID]
            .to_dict()
        )
        petab_ids = set(model_id_map.values())

        parameters_map = self._petab_problem.get_x_nominal_dict()
        parameters_map.update(zip(self.parameter_ids, self.parameters))

        condition_input_map = {
            pid: parameters_map[pid]
            for pid in petab_ids
            if pid in parameters_map
        }
        condition_input_map.update(
            {
                pid: parameters_map[target]
                for pid, target in self._parameter_mappings[
                    "hybrid_map"
                ].items()
                if pid in petab_ids and target in parameters_map
            }
        )

        nn_inputs = getattr(nn, "inputs", {})

        # Build a map from input slot index to the values at that slot.
        # Scalar inputs (e.g. "net1.inputs[1]") have a single None-keyed entry.
        # Vector inputs (e.g. "net1.inputs[0][0]", "net1.inputs[0][1]") have
        # one integer-keyed entry per element.
        # slot_idx -> {element_idx_or_None: value}
        input_slots: dict[int, dict[int | None, object]] = {}
        for model_id, petab_id in model_id_map.items():
            if petab_id in condition_input_map:
                val = condition_input_map[petab_id]
            elif (
                petab_id in nn_inputs
                and original_condition_id in nn_inputs[petab_id]
            ):
                val = nn_inputs[petab_id][original_condition_id]
            else:
                val = nn_inputs[petab_id]["0"]

            idx_strs = re.findall(r"\[(\d+)\]", model_id.split(".inputs")[1])
            slot_idx = int(idx_strs[0])
            element_idx = int(idx_strs[1]) if len(idx_strs) > 1 else None
            if slot_idx not in input_slots:
                input_slots[slot_idx] = {}
            input_slots[slot_idx][element_idx] = val

        def _slot_to_array(slot: dict[int | None, object]) -> jnp.ndarray:
            if None in slot:
                # Scalar input: wrap the single value
                return jnp.asarray(slot[None])
            # Vector input: assemble elements in index order
            return jnp.array([slot[i] for i in sorted(slot)])

        sorted_slots = sorted(input_slots.keys())
        if len(sorted_slots) == 1:
            net_input = _slot_to_array(input_slots[sorted_slots[0]])
        else:
            net_input = [_slot_to_array(input_slots[k]) for k in sorted_slots]

        return nn.forward(net_input)[ind].squeeze()

    def _dynamic_periods(
        self, experiment: petabv2.Experiment
    ) -> list[petabv2.ExperimentPeriod]:
        """Non-pre-equilibration periods of ``experiment``, in time order."""
        return [
            period
            for period in experiment.sorted_periods
            if not period.is_preequilibration
        ]

    def load_model_parameters(
        self,
        experiment: petabv2.Experiment,
        is_preeq: bool,
        period_index: int | None = None,
    ) -> jt.Float[jt.Array, "np"]:
        """
        Load parameters for an experiment.

        :param experiment:
            Experiment to load parameters for.
        :param is_preeq:
            Whether to load preequilibration or simulation parameters.
        :param period_index:
            Index (within :meth:`_dynamic_periods`) of the period to load
            simulation parameters for. Ignored if ``is_preeq`` is ``True``.
            Indices beyond the experiment's own period count are clamped to
            its last real period (used for padding periods, which are never
            actually integrated).
        :return:
            Parameters for the experiment.
        """
        if not self.model.parameter_ids:
            # a model with no free SBML parameters (e.g. only literal rate
            # constants); `jnp.stack` of an empty sequence raises.
            return jnp.array([])
        p = jnp.stack(
            [
                self._map_experiment_model_parameter_value(
                    pname, ind, experiment, is_preeq, period_index
                )
                for ind, pname in enumerate(self.model.parameter_ids)
            ]
        )
        pscale = tuple([petabv2.C.LIN for _ in self.model.parameter_ids])

        return self._unscale(p, pscale)

    def _map_experiment_model_parameter_value(
        self,
        pname: str,
        p_index: int,
        experiment: petabv2.Experiment,
        is_preeq: bool,
        period_index: int | None = None,
    ):
        """
        Get values for the given parameter `pname` from the relevant petab tables.

        :param pname: PEtab parameter id
        :param p_index: Index of the parameter in the model's parameter list
        :param experiment: PEtab experiment
        :param is_preeq: Whether to get preequilibration or simulation parameter value
        :param period_index: see :meth:`load_model_parameters`
        :return: Value of the parameter
        """
        if period_index is not None:
            dyn_periods = self._dynamic_periods(experiment)
            clamped_index = min(period_index, len(dyn_periods) - 1)
            condition_ids = dyn_periods[clamped_index].condition_ids
        else:
            # Find the first period matching the requested phase (preeq vs. sim)
            condition_ids = []
            for period in experiment.sorted_periods:
                if period.is_preequilibration == is_preeq:
                    condition_ids = period.condition_ids
                    break

        _petab_param_map = {
            param.id: param.nominal_value
            for param in self._petab_problem.parameters
        }
        if pname in self.parameter_ids:
            init_val = self.parameters[self.parameter_ids.index(pname)]
        elif pname in _petab_param_map:
            init_val = _petab_param_map[pname]
        else:
            init_val = self.model.parameters[p_index]

        targets_filtered = {
            param: value
            for condition, target in self._parameter_mappings[
                "targets_map"
            ].items()
            for param, value in target.items()
            if condition in condition_ids
        }

        if pname in targets_filtered:
            return self._resolve_parameter_reference(targets_filtered[pname])
        elif pname in self._parameter_mappings["hybrid_map"]:
            return jnp.asarray(
                self._eval_nn(
                    self._parameter_mappings["hybrid_map"][pname],
                    condition_ids[0],
                ),
                dtype=self.model.parameters.dtype,
            )
        else:
            for placeholder_attr, param_attr in (
                ("observable_placeholders", "observable_parameters"),
                ("noise_placeholders", "noise_parameters"),
            ):
                # params_list is the same for all observables; compute once
                params_list = getattr(
                    self._petab_problem.measurements[0], param_attr
                )
                for observable in self._petab_problem.observables:
                    for i, placeholder in enumerate(
                        getattr(observable, placeholder_attr)
                    ):
                        if str(placeholder) == pname:
                            return self._find_val(str(params_list[i]))
        return jnp.asarray(init_val, dtype=self.model.parameters.dtype)

    def _find_val(self, param_entry: str):
        val_float = _try_float(param_entry)
        if isinstance(val_float, float):
            return jnp.asarray(val_float, dtype=self.model.parameters.dtype)
        elif param_entry in self.parameter_ids:
            return jnp.asarray(
                self.parameters[self.parameter_ids.index(param_entry)],
                dtype=self.model.parameters.dtype,
            )
        else:
            return jnp.asarray(param_entry, dtype=self.model.parameters.dtype)

    def _first_condition_value(
        self, condition_ids: list[str], state_id: str
    ):
        """
        Find the value assigned to ``state_id`` by the first of
        ``condition_ids`` (applied simultaneously, e.g. for one experiment
        period) that actually sets it.

        Numeric values are returned as plain Python ``float``s; a
        reference to a single other parameter is returned as that
        parameter's id (``str``). Compound symbolic expressions are not
        supported.

        :param condition_ids:
            Condition IDs to check, in priority order.
        :param state_id:
            State (or parameter) id to look up.
        :return:
            The assigned value, or ``None`` if none of ``condition_ids``
            sets ``state_id``.
        """
        for condition_id in condition_ids:
            for change in self._petab_problem[condition_id].changes:
                if change.target_id != state_id:
                    continue
                return _resolve_petab_change_value(change.target_value)
        return None

    def _state_needs_reinitialisation(
        self,
        condition_ids: list[str],
        state_id: str,
    ) -> bool:
        """
        Check if a state needs reinitialisation for the given (simultaneous)
        conditions.

        :param condition_ids:
            condition ids (e.g. of one experiment period) to check
            reinitialisation for
        :param state_id:
            state id to check reinitialisation for
        :return:
            True if state needs reinitialisation, False otherwise
        """
        if state_id in self.nn_output_ids:
            return True

        if state_id in self._parameter_mappings["hybrid_map"]:
            return True

        return self._first_condition_value(condition_ids, state_id) is not None

    def _state_reinitialisation_value(
        self,
        condition_ids: list[str],
        state_id: str,
        p: jt.Float[jt.Array, "np"],
    ) -> jt.Float[jt.Scalar, ""] | float:  # noqa: F722
        """
        Get the reinitialisation value for a state.

        :param condition_ids:
            condition ids (e.g. of one experiment period) to get the
            reinitialisation value for
        :param state_id:
            state id to get reinitialisation value for
        :param p:
            parameters for the simulation condition
        :return:
            reinitialisation value for the state
        """
        if state_id in self.nn_output_ids:
            return self._eval_nn(state_id, condition_ids[0])

        if state_id in self._parameter_mappings["hybrid_map"]:
            return self._eval_nn(
                self._parameter_mappings["hybrid_map"][state_id],
                condition_ids[0],
            )

        if (
            xval := self._first_condition_value(condition_ids, state_id)
        ) is None:
            # no reinitialisation, return dummy value
            return 0.0
        if isinstance(xval, Number):
            # numerical value, return as is
            return xval
        if xval in self.model.parameter_ids:
            # model parameter, return value
            return p[self.model.parameter_ids.index(xval)]
        if xval in self.parameter_ids:
            # estimated PEtab parameter, return unscaled value. PEtab v2
            # has no parameterScale column -- all parameters are linear.
            return jax_unscale(
                self.get_petab_parameter_by_id(xval), petabv2.C.LIN
            )
        # only remaining option is nominal value for PEtab parameter
        # that is not estimated, return nominal value
        return self._petab_problem.parameter_df.loc[
            xval, petabv2.C.NOMINAL_VALUE
        ]

    def load_reinitialisation(
        self,
        condition_ids: list[str] | str,
        p: jt.Float[jt.Array, "np"],
    ) -> tuple[jt.Bool[jt.Array, "nx"], jt.Float[jt.Array, "nx"]]:  # noqa: F821
        """
        Load reinitialisation values and mask for the state vector for the
        given (simultaneous) conditions.

        :param condition_ids:
            Condition id(s) (e.g. of one experiment period) to load
            reinitialisation for. A bare string is treated as a
            single-element list. An empty list means no reinitialisation
            (e.g. for a padding period).
        :param p:
            Parameters for the simulation condition.
        :return:
            Tuple of reinitialisation masm and value for states.
        """
        if isinstance(condition_ids, str):
            condition_ids = [condition_ids]

        has_reinitialisable_states = any(
            x_id in self._all_condition_targets
            or hasattr(self, "nn_output_ids")
            and x_id in self._parameter_mappings["hybrid_map"]
            for x_id in self.model.state_ids
        )
        if not has_reinitialisable_states:
            return jnp.array([]), jnp.array([])

        if not condition_ids:
            # padding period: no reinitialisation, but still return
            # full-shaped, all-False/all-zero arrays for stacking
            nx = len(self.model.state_ids)
            return jnp.zeros(nx, dtype=bool), jnp.zeros(nx)

        mask = jnp.array(
            [
                self._state_needs_reinitialisation(condition_ids, x_id)
                for x_id in self.model.state_ids
            ]
        )
        reinit_x = jnp.array(
            [
                self._state_reinitialisation_value(condition_ids, x_id, p)
                for x_id in self.model.state_ids
            ]
        )
        return mask, reinit_x

    def update_parameters(self, p: jt.Float[jt.Array, "np"]) -> "JAXProblem":
        """
        Update parameters for the model.

        :param p:
            New problem instance with updated parameters.
        """
        return eqx.tree_at(lambda p: p.parameters, self, p)

    def _prepare_experiments(
        self,
        experiments: list[petabv2.Experiment],
        is_preeq: bool,
        op_numeric: np.ndarray | None = None,
        op_mask: np.ndarray | None = None,
        op_indices: np.ndarray | None = None,
        np_numeric: np.ndarray | None = None,
        np_mask: np.ndarray | None = None,
        np_indices: np.ndarray | None = None,
    ) -> tuple[
        jt.Float[jt.Array, "nc ... np"],  # noqa: F821, F722
        jt.Bool[jt.Array, "nx"],  # noqa: F821
        jt.Float[jt.Array, "nx"],  # noqa: F821
        jt.Float[jt.Array, "nc ... nt nop"],  # noqa: F821, F722
        jt.Float[jt.Array, "nc ... nt nnp"],  # noqa: F821, F722
    ]:
        """
        Prepare experiments for simulation.

        For the main simulation (``is_preeq=False``), all returned arrays
        (except ``h_mask``) gain a period axis right after the experiment
        axis, of static size :attr:`_max_periods`, so that
        :meth:`JAXModel.simulate_condition` can chain one ODE integration
        per period. Periods beyond an experiment's own period count are
        padding periods: they are never reinitialised and never contribute
        to the log-likelihood (see :meth:`_get_measurements`), but must
        still resolve to *some* valid parameter vector, so they simply
        reuse the experiment's own last real period.

        :param experiments:
            Experiments to prepare simulation arrays for.
        :param is_preeq:
            Whether to load preequilibration or simulation parameters.
        :param op_numeric:
            Numeric values for observable parameter overrides. If None, no overrides are used.
        :param op_mask:
            Mask for observable parameter overrides. True for free parameter overrides, False for numeric values.
        :param op_indices:
            Free parameter indices (wrt. `self.parameters`) for observable parameter overrides.
        :param np_numeric:
            Numeric values for noise parameter overrides. If None, no overrides are used.
        :param np_mask:
            Mask for noise parameter overrides. True for free parameter overrides, False for numeric values.
        :param np_indices:
            Free parameter indices (wrt. `self.parameters`) for noise parameter overrides.
        :return:
            Tuple of parameter arrays, reinitialisation masks and reinitialisation values, observable parameters and
            noise parameters.
        """
        if is_preeq:
            p_array = jnp.stack(
                [
                    self.load_model_parameters(exp, is_preeq=True)
                    for exp in experiments
                ]
            )
            t_zeros = jnp.stack(
                [
                    exp.periods[0].time if exp.periods[0].time >= 0.0 else 0.0
                    for exp in experiments
                ]
            )

            def preeq_condition_ids(exp: petabv2.Experiment) -> list[str]:
                for period in exp.sorted_periods:
                    if period.is_preequilibration:
                        return period.condition_ids
                return []

            reinit_condition_ids = [
                preeq_condition_ids(exp) for exp in experiments
            ]
        else:
            p_array = jnp.stack(
                [
                    jnp.stack(
                        [
                            self.load_model_parameters(
                                exp, is_preeq=False, period_index=i
                            )
                            for i in range(self._max_periods)
                        ]
                    )
                    for exp in experiments
                ]
            )

            def period_start_times(exp: petabv2.Experiment) -> jnp.ndarray:
                dyn_periods = self._dynamic_periods(exp)
                last_time = dyn_periods[-1].time if dyn_periods else 0.0
                return jnp.array(
                    [
                        dyn_periods[i].time
                        if i < len(dyn_periods)
                        else last_time
                        for i in range(self._max_periods)
                    ]
                )

            t_zeros = jnp.stack(
                [period_start_times(exp) for exp in experiments]
            )

            def period_condition_ids(
                exp: petabv2.Experiment, i: int
            ) -> list[str]:
                dyn_periods = self._dynamic_periods(exp)
                return dyn_periods[i].condition_ids if i < len(dyn_periods) else []

            reinit_condition_ids = [
                [period_condition_ids(exp, i) for i in range(self._max_periods)]
                for exp in experiments
            ]

        exp_ids = [exp.id for exp in experiments]
        all_exp_ids = [exp.id for exp in self._petab_problem.experiments]

        h_mask = jnp.stack(
            [
                jnp.ones(self.model.n_events)
                if (exp_id in exp_ids)
                else jnp.zeros(self.model.n_events)
                for exp_id in all_exp_ids
            ]
        )

        if self.parameters.size:
            if isinstance(self._petab_problem, petabv2.Problem):
                unscaled_parameters = jnp.stack(
                    [
                        self.parameters[ip]
                        for ip, p_id in enumerate(self.parameter_ids)
                    ]
                )
            else:
                unscaled_parameters = jnp.stack(
                    [
                        jax_unscale(
                            self.parameters[ip],
                            self._petab_problem.parameter_df.loc[
                                p_id, petabv2.C.PARAMETER_SCALE
                            ],
                        )
                        for ip, p_id in enumerate(self.parameter_ids)
                    ]
                )
        else:
            # No free parameters to gather from. `op_indices`/`np_indices`
            # may still contain the dummy index 0 for masked-out (i.e. not
            # actually free-parameter-referencing) entries, so keep a
            # single dummy slot to gather from -- it is never selected by
            # `jnp.where` since the corresponding mask is always False.
            unscaled_parameters = jnp.zeros((1,))

        def gather_unscaled(indices: np.ndarray) -> jnp.ndarray:
            fn = lambda ip: unscaled_parameters[ip]  # noqa: E731
            for _ in range(indices.ndim):
                fn = jax.vmap(fn)
            return fn(indices)

        # placeholder values from sundials code may be needed here
        if op_numeric is not None and op_numeric.size:
            op_array = jnp.where(
                op_mask,
                gather_unscaled(op_indices),
                op_numeric,
            )
        else:
            op_array = jnp.zeros((*self._ts_masks.shape, 0))

        if np_numeric is not None and np_numeric.size:
            np_array = jnp.where(
                np_mask,
                gather_unscaled(np_indices),
                np_numeric,
            )
        else:
            np_array = jnp.zeros((*self._ts_masks.shape, 0))

        if is_preeq:
            mask_reinit_array = jnp.stack(
                [
                    self.load_reinitialisation(cids, p)[0]
                    for cids, p in zip(reinit_condition_ids, p_array)
                ]
            )
            x_reinit_array = jnp.stack(
                [
                    self.load_reinitialisation(cids, p)[1]
                    for cids, p in zip(reinit_condition_ids, p_array)
                ]
            )
        else:
            mask_reinit_array = jnp.stack(
                [
                    jnp.stack(
                        [
                            self.load_reinitialisation(cids_i, p_i)[0]
                            for cids_i, p_i in zip(cids_per_period, p_per_period)
                        ]
                    )
                    for cids_per_period, p_per_period in zip(
                        reinit_condition_ids, p_array
                    )
                ]
            )
            x_reinit_array = jnp.stack(
                [
                    jnp.stack(
                        [
                            self.load_reinitialisation(cids_i, p_i)[1]
                            for cids_i, p_i in zip(cids_per_period, p_per_period)
                        ]
                    )
                    for cids_per_period, p_per_period in zip(
                        reinit_condition_ids, p_array
                    )
                ]
            )
        return (
            p_array,
            mask_reinit_array,
            x_reinit_array,
            op_array,
            np_array,
            h_mask,
            t_zeros,
        )

    @eqx.filter_vmap(
        in_axes={
            "max_steps": None,
            "self": None,
        },  # only list arguments here where eqx.is_array(0) is not the right thing
    )
    def run_simulation(
        self,
        p: jt.Float[jt.Array, "np"],  # noqa: F821, F722
        ts_dyn: np.ndarray,
        ts_posteq: np.ndarray,
        my: np.ndarray,
        iys: np.ndarray,
        iy_trafos: np.ndarray,
        ops: jt.Float[jt.Array, "nt *nop"],  # noqa: F821, F722
        nps: jt.Float[jt.Array, "nt *nnp"],  # noqa: F821, F722
        mask_reinit: jt.Bool[jt.Array, "nx"],  # noqa: F821, F722
        x_reinit: jt.Float[jt.Array, "nx"],  # noqa: F821, F722
        init_override: jt.Float[jt.Array, "nx"],  # noqa: F821, F722
        init_override_mask: jt.Bool[jt.Array, "nx"],  # noqa: F821, F722
        h_mask: jt.Bool[jt.Array, "nx"],  # noqa: F821, F722
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        root_finder: AbstractRootFinder,
        steady_state_event: Callable[
            ..., diffrax._custom_types.BoolScalarLike
        ],
        max_steps: jnp.int_,
        x_preeq: jt.Float[jt.Array, "*nx"] = jnp.array([]),  # noqa: F821, F722
        h_preeq: jt.Bool[jt.Array, "*ne"] = jnp.array([]),  # noqa: F821, F722
        ts_mask: np.ndarray = np.array([]),
        t_zeros: jnp.float_ = 0.0,
        experiment_index: jnp.int32 = jnp.int32(0),
        ret: ReturnValue = ReturnValue.llh,
    ) -> tuple[jnp.float_, dict]:
        """
        Run a simulation for a given simulation experiment.

        :param p:
            Parameters for the simulation experiment
        :param ts_dyn:
            (Padded) dynamic time points
        :param ts_posteq:
            (Padded) post-equilibrium time points
        :param my:
            (Padded) measurements
        :param iys:
            (Padded) observable indices
        :param iy_trafos:
            (Padded) observable transformations indices
        :param ops:
            (Padded) observable parameters
        :param nps:
            (Padded) noise parameters
        :param mask_reinit:
            Mask for states that need reinitialisation
        :param x_reinit:
            Reinitialisation values for states
        :param h_mask:
            Mask for the events that are part of the current experiment
        :param solver:
            ODE solver to use for simulation
        :param controller:
            Step size controller to use for simulation
        :param steady_state_event:
            Steady state event function to use for post-equilibration. Allows customisation of the steady state
            condition, see :func:`diffrax.steady_state_event` for details.
        :param max_steps:
            Maximum number of steps to take during simulation
        :param x_preeq:
            Pre-equilibration state. Can be empty if no pre-equilibration is available, in which case the states will
            be initialised to the model default values.
        :param h_preeq:
            Pre-equilibration event mask. Can be empty if no pre-equilibration is available
        :param ts_mask:
            padding mask, see :meth:`JAXModel.simulate_condition` for details.
        :param t_zeros:
            simulation start time for the current experiment.
        :param ret:
            which output to return. See :class:`ReturnValue` for available options.
        :return:
            Tuple of output value and simulation statistics
        """
        model = eqx.tree_at(
            lambda m: m._array_input_index, self.model, experiment_index
        )
        return model.simulate_condition(
            p=p,
            ts_dyn=jax.lax.stop_gradient(jnp.array(ts_dyn)),
            ts_posteq=jax.lax.stop_gradient(jnp.array(ts_posteq)),
            my=jax.lax.stop_gradient(jnp.array(my)),
            iys=jax.lax.stop_gradient(jnp.array(iys)),
            iy_trafos=jax.lax.stop_gradient(jnp.array(iy_trafos)),
            nps=nps,
            ops=ops,
            x_preeq=x_preeq,
            h_preeq=h_preeq,
            mask_reinit=jax.lax.stop_gradient(mask_reinit),
            x_reinit=x_reinit,
            init_override=init_override,
            init_override_mask=jax.lax.stop_gradient(init_override_mask),
            ts_mask=jax.lax.stop_gradient(jnp.array(ts_mask)),
            h_mask=jax.lax.stop_gradient(jnp.array(h_mask)),
            t_zero=t_zeros,
            solver=solver,
            controller=controller,
            root_finder=root_finder,
            max_steps=max_steps,
            steady_state_event=steady_state_event,
            adjoint=diffrax.RecursiveCheckpointAdjoint()
            if ret in (ReturnValue.llh, ReturnValue.chi2)
            else diffrax.DirectAdjoint(),
            ret=ret,
        )

    def run_simulations(
        self,
        experiments: list[petabv2.Experiment],
        preeq_array: jt.Float[jt.Array, "ncond *nx"],  # noqa: F821, F722
        h_preeqs: jt.Bool[jt.Array, "ncond *ne"],  # noqa: F821
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        root_finder: AbstractRootFinder,
        steady_state_event: Callable[
            ..., diffrax._custom_types.BoolScalarLike
        ],
        max_steps: jnp.int_,
        ret: ReturnValue = ReturnValue.llh,
    ):
        """
        Run simulations for a list of simulation experiments.

        :param experiments:
            Experiments to run simulations for.
        :param preeq_array:
            Matrix of pre-equilibrated states for the simulation conditions. Ordering must match the simulation
            conditions. If no pre-equilibration is available for a condition, the corresponding row must be empty.
        :param h_preeqs:
            Matrix of pre-equilibration event heaviside variables indicating whether an event condition is false or
            true after preequilibration.
        :param solver:
            ODE solver to use for simulation.
        :param controller:
            Step size controller to use for simulation.
        :param steady_state_event:
            Steady state event function to use for post-equilibration. Allows customisation of the steady state
            condition, see :func:`diffrax.steady_state_event` for details.
        :param max_steps:
            Maximum number of steps to take during simulation.
        :param ret:
            which output to return. See :class:`ReturnValue` for available options.
        :return:
            Output value and condition specific results and statistics. Results and statistics are returned as a dict
            with arrays with the leading dimension corresponding to the simulation conditions.
        """
        (
            p_array,
            mask_reinit_array,
            x_reinit_array,
            op_array,
            np_array,
            h_mask,
            t_zeros,
        ) = self._prepare_experiments(
            experiments,
            False,
            self._op_numeric,
            self._op_mask,
            self._op_indices,
            self._np_numeric,
            self._np_mask,
            self._np_indices,
        )

        init_override_mask = jnp.stack(
            [
                jnp.array(
                    [
                        p in set(self.model.parameter_ids)
                        for p in self.model.state_ids
                    ]
                )
                for _ in experiments
            ]
        )
        init_override = jnp.stack(
            [
                jnp.array(
                    [
                        self._eval_nn(
                            p, exp.sorted_periods[-1].condition_ids[0]
                        )  # TODO: Add mapping of p to eval_nn?
                        if p in set(self.model.parameter_ids)
                        else 1.0
                        for p in self.model.state_ids
                    ]
                )
                for exp in experiments
            ]
        )

        return self.run_simulation(
            p_array,
            self._ts_dyn,
            self._ts_posteq,
            self._my,
            self._iys,
            self._iy_trafos,
            op_array,
            np_array,
            mask_reinit_array,
            x_reinit_array,
            init_override,
            init_override_mask,
            h_mask,
            solver,
            controller,
            root_finder,
            steady_state_event,
            max_steps,
            preeq_array,
            h_preeqs,
            self._ts_masks,
            t_zeros,
            jnp.arange(len(experiments)),
            ret,
        )

    @eqx.filter_vmap(
        in_axes={
            "max_steps": None,
            "self": None,
        },  # only list arguments here where eqx.is_array(0) is not the right thing
    )
    def run_preequilibration(
        self,
        p: jt.Float[jt.Array, "np"],  # noqa: F821, F722
        mask_reinit: jt.Bool[jt.Array, "nx"],  # noqa: F821, F722
        x_reinit: jt.Float[jt.Array, "nx"],  # noqa: F821, F722
        h_mask: jt.Bool[jt.Array, "ne"],  # noqa: F821, F722
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        root_finder: AbstractRootFinder,
        steady_state_event: Callable[
            ..., diffrax._custom_types.BoolScalarLike
        ],
        max_steps: jnp.int_,
    ) -> tuple[jt.Float[jt.Array, "nx"], dict]:  # noqa: F821
        """
        Run a pre-equilibration simulation for a given simulation experiment.

        :param p:
            Parameters for the simulation experiment
        :param mask_reinit:
            Mask for states that need reinitialisation
        :param x_reinit:
            Reinitialisation values for states
        :param h_mask:
            Mask for the events that are part of the current experiment
        :param solver:
            ODE solver to use for simulation
        :param controller:
            Step size controller to use for simulation
        :param steady_state_event:
            Steady state event function to use for pre-equilibration. Allows customisation of the steady state
            condition, see :func:`diffrax.steady_state_event` for details.
        :param max_steps:
            Maximum number of steps to take during simulation
        :return:
            Pre-equilibration state
        """
        return self.model.preequilibrate_condition(
            p=p,
            mask_reinit=mask_reinit,
            x_reinit=x_reinit,
            h_mask=h_mask,
            solver=solver,
            controller=controller,
            root_finder=root_finder,
            max_steps=max_steps,
            steady_state_event=steady_state_event,
        )

    def run_preequilibrations(
        self,
        experiments: list[petabv2.Experiment],
        solver: diffrax.AbstractSolver,
        controller: diffrax.AbstractStepSizeController,
        root_finder: AbstractRootFinder,
        steady_state_event: Callable[
            ..., diffrax._custom_types.BoolScalarLike
        ],
        max_steps: jnp.int_,
    ):
        p_array, mask_reinit_array, x_reinit_array, _, _, h_mask, _ = (
            self._prepare_experiments(experiments, True, None, None)
        )
        return self.run_preequilibration(
            p_array,
            mask_reinit_array,
            x_reinit_array,
            h_mask,
            solver,
            controller,
            root_finder,
            steady_state_event,
            max_steps,
        )


def run_simulations(
    problem: JAXProblem,
    simulation_experiments: Iterable[str] | None = None,
    solver: diffrax.AbstractSolver = diffrax.Kvaerno5(),
    controller: diffrax.AbstractStepSizeController = diffrax.PIDController(
        **DEFAULT_CONTROLLER_SETTINGS
    ),
    root_finder: AbstractRootFinder = optimistix.Newton(
        **DEFAULT_ROOT_FINDER_SETTINGS
    ),
    steady_state_event: Callable[
        ..., diffrax._custom_types.BoolScalarLike
    ] = diffrax.steady_state_event(),
    max_steps: int = 2**13,
    ret: ReturnValue | str = ReturnValue.llh,
):
    """
    Run simulations for a problem.

    :param problem:
        Problem to run simulations for.
    :param simulation_experiments:
        Simulation experiments to run simulations for. This is an iterable of experiment ids.
        Default is to run simulations for all experiments.
    :param solver:
        ODE solver to use for simulation.
    :param controller:
        Step size controller to use for simulation.
    :param root_finder:
        Root finder to use for event detection.
    :param steady_state_event:
        Steady state event function to use for pre-/post-equilibration. Allows customisation of the steady state
        condition, see :func:`diffrax.steady_state_event` for details.
    :param max_steps:
        Maximum number of steps to take during simulation.
    :param ret:
        which output to return. See :class:`ReturnValue` for available options.
    :return:
        Overall output value and condition specific results and statistics.
    """
    if isinstance(problem._petab_problem, petabv1.Problem):
        raise TypeError(
            "run_simulations does not support PEtab v1 problems. Upgrade the problem to PEtab v2."
        )

    if isinstance(ret, str):
        ret = ReturnValue[ret]

    if simulation_experiments is None:
        experiments = problem._petab_problem.experiments
    else:
        experiments = [
            exp
            for exp in problem._petab_problem.experiments
            if exp.id in simulation_experiments
        ]

    preeq_condition_ids = _get_preequilibration_condition_ids(experiments)
    dynamic_conditions = [
        _dynamic_period_label(exp, period, preeq_condition_ids)
        for exp in experiments
        for period in exp.sorted_periods
        if not period.is_preequilibration
    ]
    dynamic_conditions = list(dict.fromkeys(dynamic_conditions))
    conditions = {
        "dynamic_conditions": dynamic_conditions,
    }

    if has_preeq := any(exp.has_preequilibration for exp in experiments):
        preeqs, preresults, h_preeqs = problem.run_preequilibrations(
            experiments,
            solver,
            controller,
            root_finder,
            steady_state_event,
            max_steps,
        )
        preeqs_array = preeqs
    else:
        preresults = {
            "stats_preeq": None,
        }
        preeqs_array = jnp.stack([jnp.array([]) for _ in experiments])
        h_preeqs = jnp.stack([jnp.array([]) for _ in experiments])

    output, results = problem.run_simulations(
        experiments,
        preeqs_array,
        h_preeqs,
        solver,
        controller,
        root_finder,
        steady_state_event,
        max_steps,
        ret,
    )

    if ret in (ReturnValue.llh, ReturnValue.chi2):
        if os.getenv("JAX_DEBUG") == "1":
            jax.debug.print(
                "ret: {}",
                ret,
            )
        output = jnp.sum(output)

    return output, results | preresults | conditions


def petab_simulate(
    problem: JAXProblem,
    solver: diffrax.AbstractSolver = diffrax.Kvaerno5(),
    controller: diffrax.AbstractStepSizeController = diffrax.PIDController(
        **DEFAULT_CONTROLLER_SETTINGS
    ),
    steady_state_event: Callable[
        ..., diffrax._custom_types.BoolScalarLike
    ] = diffrax.steady_state_event(),
    max_steps: int = 2**13,
):
    """
    Run simulations for a problem and return the results as a petab simulation dataframe.

    :param problem:
        Problem to run simulations for.
    :param solver:
        ODE solver to use for simulation.
    :param controller:
        Step size controller to use for simulation.
    :param max_steps:
        Maximum number of steps to take during simulation.
    :param steady_state_event:
        Steady state event function to use for pre-/post-equilibration. Allows customisation of the steady state
        condition, see :func:`diffrax.steady_state_event` for details.
    :return:
        petab simulation dataframe.
    """
    y, r = run_simulations(
        problem,
        solver=solver,
        controller=controller,
        steady_state_event=steady_state_event,
        max_steps=max_steps,
        ret=ReturnValue.y,
    )
    if isinstance(problem._petab_problem, petabv2.Problem):
        return _build_simulation_df_v2(problem, y, r["dynamic_conditions"])
    else:
        dfs = []
        for ic, sc in enumerate(r["dynamic_conditions"]):
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
                    petabv2.C.SIMULATION: y[ic, problem._ts_masks[ic, :]],
                    petabv2.C.TIME: t[problem._ts_masks[ic, :]],
                    petabv2.C.OBSERVABLE_ID: obs,
                    petabv2.C.CONDITION_ID: [sc] * len(t),
                },
                index=problem._petab_measurement_indices[ic, :],
            )
            if (
                petabv2.C.OBSERVABLE_PARAMETERS
                in problem._petab_problem.measurement_df
            ):
                df_sc[petabv2.C.OBSERVABLE_PARAMETERS] = (
                    problem._petab_problem.measurement_df.query(
                        f"{petabv2.C.CONDITION_ID} == '{sc}'"
                    )[petabv2.C.OBSERVABLE_PARAMETERS]
                )
            if (
                petabv2.C.NOISE_PARAMETERS
                in problem._petab_problem.measurement_df
            ):
                df_sc[petabv2.C.NOISE_PARAMETERS] = (
                    problem._petab_problem.measurement_df.query(
                        f"{petabv2.C.CONDITION_ID} == '{sc}'"
                    )[petabv2.C.NOISE_PARAMETERS]
                )
            if (
                petabv2.C.PREEQUILIBRATION_CONDITION_ID
                in problem._petab_problem.measurement_df
            ):
                df_sc[petabv2.C.PREEQUILIBRATION_CONDITION_ID] = (
                    problem._petab_problem.measurement_df.query(
                        f"{petabv2.C.CONDITION_ID} == '{sc}'"
                    )[petabv2.C.PREEQUILIBRATION_CONDITION_ID]
                )
            dfs.append(df_sc)
        return pd.concat(dfs).sort_index()


def add_default_experiment_names_to_v2_problem(petab_problem: petabv2.Problem):
    """Add default experiment names to PEtab v2 problem.

    Args:
        petab_problem: PEtab v2 problem to modify.
    """
    petab_problem.visualization_df = None

    if petab_problem.condition_df is None:
        default_condition = petabv2.core.Condition(
            id="__default__", changes=[], conditionId="__default__"
        )
        petab_problem.condition_tables[0].elements = [default_condition]

    if (
        petab_problem.experiment_df is None
        or petab_problem.experiment_df.empty
    ):
        # NOTE: `condition_df` is long-format (one row per condition
        #  *change*), so a condition with no changes (e.g. the just-created
        #  default condition) would not appear in it at all. Read condition
        #  ids from the `Condition` objects directly instead.
        condition_ids = [
            c.id
            for c in petab_problem.conditions
            if "preequilibration" not in c.id
        ]
        default_experiment = petabv2.core.Experiment(
            id="__default__",
            periods=[
                petabv2.core.ExperimentPeriod(
                    time=0.0, condition_ids=condition_ids
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

    experiment_df = experiment_df.drop(columns=[petabv2.C.TIME])
    return experiment_df


def _dynamic_period_label(
    experiment: petabv2.Experiment,
    period: petabv2.ExperimentPeriod,
    preeq_condition_ids: set[str],
) -> str:
    """Label identifying a non-pre-equilibration ``period`` for the purposes
    of building/looking up the ``dynamic_conditions`` list used by
    :func:`run_simulations`/:func:`_build_simulation_df_v2`.

    Each period contributes exactly one label -- even if it has several
    simultaneous condition ids (e.g. an indicator-condition encoding, where
    a period may carry both an experiment-indicator and a
    pre-equilibration-toggle condition id) -- since a period is always one
    simulation leg/row-group, regardless of how many condition ids describe
    it. The first non-pre-equilibration condition id is used; periods
    without any condition table changes at all (e.g. a period that only
    fixes a non-zero start time, with all parameters/states left at their
    default) have no condition id to key off, so synthesize one unique to
    this ``(experiment, period)`` pair.
    """
    cids = [
        cid for cid in period.condition_ids if cid not in preeq_condition_ids
    ]
    if cids:
        return cids[0]
    return f"__no_condition__{experiment.id}__{period.time}"


def _dynamic_condition_index_map(
    experiments: list[petabv2.Experiment],
) -> dict[str, tuple[int, int]]:
    """Map each non-pre-equilibration period's label (see
    :func:`_dynamic_period_label`) to the ``(experiment_index,
    period_index)`` position of the period it belongs to, matching the
    traversal/de-duplication order used to build ``dynamic_conditions`` in
    :func:`run_simulations`.

    :param experiments:
        Experiments to build the mapping for, in the same order as used to
        build :attr:`JAXProblem._ts_dyn` etc. (i.e.
        ``problem._petab_problem.experiments``).
    """
    preeq_ids = _get_preequilibration_condition_ids(experiments)
    positions: dict[str, tuple[int, int]] = {}
    for exp_idx, exp in enumerate(experiments):
        period_idx = 0
        for period in exp.sorted_periods:
            if period.is_preequilibration:
                continue
            label = _dynamic_period_label(exp, period, preeq_ids)
            positions.setdefault(label, (exp_idx, period_idx))
            period_idx += 1
    return positions


def _build_simulation_df_v2(problem, y, dyn_conditions):
    """Build petab simulation DataFrame of similation results from a PEtab v2 problem."""
    experiments = problem._petab_problem.experiments
    position_map = _dynamic_condition_index_map(experiments)
    nt_per_period = problem._ts_masks.shape[-1]

    dfs = []
    for sc in dyn_conditions:
        exp_idx, period_idx = position_map[sc]
        experiment_id = experiments[exp_idx].id

        if experiment_id == "__default__":
            experiment_id = jnp.nan

        mask = problem._ts_masks[exp_idx, period_idx, :]
        obs = [
            problem.model.observable_ids[io]
            for io in problem._iys[exp_idx, period_idx, mask]
        ]
        t = jnp.concatenate(
            (
                problem._ts_dyn[exp_idx, period_idx, :],
                problem._ts_posteq[exp_idx, period_idx, :],
            )
        )
        y_period = y[exp_idx].reshape(-1, nt_per_period)[period_idx, :]
        n_real = int(mask.sum())
        df_sc = pd.DataFrame(
            {
                petabv2.C.MODEL_ID: [float("nan")] * n_real,
                petabv2.C.OBSERVABLE_ID: obs,
                petabv2.C.EXPERIMENT_ID: [experiment_id] * n_real,
                petabv2.C.TIME: t[mask],
                petabv2.C.SIMULATION: y_period[mask],
            },
            index=problem._petab_measurement_indices[exp_idx, period_idx, mask],
        )
        measurement_df = problem._petab_problem.measurement_df
        # `experiment_id` may be `jnp.nan` (the "__default__" experiment
        # sentinel): a string `.query()` (`== 'nan'`) does not match actual
        # NaN values in the column, silently selecting zero rows and
        # leaving the assigned column all-NaN below. Select via a boolean
        # mask instead, which handles both the NaN and real-id cases.
        if isinstance(experiment_id, float):  # NaN sentinel
            exp_rows = measurement_df[
                measurement_df[petabv2.C.EXPERIMENT_ID].isna()
            ]
        else:
            exp_rows = measurement_df[
                measurement_df[petabv2.C.EXPERIMENT_ID] == experiment_id
            ]
        if petabv2.C.OBSERVABLE_PARAMETERS in measurement_df:
            df_sc[petabv2.C.OBSERVABLE_PARAMETERS] = exp_rows[
                petabv2.C.OBSERVABLE_PARAMETERS
            ]
        if petabv2.C.NOISE_PARAMETERS in measurement_df:
            df_sc[petabv2.C.NOISE_PARAMETERS] = exp_rows[
                petabv2.C.NOISE_PARAMETERS
            ]
        dfs.append(df_sc)
    return pd.concat(dfs).sort_index()


def _conditions_to_experiment_map(
    experiment_df: pd.DataFrame,
) -> dict[str, str]:
    return {
        row.conditionId: row.experimentId
        for row in experiment_df.itertuples()
    }


def _get_preequilibration_condition_ids(
    experiments: Iterable[petabv2.Experiment],
) -> set[str]:
    """Get the condition IDs used by pre-equilibration periods.

    Determined from :attr:`petabv2.ExperimentPeriod.is_preequilibration`
    rather than by pattern-matching condition IDs, since condition IDs are no
    longer guaranteed to follow the naming convention introduced by
    ``ExperimentsToSbmlConverter`` (which is not applied for the JAX
    backend).
    """
    return {
        cid
        for experiment in experiments
        for period in experiment.sorted_periods
        if period.is_preequilibration
        for cid in period.condition_ids
    }


def _parse_model_entity_id(
    model_entity_id: str, nn: dict
) -> list[tuple[str, str]]:
    """Parse a PEtab SciML model entity ID to find which NN layers/attributes to set.

    Handles ``net.parameters[layer]`` (all attributes of that layer),
    ``net.parameters[layer].weight`` (only that attribute), and
    ``net.parameters`` (all layers, fallback).

    :param model_entity_id:
        Model entity ID from the mapping table.
    :param nn:
        NN parameter dict for the net (``model_pars[net_id]``).
    :return:
        List of ``(layer_name, attribute)`` tuples to set.
    """
    to_set: list[tuple[str, str]] = []
    match = re.search(r"\[([^\]]+)\]", model_entity_id)
    if not match:
        return [
            (layer_name, attr)
            for layer_name, layer in nn.items()
            for attr in layer.keys()
        ]
    layer_name = match.group(1)
    if layer_name not in nn:
        return to_set
    layer = nn[layer_name]
    after_bracket = model_entity_id[match.end() :]
    if after_bracket.startswith("."):
        attribute = after_bracket[1:]
        if attribute in layer:
            to_set.append((layer_name, attribute))
        else:
            logger.warning(
                f"Attribute '{attribute}' not found in layer '{layer_name}' "
                f"while parsing model entity ID '{model_entity_id}'"
            )
    else:
        to_set.extend([(layer_name, attr) for attr in layer.keys()])
    return to_set


def _try_float(value):
    try:
        return jnp.asarray(float(value))
    except Exception as e:
        msg = str(e).lower()
        if isinstance(e, ValueError) and "could not convert" in msg:
            return value
        raise


def _resolve_petab_change_value(target_value) -> float | str:
    """
    Resolve a :class:`petabv2.Change.target_value` (a sympy expression, or
    already a plain number) to either a numeric literal or a single
    parameter id.

    Only numeric literals and references to a single other parameter are
    supported; compound symbolic expressions (e.g. ``"k1 + k2"``) are not.

    :param target_value:
        Value to resolve.
    :return:
        A ``float``, or a ``str`` naming a PEtab/model parameter id.
    """
    if getattr(target_value, "is_number", True):
        return float(target_value)
    if getattr(target_value, "is_Symbol", False):
        return str(target_value)
    raise NotImplementedError(
        "Condition table changes with compound symbolic expressions are "
        f"not supported, got {target_value!r}."
    )
