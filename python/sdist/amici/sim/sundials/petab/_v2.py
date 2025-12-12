"""PEtab v2 simulation."""

from __future__ import annotations

import logging
import numbers
from collections import Counter
from collections.abc import Sequence
from typing import Any

import numpy as np
import sympy as sp
from petab import v2 as v2
from petab.v2 import ExperimentPeriod

import amici
from amici.logging import get_logger
from amici.sim.sundials import SensitivityOrder
from amici.sim.sundials.petab.v1._simulations import (
    EDATAS,
    LLH,
    RDATAS,
    RES,
    S2LLH,
    SLLH,
    SRES,
)

logger = get_logger(__name__, log_level=logging.INFO)

__all__ = [
    "PetabSimulator",
    "ExperimentManager",
]


class ExperimentManager:
    # TODO: support for pscale?
    """
    Handles the creation of :class:`ExpData` objects for a given model and
    PEtab problem.

    The assumption is that we have a set of :class:`ExpData` objects,
    one for each PEtab experiment.
    Those are updated based on a set of global parameters (PEtab
    problem parameters, as opposed to model parameters for a single experiment
    period).
    """

    # TODO debug, remove
    _debug = False

    def __init__(
        self,
        model: amici.sim.sundials.Model,
        petab_problem: v2.Problem,
    ):
        """
        Initialize the `ExperimentManager`.

        :param model: The AMICI model to use.
        :param petab_problem: The PEtab problem to use.
            This is expected to be the output of
            :class:`petab.v2.ExperimentsToSbmlConverter` or an equivalent problem.
            This object must not be modified after the creation of this
            :class:`ExperimentManager` instance.
        """
        self._model: amici.sim.sundials.Model = model
        self._petab_problem: v2.Problem = petab_problem
        self._state_ids: tuple[str, ...] = tuple(self._model.get_state_ids())
        self._free_parameter_ids: tuple[str, ...] = tuple(
            self._model.get_free_parameter_ids()
        )
        self._fixed_parameter_ids: tuple[str, ...] = tuple(
            self._model.get_fixed_parameter_ids()
        )
        # maps parameter IDs to parameter indices in the model
        self._pid_to_idx: dict[str, int] = {
            id_: i for i, id_ in enumerate(self._free_parameter_ids)
        }
        self._fixed_pid_to_idx: dict[str, int] = {
            id_: i for i, id_ in enumerate(self._fixed_parameter_ids)
        }
        # maps PEtab observable IDs to petab Observable instances
        self._petab_id_to_obs: dict[str, v2.Observable] = {
            obs.id: obs for obs in self._petab_problem.observables
        }

        # create a new model instance from the model module from which
        #  we can get the default parameters
        model0 = model.module.get_model()
        self._original_p = np.array(model0.get_free_parameters())
        self._original_k = np.array(model0.get_fixed_parameters())

    def create_edatas(self) -> list[amici.sim.sundials.ExpData]:
        """Create ExpData objects for all experiments."""
        return [
            self.create_edata(experiment)
            for experiment in self._petab_problem.experiments
        ]

    def create_edata(
        self,
        experiment: v2.core.Experiment | str | None,
        problem_parameters: dict[str, float] | None = None,
    ) -> amici.sim.sundials.ExpData:
        """Create an ExpData object for a single experiment.

        Sets timepoints, measurements, initial conditions, ... based on the
        given experiment and the nominal parameters of the PEtab problem.

        :param experiment:
            The experiment or experiment ID to create the `ExpData` for.
        :param problem_parameters:
            Optional dictionary of problem parameters to apply to the
            `ExpData`. If `None`, the nominal parameters of the PEtab problem
            are used.
        :return:
            The created `ExpData` object for the given experiment.
        """
        from amici.sim.sundials import ExpData

        if isinstance(experiment, str):
            experiment = self._petab_problem[experiment]

        if len(experiment.periods) > 2:
            raise AssertionError(
                f"Expected <= 2 periods, got {len(experiment.periods)} "
                f"for experiment {experiment.id}."
            )

        edata = ExpData(self._model)
        edata.id = experiment.id

        self._set_constants(edata, experiment)
        self._set_timepoints_and_measurements(edata, experiment)

        if self._debug:
            logger.debug(
                f"Created ExpData id={edata.id}, "
                f"k_preeq={edata.fixed_parameters_pre_equilibration}, "
                f"k={edata.fixed_parameters}"
            )

        if problem_parameters is None:
            problem_parameters = self._petab_problem.get_x_nominal_dict()
        self.apply_parameters(edata, problem_parameters=problem_parameters)

        return edata

    def _set_constants(
        self, edata: amici.sim.sundials.ExpData, experiment: v2.core.Experiment
    ) -> None:
        """
        Set constant parameters for the given experiment.

        :param edata:
            The ExpData instance to set the constants for.
        :param experiment:
            The PEtab experiment to set the constants for.
        """
        # After converting experiments to events, all remaining
        # condition parameters are constants.
        if not experiment.periods:
            # No periods, no changes to apply.
            # Use the original fixed parameters that are encoded in the model.
            return

        def get_k(period: ExperimentPeriod):
            """Get the fixed parameters for the period."""
            changes = self._petab_problem.get_changes_for_period(period)
            fixed_pars_vals = self._original_k.copy()
            for change in changes:
                pid = self._fixed_pid_to_idx[change.target_id]
                # those are only indicator variables that are always number
                #  literals
                fixed_pars_vals[pid] = change.target_value
            return fixed_pars_vals

        if experiment.sorted_periods[0].time == -np.inf:
            # pre-equilibration period
            edata.fixed_parameters_pre_equilibration = get_k(
                experiment.sorted_periods[0]
            )
            # In PEtab, pre-equilibration always starts at t=0, since SBML
            # does not support specifying a different start time (yet).
            edata.t_start_preeq = 0

        if len(experiment.periods) >= int(1 + experiment.has_preequilibration):
            # simulation period
            main_period = experiment.sorted_periods[
                int(experiment.has_preequilibration)
            ]
            edata.fixed_parameters = get_k(main_period)
            edata.t_start = main_period.time

    def _set_timepoints_and_measurements(
        self, edata: amici.sim.sundials.ExpData, experiment: v2.core.Experiment
    ) -> None:
        """
        Set timepoints and measurements for the given experiment.

        :param edata:
            The `ExpData` instance to update.
        :param experiment:
            The PEtab experiment to set the timepoints and measurements for.
        """
        # Get the required time points: this is the superset of timepoints
        #  of the measurements of all observables, including the different
        #  replicates
        measurements = self._petab_problem.get_measurements_for_experiment(
            experiment
        )
        t_counters = {o.id: Counter() for o in self._petab_problem.observables}
        unique_t = set()
        for m in measurements:
            t_counters[m.observable_id].update([m.time])
            unique_t.add(m.time)

        max_counter = Counter()
        for t in unique_t:
            for counter in t_counters.values():
                max_counter[t] = max(max_counter[t], counter[t])

        timepoints_w_reps = sorted(max_counter.elements())

        edata.set_timepoints(timepoints_w_reps)

        # measurements and sigmas
        y, sigma_y = self._get_measurements_and_sigmas(
            measurements=measurements,
            timepoints_w_reps=timepoints_w_reps,
            observable_ids=self._model.get_observable_ids(),
        )
        edata.set_measurements(y.flatten())
        edata.set_noise_scales(sigma_y.flatten())

    def _get_measurements_and_sigmas(
        self,
        measurements: list[v2.Measurement],
        timepoints_w_reps: Sequence[numbers.Number],
        observable_ids: Sequence[str],
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Get measurements and sigmas

        Generate arrays with measured values and sigmas in AMICI format from a
        list of PEtab measurements.

        :param measurements:
            Subset of PEtab measurement table for one experiment
        :param timepoints_w_reps:
            Timepoints for which there exist measurements, including
            replicates.
        :param observable_ids:
            List of observable IDs for mapping IDs to indices.
        :return:
            arrays for measurement and sigmas
        """
        # prepare measurement matrix
        y = np.full(
            shape=(len(timepoints_w_reps), len(observable_ids)),
            fill_value=np.nan,
        )
        # prepare sigma matrix
        sigma_y = y.copy()

        time_to_meas = {}
        for m in measurements:
            time_to_meas.setdefault(m.time, []).append(m)

        for time in sorted(time_to_meas):
            # subselect for time
            time_ix_0 = timepoints_w_reps.index(time)

            # remember used time indices for each observable
            time_ix_for_obs_ix = {}

            # iterate over measurements
            for m in time_to_meas[time]:
                # extract observable index
                observable_ix = observable_ids.index(m.observable_id)

                # update time index for observable
                if observable_ix in time_ix_for_obs_ix:
                    time_ix_for_obs_ix[observable_ix] += 1
                else:
                    time_ix_for_obs_ix[observable_ix] = time_ix_0

                # fill observable and possibly noise parameter
                y[time_ix_for_obs_ix[observable_ix], observable_ix] = (
                    m.measurement
                )
                if (
                    len(m.noise_parameters) == 1
                    and m.noise_parameters[0].is_Number
                ):
                    sigma_y[
                        time_ix_for_obs_ix[observable_ix], observable_ix
                    ] = float(m.noise_parameters[0])
        return y, sigma_y

    def apply_parameters(
        self,
        edata: amici.sim.sundials.ExpData,
        problem_parameters: dict[str, float],
    ) -> None:
        """Apply problem parameters.

        Update the parameter-dependent values of the given `ExpData` instance
        according to the provided problem parameters (i.e., values of the
        parameters in the PEtab parameter table).

        This assumes that:

        * the `ExpData` instance was created by :meth:`create_edata`,
        * no other changes except for calls to :meth:`apply_parameters`
          were made,
        * and the PEtab problem was not modified since the creation of this
          :class:`ExperimentManager` instance.

        :param edata: The :class:`ExpData` instance to be updated.
            In case of errors, the state of `edata` is undefined.
        :param problem_parameters: Problem parameters to be applied.
        """
        # TODO: support ndarray in addition to dict?

        # check parameter IDs
        if set(problem_parameters) != set(self._petab_problem.x_ids):
            missing = set(self._petab_problem.x_ids) - set(problem_parameters)
            extra = set(problem_parameters) - set(self._petab_problem.x_ids)
            msg_parts = []
            if missing:
                # TODO: support a subset of parameters?
                #  if so, update only those parameters and leave the rest as is
                #  or use nominal values for missing ones?
                msg_parts.append(f"missing parameters: {missing}")
            if extra:
                msg_parts.append(f"unknown parameters: {extra}")
            raise ValueError(
                "Provided problem parameters do not match "
                "PEtab problem parameters: " + "; ".join(msg_parts)
            )

        # get the original parameter values
        #  (model parameters set during model creation)
        par_vals = np.array(self._original_p)
        pid_to_idx = self._pid_to_idx
        experiment_id = edata.id
        experiment = self._petab_problem[experiment_id]

        # plist -- estimated parameters + those mapped via placeholders
        # TODO sufficient to set them during creation of edata
        #   or allow dynamic fixing of parameters?
        plist = []
        placeholder_mappings = self._get_placeholder_mapping(experiment)
        estimated_par_ids = self._petab_problem.x_free_ids
        for model_par_idx, model_par_id in enumerate(
            self._model.get_free_parameter_ids()
        ):
            if model_par_id in estimated_par_ids or (
                (maps_to := placeholder_mappings.get(model_par_id)) is not None
                and maps_to in estimated_par_ids
            ):
                plist.append(model_par_idx)
        edata.plist = plist

        # Update fixed parameters in case they are affected by problem
        #  parameters (i.e., parameter table parameters)
        fixed_par_vals = np.asarray(edata.fixed_parameters)
        for p_id, p_val in problem_parameters.items():
            if (p_idx := self._fixed_pid_to_idx.get(p_id)) is not None:
                fixed_par_vals[p_idx] = p_val
        edata.fixed_parameters = fixed_par_vals

        if edata.fixed_parameters_pre_equilibration:
            fixed_par_vals = np.array(edata.fixed_parameters_pre_equilibration)
            for p_id, p_val in problem_parameters.items():
                if (p_idx := self._fixed_pid_to_idx.get(p_id)) is not None:
                    fixed_par_vals[p_idx] = p_val
            edata.fixed_parameters_pre_equilibration = fixed_par_vals

        # Apply problem parameter values to identical model parameters.
        #  Any other parameter mapping, except for output parameter
        #  placeholders, is handled by events.
        for k, v in problem_parameters.items():
            if (idx := pid_to_idx.get(k)) is not None:
                par_vals[idx] = v

        # Handle measurement-specific mappings to placeholders
        measurements = self._petab_problem.get_measurements_for_experiment(
            experiment
        )

        def apply_override(placeholder: str, override: sp.Basic):
            """Apply a single placeholder override."""
            if (idx := pid_to_idx.get(placeholder)) is not None:
                if override.is_Number:
                    par_vals[idx] = float(override)
                elif override.is_Symbol:
                    par_vals[idx] = problem_parameters[str(override)]
                else:
                    raise AssertionError(
                        f"Unexpected override type: {override} for {placeholder} in experiment {experiment_id}"
                    )
            else:
                raise NotImplementedError(
                    f"Cannot handle override `{placeholder}' => '{override}'"
                )

        # tracks encountered placeholders and their overrides
        # (across all observables -- placeholders IDs are globally unique)
        #  and check that all periods use the same overrides
        #  (except for numeric sigmas)
        # TODO: this can be simplified. we only need to process overrides
        #  that are parameters. the rest was handled during create_edata
        overrides = {}
        for m in measurements:
            obs = self._petab_id_to_obs[m.observable_id]
            if obs.observable_placeholders:
                for placeholder, override in zip(
                    obs.observable_placeholders,
                    m.observable_parameters,
                    strict=True,
                ):
                    placeholder = str(placeholder)
                    if (
                        prev_override := overrides.get(placeholder)
                    ) is not None and prev_override != override:
                        raise NotImplementedError(
                            "Timepoint-specific observable placeholder "
                            "overrides are not supported"
                        )
                    apply_override(placeholder, override)
            if obs.noise_placeholders:
                for placeholder, override in zip(
                    obs.noise_placeholders, m.noise_parameters, strict=True
                ):
                    placeholder = str(placeholder)
                    if (
                        prev_override := overrides.get(placeholder)
                    ) is not None and prev_override != override:
                        # TODO: this might have been handled
                        #  via .sigmay if numeric
                        raise NotImplementedError(
                            "Timepoint-specific observable placeholder "
                            "overrides are not supported"
                        )
                    apply_override(placeholder, override)

        # TODO: set all unused placeholders to NaN to make it easier to spot problems?
        edata.free_parameters = par_vals

        if self._debug:
            logger.debug("ExperimentManager.apply_parameters:")
            logger.debug(
                f"Parameters: "
                f"{dict(zip(self._free_parameter_ids, map(float, par_vals)))}"
            )

    @property
    def petab_problem(self) -> v2.Problem:
        """The PEtab problem used by this ExperimentManager.

        This must not be modified.
        """
        return self._petab_problem

    @property
    def model(self) -> amici.sim.sundials.Model:
        """The AMICI model used by this ExperimentManager."""
        return self._model

    def _get_placeholder_mapping(
        self, experiment: v2.Experiment
    ) -> dict[str, str]:
        """Get the mapping from model parameter IDs (= PEtab placeholder ID)
        to problem parameter IDs for placeholders in the given experiment.

        Because AMICI does not support timepoint-specific overrides,
        this mapping is unique.

        :param experiment: The experiment to get the mapping for.
        :return: The mapping from model parameter IDs to problem parameter IDs.
        """
        mapping = {}
        for measurement in self._petab_problem.get_measurements_for_experiment(
            experiment
        ):
            observable = self._petab_problem[measurement.observable_id]
            for placeholder, override in zip(
                observable.observable_placeholders,
                measurement.observable_parameters,
                strict=True,
            ):
                # we don't care about numeric overrides here
                if isinstance(override, sp.Symbol):
                    mapping[str(placeholder)] = str(override)

            for placeholder, override in zip(
                observable.noise_placeholders,
                measurement.noise_parameters,
                strict=True,
            ):
                # we don't care about numeric overrides here
                if isinstance(override, sp.Symbol):
                    mapping[str(placeholder)] = str(override)
        return mapping


class PetabSimulator:
    """
    Simulator for PEtab2 problems.

    This class is used to simulate all experiments of a given PEtab problem
    using a given AMICI model and solver, and to aggregate the results.
    """

    def __init__(
        self,
        em: ExperimentManager,
        *,
        solver: amici.sim.sundials.Solver | None = None,
        num_threads: int = 1,
        # TODO: allow selecting specific experiments?
        # TODO: store_edatas: bool
    ):
        """
        Initialize the simulator.

        :param em:
            The :class:`ExperimentManager` to generate the :class:`amici.ExpData`
            objects.
        :param solver: The AMICI solver to use for the simulations.
            If not provided, a new solver with default settings will be used.
        :param num_threads:
            The number of threads to use for parallel simulation of experiments.
            Only relevant if multiple experiments are present in the PEtab problem,
            and if AMICI was compiled with OpenMP support.
        """
        self._petab_problem: v2.Problem = em.petab_problem
        self._model = em.model
        self._solver = (
            solver if solver is not None else self._model.create_solver()
        )
        self._exp_man: ExperimentManager = em
        self.num_threads = num_threads

    @property
    def model(self) -> amici.sim.sundials.Model:
        """The AMICI model used by this simulator."""
        return self._model

    @property
    def solver(self) -> amici.sim.sundials.Solver:
        """The AMICI solver used by this simulator."""
        return self._solver

    @property
    def exp_man(self) -> ExperimentManager:
        """The ExperimentManager used by this simulator."""
        return self._exp_man

    def simulate(
        self, problem_parameters: dict[str, float] = None
    ) -> dict[str, Any]:
        # TODO params: dict|np.ndarray|None?
        """Simulate all experiments of the given PEtab problem.

        :return:
            Dictionary of

            * the summed cost function value (``LLH``),
            * list of :class:`amici.sim.sundials.ReturnData` (``RDATAS``)
              for each experiment,
            * list of :class:`amici.sim.sundials.ExpData` (``EDATAS``)
              for each experiment

           Note that the returned :class:`amici.amiciamici.sim.sundials.ExpData` instances
           may be changed by subsequent calls to this function.
           Create a copy if needed.
        """
        if problem_parameters is None:
            # use default parameters, i.e., nominal values for all parameters
            # TODO: Nominal parameters, or previously used parameters?
            problem_parameters = {}

        # use nominal values for all unspecified parameters
        problem_parameters_default = self._petab_problem.get_x_nominal_dict()
        problem_parameters = problem_parameters_default | problem_parameters

        # TODO cache edatas
        edatas = self._exp_man.create_edatas()

        for edata in edatas:
            self._exp_man.apply_parameters(
                edata=edata, problem_parameters=problem_parameters
            )

        rdatas = amici.sim.sundials.run_simulations(
            self._model, self._solver, edatas, num_threads=self.num_threads
        )

        return {
            EDATAS: edatas,
            RDATAS: rdatas,
            LLH: sum(rdata.llh for rdata in rdatas),
            SLLH: self._aggregate_sllh(rdatas),
            S2LLH: self._aggregate_s2llh(rdatas, use_fim=True),
            RES: np.hstack([rdata.res for rdata in rdatas]),
            # TODO: implement residual sensitivity aggregation
            SRES: None,
        }

    def _aggregate_sllh(
        self, rdatas: Sequence[amici.sim.sundials.ReturnDataView]
    ) -> dict[str, float] | None:
        """Aggregate the sensitivities of the log-likelihoods.

        :param rdatas:
            The ReturnData objects to aggregate the sensitivities from.
        :return:
            The aggregated sensitivities (parameter ID -> sensitivity value).
        """
        if self._solver.get_sensitivity_order() < SensitivityOrder.first:
            return None

        sllh_total: dict[str, float] = {}

        # Check for issues in all condition simulation results.
        for rdata in rdatas:
            # Condition failed during simulation.
            if rdata.status != amici.sim.sundials.AMICI_SUCCESS:
                return None
            # Condition simulation result does not provide SLLH.
            if rdata.sllh is None:
                raise ValueError(
                    f"The sensitivities of the likelihood for experiment "
                    f"{rdata.id} were not computed."
                )

        free_parameter_ids = self._model.get_free_parameter_ids()

        # still needs parameter mapping for placeholders
        for rdata in rdatas:
            experiment = self._petab_problem[rdata.id]
            placeholder_mappings = self._exp_man._get_placeholder_mapping(
                experiment
            )
            for model_par_idx, sllh in zip(
                rdata.plist, rdata.sllh, strict=True
            ):
                model_par_id = problem_par_id = free_parameter_ids[
                    model_par_idx
                ]
                if maps_to := placeholder_mappings.get(model_par_id):
                    problem_par_id = maps_to

                sllh_total[problem_par_id] = (
                    sllh_total.get(problem_par_id, 0.0) + sllh
                )
        return sllh_total

    def _aggregate_s2llh(
        self,
        rdatas: Sequence[amici.sim.sundials.ReturnDataView],
        use_fim: bool = True,
    ) -> np.ndarray | None:
        """Aggregate the Hessians from individual experiments.

        Compute the total second-order sensitivities of the log-likelihoods
        w.r.t. estimated PEtab problem parameters.

        :param rdatas:
            The ReturnData objects to aggregate the sensitivities from.
        :param use_fim:
            Whether to use the Fisher Information Matrix (FIM) to compute
            the 2nd order sensitivities. Only ``True`` is currently supported.
        :return:
            The aggregated 2nd order sensitivities as a 2D numpy array
            in the order of the estimated PEtab problem parameters
            (``Problem.x_free_ids``),
            or `None` if sensitivities were not computed.
        """
        # TODO: add tests
        if self._solver.get_sensitivity_order() < SensitivityOrder.first:
            return None

        if not use_fim:
            raise NotImplementedError(
                "Computation of 2nd order sensitivities without FIM is not "
                "implemented yet."
            )

        # Check for issues in all condition simulation results.
        for rdata in rdatas:
            # Condition failed during simulation.
            if rdata.status != amici.sim.sundials.AMICI_SUCCESS:
                return None
            # Condition simulation result does not provide FIM.
            if rdata.FIM is None:
                raise ValueError(
                    f"The FIM was not computed for experiment {rdata.id!r}."
                )

        # Model parameter index to problem parameter index map for estimated
        #  parameters except placeholders.
        #  This is the same for all experiments.
        global_ix_map: dict[int, int] = {
            model_ix: self._petab_problem.x_free_ids.index(model_pid)
            for model_ix, model_pid in enumerate(
                self._model.get_free_parameter_ids()
            )
            if model_pid in self._petab_problem.x_free_ids
        }
        s2llh_total = np.zeros(
            shape=(
                self._petab_problem.n_estimated,
                self._petab_problem.n_estimated,
            ),
            dtype=float,
        )

        for rdata in rdatas:
            ix_map = global_ix_map.copy()
            # still needs experiment-specific parameter mapping for
            # placeholders
            experiment = self._petab_problem[rdata.id]
            placeholder_mappings = self._exp_man._get_placeholder_mapping(
                experiment
            )
            for model_pid, problem_pid in placeholder_mappings.items():
                try:
                    ix_map[
                        self.model.get_free_parameter_ids().index(model_pid)
                    ] = self._petab_problem.x_free_ids.index(problem_pid)
                except ValueError:
                    # mapped-to parameter is not estimated
                    pass

            # translate model parameter index to plist index
            ix_map: dict[int, int] = {
                tuple(rdata.plist).index(model_par_ix): problem_par_ix
                for model_par_ix, problem_par_ix in ix_map.items()
            }
            if use_fim:
                model_s2llh = rdata.FIM
            else:
                raise NotImplementedError()

            model_par_slice = np.fromiter(ix_map.keys(), dtype=int)
            problem_par_slice = np.fromiter(ix_map.values(), dtype=int)

            # handle possible non-unique indices in problem_par_slice
            # (i.e. multiple model parameters mapping to the same problem
            # parameter)
            problem_par_slice_unique, unique_index = np.unique(
                problem_par_slice, return_index=True
            )
            # handle unique mappings
            s2llh_total[
                np.ix_(problem_par_slice_unique, problem_par_slice_unique)
            ] += model_s2llh[
                np.ix_(
                    model_par_slice[unique_index],
                    model_par_slice[unique_index],
                )
            ]
            # handle non-unique mappings if any
            if problem_par_slice_unique.size < problem_par_slice.size:
                # index in the mapping arrays of non-unique entries
                non_unique_indices = [
                    idx
                    for idx in range(len(problem_par_slice))
                    if idx not in unique_index
                ]
                for idx in non_unique_indices:
                    s2llh_total[
                        problem_par_slice[idx], problem_par_slice_unique
                    ] += model_s2llh[
                        model_par_slice[idx], model_par_slice[unique_index]
                    ]
                    s2llh_total[
                        problem_par_slice_unique, problem_par_slice[idx]
                    ] += model_s2llh[
                        model_par_slice[unique_index], model_par_slice[idx]
                    ]
                    for jdx in non_unique_indices:
                        s2llh_total[
                            problem_par_slice[idx], problem_par_slice[jdx]
                        ] += model_s2llh[
                            model_par_slice[idx], model_par_slice[jdx]
                        ]

        return s2llh_total
