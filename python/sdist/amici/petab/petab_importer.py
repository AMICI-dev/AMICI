"""PEtab v2 handling.

Functionality for importing and simulating
`PEtab v2 <https://petab.readthedocs.io/en/latest/v2/documentation_data_format.html#>`__
problems.
"""

from __future__ import annotations

import copy
import logging
import numbers
from collections import Counter
from collections.abc import Sequence
from pathlib import Path
from pprint import pprint
from typing import Any

import numpy as np
import pandas as pd
import sympy as sp
from petab import v1 as v1
from petab import v2 as v2
from petab.v2 import ExperimentPeriod, Observable
from petab.v2.converters import ExperimentsToSbmlConverter
from petab.v2.models import MODEL_TYPE_SBML

import amici
from amici import MeasurementChannel, SensitivityOrder

from ..de_model import DEModel
from ..logging import get_logger
from .sbml_import import _add_global_parameter

__all__ = [
    "PetabImporter",
    "ExperimentManager",
    "PetabSimulator",
    "rdatas_to_measurement_df",
    "flatten_timepoint_specific_output_overrides",
    "unflatten_simulation_df",
    "has_timepoint_specific_overrides",
]
logger = get_logger(__name__, log_level=logging.DEBUG)


#: Default experiment ID to be used for measurements without an experiment ID.
_DEFAULT_EXPERIMENT_ID = "__default__"

#: PEtab measurement table columns to consider for detecting timepoint-specific
#:  parameter overrides
_POSSIBLE_GROUPVARS_FLATTENED_PROBLEM = [
    v2.C.MODEL_ID,
    v2.C.EXPERIMENT_ID,
    v2.C.OBSERVABLE_ID,
    v2.C.OBSERVABLE_PARAMETERS,
    v2.C.NOISE_PARAMETERS,
]
# TODO: how to handle SBML vs PySB, jax vs sundials?
#  -> separate importers or subclasses?
# TODO: How to handle multi-model-problems?
# TODO: test simulation of up-converted benchmark problems


class PetabImporter:
    """
    Importer for PEtab problems.

    This class is used to create an AMICI model from a PEtab problem.

    The underlying SBML model will be modified to encode the experiments
    defined in the PEtab problem as events.

    Be careful when using the imported model for anything other than the
    PEtab-encoded experiments.

    All PEtab experiments will be encoded in the model, independent of
    whether they have measurements. This is to make it easier to simulate
    the respective experiments with the resulting AMICI model.
    This may make the resulting model more bloated. If this is not desired,
    the problem should be simplified before import.

    :param petab_problem:
        The PEtab problem to import. The problem must not be changed after
        construction of the importer.
    """

    def __init__(
        self,
        petab_problem: v2.Problem | v1.Problem,
        compile_: bool = None,
        validate: bool = True,
        # TODO: override the PEtab model ID vs selecting one of multiple models
        model_id: str = None,
        outdir: str | Path = None,
        jax: bool = False,
    ):
        """
        Create a new PetabImporter instance.

        :param petab_problem: The PEtab problem to import.
        :param compile_: Whether to compile the model extension after import.
        :param validate: Whether to validate the PEtab problem before import.
        :param model_id: TODO
        :param outdir:
            The output directory where the model files are written to.
        :param jax: TODO
        """
        self.petab_problem: v2.Problem = self._upgrade_if_needed(petab_problem)

        if len(self.petab_problem.models) > 1:
            raise NotImplementedError(
                "PEtab v2 importer currently only supports single-model "
                "problems."
            )

        if self.petab_problem.model.type_id != MODEL_TYPE_SBML:
            raise NotImplementedError(
                "PEtab v2 importer currently only supports SBML models. "
                f"Got {self.petab_problem.model.type_id}."
            )
        if jax:
            raise NotImplementedError(
                "PEtab v2 importer currently does not support JAX. "
            )

        pprint(self.petab_problem.model_dump())
        print(self.petab_problem.model.to_antimony())

        self._check_support(self.petab_problem)

        self._compile = compile_
        self._sym_model: DEModel | None = None
        self._model_id = model_id
        self._outdir: Path | None = (
            None if outdir is None else Path(outdir).absolute()
        )
        self._jax = jax
        self._verbose = logging.DEBUG

        if validate:
            logger.info("Validating PEtab problem ...")
            validation_result = petab_problem.validate()
            if validation_result:
                validation_result.log()

            if validation_result.has_errors():
                raise ValueError(
                    "PEtab problem is not valid, see log messages for details."
                )

        # ensure each measurement has an experimentId
        _set_default_experiment(self.petab_problem)

        # Convert petab experiments to events, because so far,
        #  AMICI only supports preequilibration/presimulation/simulation, but
        #  no arbitrary list of periods.
        self._exp_event_conv = ExperimentsToSbmlConverter(self.petab_problem)
        # This will always create a copy of the problem.
        self.petab_problem = self._exp_event_conv.convert()
        for experiment in self.petab_problem.experiments:
            if len(experiment.periods) > 2:
                raise NotImplementedError(
                    "AMICI currently does not support more than two periods."
                )

        # TODO remove dbg
        pprint(self.petab_problem.model_dump())
        print(self.petab_problem.model.to_antimony())

    def _upgrade_if_needed(
        self, problem: v1.Problem | v2.Problem
    ) -> v2.Problem:
        """Upgrade the problem to petab v2 if necessary."""
        if isinstance(problem, v2.Problem):
            return problem

        # TODO: So far, PEtab can only upgrade file-based problems,
        #  not petab.v1.Problem objects.
        raise NotImplementedError("Only petab.v2.Problem is supported.")

    @classmethod
    def _check_support(cls, petab_problem: v2.Problem):
        """Check if the PEtab problem requires unsupported features."""

        # check support for mapping tables
        relevant_mappings = [
            m
            for m in petab_problem.mappings
            # we can ignore annotation-only entries
            if m.model_id is not None
            # we can ignore identity mappings
            and m.petab_id != m.model_id
        ]
        if relevant_mappings:
            # It's partially supported. Remove at your own risk...
            raise NotImplementedError(
                "PEtab v2.0.0 mapping tables are not yet supported."
            )

    @property
    def model_id(self) -> str:
        """The model ID."""
        if self._model_id is None:
            self._model_id = self.petab_problem.model.model_id

        return self._model_id

    @property
    def outdir(self) -> Path:
        """The output directory where the model files are written to."""
        if self._outdir is None:
            self._outdir = Path(
                f"{self.model_id}-amici{amici.__version__}"
            ).absolute()
        return self._outdir

    def _do_import(self, non_estimated_parameters_as_constants: bool = True):
        """Import the model.

        Generate the symbolic model according to the given PEtab problem and
        generate the corresponding Python module.

        1. Encode all PEtab experiments as events in the SBML model.
           This leaves only (maybe) a pre-equilibration and a single
           simulation period.
        2. Add the observable parameters to the SBML model.

        :param non_estimated_parameters_as_constants:
            Whether parameters marked as non-estimated in PEtab should be
            considered constant in AMICI. Setting this to ``True`` will reduce
            model size and simulation times. If sensitivities with respect to
            those parameters are required, this should be set to ``False``.
        """
        # TODO split into DEModel creation, code generation and compilation
        #   allow retrieving DEModel without compilation

        from petab.v2.models.sbml_model import SbmlModel

        if not isinstance(self.petab_problem.model, SbmlModel):
            raise ValueError("The PEtab problem must contain an SBML model.")

        # set_log_level(logger, verbose)

        logger.info("Importing model ...")

        if not self.petab_problem.observables:
            raise NotImplementedError(
                "PEtab import without observables table "
                "is currently not supported."
            )

        if self.model_id is None:
            raise ValueError(
                "No `model_id` was provided and no model "
                "ID was specified in the SBML model."
            )

        logger.info(
            f"Model ID is '{self.model_id}'.\n"
            f"Writing model code to '{self.outdir}'."
        )

        observation_model = self._get_observation_model()

        logger.info(f"#Observables: {len(observation_model)}")
        logger.debug(f"Observables: {observation_model}")

        output_parameter_defaults = {}
        self._workaround_observable_parameters(
            output_parameter_defaults=output_parameter_defaults,
        )

        sbml_model = self.petab_problem.model.sbml_model

        # All indicator variables, i.e., all remaining targets after
        #  experiments-to-event in the PEtab problem must be converted
        #  to fixed parameters
        fixed_parameters = {
            change.target_id
            for experiment in self.petab_problem.experiments
            for period in experiment.periods
            for condition_id in period.condition_ids
            for change in self.petab_problem[condition_id].changes
        }

        # show_model_info(sbml_model)
        # TODO spline stuff, to __init__
        discard_sbml_annotations = False
        sbml_importer = amici.SbmlImporter(
            sbml_model,
            discard_annotations=discard_sbml_annotations,
        )
        sbml_model = sbml_importer.sbml

        # check for time-point-specific placeholders
        #  for now, we only support:
        #  * observable placeholders that are replaced by the same expression
        #    for all measurements for a given experiment
        #  * noise placeholders that are replaced by the same expression
        #    for all measurements for a given experiment
        #  * noise placeholders if there is only a single placeholder which
        #    is replaced by literals for all measurements for a given
        #    experiment
        for experiment in self.petab_problem.experiments:
            measurements = self.petab_problem.get_measurements_for_experiment(
                experiment
            )
            observable_overrides = {}
            noise_overrides = {}
            for measurement in measurements:
                observable_overrides.setdefault(
                    measurement.observable_id, set()
                ).add(tuple(measurement.observable_parameters))
                noise_overrides.setdefault(
                    measurement.observable_id, set()
                ).add(tuple(measurement.noise_parameters))

            for observable_id, overrides in observable_overrides.items():
                if len(overrides) > 1:
                    raise NotImplementedError(
                        f"Observable {observable_id} has multiple "
                        "timepoint-specific mappings for observable "
                        "parameters. "
                        "This is not supported by AMICI."
                    )
            for observable_id, overrides in noise_overrides.items():
                if len(overrides) > 1:
                    if len(next(iter(overrides))) == 1 and all(
                        isinstance(p[0], numbers.Number) for p in overrides
                    ):
                        continue

                    raise NotImplementedError(
                        f"Observable {observable_id} has multiple "
                        "timepoint-specific mappings for noise parameters. "
                        "This is not supported by AMICI."
                    )
                if len(overrides) == 1 and next(iter(overrides)) == ():
                    # this is a single literal, which is fine
                    continue
            print("TODO")
            print(observable_overrides)
            print(noise_overrides)

        # TODO
        # allow_n_noise_pars = (
        #     not petab.lint.observable_table_has_nontrivial_noise_formula(
        #         petab_problem.observable_df
        #     )
        # )
        # if (
        #         not jax
        #         and petab_problem.measurement_df is not None
        #         and petab.lint.measurement_table_has_timepoint_specific_mappings(
        #     petab_problem.measurement_df,
        #     allow_scalar_numeric_noise_parameters=allow_n_noise_pars,
        # )
        # ):
        #     raise ValueError(
        #         "AMICI does not support importing models with timepoint specific "
        #         "mappings for noise or observable parameters. Please flatten "
        #         "the problem and try again."
        #     )

        fixed_parameters |= _get_fixed_parameters_sbml(
            petab_problem=self.petab_problem,
            non_estimated_parameters_as_constants=non_estimated_parameters_as_constants,
        )

        fixed_parameters = list(sorted(fixed_parameters))

        logger.debug(f"Fixed parameters are {fixed_parameters}")
        logger.info(f"Overall fixed parameters: {len(fixed_parameters)}")
        logger.info(
            "Variable parameters: "
            + str(
                len(sbml_model.getListOfParameters()) - len(fixed_parameters)
            )
        )

        # Create Python module from SBML model
        if self._jax:
            sbml_importer.sbml2jax(
                model_name=self.model_id,
                output_dir=self.outdir,
                observation_model=observation_model,
                verbose=self._verbose,
                # **kwargs,
            )
            return sbml_importer
        else:
            # TODO:
            allow_reinit_fixpar_initcond = True
            sbml_importer.sbml2amici(
                model_name=self.model_id,
                output_dir=self.outdir,
                observation_model=observation_model,
                constant_parameters=fixed_parameters,
                allow_reinit_fixpar_initcond=allow_reinit_fixpar_initcond,
                verbose=self._verbose,
                # FIXME: simplification takes ages for Smith_BMCSystBiol2013
                #  due to nested piecewises / Heavisides?!
                simplify=None,
                # **kwargs,
            )
        # TODO: ensure that all estimated parameters are present as
        #  (non-constant) parameters in the model

        if self._compile:
            # check that the model extension was compiled successfully
            _ = self.import_module()
            # model = model_module.getModel()
            # TODO check_model(amici_model=model, petab_problem=petab_problem)

        return sbml_importer

    def _workaround_observable_parameters(
        self, output_parameter_defaults: dict[str, float] = None
    ) -> None:
        """
        Add any output parameters that are introduced via PEtab to the model.

        This can be placeholder parameters or any other parameters that are
        introduced in observableFormula or noiseFormula in the observable
        table, or in observableParameters or noiseParameters in the measurement
        table.
        """
        problem = self.petab_problem
        output_parameters = problem.get_output_parameters()

        logger.debug(
            "Adding output parameters to model: "
            f"{list(sorted(output_parameters))}"
        )
        output_parameter_defaults = output_parameter_defaults or {}
        if extra_pars := (
            set(output_parameter_defaults) - set(output_parameters)
        ):
            raise ValueError(
                "Default output parameter values were given for "
                f"{extra_pars}, but they those are not output parameters."
            )

        for par in sorted(output_parameters):
            _add_global_parameter(
                sbml_model=problem.model.sbml_model,
                parameter_id=par,
                value=output_parameter_defaults.get(par, 0.0),
            )

    def _get_observation_model(
        self,
    ) -> list[MeasurementChannel]:
        """Get the observation model from the PEtab problem."""
        return [
            MeasurementChannel(
                id_=observable.id,
                name=observable.name or observable.id,
                formula=observable.formula,
                noise_distribution=observable.noise_distribution,
                sigma=observable.noise_formula,
            )
            for observable in self.petab_problem.observables or []
        ]

    def import_module(self, force_import: bool = False) -> amici.ModelModule:
        """Import the generated model module.

        :param force_import:
            Whether to force re-import even if the model module already exists.
        :return: The imported model module.
        """
        if not self.outdir.is_dir() or force_import:
            self._do_import()

        return amici.import_model_module(
            self.model_id,
            self.outdir,
        )

    # def get_model(self):
    #     """Create the model."""
    #     if self._sym_model is None:
    #         self._sym_model = self.get_sym_model()
    #
    #     return self._sym_model

    def create_simulator(self, force_import: bool = False) -> PetabSimulator:
        """
        Create a PEtab simulator for the imported model.

        :param force_import:
            Whether to force re-import even if the model module already exists.
        :return: The created PEtab simulator.
        """
        model = self.import_module(force_import=force_import).get_model()
        em = ExperimentManager(model=model, petab_problem=self.petab_problem)
        return PetabSimulator(em=em)


class ExperimentManager:
    # TODO: support for pscale
    """
    Handles the creation of ExpData objects for a given model and PEtab
    problem.

    The assumption is that we have a set of :class:`amici.ExpData` objects,
    one for each PEtab experiment.
    Those are updated based on a set of global parameters (PEtab
    problem parameters, as opposed to model parameters for a single experiment
    period).

    :param model: The AMICI model to use.
    :param petab_problem: The PEtab problem to use.
        This is expected to be the output of
        `petab.v2.ExperimentsToSbmlConverter` or an equivalent problem.
        This object must not be modified after the creation of this
        `ExperimentManager` instance.
    """

    def __init__(
        self,
        model: amici.Model,
        petab_problem: v2.Problem,
    ):
        self._model: amici.Model = model
        self._petab_problem: v2.Problem = petab_problem
        self._state_ids: tuple[str] = self._model.get_state_ids()
        self._parameter_ids: tuple[str] = self._model.get_parameter_ids()
        self._fixed_parameter_ids: tuple[str] = (
            self._model.get_fixed_parameter_ids()
        )
        # maps parameter IDs to parameter indices in the model
        self._pid_to_idx: dict[str, int] = {
            id_: i for i, id_ in enumerate(self._parameter_ids)
        }
        self._fixed_pid_to_idx: dict[str, int] = {
            id_: i for i, id_ in enumerate(self._fixed_parameter_ids)
        }

        # create a new model instance from the model module from which
        #  we can get the default parameters
        model0 = model.module.get_model()
        self._original_p = np.array(model0.get_parameters())
        self._original_k = np.array(model0.get_fixed_parameters())

    def create_edatas(self) -> list[amici.ExpData]:
        """Create ExpData objects for all experiments."""
        # TODO: only those with measurements?
        # TODO: yield?

        edatas = []
        for experiment in self._petab_problem.experiments:
            edata = self.create_edata(experiment)
            edatas.append(edata)

        return edatas

    def create_edata(
        self, experiment: v2.core.Experiment | str | None
    ) -> amici.ExpData:
        """Create an ExpData object for a single experiment.

        Sets only parameter-independent values (timepoints, measurements,
        and constant noise). No parameters or initial conditions.

        :param experiment: The experiment or experiment ID to create the
            ExpData for.
        """
        if isinstance(experiment, str):
            experiment = self._petab_problem[experiment]

        edata = amici.ExpData(self._model)
        edata.id = experiment.id

        if len(experiment.periods) > 2:
            raise AssertionError(
                f"Expected <= 2 periods, got {len(experiment.periods)} "
                f"for experiment {experiment.id}."
            )

        # Set fixed parameters.
        # After converting experiments to events, all remaining
        # condition parameters are constants.
        if experiment.periods:

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
            if len(experiment.periods) >= int(
                1 + experiment.has_preequilibration
            ):
                # simulation period
                main_period = experiment.sorted_periods[
                    int(experiment.has_preequilibration)
                ]
                edata.fixed_parameters = get_k(main_period)
                edata.t_start = main_period.time

        ##########################################################################
        # timepoints

        # get the required time points: this is the superset of timepoints
        #  of the measurements of all observables, including the different
        #  replicates
        # TODO extract function
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

        ##########################################################################
        # measurements and sigmas
        y, sigma_y = self._get_measurements_and_sigmas(
            measurements=measurements,
            timepoints_w_reps=timepoints_w_reps,
            observable_ids=self._model.get_observable_ids(),
        )
        edata.set_observed_data(y.flatten())
        edata.set_observed_data_std_dev(sigma_y.flatten())

        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(
                f"Created ExpData id={edata.id}, "
                f"k_preeq={edata.fixed_parameters_pre_equilibration}, "
                f"k={edata.fixed_parameters}"
            )

        return edata

    def _get_measurements_and_sigmas(
        self,
        measurements: list[v2.Measurement],
        timepoints_w_reps: Sequence[numbers.Number],
        observable_ids: Sequence[str],
    ) -> tuple[np.array, np.array]:
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
        self, edata: amici.ExpData, problem_parameters: dict[str, float]
    ) -> None:
        """Apply problem parameters.

        Update the parameter-dependent values of the given ExpData instance
        according to the provided problem parameters.

        This assumes that:

        * the ExpData instance was created by `create_edata`,
        * no other changes except for calls to `apply_parameters` were made,
        * and the PEtab problem was not modified since the creation of this
          `ExperimentManager` instance.

        :param edata: The ExpData instance to be updated.
            In case of errors, the state of `edata` is undefined.
        :param problem_parameters: Problem parameters to be applied.
        """
        # TODO: support ndarray in addition to dict

        # TODO: must handle output overrides here, or add them to the events
        par_vals = np.array(self._original_p)
        pid_to_idx = self._pid_to_idx
        experiment_id = edata.id
        experiment = self._petab_problem[experiment_id]

        # plist -- estimated parameters + those mapped via placeholders
        # TODO sufficient to set them during creation of edata or allow dynamic fixing of parameters?
        #  store list of sensitivity parameter in class instead of using x_free_ids or estimate=True
        plist = []
        placeholder_mappings = self._get_placeholder_mapping(experiment)
        estimated_par_ids = self._petab_problem.x_free_ids
        for model_par_idx, model_par_id in enumerate(
            self._model.get_parameter_ids()
        ):
            if model_par_id in estimated_par_ids or (
                (maps_to := placeholder_mappings.get(model_par_id)) is not None
                and maps_to in estimated_par_ids
            ):
                plist.append(model_par_idx)
        edata.plist = plist

        # Update fixed parameters in case they are affected by problem
        #  parameters (i.e., parameter table parameters)
        fixed_par_vals = np.array(edata.fixed_parameters)
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

        # TODO: remove; this should be handled on construction of the ExpData
        # experiment set indicator (see petab's `experiments_to_events`)
        ind_id = ExperimentsToSbmlConverter.get_experiment_indicator(
            experiment_id
        )
        if (idx := pid_to_idx.get(ind_id)) is not None:
            par_vals[idx] = 1

        # apply problem parameter values to identical model parameters
        #  any other parameter mapping is handled by events
        for k, v in problem_parameters.items():
            if (idx := pid_to_idx.get(k)) is not None:
                par_vals[idx] = v

        # TODO handle placeholders
        #  check that all periods use the same overrides (except for numeric sigmas)
        #  see do_import() for details
        # TODO extract function
        measurements = self._petab_problem.get_measurements_for_experiment(
            experiment
        )
        # encountered placeholders and their overrides
        # (across all observables -- placeholders IDs are globally unique)
        overrides = {}
        for m in measurements:
            obs: Observable = self._petab_problem[m.observable_id]
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

            # TODO: set sigmas via parameters or .sigmay
            if obs.noise_placeholders:
                for placeholder, override in zip(
                    obs.noise_placeholders, m.noise_parameters, strict=True
                ):
                    placeholder = str(placeholder)
                    if (
                        prev_override := overrides.get(placeholder)
                    ) is not None and prev_override != override:
                        # TODO: via .sigmay if numeric
                        raise NotImplementedError(
                            "Timepoint-specific observable placeholder "
                            "overrides are not supported"
                        )
                    if (idx := pid_to_idx.get(placeholder)) is not None:
                        if override.is_Number:
                            par_vals[idx] = float(override)
                        else:
                            par_vals[idx] = problem_parameters[str(override)]
                    else:
                        raise NotImplementedError(
                            f"Cannot handle override `{placeholder}' => '{override}'"
                        )
            # print("ExperimentManager.apply_parameters:")
            # print(m)
            # print(dict(zip(self._parameter_ids, map(float, par_vals))))
        # TODO: set all unused placeholders to NaN to make it easier to spot problems?
        edata.parameters = par_vals

        # TODO debug, remove
        # print(
        #     f"Parameters: {dict(zip(self._parameter_ids, map(float, par_vals)))}"
        # )

    @property
    def petab_problem(self) -> v2.Problem:
        """The PEtab problem used by this ExperimentManager.

        This must not be modified.
        """
        return self._petab_problem

    @property
    def model(self) -> amici.Model:
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
    Simulator for PEtab problems.

    This class is used to simulate all experiments of a given PEtab problem
    using a given AMICI model and solver, and to aggregate the results.

    :param em: The ExperimentManager to generate the ExpData objects.
    :param solver: The AMICI solver to use for the simulations.
        If not provided, a new solver with default settings will be used.
    """

    def __init__(
        self, em: ExperimentManager, solver: amici.Solver | None = None
    ):
        self._petab_problem = em.petab_problem
        self._model = em.model
        self._solver = (
            solver if solver is not None else self._model.create_solver()
        )
        self._exp_man: ExperimentManager = em

    def simulate(
        self, problem_parameters: dict[str, float] = None
    ) -> dict[str, Any]:
        # TODO params: dict|np.ndarray|None?
        """Simulate all experiments of the given PEtab problem.

        :return:
            Dictionary of

            * the summed cost function value (``LLH``),
            * list of :class:`amici.amici.ReturnData` (``RDATAS``)
              for each experiment,
            * list of :class:`amici.amici.ExpData` (``EDATAS``)
              for each experiment

           Note that the returned :class:`amici.amici.ExpData` instances
           may be changed by subsequent calls to this function.
           Create a copy if needed.
        """
        if problem_parameters is None:
            # use default parameters, i.e., nominal values for all parameters
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

        rdatas = amici.run_simulations(self._model, self._solver, edatas)
        from . import EDATAS, LLH, RDATAS, SLLH

        return {
            RDATAS: rdatas,
            LLH: sum(rdata.llh for rdata in rdatas),
            SLLH: self._aggregate_sllh(rdatas),
            EDATAS: edatas,
        }

    def _aggregate_sllh(
        self, rdatas: Sequence[amici.ReturnDataView]
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
            if rdata.status != amici.AMICI_SUCCESS:
                return None
            # Condition simulation result does not provide SLLH.
            if rdata.sllh is None:
                raise ValueError(
                    f"The sensitivities of the likelihood for experiment "
                    f"{rdata.id} were not computed."
                )

        parameter_ids = self._model.get_parameter_ids()

        # still needs parameter mapping for placeholders
        for rdata in rdatas:
            experiment = self._petab_problem[rdata.id]
            placeholder_mappings = self._exp_man._get_placeholder_mapping(
                experiment
            )
            for model_par_idx, sllh in zip(
                rdata.plist, rdata.sllh, strict=True
            ):
                model_par_id = problem_par_id = parameter_ids[model_par_idx]
                if maps_to := placeholder_mappings.get(model_par_id):
                    problem_par_id = maps_to

                sllh_total[problem_par_id] = (
                    sllh_total.get(problem_par_id, 0.0) + sllh
                )
        return sllh_total


def _set_default_experiment(
    problem: v2.Problem, id_: str = _DEFAULT_EXPERIMENT_ID
) -> None:
    """Replace any empty experiment ID in the measurement table by
    a new dummy experiment with ID ``id_``.

    :param problem: The PEtab problem. This will be modified in place.
    """
    if not any(m.experiment_id is None for m in problem.measurements):
        return

    # create dummy experiment
    problem += v2.core.Experiment(
        id=id_,
    )

    for m in problem.measurements:
        if m.experiment_id is None:
            m.experiment_id = id_


def rdatas_to_measurement_df(
    rdatas: Sequence[amici.ReturnData],
    model: amici.AmiciModel,
    petab_problem: v2.Problem,
) -> pd.DataFrame:
    """
    Create a measurement dataframe in the PEtab format from the passed
    ``rdatas`` and own information.

    :param rdatas:
        A sequence of rdatas with the ordering of
        :func:`petab.get_simulation_conditions`.

    :param model:
        AMICI model used to generate ``rdatas``.

    :param petab_problem:
        The PEtab problem used to generate ``rdatas``.

    :return:
        A dataframe built from the rdatas in the format of ``measurement_df``.
    """

    measurement_df = petab_problem.measurement_df
    observable_ids = model.get_observable_ids()
    rows = []
    # iterate over conditions
    for rdata in rdatas:
        experiment_id = rdata.id

        # current simulation matrix
        y = rdata.y
        # time array used in rdata
        t = list(rdata.ts)

        # extract rows for condition
        cur_measurement_df = measurement_df[
            measurement_df[v2.C.EXPERIMENT_ID] == experiment_id
        ]

        # iterate over entries for the given condition
        # note: this way we only generate a dataframe entry for every
        # row that existed in the original dataframe. if we want to
        # e.g. have also timepoints non-existent in the original file,
        # we need to instead iterate over the rdata['y'] entries
        for _, row in cur_measurement_df.iterrows():
            # copy row
            row_sim = copy.deepcopy(row)

            # extract simulated measurement value
            timepoint_idx = t.index(row[v2.C.TIME])
            observable_idx = observable_ids.index(row[v2.C.OBSERVABLE_ID])
            measurement_sim = y[timepoint_idx, observable_idx]

            # change measurement entry
            row_sim[v2.C.MEASUREMENT] = measurement_sim

            rows.append(row_sim)

    return pd.DataFrame(rows)


def rdatas_to_simulation_df(
    rdatas: Sequence[amici.ReturnData],
    model: amici.AmiciModel,
    petab_problem: v2.Problem,
) -> pd.DataFrame:
    """
    Create a simulation dataframe in the PEtab format from the passed
    ``rdatas`` and own information.

    :param rdatas:
        A sequence of rdatas with the ordering of
        :func:`petab.get_simulation_conditions`.

    :param model:
        AMICI model used to generate ``rdatas``.

    :param petab_problem:
        The PEtab problem used to generate ``rdatas``.

    :return:
        A dataframe built from the rdatas in the format of
        ``petab_problem.measurement_df``.
    """
    measurement_df = rdatas_to_measurement_df(rdatas, model, petab_problem)

    simulation_df = measurement_df.rename(
        columns={v2.C.MEASUREMENT: v2.C.SIMULATION}
    )

    # revert setting default experiment Id
    simulation_df.loc[
        simulation_df[v2.C.EXPERIMENT_ID] == _DEFAULT_EXPERIMENT_ID,
        v2.C.EXPERIMENT_ID,
    ] = np.nan

    return simulation_df


def has_timepoint_specific_overrides(
    problem: v2.Problem,
    ignore_scalar_numeric_noise_parameters: bool = False,
    ignore_scalar_numeric_observable_parameters: bool = False,
) -> bool:
    """Check if the measurements have timepoint-specific observable or
    noise parameter overrides.

    :param ignore_scalar_numeric_noise_parameters:
        ignore scalar numeric assignments to noiseParameter placeholders

    :param ignore_scalar_numeric_observable_parameters:
        ignore scalar numeric assignments to observableParameter
        placeholders

    :return: True if the problem has timepoint-specific overrides, False
        otherwise.
    """
    if not problem.measurements:
        return False

    from petab.v1.core import get_notnull_columns
    from petab.v1.lint import is_scalar_float

    measurement_df = problem.measurement_df

    # mask numeric values
    for col, allow_scalar_numeric in [
        (
            v2.C.OBSERVABLE_PARAMETERS,
            ignore_scalar_numeric_observable_parameters,
        ),
        (v2.C.NOISE_PARAMETERS, ignore_scalar_numeric_noise_parameters),
    ]:
        if col not in measurement_df:
            continue

        measurement_df[col] = measurement_df[col].apply(str)

        if allow_scalar_numeric:
            measurement_df.loc[
                measurement_df[col].apply(is_scalar_float), col
            ] = ""

    grouping_cols = get_notnull_columns(
        measurement_df,
        _POSSIBLE_GROUPVARS_FLATTENED_PROBLEM,
    )
    grouped_df = measurement_df.groupby(grouping_cols, dropna=False)

    grouping_cols = get_notnull_columns(
        measurement_df,
        [
            v2.C.MODEL_ID,
            v2.C.OBSERVABLE_ID,
            v2.C.EXPERIMENT_ID,
        ],
    )
    grouped_df2 = measurement_df.groupby(grouping_cols)

    # data frame has timepoint specific overrides if grouping by noise
    # parameters and observable parameters in addition to observable and
    # experiment id yields more groups
    return len(grouped_df) != len(grouped_df2)


def _get_flattened_id_mappings(
    petab_problem: v2.Problem,
) -> dict[str, str]:
    """Get mapping from flattened to unflattened observable IDs.

    :param petab_problem:
        The unflattened PEtab problem.
    :returns:
        A mapping from flattened ID to original observable ID.
    """
    from petab.v1.core import (
        get_notnull_columns,
        get_observable_replacement_id,
    )

    groupvars = get_notnull_columns(
        petab_problem.measurement_df, _POSSIBLE_GROUPVARS_FLATTENED_PROBLEM
    )
    mappings: dict[str, str] = {}

    old_observable_ids = {obs.id for obs in petab_problem.observables}
    for groupvar, _ in petab_problem.measurement_df.groupby(
        groupvars, dropna=False
    ):
        observable_id = groupvar[groupvars.index(v2.C.OBSERVABLE_ID)]
        observable_replacement_id = get_observable_replacement_id(
            groupvars, groupvar
        )

        logger.debug(f"Creating synthetic observable {observable_id}")
        if (
            observable_id != observable_replacement_id
            and observable_replacement_id in old_observable_ids
        ):
            raise RuntimeError(
                "could not create synthetic observables "
                f"since {observable_replacement_id} was "
                "already present in observable table"
            )

        mappings[observable_replacement_id] = observable_id

    return mappings


def flatten_timepoint_specific_output_overrides(
    petab_problem: v2.Problem,
) -> None:
    """Flatten timepoint-specific output parameter overrides.

    If the PEtab problem definition has timepoint-specific
    `observableParameters` or `noiseParameters` for the same observable,
    replace those by replicating the respective observable.

    This is a helper function for some tools which may not support such
    timepoint-specific mappings. The observable table and measurement table
    are modified in place.

    :param petab_problem:
        PEtab problem to work on. Modified in place.
    """
    from petab.v1.core import (
        get_notnull_columns,
        get_observable_replacement_id,
    )

    # Update observables
    def create_new_observable(old_id, new_id) -> Observable:
        if old_id not in petab_problem.observable_df.index:
            raise ValueError(
                f"Observable {old_id} not found in observable table."
            )

        # copy original observable and update ID
        observable: Observable = copy.deepcopy(petab_problem[old_id])
        observable.id = new_id

        # update placeholders
        old_obs_placeholders = observable.observable_placeholders
        old_noise_placeholders = observable.noise_placeholders
        suffix = new_id.removeprefix(old_id)
        observable.observable_placeholders = [
            f"{sym.name}{suffix}" for sym in observable.observable_placeholders
        ]
        observable.noise_placeholders = [
            f"{sym.name}{suffix}" for sym in observable.noise_placeholders
        ]

        # placeholders in formulas
        subs = dict(
            zip(
                old_obs_placeholders,
                observable.observable_placeholders,
                strict=True,
            )
        )
        observable.formula = observable.formula.subs(subs)
        subs |= dict(
            zip(
                old_noise_placeholders,
                observable.noise_placeholders,
                strict=True,
            )
        )
        observable.noise_formula = observable.noise_formula.subs(subs)

        return observable

    mappings = _get_flattened_id_mappings(petab_problem)

    petab_problem.observable_tables = [
        v2.ObservableTable(
            [
                create_new_observable(old_id, new_id)
                for new_id, old_id in mappings.items()
            ]
        )
    ]

    # Update measurements
    groupvars = get_notnull_columns(
        petab_problem.measurement_df, _POSSIBLE_GROUPVARS_FLATTENED_PROBLEM
    )
    for measurement_table in petab_problem.measurement_tables:
        for measurement in measurement_table.measurements:
            # TODO: inefficient, but ok for a start
            group_vals = (
                v2.MeasurementTable([measurement])
                .to_df()
                .iloc[0][groupvars]
                .tolist()
            )
            new_obs_id = get_observable_replacement_id(groupvars, group_vals)
            measurement.observable_id = new_obs_id


def unflatten_simulation_df(
    simulation_df: pd.DataFrame,
    petab_problem: v2.Problem,
) -> pd.DataFrame:
    """Unflatten simulations from a flattened PEtab problem.

    A flattened PEtab problem is the output of applying
    :func:`flatten_timepoint_specific_output_overrides` to a PEtab problem.

    :param simulation_df:
        The simulation dataframe. A dataframe in the same format as a PEtab
        measurements table, but with the ``measurement`` column switched
        with a ``simulation`` column.
    :param petab_problem:
        The unflattened PEtab problem.
    :returns:
        The simulation dataframe for the unflattened PEtab problem.
    """
    mappings = _get_flattened_id_mappings(petab_problem)
    original_observable_ids = simulation_df[v2.C.OBSERVABLE_ID].replace(
        mappings
    )
    unflattened_simulation_df = simulation_df.assign(
        **{
            v2.C.OBSERVABLE_ID: original_observable_ids,
        }
    )
    return unflattened_simulation_df


def _get_fixed_parameters_sbml(
    petab_problem: v2.Problem,
    non_estimated_parameters_as_constants=True,
) -> set[str]:
    """
    Determine, set and return fixed model parameters.

    :param petab_problem:
        The PEtab problem instance

    :param non_estimated_parameters_as_constants:
        Whether parameters marked as non-estimated in PEtab should be
        considered constant in AMICI. Setting this to ``True`` will reduce
        model size and simulation times. If sensitivities with respect to those
        parameters are required, this should be set to ``False``.

    :return:
        List of IDs of (AMICI) parameters that are not estimated.
    """
    if not petab_problem.model.type_id == MODEL_TYPE_SBML:
        raise ValueError("Not an SBML model.")

    if not non_estimated_parameters_as_constants:
        raise NotImplementedError(
            "Only non_estimated_parameters_as_constants=True is supported currently."
        )

    # everything
    # 1) that is a parameter in AMICI
    # and
    # 2) that is not flagged as estimated in PEtab
    # and
    # 3) for which there is no condition, where this parameter occurs as a
    #    targetId where the targetValue expression contains any estimated
    #    parameters
    #    TODO: if we assume that condition table changes have been converted to
    #     events, we can skip this check. indicator variables are always
    #     literal numbers. right?

    sbml_model = petab_problem.model.sbml_model

    # What will be implemented as a parameter in the amici model?
    # TODO: constant SBML parameters - what else?
    amici_parameters = {
        p.getId()
        for p in sbml_model.getListOfParameters()
        if p.getConstant() is True
        # TODO: literal is okay?
        # TODO: collect IAs once
        and sbml_model.getInitialAssignmentBySymbol(p.getId()) is None
    }

    estimated_parameters = set(petab_problem.x_free_ids)

    return amici_parameters - estimated_parameters
