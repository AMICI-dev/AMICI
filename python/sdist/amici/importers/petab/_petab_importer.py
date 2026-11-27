"""PEtab v2 import."""

from __future__ import annotations

import copy
import difflib
import inspect
import logging
import numbers
import re
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from pprint import pprint
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
import sympy as sp
from petab import v1 as v1
from petab import v2 as v2
from petab.v2 import ExperimentPeriod, Observable
from petab.v2.converters import ExperimentsToSbmlConverter
from petab.v2.models import MODEL_TYPE_PYSB, MODEL_TYPE_SBML

import amici
from amici import get_model_dir
from amici._model_dir import is_model_dir_complete
from amici._symbolic import DEModel, Event
from amici.importers.utils import (
    RESERVED_SYMBOLS,
    MeasurementChannel,
    amici_time_symbol,
    symbol_with_assumptions,
)
from amici.logging import get_logger

from .v1._sbml_import import _add_global_parameter

if TYPE_CHECKING:
    import pysb

    from amici.sim._parameter_scales import ParameterScale

__all__ = [
    "PetabImporter",
    "rdatas_to_measurement_df",
    "rdatas_to_simulation_df",
    "flatten_timepoint_specific_output_overrides",
    "unflatten_simulation_df",
    "has_timepoint_specific_overrides",
]
logger = get_logger(__name__, log_level=logging.INFO)


#: Arguments of the model import functions that are set by
#: :class:`PetabImporter` and can't be passed via ``import_kwargs``.
_RESERVED_IMPORT_KWARGS = frozenset(
    {
        "self",
        "model",
        "model_name",
        "output_dir",
        "observation_model",
        "fixed_parameters",
        "verbose",
        "compile",
        "skip_if_unchanged",
        "hybridization",
        "pysb_model_has_obs_and_noise",
        "_events",
        # required for PEtab semantics
        "allow_reinit_fixpar_initcond",
        "hardcode_symbols",
    }
)

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


def _all_placeholder_parameters(observables: list[v2.Observable]) -> set[str]:
    """Return all observable/noise placeholder parameter names declared
    across the given observables (via ``observablePlaceholders`` /
    ``noisePlaceholders``)."""
    return {
        str(placeholder)
        for obs in observables
        for placeholder in (
            *obs.observable_placeholders,
            *obs.noise_placeholders,
        )
    }


def _has_bare_noise_placeholder(
    petab_problem: v2.Problem, observable: v2.Observable
) -> bool:
    """Check whether the noise formula is just the single noise placeholder.

    Only then can numeric noise parameter overrides be used as sigmas
    directly. For PySB models, this resolves the ``{observable_id}_sigma``
    expression introduced by :func:`_add_observation_model_pysb`.
    """
    if len(observable.noise_placeholders) != 1:
        return False

    formula = observable.noise_formula
    if petab_problem.model.type_id == MODEL_TYPE_PYSB and isinstance(
        formula, sp.Symbol
    ):
        expressions = petab_problem.model.model.expressions
        if formula.name == f"{observable.id}_sigma" and (
            sigma_expr := expressions.get(formula.name)
        ):
            formula = sigma_expr.expr

    return isinstance(formula, sp.Symbol) and formula.name == str(
        observable.noise_placeholders[0]
    )


def _substitute_jax_placeholder_parameters(
    observation_model: list[MeasurementChannel],
    observables: list[v2.Observable],
) -> list[MeasurementChannel]:
    """Substitute PEtab's observable/noise placeholder parameters with
    their generic, index-based form.

    PEtab declares observable/noise placeholders explicitly per observable
    via ``observablePlaceholders`` / ``noisePlaceholders`` -- these need not
    follow any particular naming convention (e.g. the PEtab v1-style
    ``noiseParameter{n}_{observableId}`` is common, but not required; see
    ``petab.v2.lint.CheckOverridesMatchPlaceholders``, which matches
    placeholders to measurement overrides purely by count and position).
    The JAX backend maps observable/noise parameter overrides to the
    generated model by INDEX (see :mod:`amici.sim.jax.petab`), not by
    placeholder name or observable identity, so all observables must share
    the same ``observableParameter{n}``/``noiseParameter{n}`` symbols for a
    given index; substitute each declared placeholder (in declaration
    order) with the corresponding generic symbol.

    This mirrors the substitution already performed for the legacy PEtab v1
    JAX import, see
    :func:`amici.importers.petab.v1._sbml_import._workaround_observable_parameters`.

    Returns new :class:`MeasurementChannel` instances; the input is not
    modified. Only affects the (ephemeral) importer-internal observation
    model used for JAX code generation -- the PEtab problem itself is left
    untouched.

    :return:
        The new measurement channels. Each declares the generic placeholders
        corresponding to the placeholders of its observable.
    """
    observables_by_id = {obs.id: obs for obs in observables}
    result = []
    for channel in observation_model:
        channel = copy.copy(channel)
        channel.noise_placeholders = list(channel.noise_placeholders)
        channel.observable_placeholders = list(channel.observable_placeholders)
        obs = observables_by_id.get(channel.id)
        if obs is not None:
            for attr, placeholders, generic_prefix, declared_attr in (
                (
                    "sigma",
                    obs.noise_placeholders,
                    "noiseParameter",
                    "noise_placeholders",
                ),
                (
                    "formula",
                    obs.observable_placeholders,
                    "observableParameter",
                    "observable_placeholders",
                ),
            ):
                expr = getattr(channel, attr)
                if expr is None or not placeholders:
                    continue
                expr = sp.sympify(expr)
                # declare all placeholders, also those that don't occur in
                # the formula: values are bound by position
                setattr(
                    channel,
                    declared_attr,
                    [
                        f"{generic_prefix}{i}"
                        for i in range(1, len(placeholders) + 1)
                    ],
                )
                for i, placeholder in enumerate(placeholders, start=1):
                    # substitute using the placeholder symbol itself, not a
                    # freshly created one: petab.v2 symbols carry assumptions
                    # (e.g. real=True) that a plain sp.Symbol(name) lacks, so
                    # they would not compare equal for .subs()
                    expr = expr.subs(
                        placeholder, sp.Symbol(f"{generic_prefix}{i}")
                    )
                setattr(channel, attr, expr)
        result.append(channel)
    return result


# TODO: how to handle SBML vs PySB, jax vs sundials?
#  -> separate importers or subclasses?
# TODO: How to handle multi-model-problems?


class PetabImporter:
    """
    Importer for PEtab2 problems.

    This class is used to create an AMICI model from a PEtab problem.

    The underlying SBML or PySB model will be modified to encode the
    experiments defined in the PEtab problem as events or initial conditions.

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

    # TODO remove: extra debug output
    _debug = False

    def __init__(
        self,
        petab_problem: v2.Problem | v1.Problem,
        *,
        compile_: bool = None,
        validate: bool = True,
        module_name: str = None,
        # TODO: model_id for selecting the model in multi-model problems
        # model_id: str = None,
        output_dir: str | Path = None,
        jax: bool = False,
        output_parameter_defaults: dict[str, float] | None = None,
        verbose: int | bool = logging.INFO,
        non_estimated_parameters_as_constants: bool = True,
        import_kwargs: Mapping[str, Any] | None = None,
    ):
        """
        Create a new PetabImporter instance.

        :param petab_problem: The PEtab problem to import.
        :param compile_: Whether to compile the model extension after import.
        :param validate: Whether to validate the PEtab problem before import.
        :param module_name: The name of model module to generate.
        :param output_dir:
            The output directory where the model files are written to.
            Any previous model import in this directory is removed.
            Directories containing files not created by AMICI are rejected.
        :param jax: Whether to generate a JAX model instead of a
            SUNDIALS model. Currently, only ``False`` is supported.
        :param output_parameter_defaults:
            Optional default parameter values for output parameters introduced
            in the PEtab observables table, in particular for placeholder
            parameters. A dictionary mapping parameter IDs to default values.
        :param verbose:
            The verbosity level. If ``True``, set to ``logging.INFO``.
            If ``False``, set to ``logging.WARNING``. Otherwise, use the given
            logging level.
        :param non_estimated_parameters_as_constants:
            Whether parameters marked as non-estimated in PEtab should be
            considered constant in AMICI. Setting this to ``True`` will reduce
            model size and simulation times. If sensitivities with respect to
            those parameters are required, this should be set to ``False``.
        :param import_kwargs:
            Additional keyword arguments for the model import function,
            e.g., ``{"compute_conservation_laws": False}``.
            Depending on the model type and on ``jax``, this is
            :meth:`amici.importers.sbml.SbmlImporter.sbml2amici`,
            :meth:`amici.importers.sbml.SbmlImporter.sbml2jax`,
            :func:`amici.importers.pysb.pysb2amici`, or
            :func:`amici.importers.pysb.pysb2jax`.
            Arguments that are set by the PetabImporter, such as
            ``model_name`` or ``fixed_parameters``, are not allowed.
        """
        self.petab_problem: v2.Problem = self._upgrade_or_copy_if_needed(
            petab_problem
        )
        self._verbose = (
            logging.INFO
            if verbose is True
            else logging.WARNING
            if verbose is False
            else verbose
        )
        self._output_parameter_defaults = output_parameter_defaults

        if len(self.petab_problem.models) > 1:
            raise NotImplementedError(
                "PEtab v2 importer currently only supports single-model "
                "problems."
            )

        if self.petab_problem.model.type_id not in (
            MODEL_TYPE_SBML,
            MODEL_TYPE_PYSB,
        ):
            raise NotImplementedError(
                "PEtab v2 importer currently only supports SBML and PySB "
                f"models. Got {self.petab_problem.model.type_id!r}."
            )

        if self._debug:
            print("PetabImpoter.__init__: petab_problem:")
            pprint(self.petab_problem.model_dump())
            if self.petab_problem.model.type_id == MODEL_TYPE_SBML:
                print(self.petab_problem.model.to_antimony())
            elif self.petab_problem.model.type_id == MODEL_TYPE_PYSB:
                print(self.petab_problem.model.to_str())

        self._compile = not jax if compile_ is None else compile_
        # whether the model in the output directory was checked to match the
        #  current inputs, or was imported from them
        self._inputs_checked = False
        self._sym_model: DEModel | None = None
        self._model_id = self.petab_problem.model.model_id
        self._module_name = module_name or (
            f"{self.petab_problem.id}_{self.model_id}"
            if self.petab_problem.id
            else self.model_id
        )
        if self._module_name is None:
            raise ValueError(
                "No `module_name` was provided and no model ID "
                "was specified in the PEtab problem."
            )

        self._output_dir: Path | None = (
            None if output_dir is None else Path(output_dir).absolute()
        )
        self._jax = jax
        self._non_estimated_parameters_as_constants: bool = (
            non_estimated_parameters_as_constants
        )
        self._import_kwargs = self._check_import_kwargs(import_kwargs)

        if validate:
            logger.info("Validating PEtab problem ...")
            validation_result = petab_problem.validate()
            if validation_result:
                validation_result.log()

            if validation_result.has_errors():
                raise ValueError(
                    "PEtab problem is not valid, see log messages for details."
                )

        # ensure each measurement has an experimentId to simplify processing
        _set_default_experiment(self.petab_problem)

        if self.petab_problem.model.type_id == MODEL_TYPE_SBML:
            self._preprocess_sbml()
        elif self.petab_problem.model.type_id == MODEL_TYPE_PYSB:
            self._preprocess_pysb()
        else:
            raise AssertionError()

    def _preprocess_sbml(self):
        """Pre-process the SBML-based PEtab problem to make it
        amici-compatible."""
        from petab.v2.models.sbml_model import SbmlModel

        if not isinstance(self.petab_problem.model, SbmlModel):
            raise ValueError("The PEtab problem must contain an SBML model.")

        # Convert petab experiments to events, because so far,
        #  AMICI only supports preequilibration/presimulation/simulation, but
        #  no arbitrary list of periods.
        exp_event_conv = ExperimentsToSbmlConverter(self.petab_problem)
        # This will always create a copy of the problem.
        if self._jax:
            self._unconverted_problem = exp_event_conv._original_problem
            condition_targets = {
                change.target_id
                for condition in self.petab_problem.conditions
                for change in condition.changes
            }
            mapping_petab_ids = {
                mapping.petab_id for mapping in self.petab_problem.mappings
            }
            if condition_targets.intersection(mapping_petab_ids):
                raise NotImplementedError(
                    "The JAX backend does not currently support PEtab problems where network "
                    "parameters appear in the conditions table. "
                )
        self.petab_problem = exp_event_conv.convert()
        for experiment in self.petab_problem.experiments:
            if len(experiment.periods) > 2:
                # This should never happen due to the conversion above
                raise NotImplementedError(
                    "AMICI currently does not support more than two periods."
                )

        if self._debug:
            print("PetabImpoter._preprocess_sbml: petab_problem:")
            pprint(self.petab_problem.model_dump())
            print(self.petab_problem.model.to_antimony())

    def _preprocess_pysb(self):
        """Pre-process the PySB-based PEtab problem to make it
        amici-compatible."""
        import pysb
        from petab.v2.models.pysb_model import PySBModel

        if not isinstance(self.petab_problem.model, PySBModel):
            raise ValueError("The PEtab problem must contain a PySB model.")

        _add_observation_model_pysb(self.petab_problem, self._jax)

        # Species targets are resolved against the reaction network
        pysb.bng.generate_equations(self.petab_problem.model.model)

        # Convert PEtab v2 experiments/conditions to events
        converter = ExperimentsToPySBConverter(self.petab_problem)
        self.petab_problem, self._events = converter.convert()
        if self._jax and self._events:
            raise NotImplementedError(
                "The JAX backend does not support condition changes after "
                "the first period for PySB models."
            )

    def _upgrade_or_copy_if_needed(
        self, problem: v1.Problem | v2.Problem
    ) -> v2.Problem:
        """Upgrade the problem to petab v2 if necessary.
        Otherwise, create a deep copy of the problem."""
        if isinstance(problem, v2.Problem):
            return copy.deepcopy(problem)

        raise NotImplementedError(
            "'petab_problem' must be a `petab.v2.Problem`. "
            "`petab.v1.Problem` is not directly supported, but "
            "file-based PEtab v1 problems can be upgraded via "
            "`petab.v2.Problem.from_yaml(petab_v1_yaml_file)`."
        )

    @property
    def model_id(self) -> str:
        """The model ID."""
        if self._model_id is None:
            self._model_id = self.petab_problem.model.model_id

        return self._model_id

    @property
    def output_dir(self) -> Path:
        """The output directory where the model files are written to."""
        if self._output_dir is None:
            self._output_dir = get_model_dir(self._module_name, jax=self._jax)
        return self._output_dir

    def _import_function(self) -> Callable:
        """Get the model import function for the model type and backend."""
        if self.petab_problem.model.type_id == MODEL_TYPE_SBML:
            return (
                amici.SbmlImporter.sbml2jax
                if self._jax
                else amici.SbmlImporter.sbml2amici
            )

        from amici.importers.pysb import pysb2amici, pysb2jax

        return pysb2jax if self._jax else pysb2amici

    def _check_import_kwargs(
        self, import_kwargs: Mapping[str, Any] | None
    ) -> dict[str, Any]:
        """Check the user-provided arguments for the model import function.

        :return: A copy of ``import_kwargs``.
        :raises ValueError:
            For arguments that are reserved or not accepted by the import
            function.
        """
        import_kwargs = dict(import_kwargs or {})
        if reserved := sorted(import_kwargs.keys() & _RESERVED_IMPORT_KWARGS):
            raise ValueError(
                f"The import arguments {reserved} are set by the "
                "PetabImporter and can't be passed via `import_kwargs`."
            )

        import_function = self._import_function()
        params = inspect.signature(import_function).parameters.values()
        if any(p.kind == p.VAR_KEYWORD for p in params):
            # validated by the import function
            return import_kwargs
        allowed = {
            p.name for p in params if p.kind != p.VAR_POSITIONAL
        } - _RESERVED_IMPORT_KWARGS
        if unknown := sorted(import_kwargs.keys() - allowed):
            hints = {
                key: match[0]
                for key in unknown
                if (match := difflib.get_close_matches(key, allowed, n=1))
            }
            raise ValueError(
                f"Unknown arguments for `{import_function.__name__}` in "
                f"`import_kwargs`: {unknown}."
                + (f" Did you mean: {hints}?" if hints else "")
            )
        return import_kwargs

    def _do_import_sbml(self, skip_if_unchanged: bool = False):
        """Import the model.

        Generate the symbolic model according to the given PEtab problem and
        generate the corresponding Python module.

        1. Encode all PEtab experiments as events and initial assignments
           in the SBML model.
           This leaves only (maybe) a pre-equilibration and a single
           simulation period.
        2. Add the observable parameters to the SBML model.
        """
        logger.info(f"Importing model {self.model_id!r}...")

        if not self.petab_problem.observables:
            raise NotImplementedError(
                "PEtab import without observables table "
                "is currently not supported."
            )

        logger.info(
            f"Module name is '{self._module_name}'.\n"
            f"Writing model code to '{self.output_dir}'."
        )

        observation_model = self._get_observation_model()
        if self._jax:
            observation_model = _substitute_jax_placeholder_parameters(
                observation_model, self.petab_problem.observables
            )

        logger.info(f"#Observables: {len(observation_model)}")
        logger.debug(f"Observables: {observation_model}")

        self._workaround_observable_parameters_sbml(
            output_parameter_defaults=self._output_parameter_defaults,
        )

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

        from .v1._sbml_import import show_model_info

        show_model_info(self.petab_problem.model.sbml_model)
        sbml_importer = amici.SbmlImporter(
            self.petab_problem.model.sbml_model, jax=self._jax
        )

        self._check_placeholders()

        fixed_parameters |= _get_fixed_parameters_sbml(
            petab_problem=self.petab_problem,
            non_estimated_parameters_as_constants=self._non_estimated_parameters_as_constants,
        )

        fixed_parameters = list(sorted(fixed_parameters))
        logger.info(f"Number of fixed parameters: {len(fixed_parameters)}")
        logger.debug(f"Fixed parameters are {fixed_parameters}")

        # Create Python module from SBML model
        if self._jax:
            if self.petab_problem.config.extensions.get("sciml"):
                hybridization = self._build_hybridization()
            else:
                hybridization = None
            sbml_importer.sbml2jax(
                model_name=self._module_name,
                output_dir=self.output_dir,
                observation_model=observation_model,
                fixed_parameters=fixed_parameters,
                verbose=self._verbose,
                hybridization=hybridization,
                skip_if_unchanged=skip_if_unchanged,
                **self._import_kwargs,
            )
            return sbml_importer
        else:
            sbml_importer.sbml2amici(
                model_name=self._module_name,
                output_dir=self.output_dir,
                observation_model=observation_model,
                fixed_parameters=fixed_parameters,
                allow_reinit_fixpar_initcond=True,
                verbose=self._verbose,
                compile=self._compile,
                skip_if_unchanged=skip_if_unchanged,
                **self._import_kwargs,
            )
        # TODO check_model(amici_model=model, petab_problem=petab_problem)

        return sbml_importer

    def _do_import_pysb(self, skip_if_unchanged: bool = False) -> None:
        """Import the PySB model.

        Generate the symbolic model according to the given PEtab problem and
        generate the corresponding Python module.
        """
        logger.info(f"Importing PySB model {self.model_id!r}...")

        if not self.petab_problem.observables:
            raise NotImplementedError(
                "PEtab import without observables table "
                "is currently not supported."
            )

        logger.info(
            f"Module name is '{self._module_name}'.\n"
            f"Writing model code to '{self.output_dir}'."
        )

        observation_model = self._get_observation_model()

        logger.info(f"#Observables: {len(observation_model)}")
        logger.debug(f"Observables: {observation_model}")

        pysb_model = self.petab_problem.model.model

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
        # TODO: handle self._non_estimated_parameters_as_constants

        self._check_placeholders()
        fixed_parameters = list(sorted(fixed_parameters))

        logger.info(f"Number of fixed parameters: {len(fixed_parameters)}")
        logger.debug(f"Fixed parameters are {fixed_parameters}")

        from amici.importers.pysb import pysb2amici, pysb2jax

        # Create Python module from PySB model
        if self._jax:
            pysb2jax(
                model=pysb_model,
                model_name=self._module_name,
                output_dir=self.output_dir,
                observation_model=observation_model,
                verbose=self._verbose,
                pysb_model_has_obs_and_noise=True,
                skip_if_unchanged=skip_if_unchanged,
                **self._import_kwargs,
            )
            return
        else:
            pysb2amici(
                model=pysb_model,
                model_name=self._module_name,
                output_dir=self.output_dir,
                verbose=True,
                fixed_parameters=fixed_parameters,
                observation_model=observation_model,
                pysb_model_has_obs_and_noise=True,
                compile=self._compile,
                _events=self._events,
                skip_if_unchanged=skip_if_unchanged,
                **self._import_kwargs,
            )
        # TODO check_model(amici_model=model, petab_problem=petab_problem)

    def _check_placeholders(self):
        # check for time-point-specific placeholders
        #  for now, we only support:
        #  * observable placeholders that are replaced by the same expression
        #    for all measurements for a given experiment
        #  * noise placeholders that are replaced by the same expression
        #    for all measurements for a given experiment
        #  * noise placeholders if there is only a single placeholder which
        #    is replaced by literals for all measurements for a given
        #    experiment, and the noise formula is just that placeholder
        #    (the literals are then passed as sigmas via ExpData)
        id_to_observable = {o.id: o for o in self.petab_problem.observables}
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
                        "This is not supported by AMICI. "
                        "Use `flatten_timepoint_specific_output_overrides` "
                        "on the PEtab problem before importing into amici."
                    )
            for observable_id, overrides in noise_overrides.items():
                if len(overrides) > 1:
                    if all(
                        len(p) == 1 and isinstance(p[0], numbers.Number)
                        for p in overrides
                    ) and _has_bare_noise_placeholder(
                        self.petab_problem, id_to_observable[observable_id]
                    ):
                        continue

                    raise NotImplementedError(
                        f"Observable {observable_id} has multiple "
                        "timepoint-specific mappings for noise parameters. "
                        "This is not supported by AMICI. "
                        "Use `flatten_timepoint_specific_output_overrides` "
                        "on the PEtab problem before importing into amici."
                    )
                if len(overrides) == 1 and next(iter(overrides)) == ():
                    # this is a single literal, which is fine
                    continue
            if self._debug:
                print(experiment.id)
                print(observable_overrides)
                print(noise_overrides)

    def _workaround_observable_parameters_sbml(
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

        if self._jax:
            # Declared observable/noise placeholder parameters are not
            # declared as SBML parameters for the JAX backend; they are
            # recognized as placeholders via the (substituted) observation
            # model instead, see `_substitute_jax_placeholder_parameters`.
            declared_placeholders = _all_placeholder_parameters(
                problem.observables
            )
            output_parameters = [
                par
                for par in output_parameters
                if par not in declared_placeholders
            ]

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

    def _get_net_mappings(self, net_id: str):
        """Get (petab_id, model_id) pairs for all mappings belonging to ``net_id``."""
        return (
            self.petab_problem.mapping_df.loc[
                self.petab_problem.mapping_df[v2.C.MODEL_ENTITY_ID]
                .str.split(".")
                .str[0]
                == net_id,
                v2.C.MODEL_ENTITY_ID,
            ]
            .to_dict()
            .items()
        )

    def _build_hybridization(self) -> dict[str, dict]:
        if "sciml" not in self.petab_problem.config.extensions:
            return None

        if (
            not self._jax
            or self.petab_problem.model.type_id != MODEL_TYPE_SBML
        ):
            raise NotImplementedError(
                "petab_sciml extension is currently only supported for sbml models"
            )

        from petab_sciml.standard import NNModelStandard

        config = self.petab_problem.config.extensions["sciml"]
        # TODO: only accept YAML format for now
        hybridizations = [
            pd.read_csv(hf, sep="\t") for hf in config.hybridization_files
        ]
        hybridization_table = pd.concat(hybridizations)

        input_hybridization = dict(
            zip(
                hybridization_table["targetId"],
                hybridization_table["targetValue"],
            )
        )
        output_hybridization = dict(
            zip(
                hybridization_table["targetValue"],
                hybridization_table["targetId"],
            )
        )
        observable_mapping = dict(
            zip(
                self.petab_problem.observable_df["observableFormula"],
                self.petab_problem.observable_df.index,
            )
        )
        # Build a mapping from petab entity IDs to HDF5 files for array inputs
        array_files = config.array_files
        array_input_files = {}
        for file_spec in array_files:
            import h5py

            with h5py.File(file_spec, "r") as hf:
                if "inputs" in hf:
                    for petab_id in hf["inputs"]:
                        array_input_files[petab_id] = Path(file_spec).resolve()

        condition_id_mapping = {}
        if self._unconverted_problem is not None:
            for orig_exp, conv_exp in zip(
                self._unconverted_problem.experiments,
                self.petab_problem.experiments,
            ):
                condition_id_mapping[orig_exp.id] = {
                    "original": [
                        cid
                        for p in orig_exp.periods
                        for cid in p.condition_ids
                    ],
                    "converted": [
                        cid
                        for p in conv_exp.periods
                        for cid in p.condition_ids
                    ],
                }

        non_estimated_param_ids = set(self.petab_problem.x_fixed_ids)

        hybridization = {}
        for net_id, net_config in config.neural_networks.items():
            net_mappings = list(self._get_net_mappings(net_id))
            input_mappings = [
                (pid, mid)
                for pid, mid in net_mappings
                if mid.split(".")[1].startswith("input")
            ]
            output_mappings = [
                (pid, mid)
                for pid, mid in net_mappings
                if mid.split(".")[1].startswith("output")
            ]

            model = (
                Path(net_config.location).resolve()
                if net_config.format == "equinox"
                else NNModelStandard.load_data(Path(net_config.location))
            )

            hybridization[net_id] = {
                "model": model,
                "input_vars": [
                    input_hybridization[petab_id]
                    for petab_id, _ in input_mappings
                    if petab_id in input_hybridization
                ],
                "array_inputs": {
                    petab_id: array_input_files[petab_id]
                    for petab_id, _ in input_mappings
                    if petab_id in input_hybridization
                    and input_hybridization[petab_id] == "array"
                    and petab_id in array_input_files
                },
                "output_vars": {
                    output_hybridization[petab_id]: _get_net_index(model_id)
                    for petab_id, model_id in output_mappings
                    if petab_id in output_hybridization
                },
                "observable_vars": {
                    obs_id: {
                        "index": _get_net_index(model_id),
                        "formula": formula,
                        "petab_id": petab_id,
                    }
                    for petab_id, model_id in output_mappings
                    for formula, obs_id in observable_mapping.items()
                    if re.search(r"\b" + re.escape(petab_id) + r"\b", formula)
                },
                "frozen_layers": dict(
                    _get_frozen_layers(model_id)
                    for petab_id, model_id in net_mappings
                    if petab_id in non_estimated_param_ids
                ),
                "condition_id_mapping": condition_id_mapping,
                **dict(net_config),
            }

        return hybridization

    def import_module(self, force_import: bool = False) -> amici.ModelModule:
        """Import the generated model module.

        The model is (re-)imported, unless the output directory contains
        a complete model that was imported from the same inputs
        (see :ref:`amici_model_reuse`).

        :param force_import:
            Whether to force re-import even if the model module already exists.
        :return: The imported model module.
        """
        # the inputs only need to be checked once, since the problem must
        #  not be changed after construction of the importer
        if (
            force_import
            or not self._inputs_checked
            or not is_model_dir_complete(
                self.output_dir,
                self._module_name,
                jax=self._jax,
                compiled=self._compile,
            )
        ):
            if self.petab_problem.model.type_id == MODEL_TYPE_SBML:
                self._do_import_sbml(skip_if_unchanged=not force_import)
            else:
                self._do_import_pysb(skip_if_unchanged=not force_import)
            self._inputs_checked = True

        if self._jax:
            return amici.import_model_module(
                Path(self.output_dir).stem, Path(self.output_dir).parent
            )

        return amici.import_model_module(
            self._module_name,
            self.output_dir,
        )

    def create_model(self) -> amici.sim.sundials.Model:
        """
        Create a :class:`amici.sim.sundials.Model` instance from the imported
        model.
        """
        return self.import_module().get_model()

    def create_simulator(
        self,
        force_import: bool = False,
        *,
        parameter_scales: Mapping[str, ParameterScale] | None = None,
    ) -> amici.sim.sundials.petab.PetabSimulator | amici.sim.jax.JAXProblem:
        """
        Create a PEtab simulator for the imported model.

        :param force_import:
            Whether to force re-import even if the model module already exists.
        :param parameter_scales:
            The scales of the PEtab problem parameters, as a mapping from
            problem parameter IDs to ``"lin"``, ``"log"`` (natural
            logarithm), or ``"log10"``. Parameters not included are on
            linear scale. Parameter values passed to the simulator, and the
            sensitivities it returns, are on these scales.
            See :class:`amici.sim.sundials.petab.ExperimentManager` and
            :class:`amici.sim.jax.JAXProblem`.
        :return: The created PEtab simulator.
        """
        from amici.sim._parameter_scales import get_parameter_scales
        from amici.sim.sundials.petab import ExperimentManager, PetabSimulator

        # validate before the possibly expensive model import
        parameter_scales = get_parameter_scales(
            self.petab_problem, parameter_scales
        )

        if self._jax:
            model_module = self.import_module(force_import=force_import)
            model = model_module.Model()

            from amici.sim.jax.petab import JAXProblem

            return JAXProblem(
                model,
                self.petab_problem,
                unconverted_problem=getattr(
                    self, "_unconverted_problem", None
                ),
                parameter_scales=parameter_scales,
            )

        model = self.import_module(force_import=force_import).get_model()
        em = ExperimentManager(
            model=model,
            petab_problem=self.petab_problem,
            parameter_scales=parameter_scales,
        )
        return PetabSimulator(em=em)


def _get_frozen_layers(model_id):
    layers = re.findall(r"\[(.*?)\]", model_id)
    array_attr = model_id.split(".")[-1]
    layer_id = layers[0] if len(layers) else None
    array_attr = array_attr if array_attr in ("weight", "bias") else None
    return layer_id, array_attr


def _get_net_index(model_id: str):
    matches = re.findall(r"\[(\d+)\]", model_id)
    if matches:
        return int(matches[-1])


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
    rdatas: Sequence[amici.sim.sundials.ReturnData],
    model: amici.sim.sundials.AmiciModel,
    petab_problem: v2.Problem,
) -> pd.DataFrame:
    """
    Create a measurement dataframe in the PEtab format from the passed
    ``rdatas`` and own information.

    :param rdatas:
        A sequence of :class:`amici.sim.sundials.ReturnData`.
    :param model:
        AMICI model used to generate ``rdatas``.
    :param petab_problem:
        The PEtab problem used to generate ``rdatas``.
    :return:
        A dataframe built from simulation results in `rdatas` in the format
        of the PEtab measurement table.
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
    rdatas: Sequence[amici.sim.sundials.ReturnData],
    model: amici.sim.sundials.AmiciModel,
    petab_problem: v2.Problem,
) -> pd.DataFrame:
    """
    Create a simulation dataframe in the PEtab format from the passed
    ``rdatas`` and own information.

    :param rdatas:
        A sequence of :class:`amici.sim.sundials.ReturnData`.
    :param model:
        AMICI model used to generate ``rdatas``.
    :param petab_problem:
        The PEtab problem used to generate ``rdatas``.
    :return:
        A dataframe built from simulation results in `rdatas` in the format
        of the PEtab simulation table.
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
    petab_problem: v2.Problem,
    ignore_scalar_numeric_noise_parameters: bool = False,
    ignore_scalar_numeric_observable_parameters: bool = False,
) -> bool:
    """Check if the measurements have timepoint-specific observable or
    noise parameter overrides.

    :param petab_problem:
        PEtab problem to check.
    :param ignore_scalar_numeric_noise_parameters:
        ignore scalar numeric assignments to noiseParameter placeholders
    :param ignore_scalar_numeric_observable_parameters:
        ignore scalar numeric assignments to observableParameter
        placeholders
    :return: `True` if the problem has timepoint-specific overrides, `False`
        otherwise.
    """
    if not petab_problem.measurements:
        return False

    from petab.v1.core import get_notnull_columns
    from petab.v1.lint import is_scalar_float

    measurement_df = petab_problem.measurement_df

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
        subs[old_id] = new_id
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
            "Only non_estimated_parameters_as_constants=True is supported."
        )

    # For amici constants we select everything
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
    amici_parameters = {
        p.getId()
        for p in sbml_model.getListOfParameters()
        if p.getConstant() is True
        # TODO: IAs with literals can be ignored
        # TODO(performance): collect IAs once
        and sbml_model.getInitialAssignmentBySymbol(p.getId()) is None
    }

    estimated_parameters = set(petab_problem.x_free_ids)

    return amici_parameters - estimated_parameters


def _add_observation_model_pysb(petab_problem: v2.Problem, jax: bool = False):
    """Extend PySB model by observation model as defined in the PEtab
    observables table"""
    import pysb

    pysb_model: pysb.Model = petab_problem.model.model

    # add any required output parameters
    local_syms = {
        sp.Symbol(sp.Symbol.__str__(comp), real=True): comp
        for comp in pysb_model.components
        if isinstance(comp, sp.Symbol)
    }
    # PEtab's `time` is the simulation time
    local_syms[sp.Symbol("time", real=True)] = pysb.core.time

    def process_formula(sym: sp.Basic):
        changed_formula = False
        sym = sym.subs(local_syms)
        for s in sym.free_symbols:
            if not isinstance(s, pysb.Component | pysb.core.SpecialSymbol):
                if not isinstance(s, sp.Symbol):
                    raise AssertionError(
                        f"Unexpected symbol type in observable formula: {s}, "
                        f"{type(s)}"
                    )
                name = str(s)
                p = pysb.Parameter(name, 1.0, _export=False)
                pysb_model.add_component(p)

                # placeholders for multiple observables are mapped to the
                # same symbol, so only add to local_syms when necessary
                if name not in local_syms:
                    local_syms[sp.Symbol(name, real=True)] = p

                # replace placeholder with parameter
                if jax and name != str(s):
                    changed_formula = True
                    sym = sym.subs(s, local_syms[name])
        return sym, changed_formula

    for observable in petab_problem.observables:
        sym, changed_formula = process_formula(observable.formula)
        observable.formula = sym
        sym, changed_formula = process_formula(observable.noise_formula)
        observable.noise_formula = sym

    # add observables and sigmas to pysb model
    for observable in petab_problem.observables:
        # obs_symbol = sp.sympify(observable_formula, locals=local_syms)
        if observable.id in pysb_model.expressions.keys():
            obs_expr = pysb_model.expressions[observable.id]
        else:
            obs_expr = pysb.Expression(
                observable.id, observable.formula, _export=False
            )
            pysb_model.add_component(obs_expr)
        local_syms[sp.Symbol(observable.id, real=True)] = obs_expr

        sigma_id = f"{observable.id}_sigma"
        sigma_expr = pysb.Expression(
            sigma_id, observable.noise_formula.subs(local_syms), _export=False
        )
        observable.noise_formula = sp.Symbol(sigma_id, real=True)
        pysb_model.add_component(sigma_expr)
        local_syms[sp.Symbol(sigma_id, real=True)] = sigma_expr


class ExperimentsToPySBConverter:
    """
    Convert PEtab experiments to amici events and PySB initials.

    See :meth:`convert` for details.
    """

    #: ID of the parameter that indicates whether the model is in
    #  the pre-equilibration phase (1) or not (0).
    PREEQ_INDICATOR = "_petab_preequilibration_indicator"

    #: The condition ID of the condition that sets the
    #: pre-equilibration indicator to 1.
    CONDITION_ID_PREEQ_ON = "_petab_preequilibration_on"

    #: The condition ID of the condition that sets the
    #: pre-equilibration indicator to 0.
    CONDITION_ID_PREEQ_OFF = "_petab_preequilibration_off"

    def __init__(self, petab_problem: v2.Problem):
        """Initialize the converter.

        :param petab_problem:
            The PEtab problem to convert.
            This will not be modified by this class.
        """
        from petab.v2.models.pysb_model import PySBModel

        from amici.importers.pysb import _compartment_size_parameters

        if len(petab_problem.models) > 1:
            #  https://github.com/PEtab-dev/libpetab-python/issues/392
            raise NotImplementedError(
                "Only single-model PEtab problems are supported."
            )
        if not isinstance(petab_problem.model, PySBModel):
            raise ValueError("Only PySB models are supported.")
        model = petab_problem.model.model
        # BNG evaluates compartment sizes during network generation.
        #  Changing those values later on would lead to incorrect results.
        compartment_ids = {c.name for c in model.compartments} | (
            _compartment_size_parameters(model)
        )
        petab_to_model = {
            mapping.petab_id: mapping.model_id
            for mapping in petab_problem.mappings
        }
        for cond in petab_problem.conditions:
            for change in cond.changes:
                model_id = (
                    petab_to_model.get(change.target_id) or change.target_id
                )
                if model_id in compartment_ids:
                    raise NotImplementedError(
                        "Changes to compartment sizes are not supported for "
                        f"PySB models. Got a change of {change.target_id!r}."
                    )

        #: The PEtab problem to convert. Not modified by this class.
        self._petab_problem = petab_problem
        self._events: list[Event] = []
        self._new_problem: v2.Problem | None = None

    @staticmethod
    def _get_experiment_indicator_condition_id(experiment_id: str) -> str:
        """Get the condition ID for the experiment indicator parameter."""
        return f"_petab_experiment_condition_{experiment_id}"

    def convert(
        self,
    ) -> tuple[v2.Problem, list[Event]]:
        """Convert PEtab experiments to amici events and pysb initials.

        Generate events and initial values that implement the changes
        encoded in the PEtab v2 experiment / condition table.
        This adds indicator variables to the PEtab problem and removes all
        condition changes that are implemented as events or initial values.

        Changes of the first period are applied to the initial values of
        species.
        Changes of later periods are implemented as events.
        The model must already contain the generated reaction network
        (see :func:`pysb.bng.generate_equations`), since species targets
        are resolved against it.

        :returns:
            A PEtab problem with only indicator parameters left in the
            condition table and a maximum of two periods per experiment
            (pre-equilibration and main simulation),
            and a list of events to be passed to `pysb2amici`.
        """
        self._new_problem = copy.deepcopy(self._petab_problem)
        self._events: list[Event] = []
        # target ID -> species index
        self._species_indices: dict[str, int] = {}
        model = self._new_problem.model.model
        # mapping table mappings
        self.map_petab_to_pysb: dict[str, str] = {
            mapping.petab_id: mapping.model_id
            for mapping in self._petab_problem.mappings
            if mapping.petab_id is not None and mapping.model_id is not None
        }
        # first-period changes across all experiments
        # target_id -> [(experiment_indicator, target_value), ...]
        period0_assignments: dict[str, list[tuple[str, sp.Basic]]] = {}

        self._add_preequilibration_indicator()

        for experiment in self._new_problem.experiments:
            self._convert_experiment(experiment, period0_assignments)

        # the target values may depend on parameters that are only
        #  introduced in the PEtab parameter table - those need
        #  to be added to the model
        component_names = {c.name for c in model.components}
        for sym in sorted(
            {
                sym
                for _, changes in period0_assignments.items()
                for _, target_value in changes
                for sym in target_value.free_symbols
            },
            key=str,
        ):
            if sym.name not in component_names:
                self._add_parameter(sym.name, 0)

        for target_id, changes in period0_assignments.items():
            self._apply_first_period_changes(target_id, changes)

        self._add_indicators_to_conditions()

        validation_results = self._new_problem.validate()
        validation_results.log()

        return self._new_problem, self._events

    def _convert_experiment(
        self,
        experiment: v2.Experiment,
        period0_assignments: dict[str, list[tuple[str, sp.Basic]]],
    ) -> None:
        """
        Convert a single experiment to events or initial values.

        :param experiment: The experiment to convert.
        :param period0_assignments: The changes of the first period are not
            applied directly, but are added to this dictionary
            (``target_id -> [(experiment_indicator, target_value), ...]``).
            They are combined across experiments by
            :meth:`_apply_first_period_changes`.
        """
        model = self._new_problem.model.model
        experiment.sort_periods()
        has_preequilibration = experiment.has_preequilibration

        # add experiment indicator
        exp_ind_id = self.get_experiment_indicator(experiment.id)
        if exp_ind_id in map(str, model.components):
            raise ValueError(
                f"The model has entity with ID `{exp_ind_id}`. "
                "IDs starting with `petab_` are reserved for "
                f"{self.__class__.__name__} and should not be used in the "
                "model."
            )
        self._add_parameter(exp_ind_id, 0)
        kept_periods: list[ExperimentPeriod] = []

        for i_period, period in enumerate(experiment.sorted_periods):
            if period.is_preequilibration:
                # pre-equilibration cannot be encoded as event,
                #  so we need to keep this period in the Problem.
                kept_periods.append(period)
            elif i_period == int(has_preequilibration):
                # we always keep the first non-pre-equilibration period
                #  to set the indicator parameters
                kept_periods.append(period)
            elif not period.condition_ids:
                # no condition, no changes, no need for an event,
                #  no need to keep the period unless it's the pre-equilibration
                #  or the only non-equilibration period (handled above)
                continue

            # Encode the period changes as events that trigger at the start
            #  of the period or, for the first period, as initial values.
            #  Initial values are required for the first period,
            #  because other initial values may depend on the changed values.
            if i_period == 0:
                for change in self._new_problem.get_changes_for_period(period):
                    target_value = self._check_first_period_change(
                        change, period
                    )
                    period0_assignments.setdefault(
                        change.target_id, []
                    ).append((exp_ind_id, target_value))
            else:
                self._create_period_start_event(
                    experiment=experiment,
                    i_period=i_period,
                    period=period,
                )

        if len(kept_periods) > 2:
            raise AssertionError("Expected at most two periods to be kept.")

        # add conditions that set the indicator parameters
        for period in kept_periods:
            period.condition_ids = [
                self._get_experiment_indicator_condition_id(experiment.id),
                self.CONDITION_ID_PREEQ_ON
                if period.is_preequilibration
                else self.CONDITION_ID_PREEQ_OFF,
            ]

        experiment.periods = kept_periods

    def _check_first_period_change(
        self, change: v2.Change, period: ExperimentPeriod
    ) -> sp.Basic:
        """Check that the target value of a first-period change can be
        evaluated before model initialization.

        :returns: The target value with `time` replaced by the period's
            start time.
        """
        model = self._new_problem.model.model
        allowed = set(self._new_problem.x_ids) | {
            p.name for p in model.parameters
        }
        time_syms = {
            s
            for s in change.target_value.free_symbols
            if s.name == v2.C.TIME_SYMBOL
        }
        if (
            invalid := {str(s) for s in change.target_value.free_symbols}
            - allowed
            - {v2.C.TIME_SYMBOL}
        ):
            raise ValueError(
                "Target values of changes in the first period must only "
                f"refer to parameters or `{v2.C.TIME_SYMBOL}`. Got "
                f"{str(change.target_value)!r} for target "
                f"{change.target_id!r}, which refers to {invalid}."
            )
        if not time_syms:
            return change.target_value
        if not np.isfinite(period.time):
            raise NotImplementedError(
                f"`{v2.C.TIME_SYMBOL}` in target values of "
                "pre-equilibration changes is "
                f"not supported. Got {str(change.target_value)!r} for "
                f"target {change.target_id!r}."
            )
        return change.target_value.xreplace(
            {s: sp.Float(period.time) for s in time_syms}
        )

    def _get_target_species_idx(self, target_id: str) -> int:
        """Get the index of the species that a condition target refers to.

        :returns:
            The index of the species in ``model.species``.
        """
        import pysb
        from petab.v2.models.pysb_model import pattern_from_string

        if (ix := self._species_indices.get(target_id)) is not None:
            return ix

        model = self._new_problem.model.model
        model_id = self.map_petab_to_pysb.get(target_id, target_id)

        if (component := model.components.get(model_id)) is not None:
            if isinstance(component, pysb.Parameter):
                raise NotImplementedError(
                    f"Changes of pysb.Parameters ({target_id!r}) are not "
                    "supported yet."
                )
            # Expressions, observables, ... are algebraic entities
            raise ValueError(
                f"Invalid condition target {target_id!r}: "
                f"{type(component).__name__}s cannot be targets of "
                "condition changes. Only species (via the mapping table) "
                "can be changed."
            )

        try:
            pattern = pattern_from_string(model_id, model)
        except Exception as e:
            raise ValueError(
                f"Invalid condition target {target_id!r} ({model_id!r}). "
                "Expected a fully specified species."
            ) from e
        if not pattern.is_concrete():
            raise ValueError(
                f"Invalid condition target {target_id!r} ({model_id!r}). "
                "Species targets must be fully specified (all sites and "
                "compartments)."
            )
        for ix, species in enumerate(model.species):
            if pysb.pattern.match_complex_pattern(
                pattern, species, exact=True
            ):
                self._species_indices[target_id] = ix
                return ix
        raise ValueError(
            f"Species {model_id!r} (condition target {target_id!r}) is not "
            "part of the generated reaction network."
        )

    def _apply_first_period_changes(
        self, target_id: str, changes: list[tuple[str, sp.Basic]]
    ) -> None:
        """Apply the first-period changes of all experiments to the initial
        value of the given target.

        :param target_id: The target ID of the changes.
        :param changes: ``[(experiment_indicator, target_value), ...]``
        """
        import pysb

        model = self._new_problem.model.model
        species_idx = self._get_target_species_idx(target_id)

        # update the matching pysb.Initial or create a new one
        species = model.species[species_idx]
        initial = next(
            (
                initial
                for initial in model.initials
                if pysb.pattern.match_complex_pattern(
                    initial.pattern, species, exact=True
                )
            ),
            None,
        )
        default = initial.value if initial is not None else sp.Float(0.0)
        expr = self._piecewise_by_experiment(changes, default)
        if expr == default:
            return
        # Initial.value needs to be a parameter or an expression
        expr_expr = pysb.Expression(
            f"_petab_initial_{target_id}", expr, _export=False
        )
        model.add_component(expr_expr)
        if initial is not None:
            initial.value = expr_expr
        else:
            # pysb.Initial only accepts constant expressions, which excludes
            #  piecewise functions. Set the actual value afterwards.
            expr_expr.expr = sp.Float(0.0)
            model.add_initial(pysb.Initial(species, expr_expr, _export=False))
            expr_expr.expr = expr

    def _piecewise_by_experiment(
        self, changes: list[tuple[str, sp.Basic]], default: sp.Basic
    ) -> sp.Basic:
        """Combine the first-period changes of different experiments into
        a single expression, conditional on the experiment indicators."""
        parameters = self._new_problem.model.model.parameters
        if expr_cond_pairs := [
            (target_value, parameters[exp_ind] > 0.5)
            for exp_ind, target_value in changes
            if target_value != default
        ]:
            return sp.Piecewise(*expr_cond_pairs, (default, True))
        return default

    def _create_period_start_event(
        self,
        experiment: v2.Experiment,
        i_period: int,
        period: ExperimentPeriod,
    ) -> None:
        """Create an event that triggers at the start of a period."""
        import pysb

        if not (changes := self._new_problem.get_changes_for_period(period)):
            return

        exp_ind_sym = symbol_with_assumptions(
            self.get_experiment_indicator(experiment.id)
        )
        preeq_ind_sym = symbol_with_assumptions(self.PREEQ_INDICATOR)

        # Create trigger expressions
        # Since handling of == and !=, and distinguishing < and <=
        # (and > and >=), is a bit tricky in terms of root-finding,
        # we use these slightly more convoluted expressions.
        # (assuming that the indicator parameters are {0, 1})
        if period.is_preequilibration:
            root_fun = sp.Min(exp_ind_sym - 0.5, preeq_ind_sym - 0.5)
        else:
            root_fun = sp.Min(
                exp_ind_sym - 0.5,
                0.5 - preeq_ind_sym,
                amici_time_symbol - period.time,
            )

        event_id = f"_petab_event_{experiment.id}_{i_period}"
        model = self._new_problem.model.model
        component_names = {c.name for c in model.components}
        # symbols in target values: species (via the mapping table),
        #  model components, `time`, or parameters that are only defined in
        #  the parameter table
        # PEtab ID -> AMICI symbol
        symbols: dict[str, sp.Symbol] = {}
        for petab_id in sorted(
            {
                *(change.target_id for change in changes),
                *(
                    str(sym)
                    for change in changes
                    for sym in change.target_value.free_symbols
                ),
            }
        ):
            model_id = self.map_petab_to_pysb.get(petab_id, petab_id)
            if model_id in RESERVED_SYMBOLS:
                # The PySB import renames components with reserved names
                #  (`resolve_reserved_symbol_renames`), but not in events.
                #  The first period isn't affected: those changes end up in
                #  PySB initials and expressions, which are renamed.
                raise NotImplementedError(
                    f"Entities named {model_id!r} can't be used in "
                    "condition changes after the first period."
                )
            if model_id in component_names:
                component = model.components[model_id]
                if not isinstance(
                    component,
                    pysb.Parameter | pysb.Expression | pysb.Observable,
                ):
                    raise ValueError(
                        f"Invalid symbol {petab_id!r} in a condition change: "
                        f"{type(component).__name__}s can't be used in "
                        "target values."
                    )
                sym = symbol_with_assumptions(model_id)
            elif petab_id == v2.C.TIME_SYMBOL:
                sym = amici_time_symbol
            elif petab_id in self.map_petab_to_pysb:
                sym = symbol_with_assumptions(
                    f"__s{self._get_target_species_idx(petab_id)}"
                )
            else:
                self._add_parameter(model_id, 0)
                component_names.add(model_id)
                sym = symbol_with_assumptions(model_id)
            symbols[petab_id] = sym

        assignments: dict[sp.Symbol, sp.Expr] = {}
        for change in changes:
            # validates the target type
            self._get_target_species_idx(change.target_id)
            # match by name: depending on how the target value was created,
            #  its symbols may have different assumptions
            assignments[symbols[change.target_id]] = (
                change.target_value.xreplace(
                    {
                        sym: symbols[str(sym)]
                        for sym in change.target_value.free_symbols
                    }
                )
            )

        event = Event(
            symbol=symbol_with_assumptions(event_id),
            name=event_id,
            value=root_fun,
            assignments=assignments,
            initial_value=False,
            use_values_from_trigger_time=False,
        )

        self._events.append(event)

    def _add_parameter(self, par_id: str, value: float) -> None:
        """Add a parameter to the PySB model."""
        import pysb

        p = pysb.Parameter(par_id, value, _export=False)
        self._new_problem.model.model.add_component(p)

    def _add_preequilibration_indicator(
        self,
    ) -> None:
        """Add an indicator parameter for the pre-equilibration to the SBML
        model."""
        par_id = self.PREEQ_INDICATOR
        if par_id in map(str, self._new_problem.model.model.components):
            raise ValueError(
                f"Entity with ID {par_id} already exists in the model."
            )

        # add the pre-steady-state indicator parameter
        self._add_parameter(par_id, 0.0)

    @staticmethod
    def get_experiment_indicator(experiment_id: str) -> str:
        """The ID of the experiment indicator parameter.

        The experiment indicator parameter is used to identify the
        experiment in the SBML model. It is a parameter that is set
        to 1 for the current experiment and 0 for all other
        experiments. The parameter is used in the event trigger
        to determine whether the event should be triggered.

        :param experiment_id: The ID of the experiment for which to create
            the experiment indicator parameter ID.
        """
        return f"_petab_experiment_indicator_{experiment_id}"

    def _add_indicators_to_conditions(self) -> None:
        """After converting the experiments to events, add the indicator
        parameters for the pre-equilibration period and for the different
        experiments to the remaining conditions.
        Then remove all other conditions."""
        from petab.v2 import Change, Condition, ConditionTable

        problem = self._new_problem

        # create conditions for indicator parameters
        problem += Condition(
            id=self.CONDITION_ID_PREEQ_ON,
            changes=[Change(target_id=self.PREEQ_INDICATOR, target_value=1)],
        )

        problem += Condition(
            id=self.CONDITION_ID_PREEQ_OFF,
            changes=[Change(target_id=self.PREEQ_INDICATOR, target_value=0)],
        )

        # add conditions for the experiment indicators
        for experiment in problem.experiments:
            cond_id = self._get_experiment_indicator_condition_id(
                experiment.id
            )
            changes = [
                Change(
                    target_id=self.get_experiment_indicator(experiment.id),
                    target_value=1,
                )
            ]
            problem += Condition(
                id=cond_id,
                changes=changes,
            )

        #  All changes have been encoded in event assignments and can be
        #  removed. Only keep the conditions setting our indicators.
        problem.condition_tables = [
            ConditionTable(
                [
                    condition
                    for condition in problem.conditions
                    if condition.id.startswith("_petab")
                ]
            )
        ]
