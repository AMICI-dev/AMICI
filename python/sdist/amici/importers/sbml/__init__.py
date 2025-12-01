"""
SBML Import
-----------
This module provides all necessary functionality to import a model specified
in the `Systems Biology Markup Language (SBML) <https://sbml.org/>`_.
"""

from __future__ import annotations

import contextlib
import copy
import itertools as itt
import logging
import math
import os
import re
import warnings
import xml.etree.ElementTree as ET
from collections.abc import Callable, Iterable, Sequence
from itertools import chain
from pathlib import Path
from typing import (
    Any,
)

import libsbml
import numpy as np
import sympy as sp
from sbmlmath import SBMLMathMLParser, TimeSymbol, avogadro
from sympy.logic.boolalg import Boolean, BooleanFalse, BooleanTrue
from sympy.matrices.dense import MutableDenseMatrix

import amici
from amici import get_model_dir
from amici._symbolic.de_model import DEModel
from amici._symbolic.de_model_components import (
    DifferentialState,
    Event,
    Expression,
    symbol_to_type,
)
from amici._symbolic.sympy_utils import (
    _monkeypatch_sympy,
    _piecewise_to_minmax,
    smart_is_zero_matrix,
    smart_multiply,
)
from amici.constants import SymbolId
from amici.importers.sbml.splines import AbstractSpline
from amici.importers.sbml.utils import SBMLException
from amici.importers.utils import (
    RESERVED_SYMBOLS,
    MeasurementChannel,
    _check_unsupported_functions,
    _default_simplify,
    _eq_to_and,
    _get_str_symbol_identifiers,
    _ne_to_or,
    _parse_piecewise_to_heaviside,
    _xor_to_or,
    amici_time_symbol,
    annotation_namespace,
    generate_flux_symbol,
    generate_measurement_symbol,
    generate_regularization_symbol,
    noise_distribution_to_cost_function,
    noise_distribution_to_observable_transformation,
    sbml_time_symbol,
    smart_subs,
    smart_subs_dict,
    symbol_with_assumptions,
    toposort_symbols,
)
from amici.logging import get_logger, log_execution_time, set_log_level

__all__ = [
    "SbmlImporter",
    "SBMLException",
    "MeasurementChannel",
    "assignment_rules_to_observables",
]

SymbolicFormula = dict[sp.Symbol, sp.Expr]

default_symbols = {symbol: {} for symbol in SymbolId}

ConservationLaw = dict[str, str | sp.Expr]

logger = get_logger(__name__, logging.ERROR)


class SbmlImporter:
    """
    Importer for models provided in the Systems Biology Markup Language (SBML).
    """

    def __init__(
        self,
        sbml_source: str | Path | libsbml.Model,
        show_sbml_warnings: bool = False,
        from_file: bool = True,
        discard_annotations: bool = False,
        jax: bool = False,
    ) -> None:
        """
        Initialize.

        :param sbml_source:
            Either a path to SBML file where the model is specified,
            or a model string as created by sbml.sbmlWriter(
            ).writeSBMLToString() or an instance of `libsbml.Model`.

        :param show_sbml_warnings:
            Indicates whether libSBML warnings should be displayed.

        :param from_file:
            Whether `sbml_source` is a file name (True, default), or an SBML
            string

        :param discard_annotations:
            discard information contained in AMICI SBML (spline) annotations
            (debug).
        """
        # The libSBML sbml reader.
        # Warning: Not storing this may result in a segfault.
        self._sbml_reader: libsbml.SBMLReader | None = None

        self._sbml_doc: libsbml.Document = None
        if isinstance(sbml_source, libsbml.Model):
            self._sbml_doc: libsbml.Document = sbml_source.getSBMLDocument()
        else:
            self._sbml_reader: libsbml.SBMLReader = libsbml.SBMLReader()
            if from_file:
                sbml_doc = self._sbml_reader.readSBMLFromFile(str(sbml_source))
            else:
                sbml_doc = self._sbml_reader.readSBMLFromString(sbml_source)
            self._sbml_doc = sbml_doc

        # whether to show SBML warnings
        self.show_sbml_warnings: bool = show_sbml_warnings

        # process document
        self._process_document()
        #: the SBML model to import
        self.sbml_model: libsbml.Model = self._sbml_doc.getModel()

        # Long and short names for model components
        self._symbols: dict[SymbolId, dict[sp.Symbol, dict[str, Any]]] = {}

        # model symbols for sympy to consider during sympification
        # see `locals`argument in `sympy.sympify`
        self._local_symbols: dict[str, sp.Expr | sp.Function] = {}
        # dict of compartment ids and compartment volumes
        self._compartments: SymbolicFormula = {}
        # Assignment rules for compartments.
        # Key is _symbolic identifier and value is assignment value
        self._compartment_assignment_rules: SymbolicFormula = {}
        # Assignment rules for species.
        # Key is _symbolic identifier and value is assignment value
        self._species_assignment_rules: SymbolicFormula = {}
        # assignment rules for parameters, these parameters are not permissible
        # for sensitivity analysis
        self._parameter_assignment_rules: SymbolicFormula = {}
        # initial assignments for parameters, these parameters are not
        # permissible for sensitivity analysis
        self._initial_assignments: SymbolicFormula = {}
        self._splines: list[AbstractSpline] = []
        # stoichiometric matrix of the model
        self._stoichiometric_matrix: sp.Matrix | None = None
        # reaction kinetic laws
        self._flux_vector = None
        # identifiers for elements of flux_vector
        self._flux_ids = None
        self._reset_symbols()

        # sets behavior of SBML Formula parsing
        # https://sbml.org/software/libsbml/5.18.0/docs/formatted/python-api/classlibsbml_1_1_l3_parser_settings.html#ab30d7ed52ca24cbb842d0a7fed7f4bfd
        # all defaults except disable unit parsing
        self._sbml_parser_settings = libsbml.L3ParserSettings()
        self._sbml_parser_settings.setModel(self.sbml_model)
        self._sbml_parser_settings.setParseUnits(libsbml.L3P_NO_UNITS)

        self._discard_annotations: bool = discard_annotations
        self._mathml_parser = SBMLMathMLParser(
            sbml_level=self.sbml_model.getLevel(),
            sbml_version=self.sbml_model.getVersion(),
            symbol_kwargs={"real": True},
            ignore_units=True,
            evaluate=True,
        )
        self.jax = jax

    @log_execution_time("loading SBML", logger)
    def _process_document(self) -> None:
        """
        Validate and simplify document.
        """
        # Ensure we got a valid SBML model, otherwise further processing
        # might lead to undefined results
        log_execution_time("validating SBML", logger)(
            self._sbml_doc.validateSBML
        )()
        _check_lib_sbml_errors(self._sbml_doc, self.show_sbml_warnings)

        # Flatten "comp" model? Do that before any other converters are run
        if any(
            self._sbml_doc.getPlugin(i_plugin).getPackageName() == "comp"
            for i_plugin in range(self._sbml_doc.getNumPlugins())
        ):
            # see libsbml CompFlatteningConverter for options
            conversion_properties = libsbml.ConversionProperties()
            conversion_properties.addOption("flatten comp", True)
            conversion_properties.addOption("leave_ports", False)
            conversion_properties.addOption("performValidation", False)
            conversion_properties.addOption("abortIfUnflattenable", "none")
            if (
                log_execution_time("flattening hierarchical SBML", logger)(
                    self._sbml_doc.convert
                )(conversion_properties)
                != libsbml.LIBSBML_OPERATION_SUCCESS
            ):
                raise SBMLException(
                    "Required SBML comp extension is currently not supported "
                    "and flattening the model failed."
                )
                # check the flattened model is still valid
            log_execution_time("re-validating SBML", logger)(
                self._sbml_doc.validateSBML
            )()
            _check_lib_sbml_errors(self._sbml_doc, self.show_sbml_warnings)

        # apply several model simplifications that make our life substantially
        # easier
        if self._sbml_doc.getModel().getNumFunctionDefinitions():
            convert_config = libsbml.SBMLFunctionDefinitionConverter().getDefaultProperties()
            log_execution_time("converting SBML functions", logger)(
                self._sbml_doc.convert
            )(convert_config)

        convert_config = (
            libsbml.SBMLLocalParameterConverter().getDefaultProperties()
        )
        log_execution_time("converting SBML local parameters", logger)(
            self._sbml_doc.convert
        )(convert_config)

        # If any of the above calls produces an error, this will be added to
        # the SBMLError log in the sbml document. Thus, it is sufficient to
        # check the error log just once after all conversion/validation calls.
        _check_lib_sbml_errors(self._sbml_doc, self.show_sbml_warnings)

        # need to reload the converted model
        self.sbml_model = self._sbml_doc.getModel()

    def _reset_symbols(self) -> None:
        """
        Reset the symbols attribute to default values
        """
        self._symbols = copy.deepcopy(default_symbols)
        self._local_symbols = {}

    def sbml2amici(
        self,
        model_name: str,
        output_dir: str | Path = None,
        fixed_parameters: Iterable[str] = None,
        observation_model: list[MeasurementChannel] = None,
        verbose: int | bool = logging.ERROR,
        assume_pow_positivity: bool = False,
        compiler: str = None,
        allow_reinit_fixpar_initcond: bool = True,
        compile: bool = True,
        compute_conservation_laws: bool = True,
        simplify: Callable | None = _default_simplify,
        cache_simplify: bool = False,
        generate_sensitivity_code: bool = True,
        hardcode_symbols: Sequence[str] = None,
    ) -> amici.sim.sundials.Model | None:
        """
        Generate and compile AMICI C++ files for the model provided to the
        constructor.

        The resulting model can be imported as a regular Python module (if
        `compile=True`), or used from C++ as described in the documentation of
        AMICI's C++ interface.

        Note that this generates model ODEs for changes in concentrations, not
        amounts unless the `hasOnlySubstanceUnits` attribute has been
        defined in the SBML model for a particular species.

        Sensitivity analysis for local parameters is enabled by creating
        global parameters ``_{reactionId}_{localParameterName}``.

        :param model_name:
            Name of the generated model package.
            Note that in a given Python session, only one model with a given
            name can be loaded at a time.
            The generated Python extensions cannot be unloaded. Therefore,
            make sure to choose a unique name for each model.

        :param output_dir:
            Directory where the generated model package will be stored.
            Defaults to :func:`amici.get_model_dir`.

        :param fixed_parameters:
            SBML Ids to be excluded from sensitivity analysis

        :param observation_model:
            The different measurement channels that make up the observation
            model, see :class:`MeasurementChannel`.
            If ``None``, default observables will be added (for all
            state variables of the model and all species, compartments,
            and assignment rule targets) and normally distributed
            measurement noise is assumed.
            If ``[]``, no observables will be added.
        :param verbose:
            Verbosity level for logging, ``True``/``False`` defaults to
            ``logging.Error``/``logging.DEBUG``.

        :param assume_pow_positivity:
            if set to ``True``, a special pow function is
            used to avoid problems with state variables that may become
            negative due to numerical errors

        :param compiler:
            Absolute path to the compiler executable to be used to build the Python
            extension, e.g. ``/usr/bin/clang``.

        :param allow_reinit_fixpar_initcond:
            See :class:`amici.de_export.DEExporter`.

        :param compile:
            If ``True``, compile the generated Python package,
            if ``False``, just generate code.

        :param compute_conservation_laws:
            if set to ``True``, conservation laws are automatically computed
            and applied such that the state-jacobian of the ODE
            right-hand-side has full rank. This option should be set to
            ``True`` when using the Newton algorithm to compute steadystate
            sensitivities.
            Conservation laws for constant species are enabled by default.
            Support for conservation laws for non-constant species is
            experimental and may be enabled by setting an environment variable
            ``AMICI_EXPERIMENTAL_SBML_NONCONST_CLS`` to either ``demartino``
            to use the algorithm proposed by De Martino et al. (2014)
            https://doi.org/10.1371/journal.pone.0100750, or to any other value
            to use the deterministic algorithm implemented in
            ``conserved_moieties2.py``. In some cases, the ``demartino`` may
            run for a very long time. This has been observed for example in the
            case of stoichiometric coefficients with many significant digits.

        :param simplify:
            See :attr:`amici.DEModel._simplify`.

        :param cache_simplify:
            See :meth:`amici.DEModel.__init__`.

        :param generate_sensitivity_code:
            If ``False``, the code required for sensitivity computation will
            not be generated.

        :param hardcode_symbols:
            List of SBML entity IDs that are to be hardcoded in the generated model.
            Their values cannot be changed anymore after model import.
            Currently, only parameters that are not targets of rules or
            initial assignments are supported.

        :return:
            If `compile` is `True` and compilation was successful, an instance
            of the generated model class, otherwise `None`.
        """
        set_log_level(logger, verbose)

        ode_model = self._build_ode_model(
            fixed_parameters=fixed_parameters,
            observation_model=observation_model,
            verbose=verbose,
            compute_conservation_laws=compute_conservation_laws,
            simplify=simplify,
            cache_simplify=cache_simplify,
            hardcode_symbols=hardcode_symbols,
        )

        output_dir = output_dir or get_model_dir(model_name)

        from amici.exporters.sundials.de_export import (
            DEExporter,
        )

        try:
            exporter = DEExporter(
                ode_model,
                model_name=model_name,
                outdir=output_dir,
                verbose=verbose,
                assume_pow_positivity=assume_pow_positivity,
                compiler=compiler,
                allow_reinit_fixpar_initcond=allow_reinit_fixpar_initcond,
                generate_sensitivity_code=generate_sensitivity_code,
            )
            exporter.generate_model_code()
        except NotImplementedError as e:
            raise SBMLException(
                "The SBML model contains features that are currently not "
                "supported by the AMICI SBML importer."
            ) from e

        if compile:
            from amici.sim.sundials import has_clibs

            if not has_clibs:
                warnings.warn(
                    "AMICI C++ extensions have not been built. "
                    "Generated model code, but unable to compile.",
                    stacklevel=2,
                )
            exporter.compile_model()

            from amici import import_model_module

            return import_model_module(
                module_name=model_name, module_path=output_dir
            ).get_model()

        return None

    def sbml2jax(
        self,
        model_name: str,
        output_dir: str | Path = None,
        observation_model: list[MeasurementChannel] = None,
        verbose: int | bool = logging.ERROR,
        compute_conservation_laws: bool = True,
        simplify: Callable | None = _default_simplify,
        cache_simplify: bool = False,
        hybridization: dict = None,
    ) -> None:
        """
        Generate and compile AMICI jax files for the model provided to the
        constructor.

        The resulting model can be imported as a regular Python module.

        Note that this generates model ODEs for changes in concentrations, not
        amounts unless the `hasOnlySubstanceUnits` attribute has been
        defined for a particular species.

        :param model_name:
            Name of the generated model package.
            Note that in a given Python session, only one model with a given
            name can be loaded at a time.
            The generated Python extensions cannot be unloaded. Therefore,
            make sure to choose a unique name for each model.

        :param output_dir:
            Directory where the generated model package will be stored.

        :param observation_model:
            The different measurement channels that make up the observation
            model, see :class:`amici.importers.utils.MeasurementChannel`.
            Only time-resolved observables are supported here.
            If ``None``, default observables will be added (for all
            state variables of the model and all species, compartments,
            and assignment rule targets) and normally distributed
            measurement noise is assumed.
            If ``[]``, no observables will be added.

        :param verbose:
            verbosity level for logging, ``True``/``False`` default to
            ``logging.Error``/``logging.DEBUG``

        :param compute_conservation_laws:
            if set to ``True``, conservation laws are automatically computed
            and applied such that the state-jacobian of the ODE
            right-hand-side has full rank. This option should be set to
            ``True`` when using the Newton algorithm to compute steadystate
            sensitivities.
            Conservation laws for constant species are enabled by default.
            Support for conservation laws for non-constant species is
            experimental and may be enabled by setting an environment variable
            ``AMICI_EXPERIMENTAL_SBML_NONCONST_CLS`` to either ``demartino``
            to use the algorithm proposed by De Martino et al. (2014)
            https://doi.org/10.1371/journal.pone.0100750, or to any other value
            to use the deterministic algorithm implemented in
            ``conserved_moieties2.py``. In some cases, the ``demartino`` may
            run for a very long time. This has been observed for example in the
            case of stoichiometric coefficients with many significant digits.

        :param simplify:
            see :attr:`amici.DEModel._simplify`

        :param cache_simplify:
            see :meth:`amici.DEModel.__init__`

        :param hybridization:
            dict representation of the hybridization information in the PEtab YAML file, see
            https://petab-sciml.readthedocs.io/latest/format.html#problem-yaml-file
        """
        set_log_level(logger, verbose)

        ode_model = self._build_ode_model(
            observation_model=observation_model,
            verbose=verbose,
            compute_conservation_laws=compute_conservation_laws,
            simplify=simplify,
            cache_simplify=cache_simplify,
            hybridization=hybridization,
        )

        from amici.jax.ode_export import ODEExporter

        exporter = ODEExporter(
            ode_model,
            model_name=model_name,
            outdir=output_dir,
            verbose=verbose,
            hybridization=hybridization,
        )
        exporter.generate_model_code()

    @_monkeypatch_sympy
    def _build_ode_model(
        self,
        fixed_parameters: Iterable[str] = None,
        observation_model: list[MeasurementChannel] = None,
        verbose: int | bool = logging.ERROR,
        compute_conservation_laws: bool = True,
        simplify: Callable | None = _default_simplify,
        cache_simplify: bool = False,
        hardcode_symbols: Sequence[str] = None,
        hybridization: dict = None,
    ) -> DEModel:
        """Generate a DEModel from this SBML model.

        See :py:func:`sbml2amici` for parameters.
        """
        fixed_parameters = list(fixed_parameters) if fixed_parameters else []

        hardcode_symbols = set(hardcode_symbols) if hardcode_symbols else {}
        if invalid := (set(fixed_parameters) & set(hardcode_symbols)):
            raise ValueError(
                "The following parameters were selected as both constant "
                f"and hard-coded which is not allowed: {invalid}"
            )

        self._reset_symbols()
        self._process_sbml(
            fixed_parameters=fixed_parameters,
            hardcode_symbols=hardcode_symbols,
        )

        if (
            self._symbols.get(SymbolId.EVENT, False)
            or any(
                x["value"].has(sp.Heaviside, sp.Piecewise)
                for x in self._symbols[SymbolId.EXPRESSION].values()
            )
            or self._flux_vector.has(sp.Heaviside, sp.Piecewise)
        ):
            if compute_conservation_laws:
                logger.warning(
                    "Conservation laws are currently not supported for models "
                    "with events, piecewise or Heaviside functions, "
                    "and will be turned off."
                )
            compute_conservation_laws = False

        self._process_observation_model(observation_model)
        self._replace_compartments_with_volumes()

        self._clean_reserved_symbols()
        self._process_time()

        ode_model = DEModel(
            verbose=verbose,
            simplify=simplify,
            cache_simplify=cache_simplify,
        )

        ode_model._has_quadratic_nllh = all(
            llh["dist"]
            in ["normal", "lin-normal", "log-normal", "log10-normal"]
            for llh in self._symbols[SymbolId.LLHY].values()
        )

        # add splines as expressions to the model
        # saved for later substituting into the fluxes
        spline_subs = {}
        for ispl, spl in enumerate(self._splines):
            spline_expr = spl.ode_model_symbol(self)
            spline_subs[spl.sbml_id] = spline_expr
            ode_model.add_spline(spl, spline_expr)

        # assemble fluxes and add them as expressions to the model
        assert len(self._flux_ids) == len(self._flux_vector)
        fluxes = [
            generate_flux_symbol(ir, name=flux_id)
            for ir, flux_id in enumerate(self._flux_ids)
        ]

        # create dynamics without respecting conservation laws first
        dxdt = (
            smart_multiply(
                self._stoichiometric_matrix, MutableDenseMatrix(fluxes)
            )
            if len(self._stoichiometric_matrix)
            else sp.zeros(self._stoichiometric_matrix.shape[0], 1)
        )
        # dxdt has algebraic states at the end
        assert dxdt.shape[0] - len(self._symbols[SymbolId.SPECIES]) == len(
            self._symbols.get(SymbolId.ALGEBRAIC_STATE, [])
        ), (
            self._symbols.get(SymbolId.SPECIES),
            dxdt,
            self._symbols.get(SymbolId.ALGEBRAIC_STATE),
        )

        # correct time derivatives for compartment changes
        for ix, ((species_id, species), formula) in enumerate(
            zip(self._symbols[SymbolId.SPECIES].items(), dxdt, strict=False)
        ):
            # rate rules and amount species don't need to be updated
            if "dt" in species:
                continue
            if species["amount"]:
                species["dt"] = formula
            else:
                species["dt"] = self._transform_dxdt_to_concentration(
                    species_id, formula
                )

        # create all basic components of the DE model and add them.
        for symbol_name in self._symbols:
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
                args += ["transformation"]
            elif symbol_name == SymbolId.EVENT_OBSERVABLE:
                args += ["event"]

            comp_kwargs = [
                {
                    "symbol": var_id,
                    **{k: v for k, v in var.items() if k in args},
                }
                for var_id, var in self._symbols[symbol_name].items()
            ]

            for comp_kwarg in comp_kwargs:
                ode_model.add_component(
                    symbol_to_type[symbol_name](**comp_kwarg)
                )

        # add fluxes as expressions, this needs to happen after base
        # expressions from symbols have been parsed
        for flux_id, flux in zip(fluxes, self._flux_vector, strict=True):
            # replace splines inside fluxes
            flux = flux.subs(spline_subs)
            ode_model.add_component(
                Expression(symbol=flux_id, name=str(flux_id), value=flux)
            )

        if compute_conservation_laws:
            self._process_conservation_laws(ode_model)

        if hybridization:
            ode_model._process_hybridization(hybridization)

        # substitute SBML-rateOf constructs
        self._process_sbml_rate_of(ode_model)

        # fill in 'self._sym' based on prototypes and components in ode_model
        ode_model.generate_basic_variables()

        return ode_model

    @log_execution_time("importing SBML", logger)
    def _process_sbml(
        self,
        fixed_parameters: list[str] = None,
        hardcode_symbols: Sequence[str] = None,
    ) -> None:
        """
        Read parameters, species, reactions, and so on from SBML model

        :param fixed_parameters:
            SBML Ids to be excluded from sensitivity analysis
        :param hardcode_symbols:
            Parameter IDs to be replaced by their values in the generated model.
        """
        if not self._discard_annotations:
            self._process_annotations()
        self._check_support()
        self._gather_locals(hardcode_symbols=hardcode_symbols)
        self._process_parameters(
            fixed_parameter_ids=fixed_parameters,
            hardcode_symbols=hardcode_symbols,
        )
        self._process_compartments()
        self._process_species()
        self._process_reactions()
        self._process_rules()
        self._process_events()
        self._process_initial_assignments()
        self._process_species_references()

    def _check_support(self) -> None:
        """
        Check whether all required SBML features are supported.
        Also ensures that the SBML contains at least one reaction, or rate
        rule, or assignment rule, to produce change in the system over time.
        """

        # Check for required but unsupported SBML extensions
        if (
            self._sbml_doc.getLevel() != 3
            and hasattr(self.sbml_model, "all_elements_from_plugins")
            and self.sbml_model.all_elements_from_plugins.getSize()
        ):
            raise SBMLException("SBML extensions are currently not supported!")

        if self._sbml_doc.getLevel() == 3:
            # the "required" attribute is only available in SBML Level 3
            for i_plugin in range(self.sbml_model.getNumPlugins()):
                plugin = self.sbml_model.getPlugin(i_plugin)
                if (
                    self._sbml_doc.getPkgRequired(plugin.getPackageName())
                    is False
                ):
                    # if not "required", this has no impact on model
                    #  simulation, and we can safely ignore it

                    if (
                        plugin.getPackageName() == "fbc"
                        and plugin.getListOfAllElements()
                    ):
                        # fbc is labeled not-required, but in fact it is.
                        # we don't care about the extra attributes of core
                        # elements, such as fbc:chemicalFormula, but we can't
                        # do anything meaningful with fbc:objective or
                        # fbc:fluxBounds
                        raise SBMLException(
                            "The following fbc extension elements are "
                            "currently not supported: "
                            + ", ".join(
                                list(map(str, plugin.getListOfAllElements()))
                            )
                        )

                    continue

                # Check if there are extension elements. If not, we can safely
                #  ignore the enabled package
                if plugin.getListOfAllElements():
                    raise SBMLException(
                        f"Required SBML extension {plugin.getPackageName()} "
                        f"is currently not supported!"
                    )

        if any(
            rule.isRate()
            and not isinstance(
                self.sbml_model.getElementBySId(rule.getVariable()),
                libsbml.Compartment | libsbml.Species | libsbml.Parameter,
            )
            for rule in self.sbml_model.getListOfRules()
        ):
            raise SBMLException(
                "Rate rules are only supported for "
                "species, compartments, and parameters."
            )

        if any(r.getFast() for r in self.sbml_model.getListOfReactions()):
            raise SBMLException("Fast reactions are currently not supported!")

        # Check events for unsupported functionality
        self._check_event_support()

    def _check_event_support(self) -> None:
        """
        Check possible events in the model, as AMICI does currently not support

        * delays in events
        * priorities of events
        * events fired at initial time

        Furthermore, event triggers are optional (e.g., if an event is fired at
        initial time, no trigger function is necessary).
        In this case, warn that this event will have no effect.
        """
        for event in self.sbml_model.getListOfEvents():
            event_id = event.getId()
            # Check for delays in events
            delay = event.getDelay()
            if delay is not None:
                try:
                    delay_time = float(self._sympify(delay))
                    if delay_time != 0:
                        raise ValueError
                # `TypeError` would be raised in the above `float(...)`
                # if the delay is not a fixed time
                except (TypeError, ValueError):
                    raise SBMLException(
                        "Events with execution delays are "
                        "currently not supported in AMICI."
                    )
            # Check for priorities
            if (prio := event.getPriority()) is not None:
                if (prio := self._sympify(prio)) and not prio.is_Number:
                    # Computing sensitivities with respect to event priorities
                    #  is not implemented, so let's import such models at all.
                    #  We could support expressions that only depend on
                    #  constant parameters, though. But who needs that anyway?
                    raise SBMLException(
                        f"Event {event_id} has a non-numeric priority "
                        "specified. This is currently not "
                        "supported in AMICI."
                    )

            # check trigger
            trigger_sbml = event.getTrigger()
            if trigger_sbml is None:
                logger.warning(
                    f"Event {event_id} trigger has no trigger, "
                    "so will be skipped."
                )
                continue
            if trigger_sbml.getMath() is None:
                logger.warning(
                    f"Event {event_id} trigger has no trigger "
                    "expression, so a dummy trigger will be set."
                )

            if not trigger_sbml.getPersistent():
                raise SBMLException(
                    f"Event {event_id} has a non-persistent trigger."
                    "This is currently not supported in AMICI."
                )

    @log_execution_time("gathering local SBML symbols", logger)
    def _gather_locals(self, hardcode_symbols: Sequence[str] = None) -> None:
        """
        Populate self.local_symbols with all model entities.

        This is later used during sympifications to avoid sympy builtins
        shadowing model entities as well as to avoid possibly costly
        _symbolic substitutions
        """
        self._gather_base_locals(hardcode_symbols=hardcode_symbols)
        self._gather_dependent_locals()

    def _gather_base_locals(
        self, hardcode_symbols: Sequence[str] = None
    ) -> None:
        """
        Populate self.local_symbols with pure symbol definitions that do not
        depend on any other symbol.
        """

        special_symbols_and_funs = {
            # oo is sympy infinity
            "INF": sp.oo,
            "NaN": sp.nan,
            "rem": sp.Mod,
            "time": sbml_time_symbol,
            # SBML L3 explicitly defines this value, which is not equal
            # to the most recent SI definition.
            "avogadro": sp.Float(6.02214179e23),
            "exponentiale": sp.E,
            "log10": lambda x: sp.log(x, 10),
        }
        for s, v in special_symbols_and_funs.items():
            self._add_local_symbol(s, v)

        for c in itt.chain(
            self.sbml_model.getListOfSpecies(),
            self.sbml_model.getListOfParameters(),
            self.sbml_model.getListOfCompartments(),
        ):
            if not c.isSetId():
                continue
            if c.getId() in hardcode_symbols:
                if c.getConstant() is not True:
                    # disallow anything that can be changed by rules/reaction/events
                    raise ValueError(
                        f"Cannot hardcode non-constant symbol `{c.getId()}`."
                    )
                if self.sbml_model.getInitialAssignment(c.getId()):
                    raise NotImplementedError(
                        f"Cannot hardcode symbol `{c.getId()}` "
                        "that is an initial assignment target."
                    )
                self._add_local_symbol(c.getId(), sp.Float(c.getValue()))
            else:
                self._add_local_symbol(c.getId(), _get_identifier_symbol(c))

        for x_ref in _get_list_of_species_references(self.sbml_model):
            if not x_ref.isSetId():
                continue
            if (
                hasattr(x_ref, "isSetStoichiometry")
                and x_ref.isSetStoichiometry()
                and not self._is_assignment_rule_target(x_ref)
            ):
                value = sp.Float(x_ref.getStoichiometry())
            else:
                value = _get_identifier_symbol(x_ref)

            ia_sym = self._get_element_initial_assignment(x_ref.getId())
            if ia_sym is not None:
                value = ia_sym

            self._add_local_symbol(x_ref.getId(), value)

    def _gather_dependent_locals(self):
        """
        Populate self.local_symbols with symbol definitions that may depend on
        other symbol definitions.
        """
        for r in self.sbml_model.getListOfReactions():
            if not r.isSetId():
                continue
            self._add_local_symbol(
                r.getId(),
                self._sympify(r.getKineticLaw() or sp.Float(0)),
            )

    def _add_local_symbol(self, key: str, value: sp.Expr):
        """
        Add local symbols with some sanity checking for duplication which
        would indicate redefinition of internals, which SBML permits,
        but we don't.

        :param key:
            local symbol key

        :param value:
            local symbol value
        """
        if key in self._local_symbols.keys():
            raise SBMLException(
                f"AMICI tried to add a local symbol {key} with value {value}, "
                f"but {key} was already instantiated with "
                f"{self._local_symbols[key]}. This means that there "
                f"are multiple SBML elements with SId {key}, which is "
                f"invalid SBML. This can be fixed by renaming "
                f"the elements with SId {key}."
            )
        if key in {"True", "False", "true", "false", "pi"}:
            raise SBMLException(
                f"AMICI tried to add a local symbol {key} with value {value}, "
                f"but {key} is a reserved symbol in AMICI. This can be fixed "
                f"by renaming the element with SId {key}."
            )
        self._local_symbols[key] = value

    @log_execution_time("processing SBML compartments", logger)
    def _process_compartments(self) -> None:
        """
        Get compartment information, stoichiometric matrix and fluxes from
        SBML model.
        """
        compartments = self.sbml_model.getListOfCompartments()
        self._compartments = {}
        for comp in compartments:
            init = sp.Float(1.0)

            if comp.isSetVolume():
                init = self._sympify(comp.getVolume())

            ia_sym = self._get_element_initial_assignment(comp.getId())
            if ia_sym is not None:
                init = ia_sym

            self._compartments[_get_identifier_symbol(comp)] = init

    @log_execution_time("processing SBML species", logger)
    def _process_species(self) -> None:
        """
        Get species information from SBML model.
        """
        if self.sbml_model.isSetConversionFactor():
            conversion_factor = symbol_with_assumptions(
                self.sbml_model.getConversionFactor()
            )
        else:
            conversion_factor = 1

        for s in self.sbml_model.getListOfSpecies():
            if self._is_assignment_rule_target(s):
                continue
            self._symbols[SymbolId.SPECIES][_get_identifier_symbol(s)] = {
                "name": s.getName() if s.isSetName() else s.getId(),
                "compartment": _get_species_compartment_symbol(s),
                "constant": s.getConstant() or s.getBoundaryCondition(),
                "amount": s.getHasOnlySubstanceUnits(),
                "conversion_factor": symbol_with_assumptions(
                    s.getConversionFactor()
                )
                if s.isSetConversionFactor()
                else conversion_factor,
                "index": len(self._symbols[SymbolId.SPECIES]),
            }

        self._convert_event_assignment_parameter_targets_to_species()
        self._convert_event_assignment_compartment_targets_to_species()
        self._process_species_initial()
        self._process_rate_rules()

    @log_execution_time("processing SBML species initials", logger)
    def _process_species_initial(self):
        """
        Extract initial values and initial assignments from species
        """
        for species_variable in self.sbml_model.getListOfSpecies():
            initial = get_species_initial(species_variable)

            species_id = _get_identifier_symbol(species_variable)
            # If species_id is a target of an AssignmentRule, species will be
            # None, but we don't have to account for the initial definition
            # of the species itself and SBML doesn't permit AssignmentRule
            # targets to have InitialAssignments.
            species = self._symbols[SymbolId.SPECIES].get(species_id, None)

            ia_initial = self._get_element_initial_assignment(
                species_variable.getId()
            )
            if ia_initial is not None:
                initial = ia_initial
            if species:
                species["init"] = initial

        # hide rateOf-arguments from toposort and the substitution below
        all_rateof_dummies = []
        for species in self._symbols[SymbolId.SPECIES].values():
            species["init"], rateof_dummies = _rateof_to_dummy(species["init"])
            all_rateof_dummies.append(rateof_dummies)

        # don't assign this since they need to stay in order
        sorted_species = toposort_symbols(
            self._symbols[SymbolId.SPECIES], "init"
        )
        for species, rateof_dummies in zip(
            self._symbols[SymbolId.SPECIES].values(),
            all_rateof_dummies,
            strict=True,
        ):
            species["init"] = _dummy_to_rateof(
                smart_subs_dict(species["init"], sorted_species, "init"),
                rateof_dummies,
            )

    @log_execution_time("processing SBML rate rules", logger)
    def _process_rate_rules(self):
        """
        Process rate rules for species, compartments and parameters.
        Compartments and parameters with rate rules are implemented as species.
        Note that, in the case of species, rate rules may describe the change
        in amount, not concentration, of a species.
        """
        rules = self.sbml_model.getListOfRules()
        # compartments with rules are replaced with constants in the relevant
        # equations during the _replace_in_all_expressions call inside
        # _process_rules
        for rule in rules:
            if rule.getTypeCode() != libsbml.SBML_RATE_RULE:
                continue

            variable = symbol_with_assumptions(rule.getVariable())
            formula = self._sympify(rule)
            if formula is None:
                continue

            # Species rules are processed first, to avoid processing
            # compartments twice (as compartments with rate rules are
            # implemented as species).
            ia_init = self._get_element_initial_assignment(rule.getVariable())
            if variable in self._symbols[SymbolId.SPECIES]:
                init = self._symbols[SymbolId.SPECIES][variable]["init"]
                name = None

            if variable in self._compartments:
                init = self._compartments[variable]
                name = str(variable)
                del self._compartments[variable]

            elif variable in self._symbols[SymbolId.FREE_PARAMETER]:
                init = self._sympify(
                    self._symbols[SymbolId.FREE_PARAMETER][variable]["value"],
                )
                name = self._symbols[SymbolId.FREE_PARAMETER][variable]["name"]
                del self._symbols[SymbolId.FREE_PARAMETER][variable]

            # parameter with initial assignment, cannot use
            # self.initial_assignments as it is not filled at this
            # point
            elif ia_init is not None:
                init = ia_init
                par = self.sbml_model.getElementBySId(rule.getVariable())
                name = par.getName() if par.isSetName() else par.getId()

            self._add_d_dt(formula, variable, init, name)

    def _add_d_dt(
        self,
        d_dt: sp.Expr,
        variable: sp.Symbol,
        variable0: float | sp.Expr,
        name: str,
    ) -> None:
        """
        Creates or modifies species, to implement rate rules for
        compartments and species, respectively.

        :param d_dt:
            The rate rule (or, right-hand side of an ODE).

        :param variable:
            The subject of the rate rule.

        :param variable0:
            The initial value of the variable.

        :param name:
            Species name, only applicable if this function generates a new
            species
        """
        if variable in self._symbols[SymbolId.SPECIES]:
            # only update dt if species was already generated
            self._symbols[SymbolId.SPECIES][variable]["dt"] = d_dt
        else:
            # update initial values
            for species_id, species in self._symbols[SymbolId.SPECIES].items():
                variable0 = smart_subs(variable0, species_id, species["init"])

            for species in self._symbols[SymbolId.SPECIES].values():
                species["init"] = smart_subs(
                    species["init"], variable, variable0
                )

            # add compartment/parameter species
            self._symbols[SymbolId.SPECIES][variable] = {
                "name": name,
                "init": variable0,
                "amount": False,
                "conversion_factor": 1.0,
                "constant": False,
                "index": len(self._symbols[SymbolId.SPECIES]),
                "dt": d_dt,
            }

    @log_execution_time("processing SBML annotations", logger)
    def _process_annotations(self) -> None:
        """
        Process annotations that make modifications to the
        SBML model and thus have to be run before everything else
        """
        # Remove all parameters (and corresponding rules)
        # for which amici:discard is set
        parameter_ids_to_remove = []
        for p in self.sbml_model.getListOfParameters():
            annotation = p.getAnnotationString()
            assert isinstance(annotation, str)
            if len(annotation) != 0:
                annotation = ET.fromstring(annotation)
                for child in annotation:
                    if child.tag == f"{{{annotation_namespace}}}discard":
                        parameter_ids_to_remove.append(p.getIdAttribute())
        for parameter_id in parameter_ids_to_remove:
            # Remove corresponding rules
            self.sbml_model.removeRuleByVariable(parameter_id)
            # Remove parameter
            self.sbml_model.removeParameter(parameter_id)

    @log_execution_time("processing SBML parameters", logger)
    def _process_parameters(
        self,
        fixed_parameter_ids: list[str] | None = None,
        hardcode_symbols: Sequence[str] | None = None,
    ) -> None:
        """
        Get parameter information from SBML model.

        :param fixed_parameter_ids:
            SBML Ids to be excluded from sensitivity analysis

        :param hardcode_symbols:
            Parameter IDs to be replaced by their values in the generated model.
        """

        if fixed_parameter_ids is None:
            fixed_parameter_ids = []

        # Ensure specified constant parameters exist in the model
        for parameter in fixed_parameter_ids:
            if not self.sbml_model.getParameter(parameter):
                raise KeyError(
                    f"Cannot make {parameter} a constant parameter: "
                    "Parameter does not exist."
                )

        # parameter ID => initial assignment sympy expression
        par_id_to_ia = {
            par.getId(): _try_evalf(ia)
            for par in self.sbml_model.getListOfParameters()
            if (ia := self._get_element_initial_assignment(par.getId()))
            is not None
        }

        fixed_parameters = [
            parameter
            for parameter in self.sbml_model.getListOfParameters()
            if parameter.getId() in fixed_parameter_ids
        ]
        for parameter in fixed_parameters:
            ia_math = par_id_to_ia.get(parameter.getId())
            if (
                (ia_math is not None and not ia_math.is_Number)
                or self._is_assignment_rule_target(parameter)
                or self._is_rate_rule_target(parameter)
            ):
                raise SBMLException(
                    f"Cannot turn parameter {parameter.getId()} into a "
                    "constant/fixed parameter since it either has an "
                    "initial assignment or is the target of an assignment or "
                    "rate rule."
                )

        parameters = [
            parameter
            for parameter in self.sbml_model.getListOfParameters()
            if parameter.getId() not in fixed_parameter_ids
            and (
                (ia_math := par_id_to_ia.get(parameter.getId())) is None
                or ia_math.is_Number
            )
            and not self._is_assignment_rule_target(parameter)
            and parameter.getId() not in hardcode_symbols
        ]

        loop_settings = {
            SymbolId.FREE_PARAMETER: {"var": parameters, "name": "parameter"},
            SymbolId.FIXED_PARAMETER: {
                "var": fixed_parameters,
                "name": "fixed_parameter",
            },
        }

        for partype, settings in loop_settings.items():
            for par in settings["var"]:
                self._symbols[partype][_get_identifier_symbol(par)] = {
                    "name": par.getName() if par.isSetName() else par.getId(),
                    "value": par_id_to_ia.get(
                        par.getId(), sp.Float(par.getValue())
                    ),
                }

        # Parameters that need to be turned into expressions or species
        #  so far, this concerns parameters with _symbolic initial assignments
        #  (those have been skipped above) that are not rate rule targets

        # Set of symbols in initial assignments that still allows handling them
        #  via amici expressions
        syms_allowed_in_expr_ia = set(
            self._symbols[SymbolId.FREE_PARAMETER]
        ) | set(self._symbols[SymbolId.FIXED_PARAMETER])

        for par in self.sbml_model.getListOfParameters():
            if (
                (ia := par_id_to_ia.get(par.getId())) is not None
                and not ia.is_Number
                and not self._is_rate_rule_target(par)
            ):
                if not (ia.free_symbols - syms_allowed_in_expr_ia):
                    self._symbols[SymbolId.EXPRESSION][
                        _get_identifier_symbol(par)
                    ] = {
                        "name": par.getName()
                        if par.isSetName()
                        else par.getId(),
                        "value": ia,
                    }
                else:
                    # Convert parameters with initial assignments that
                    #  explicitly depend on time to species.
                    #  We can't represent that as expression, since the
                    #  initial simulation time is only known at the time of the
                    #  simulation, so we can't substitute it.
                    # Also, any parameter with an initial assignment
                    #  that expression that is implicitly time-dependent
                    #  must be converted to a species to avoid re-evaluating
                    #  the initial assignment at every time step.
                    self._symbols[SymbolId.SPECIES][
                        _get_identifier_symbol(par)
                    ] = {
                        "name": par.getName()
                        if par.isSetName()
                        else par.getId(),
                        "init": ia,
                        "dt": sp.Float(0),
                        "amount": True,
                        "constant": True,
                        "index": len(self._symbols[SymbolId.SPECIES]),
                    }

    @log_execution_time("processing SBML reactions", logger)
    def _process_reactions(self):
        """
        Get reactions from SBML model.
        """
        reactions = self.sbml_model.getListOfReactions()
        nr = len(reactions)
        nx = len(self._symbols[SymbolId.SPECIES])
        # stoichiometric matrix
        self._stoichiometric_matrix = sp.SparseMatrix(sp.zeros(nx, nr))
        self._flux_vector = sp.zeros(nr, 1)
        # Use reaction IDs as IDs for flux expressions (note that prior to SBML
        #  level 3 version 2 the ID attribute was not mandatory and may be
        #  unset)
        self._flux_ids = [
            f"flux_{reaction.getId()}"
            if reaction.isSetId()
            else f"flux_r{reaction_idx}"
            for reaction_idx, reaction in enumerate(reactions)
        ]

        reaction_ids = [
            reaction.getId() for reaction in reactions if reaction.isSetId()
        ]

        for reaction_index, reaction in enumerate(reactions):
            for element_list, sign in [
                (reaction.getListOfReactants(), -1),
                (reaction.getListOfProducts(), 1),
            ]:
                for element in element_list:
                    stoichiometry = self._get_element_stoichiometry(element)
                    sbml_species = self.sbml_model.getSpecies(
                        element.getSpecies()
                    )
                    if self._is_assignment_rule_target(sbml_species):
                        continue
                    species_id = _get_identifier_symbol(sbml_species)
                    species = self._symbols[SymbolId.SPECIES][species_id]

                    if species["constant"]:
                        continue

                    # Division by species compartment size (to find the
                    # rate of change in species concentration) now occurs
                    # in the `dx_dt` method in "de_export.py", which also
                    # accounts for possibly variable compartments.
                    self._stoichiometric_matrix[
                        species["index"], reaction_index
                    ] += sign * stoichiometry * species["conversion_factor"]

            if reaction.isSetId():
                sym_math = self._local_symbols[reaction.getId()]
            else:
                sym_math = self._sympify(
                    reaction.getKineticLaw() or sp.Float(0)
                )

            self._flux_vector[reaction_index] = sym_math

            if any(
                str(symbol) in reaction_ids
                for symbol in self._flux_vector[reaction_index].free_symbols
            ):
                raise SBMLException(
                    "Kinetic laws involving reaction ids are currently"
                    " not supported!"
                )

    @log_execution_time("processing SBML rules", logger)
    def _process_rules(self) -> None:
        """
        Process Rules defined in the SBML model.
        """
        for rule in self.sbml_model.getListOfRules():
            # rate rules are processed in _process_species
            if rule.getTypeCode() == libsbml.SBML_RATE_RULE:
                continue

            if rule.getTypeCode() == libsbml.SBML_ALGEBRAIC_RULE:
                if self._sbml_doc.getLevel() < 3:
                    # not interested in implementing level 2 boundary condition
                    # shenanigans, see test 01787 in the sbml testsuite
                    raise SBMLException(
                        "Algebraic rules are only supported in SBML L3+"
                    )
                self._process_rule_algebraic(rule)
            else:
                self._process_rule_assignment(rule)

        self._symbols[SymbolId.EXPRESSION] = toposort_symbols(
            self._symbols[SymbolId.EXPRESSION], "value"
        )

        # substitute symbols that must not occur in the definition of x0
        # allowed symbols: amici model parameters and time
        allowed_syms = (
            set(self._symbols[SymbolId.FREE_PARAMETER])
            | set(self._symbols[SymbolId.FIXED_PARAMETER])
            | {sbml_time_symbol}
        )
        for species in self._symbols[SymbolId.SPECIES].values():
            while True:
                species["init"] = species["init"].subs(self._compartments)
                sym_math, rateof_to_dummy = _rateof_to_dummy(species["init"])
                old_init = species["init"]
                if (
                    sym_math.free_symbols
                    - allowed_syms
                    - set(rateof_to_dummy.values())
                    == set()
                ):
                    break
                species["init"] = self._make_initial(
                    smart_subs_dict(
                        species["init"],
                        self._symbols[SymbolId.EXPRESSION],
                        "value",
                    )
                )
                if species["init"] == old_init:
                    raise AssertionError(
                        f"Infinite loop detected in _process_rules {species}."
                    )

    def _process_rule_algebraic(self, rule: libsbml.AlgebraicRule):
        formula = self._sympify(rule)
        if formula is None:
            return

        free_variables = set()
        # SBML L3V2 spec, p. 61:
        # "Therefore, if an algebraic rule is introduced in a model,
        # for at least one of the entities referenced in the rule’s
        # math element the value of that entity must not be
        # completely determined by other constructs in the model"
        # find those elements:
        for symbol in formula.free_symbols:
            if symbol == sbml_time_symbol:
                continue
            sbml_var = self.sbml_model.getElementBySId(str(symbol))
            if sbml_var is None:
                raise SBMLException(
                    f"Algebraic rule references unexpected symbol '{symbol}'."
                )
            # This means that at least this entity must
            # not have the attribute constant=“true”
            if sbml_var.isSetConstant() and sbml_var.getConstant():
                continue
            # and there must also not be a rate rule or assignment
            # rule for it
            if self._is_assignment_rule_target(
                sbml_var
            ) or self._is_rate_rule_target(sbml_var):
                continue
            # Furthermore, if the entity is a Species object, its value
            # must not be determined by reactions, which means that it
            # must either have the attribute boundaryCondition=“false”
            # or else not be involved in any reaction at all.
            is_species = isinstance(sbml_var, libsbml.Species)
            is_boundary_condition = (
                is_species
                and sbml_var.isSetBoundaryCondition()
                and sbml_var.getBoundaryCondition()
            )
            is_involved_in_reaction = (
                is_species
                and symbol in self._symbols[SymbolId.SPECIES]
                and not smart_is_zero_matrix(
                    self._stoichiometric_matrix[
                        list(self._symbols[SymbolId.SPECIES].keys()).index(
                            symbol
                        ),
                        :,
                    ]
                )
            )
            if (
                is_species
                and not is_boundary_condition
                and is_involved_in_reaction
            ):
                continue
            free_variables.add(symbol)

        # this should be guaranteed by sbml validation, so better check to make sure
        # we don't mess anything up
        assert len(free_variables) >= 1

        self._symbols[SymbolId.ALGEBRAIC_EQUATION][
            symbol_with_assumptions(
                f"ae{len(self._symbols[SymbolId.ALGEBRAIC_EQUATION])}"
            )
        ] = {"value": formula}
        # remove the symbol from the original definition and add to
        # algebraic symbols (if not already done)
        for var in free_variables:
            if var in self._symbols[SymbolId.FIXED_PARAMETER]:
                raise SBMLException(
                    "There are algebraic rules that specify the "
                    f"value of {var}, which is also marked as "
                    "fixed parameter. This is currently not supported! "
                    f"If {var} is supposed to be a fixed parameter, "
                    "set its SBML attribute `constant` to True."
                )

            if var in self._symbols[SymbolId.ALGEBRAIC_STATE]:
                continue
            if var in self._compartments:
                init = self._compartments[var]
                symbol = {
                    "name": str(var),
                    "value": init,
                }
                symbol_id = "compartment"
                var_ix = np.nan
                del self._compartments[var]
            else:
                symbol_id, source_symbols = next(
                    (
                        (symbol_id, self._symbols[symbol_id])
                        for symbol_id in (
                            SymbolId.FREE_PARAMETER,
                            SymbolId.SPECIES,
                        )
                        if var in self._symbols[symbol_id]
                    ),
                )
                var_ix = list(source_symbols.keys()).index(var)
                symbol = source_symbols.pop(var)
            # update symbol and adapt stoichiometric matrix
            if symbol_id != SymbolId.SPECIES:
                # parameters have numeric values so we can use Float here
                symbol["init"] = sp.Float(symbol.pop("value"))
                # if not a species, add a zeros row to the stoichiometric
                # matrix
                if (
                    isinstance(symbol["init"], float)
                    and np.isnan(symbol["init"])
                ) or (
                    isinstance(symbol["init"], sp.Number)
                    and symbol["init"] == sp.nan
                ):
                    # placeholder, needs to be determined in IC calculation
                    symbol["init"] = sp.Float(0.0)
                if len(self._stoichiometric_matrix):
                    self._stoichiometric_matrix = (
                        self._stoichiometric_matrix.row_insert(
                            self._stoichiometric_matrix.shape[0],
                            sp.SparseMatrix(
                                1, self._stoichiometric_matrix.shape[1], {}
                            ),
                        )
                    )
                else:
                    # we can't use Matrix.row_insert for (n, 0) matrices
                    #  as this will ignore the previous number of rows
                    self._stoichiometric_matrix = sp.SparseMatrix(
                        self._stoichiometric_matrix.shape[0] + 1, 0, {}
                    )
            elif var_ix != self._stoichiometric_matrix.shape[0] - 1:
                # if not the last col, move it to the end
                # as we reorder state variables
                state_ordering = list(
                    range(
                        len(self._symbols[SymbolId.SPECIES])
                        + len(self._symbols[SymbolId.ALGEBRAIC_STATE])
                        + 1
                    )
                )
                state_ordering.append(state_ordering.pop(var_ix))
                self._stoichiometric_matrix = self._stoichiometric_matrix[
                    state_ordering, :
                ]

            self._symbols[SymbolId.ALGEBRAIC_STATE][var] = symbol

    def _process_rule_assignment(self, rule: libsbml.AssignmentRule):
        sbml_var = self.sbml_model.getElementBySId(rule.getVariable())
        sym_id = symbol_with_assumptions(rule.getVariable())

        # Check whether this rule is a spline rule.
        if not self._discard_annotations:
            if rule.getTypeCode() == libsbml.SBML_ASSIGNMENT_RULE:
                annotation = AbstractSpline.get_annotation(rule)
                if annotation is not None:
                    spline = AbstractSpline.from_annotation(
                        sym_id,
                        annotation,
                        locals_=self._local_symbols,
                    )
                    if (
                        spline.evaluate_at != amici_time_symbol
                        and spline.evaluate_at != sbml_time_symbol
                    ):
                        raise NotImplementedError(
                            "AMICI at the moment does not support splines "
                            "whose evaluation point is not the model time. "
                            f"'evaluate_at' was {spline.evaluate_at} "
                            f"({type(spline.evaluate_at)}) in spline for "
                            f"symbol {sym_id}."
                        )
                    self._splines.append(spline)
                    return

        formula = self._sympify(rule)
        if formula is None:
            return

        if isinstance(sbml_var, libsbml.Species):
            self._species_assignment_rules[sym_id] = formula

        elif isinstance(sbml_var, libsbml.Compartment):
            self._compartment_assignment_rules[sym_id] = formula
            self._compartments[sym_id] = formula

        elif isinstance(sbml_var, libsbml.Parameter):
            self._parameter_assignment_rules[sym_id] = formula

        self._symbols[SymbolId.EXPRESSION][sym_id] = {
            "name": str(sym_id),
            "value": formula,
        }

    def _process_time(self) -> None:
        """
        Convert time_symbol into cpp variable.
        """
        self._replace_in_all_expressions(sbml_time_symbol, amici_time_symbol)

    def _convert_event_assignment_parameter_targets_to_species(self):
        """
        Convert parameters that are targets of event assignments to species.

        This is for the convenience of only implementing event assignments for
        "species".
        """
        parameter_targets = _collect_event_assignment_parameter_targets(
            self.sbml_model
        )
        for parameter_target in parameter_targets:
            # Parameter rate rules already exist as species.
            if parameter_target in self._symbols[SymbolId.SPECIES]:
                continue
            if parameter_target in self._parameter_assignment_rules:
                raise SBMLException(
                    "AMICI does not currently support models with SBML events "
                    "that affect parameters that are also the target of "
                    "assignment rules."
                )
            parameter_def = None
            for symbol_id in {
                SymbolId.FREE_PARAMETER,
                SymbolId.FIXED_PARAMETER,
                SymbolId.EXPRESSION,
            }:
                if parameter_target in self._symbols[symbol_id]:
                    # `parameter_target` should only exist in one of the
                    # `symbol_id` dictionaries.
                    if parameter_def is not None:
                        raise AssertionError(
                            "Unexpected error. The parameter target of an "
                            "event assignment was processed twice."
                        )
                    parameter_def = self._symbols[symbol_id].pop(
                        parameter_target
                    )
            if parameter_def is None:
                # this happens for parameters that have initial assignments
                # or are assignment rule targets
                par = self.sbml_model.getElementBySId(str(parameter_target))
                ia_init = self._get_element_initial_assignment(par.getId())
                parameter_def = {
                    "name": par.getName() if par.isSetName() else par.getId(),
                    "value": sp.Float(par.getValue())
                    if ia_init is None
                    else ia_init,
                }
            # Fixed parameters are added as species such that they can be
            # targets of events.
            self._symbols[SymbolId.SPECIES][parameter_target] = {
                "name": parameter_def["name"],
                "init": parameter_def["value"],
                # 'compartment': None,  # can ignore for amounts
                "constant": False,
                "amount": True,
                # 'conversion_factor': 1.0,  # can be ignored
                "index": len(self._symbols[SymbolId.SPECIES]),
                "dt": sp.Float(0),
            }

    def _convert_event_assignment_compartment_targets_to_species(self):
        """Find compartments that are event assignment targets and convert
        those compartments to species."""
        for event in self.sbml_model.getListOfEvents():
            for event_assignment in event.getListOfEventAssignments():
                if event_assignment.getMath() is None:
                    # Ignore event assignments with no change in value.
                    continue
                variable = symbol_with_assumptions(
                    event_assignment.getVariable()
                )
                if variable not in self._compartments:
                    continue
                if variable in self._symbols[SymbolId.SPECIES]:
                    # Compartments with rate rules are already present as
                    # species
                    continue

                self._symbols[SymbolId.SPECIES][variable] = {
                    "name": str(variable),
                    "init": self._compartments[variable],
                    # 'compartment': None,  # can ignore for amounts
                    "constant": False,
                    "amount": True,
                    # 'conversion_factor': 1.0,  # can be ignored
                    "index": len(self._symbols[SymbolId.SPECIES]),
                    "dt": sp.Float(0),
                }
                del self._compartments[variable]

    @log_execution_time("processing SBML events", logger)
    def _process_events(self) -> None:
        """Process SBML events."""
        events = self.sbml_model.getListOfEvents()

        # Used to update species concentrations when an event affects a
        # compartment.
        concentration_species_by_compartment = {
            symbol_with_assumptions(c.getId()): []
            for c in self.sbml_model.getListOfCompartments()
        }
        for species, species_def in self._symbols[SymbolId.SPECIES].items():
            if (
                # Species is a concentration
                not species_def.get("amount", True)
                and
                # Species has a compartment
                "compartment" in species_def
            ):
                concentration_species_by_compartment[
                    species_def["compartment"]
                ].append(species)

        # Currently, all event assignment targets must exist in
        # self.symbols[SymbolId.SPECIES]
        state_vector = list(self._symbols[SymbolId.SPECIES].keys())

        for ievent, event in enumerate(events):
            # get the event id (which is optional unfortunately)
            event_id = event.getId()
            if event_id is None or event_id == "":
                event_id = f"event_{ievent}"
            event_sym = sp.Symbol(event_id)

            # get and parse the trigger function
            trigger_sbml = event.getTrigger()
            trigger_sym = self._sympify(trigger_sbml, bool2num=False)
            trigger = _parse_event_trigger(trigger_sym)

            # parse the boluses / event assignments
            assignment_exprs = {}
            event_assignments = event.getListOfEventAssignments()
            compartment_event_assignments: set[tuple[sp.Symbol, sp.Expr]] = (
                set()
            )
            for event_assignment in event_assignments:
                variable_sym = symbol_with_assumptions(
                    event_assignment.getVariable()
                )
                if event_assignment.getMath() is None:
                    # Ignore event assignments with no change in value.
                    continue
                formula = self._sympify(event_assignment)
                try:
                    # Try to find the species in the state vector.
                    _ = state_vector.index(variable_sym)
                    assignment_exprs[variable_sym] = formula
                except ValueError:
                    raise SBMLException(
                        "Could not process event assignment for "
                        f"{str(variable_sym)}. AMICI currently only allows "
                        "event assignments to species; parameters; or, "
                        "compartments."
                    )
                try:
                    # Try working with the formula now to detect errors
                    # here instead of at multiple points downstream.
                    _ = formula - variable_sym
                except TypeError:
                    raise SBMLException(
                        "Could not process event assignment for "
                        f"{str(variable_sym)}. AMICI only allows _symbolic "
                        "expressions as event assignments."
                    )
                if variable_sym in concentration_species_by_compartment:
                    compartment_event_assignments.add((variable_sym, formula))

                for (
                    comp,
                    assignment,
                ) in self._compartment_assignment_rules.items():
                    if variable_sym not in assignment.free_symbols:
                        continue
                    compartment_event_assignments.add((comp, formula))

            # Update the concentration of species with concentration units
            # in compartments that were affected by the event assignments.
            for compartment_sym, formula in compartment_event_assignments:
                for species_sym in concentration_species_by_compartment[
                    compartment_sym
                ]:
                    # If the species was not affected by an event assignment,
                    # then the old value should be updated.
                    if species_sym not in assignment_exprs:
                        species_value = species_sym
                    # else the species was affected by an event assignment,
                    # hence the updated value should be updated further.
                    else:
                        species_value = assignment_exprs[species_sym]
                    # New species value is old amount / new volume.
                    assignment_exprs[species_sym] = (
                        species_value * compartment_sym / formula
                    )

            initial_value = (
                trigger_sbml.getInitialValue()
                if trigger_sbml is not None
                else True
            )
            if (
                self._symbols[SymbolId.ALGEBRAIC_EQUATION]
                and not initial_value
            ):
                # in principle this could be implemented, requires running
                # IDACalcIc (in solver->setup) before check event initialization
                # (in model->initialize), but at the time of writing this sounded
                # like something that might break stuff in all kinds of other places
                # (it might not, but this could be checked when someone actually
                # needs the feature).
                raise SBMLException(
                    "Events with initial values are not supported in models with"
                    " algebraic rules."
                )

            # Store `useValuesFromTriggerTime` attribute for checking later
            # Since we assume valid in SBML models here, this attribute is
            # either given (mandatory in L3), or defaults to True (L2)
            use_trig_val = (
                event.getUseValuesFromTriggerTime()
                if event.isSetUseValuesFromTriggerTime()
                else True
            )

            self._symbols[SymbolId.EVENT][event_sym] = {
                "name": event_id,
                "value": trigger,
                "assignments": assignment_exprs,
                "initial_value": initial_value,
                "use_values_from_trigger_time": use_trig_val,
                "priority": self._sympify(event.getPriority()),
            }

            if self.jax:
                # Add a negative event for JAX models to handle
                neg_event_id = event_id + "_negative"
                neg_event_sym = sp.Symbol(neg_event_id)
                self.symbols[SymbolId.EVENT][neg_event_sym] = {
                    "name": neg_event_id,
                    "value": -trigger,
                    "assignments": None,
                    "initial_value": not initial_value,
                    "use_values_from_trigger_time": use_trig_val,
                    "priority": self._sympify(event.getPriority()),
                }

    @log_execution_time("processing observation model", logger)
    def _process_observation_model(
        self,
        observation_model: list[MeasurementChannel] = None,
    ) -> None:
        """
        Perform _symbolic computations required for observable and objective
        function evaluation.

        Add user-provided observables or make all species, and compartments
        with assignment rules, observable.
        """
        if observation_model:
            # prepare and validate
            for channel in observation_model:
                # gather local symbols before parsing observable and sigma
                # formulas
                self._add_local_symbol(
                    channel.id, symbol_with_assumptions(channel.id)
                )

                if (
                    channel.event_id is not None
                    and sp.Symbol(channel.event_id)
                    not in self._symbols[SymbolId.EVENT]
                ):
                    raise ValueError(
                        "Could not find an event with the event identifier "
                        f"{channel.event_id} for the event observable with "
                        f"name {channel.name or channel.id}."
                    )

            # Add time-resolved observables
            self._symbols[SymbolId.OBSERVABLE] = {
                symbol_with_assumptions(channel.id): {
                    "name": channel.name or f"y{iobs}",
                    "value": self._sympify(channel.formula),
                    "transformation": noise_distribution_to_observable_transformation(
                        channel.noise_distribution
                    ),
                }
                for iobs, channel in enumerate(
                    filter(lambda c: c.is_time_resolved, observation_model)
                )
            }
            # check for nesting of observables (unsupported)
            _check_symbol_nesting(
                self._symbols[SymbolId.OBSERVABLE], "observable"
            )

            # Add event observables
            self._symbols[SymbolId.EVENT_OBSERVABLE] = {
                symbol_with_assumptions(channel.id): {
                    "name": channel.name or f"z{iobs}",
                    "value": self._sympify(channel.formula),
                    "event": sp.Symbol(channel.event_id),
                    "transformation": noise_distribution_to_observable_transformation(
                        channel.noise_distribution
                    ),
                }
                for iobs, channel in enumerate(
                    filter(lambda c: c.is_event_resolved, observation_model)
                )
            }
            # Check for time errors
            wrong_t = sp.Symbol("t")
            for eo in self._symbols[SymbolId.EVENT_OBSERVABLE].values():
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
                self._symbols[SymbolId.EVENT_OBSERVABLE], "eventObservable"
            )

            # FIXME: that name-based approach is a bit fragile
            noise_pars = list(
                {
                    name
                    for channel in observation_model
                    if channel.is_time_resolved and channel.sigma is not None
                    for symbol in self._sympify(channel.sigma).free_symbols
                    if re.match(r"noiseParameter\d+$", (name := str(symbol)))
                }
            )

            self._symbols[SymbolId.NOISE_PARAMETER] = {
                symbol_with_assumptions(noise_par): {"name": noise_par}
                for noise_par in noise_pars
            }

            observable_pars = list(
                {
                    name
                    for channel in observation_model
                    if channel.is_time_resolved and channel.sigma is not None
                    for symbol in self._sympify(channel.formula).free_symbols
                    if re.match(
                        r"observableParameter\d+$", (name := str(symbol))
                    )
                }
            )
            self._symbols[SymbolId.OBSERVABLE_PARAMETER] = {
                symbol_with_assumptions(obs_par): {"name": obs_par}
                for obs_par in observable_pars
            }

        elif observation_model is None:
            self._generate_default_observables()

        # TODO: _process_log_likelihood everything below and
        #   needs some refactoring
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
        self._process_log_likelihood(sigmas, noise_distributions)

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
        self._process_log_likelihood(sigmas, noise_distributions, events=True)
        self._process_log_likelihood(
            sigmas,
            noise_distributions,
            events=True,
            event_reg=True,
        )

    def _generate_default_observables(self):
        """
        Generate default observables from species, compartments and
        (initial) assignment rules.
        """
        self._symbols[SymbolId.OBSERVABLE] = {
            symbol_with_assumptions(f"y{state_id}"): {
                "name": state["name"],
                "value": state_id,
            }
            for state_id, state in {
                **self._symbols[SymbolId.SPECIES],
                **self._symbols[SymbolId.ALGEBRAIC_STATE],
            }.items()
        }

        for variable, formula in itt.chain(
            self._parameter_assignment_rules.items(),
            self._initial_assignments.items(),
            self._compartment_assignment_rules.items(),
            self._species_assignment_rules.items(),
            self._compartments.items(),
        ):
            symbol = symbol_with_assumptions(f"y{variable}")
            # Assignment rules take precedence over compartment volume
            # definitions, so they need to be evaluated first.
            # Species assignment rules always overwrite.
            if (
                symbol in self._symbols[SymbolId.OBSERVABLE]
                and variable not in self._species_assignment_rules
            ):
                continue
            self._symbols[SymbolId.OBSERVABLE][symbol] = {
                "name": str(variable),
                "value": formula,
            }

    def _process_log_likelihood(
        self,
        sigmas: dict[str, str | float],
        noise_distributions: dict[str, str],
        events: bool = False,
        event_reg: bool = False,
    ):
        """
        Perform _symbolic computations required for objective function
        evaluation.

        :param sigmas:
            See :py:func:`SBMLImporter._process_observables`

        :param noise_distributions:
            See :py:func:`SBMLImporter._process_observables`

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

        for obs_id, obs in self._symbols[obs_symbol].items():
            obs["measurement_symbol"] = generate_measurement_symbol(obs_id)
            if event_reg:
                obs["reg_symbol"] = generate_regularization_symbol(obs_id)

        if not event_reg:
            sigmas = {
                obs_id: self._sympify(sigma)
                for obs_id, sigma in sigmas.items()
            }
            obs_syms = set(self._symbols[obs_symbol].keys())
            for obs_id in self._symbols[obs_symbol]:
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
            self._symbols[sigma_symbol] = {
                symbol_with_assumptions(f"sigma_{obs_id}"): {
                    "name": f"sigma_{obs['name']}",
                    "value": sigmas.get(str(obs_id), sp.Float(1.0)),
                }
                for obs_id, obs in self._symbols[obs_symbol].items()
            }

        self._symbols[llh_symbol] = {}
        for (obs_id, obs), (sigma_id, sigma) in zip(
            self._symbols[obs_symbol].items(),
            self._symbols[sigma_symbol].items(),
            strict=True,
        ):
            symbol = symbol_with_assumptions(f"J{obs_id}")
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
                value = value.subs(obs_id, obs["reg_symbol"])
            self._symbols[llh_symbol][symbol] = {
                "name": f"J{obs['name']}",
                "value": value,
                "dist": dist,
            }

    @log_execution_time("processing SBML initial assignments", logger)
    def _process_initial_assignments(self):
        """
        Accounts for initial assignments of parameters and species
        references. Initial assignments for species and compartments are
        processed in :py:func:`amici.SBMLImporter._process_initial_species` and
        :py:func:`amici.SBMLImporter._process_compartments` respectively.
        """
        for ia in self.sbml_model.getListOfInitialAssignments():
            identifier = _get_identifier_symbol(ia)
            if identifier in itt.chain(
                self._symbols[SymbolId.SPECIES],
                self._compartments,
                self._symbols[SymbolId.EXPRESSION],
                self._symbols[SymbolId.FREE_PARAMETER],
                self._symbols[SymbolId.FIXED_PARAMETER],
            ):
                continue

            sym_math = self._get_element_initial_assignment(ia.getId())
            if sym_math is None:
                continue

            sym_math = self._make_initial(
                smart_subs_dict(
                    sym_math, self._symbols[SymbolId.EXPRESSION], "value"
                )
            )
            self._initial_assignments[_get_identifier_symbol(ia)] = sym_math

        # sort and flatten
        self._initial_assignments = toposort_symbols(self._initial_assignments)
        for ia_id, ia in self._initial_assignments.items():
            self._initial_assignments[ia_id] = smart_subs_dict(
                ia, self._initial_assignments
            )

        for identifier, sym_math in list(self._initial_assignments.items()):
            self._replace_in_all_expressions(identifier, sym_math)

    @log_execution_time("processing SBML species references", logger)
    def _process_species_references(self):
        """
        Replaces species references that define anything but stoichiometries.

        Species references for stoichiometries are processed in
        :py:func:`amici.SBMLImporter._process_reactions`.
        """
        # doesnt look like there is a better way to get hold of those lists:
        species_references = _get_list_of_species_references(self.sbml_model)
        for species_reference in species_references:
            if (
                hasattr(species_reference, "getStoichiometryMath")
                and species_reference.getStoichiometryMath() is not None
            ):
                raise SBMLException(
                    "StoichiometryMath is currently not "
                    "supported for species references."
                )
            if species_reference.getId() == "":
                continue

            stoich = self._get_element_stoichiometry(species_reference)
            self._replace_in_all_expressions(
                _get_identifier_symbol(species_reference),
                self._sympify(stoich),
            )

    def _make_initial(
        self, sym_math: sp.Expr | None | float
    ) -> sp.Expr | None | float:
        """
        Transforms an expression to its value at the initial time point by
        replacing species by their initial values.

        :param sym_math:
            _symbolic expression
        :return:
            transformed expression
        """
        if not isinstance(sym_math, sp.Expr):
            return sym_math

        sym_math, rateof_to_dummy = _rateof_to_dummy(sym_math)

        # we can't rely on anything else being properly initialized at this point, so we need to
        # compute all initial values from scratch, recursively
        for var in sym_math.free_symbols:
            element_id = str(var)
            # already recursive since _get_element_initial_assignment calls _make_initial
            if (
                ia := self._get_element_initial_assignment(element_id)
            ) is not None:
                sym_math = sym_math.subs(var, ia)
            elif (
                species := self.sbml_model.getSpecies(element_id)
            ) is not None:
                # recursive!
                init = self._make_initial(get_species_initial(species))
                sym_math = sym_math.subs(var, init)
            elif var in self._symbols[SymbolId.SPECIES]:
                sym_math = sym_math.subs(
                    var, self._symbols[SymbolId.SPECIES][var]["init"]
                )
            elif var in self._symbols[SymbolId.ALGEBRAIC_STATE]:
                sym_math = sym_math.subs(
                    var, self._symbols[SymbolId.ALGEBRAIC_STATE][var]["init"]
                )
            elif (
                element := self.sbml_model.getElementBySId(element_id)
            ) and self._is_rate_rule_target(element):
                # no need to recurse here, as value is numeric
                init = sp.Float(element.getValue())
                sym_math = sym_math.subs(var, init)
            elif spline := [
                spl for spl in self._splines if spl.sbml_id == var
            ]:
                # x0 must not depend on splines -- substitute
                assert len(spline) == 1
                spline = spline[0]
                if spline.evaluate_at != sbml_time_symbol:
                    raise NotImplementedError(
                        "AMICI at the moment does not support splines "
                        "whose evaluation point is not the model time."
                        f"'evaluate_at' was {spline.evaluate_at} "
                        f"({type(spline.evaluate_at)}) in spline for "
                        f"symbol {var}."
                    )
                sym_math = sym_math.subs(
                    var, spline.evaluate(sbml_time_symbol)
                )

        sym_math = _dummy_to_rateof(sym_math, rateof_to_dummy)

        return sym_math

    def _process_conservation_laws(self, ode_model: DEModel) -> None:
        """
        Find conservation laws in reactions and species.

        :param ode_model:
            :class:`DEModel` object with basic definitions
        """
        conservation_laws = []

        # Create conservation laws for constant species
        species_solver = _add_conservation_for_constant_species(
            ode_model, conservation_laws
        )
        # Non-constant species processed here
        if (
            "AMICI_EXPERIMENTAL_SBML_NONCONST_CLS" in os.environ
            or "GITHUB_ACTIONS" in os.environ
        ):
            species_solver = list(
                set(
                    self._add_conservation_for_non_constant_species(
                        ode_model, conservation_laws
                    )
                )
                & set(species_solver)
            )

        # add algebraic variables to species_solver as they were ignored above
        ndifferential = len(ode_model._differential_states)
        nalgebraic = len(ode_model._algebraic_states)
        species_solver.extend(
            list(range(ndifferential, ndifferential + nalgebraic))
        )

        # Check, whether species_solver is empty now. As currently, AMICI
        # cannot handle ODEs without species, CLs must be switched off in this
        # case
        if not len(species_solver):
            conservation_laws = []
            species_solver = list(range(ode_model.num_states_rdata()))

        # prune out species from stoichiometry and
        if len(self._stoichiometric_matrix):
            self._stoichiometric_matrix = self._stoichiometric_matrix[
                species_solver, :
            ]

        # add the found CLs to the ode_model
        for cl in conservation_laws:
            ode_model.add_conservation_law(**cl)

    def _get_conservation_laws_demartino(
        self,
        ode_model: DEModel,
    ) -> list[tuple[int, list[int], list[float]]]:
        """Identify conservation laws based on algorithm by DeMartino et al.
        (see conserved_moieties.py).

        :param ode_model: Model for which to compute conserved quantities
        :returns: List of one tuple per conservation law, each containing:
            (0) the index of the (solver-)species to eliminate,
            (1) (solver-)indices of all species engaged in the conserved
            quantity (including the eliminated one)
            (2) coefficients for the species in (1)
        """
        if not len(self._stoichiometric_matrix) or not len(
            self._symbols[SymbolId.SPECIES]
        ):
            # no species or no reactions -> no conservation laws
            return []

        sm = self._stoichiometric_matrix[
            : len(self._symbols[SymbolId.SPECIES]), :
        ]

        try:
            stoichiometric_list = [float(entry) for entry in sm.T.flat()]
        except TypeError:
            # Due to the numerical algorithm currently used to identify
            #  conserved quantities, we can't have symbols in the
            #  stoichiometric matrix
            warnings.warn(
                "Conservation laws for non-constant species in "
                "combination with parameterized stoichiometric "
                "coefficients are not currently supported "
                "and will be turned off.",
                stacklevel=1,
            )
            return []

        if not _non_const_conservation_laws_supported(self.sbml_model):
            return []

        from .conserved_quantities_demartino import (
            compute_moiety_conservation_laws,
        )

        cls_state_idxs, cls_coefficients = compute_moiety_conservation_laws(
            stoichiometric_list,
            *sm.shape,
            rng_seed=32,
            species_names=[x.get_id() for x in ode_model._differential_states],
        )

        # Sparsify conserved quantities
        #  ``compute_moiety_conservation_laws`` identifies conserved quantities
        #  with positive coefficients. The resulting system is, therefore,
        #  often non-sparse. This leads to circular dependencies in the
        #  state expressions of eliminated states. The identified conserved
        #  quantities are linearly independent. We can construct `A` as in
        #  `A * x0 = total_cl` and bring it to reduced row echelon form. The
        #  pivot species are the ones to be eliminated. The resulting state
        #  expressions are sparse and void of any circular dependencies.
        A = sp.zeros(
            len(cls_coefficients), len(ode_model._differential_states)
        )
        for i_cl, (cl, coefficients) in enumerate(
            zip(cls_state_idxs, cls_coefficients, strict=True)
        ):
            for i, c in zip(cl, coefficients, strict=True):
                A[i_cl, i] = sp.Rational(c)
        rref, pivots = A.rref()

        raw_cls = []
        for i_cl, target_state_model_idx in enumerate(pivots):
            # collect values for species engaged in the current CL
            state_idxs = [i for i, coeff in enumerate(rref[i_cl, :]) if coeff]
            coefficients = [coeff for coeff in rref[i_cl, :] if coeff]
            raw_cls.append(
                (target_state_model_idx, state_idxs, coefficients),
            )
        return raw_cls

    def _get_conservation_laws_rref(
        self,
    ) -> list[tuple[int, list[int], list[float]]]:
        """Identify conservation laws based on left nullspace of the
        stoichiometric matrix, computed through (numeric) Gaussian elimination

        :returns: List of one tuple per conservation law, each containing:
            (0) the index of the (solver-)species to eliminate,
            (1) (solver-)indices of all species engaged in the conserved
            quantity (including the eliminated one)
            (2) coefficients for the species in (1)
        """
        if not len(self._stoichiometric_matrix) or not len(
            self._symbols[SymbolId.SPECIES]
        ):
            # no species or no reactions -> no conservation laws
            return []

        import numpy as np
        from numpy.linalg import matrix_rank

        from .conserved_quantities_rref import nullspace_by_rref, rref

        try:
            S = np.asarray(
                self._stoichiometric_matrix[
                    : len(self._symbols[SymbolId.SPECIES]), :
                ],
                dtype=float,
            )
        except TypeError:
            # Due to the numerical algorithm currently used to identify
            #  conserved quantities, we can't have symbols in the
            #  stoichiometric matrix
            warnings.warn(
                "Conservation laws for non-constant species in "
                "combination with parameterized stoichiometric "
                "coefficients are not currently supported "
                "and will be turned off.",
                stacklevel=1,
            )
            return []

        if not _non_const_conservation_laws_supported(self.sbml_model):
            return []

        # Determine rank via SVD
        rank = matrix_rank(S) if S.shape[0] else 0
        if rank == S.shape[0]:
            return []
        kernel = nullspace_by_rref(S.T)
        # Check dimensions - due to numerical errors, nullspace_by_rref may
        #  fail in certain situations
        if kernel.shape[0] != S.shape[0] - rank:
            raise AssertionError(
                "Failed to determine all conserved quantities "
                f"(found {kernel.shape[0]}, expected {S.shape[0] - rank}). "
                "Try another algorithm, disable detection of conservation "
                "laws, or submit a bug report along with the model."
            )
        kernel = rref(kernel)
        raw_cls = []
        for row in kernel:
            state_idxs = [i for i, coeff in enumerate(row) if coeff]
            coefficients = [coeff for coeff in row if coeff]
            raw_cls.append(
                (state_idxs[0], state_idxs, coefficients),
            )

        return raw_cls

    def _add_conservation_for_non_constant_species(
        self, model: DEModel, conservation_laws: list[ConservationLaw]
    ) -> list[int]:
        """Add non-constant species to conservation laws

        :param model:
            DEModel object with basic definitions
        :param conservation_laws:
            List of already known conservation laws
        :returns:
            List of species indices which later remain in the DE solver
        """
        # indices of retained species
        species_solver = list(range(len(model._differential_states)))

        algorithm = os.environ.get("AMICI_EXPERIMENTAL_SBML_NONCONST_CLS", "")
        if algorithm.lower() == "demartino":
            raw_cls = self._get_conservation_laws_demartino(model)
        else:
            raw_cls = self._get_conservation_laws_rref()

        if not raw_cls:
            # no conservation laws identified
            return species_solver

        species_to_be_removed = {x[0] for x in raw_cls}

        # keep new conservations laws separate until we know everything worked
        new_conservation_laws = []
        # previously removed constant species
        eliminated_state_ids = {cl["state"] for cl in conservation_laws}

        all_state_ids = [x.get_sym() for x in model.states()]
        all_compartment_sizes = []
        for state_id in all_state_ids:
            symbol = {
                **self._symbols[SymbolId.SPECIES],
                **self._symbols[SymbolId.ALGEBRAIC_STATE],
            }[state_id]
            if "amount" not in symbol:
                continue  # not a species
            if symbol["amount"]:
                compartment_size = sp.Integer(1)
            else:
                compartment_size = self._compartments[symbol["compartment"]]
            all_compartment_sizes.append(compartment_size)

        # iterate over list of conservation laws, create _symbolic expressions,
        for target_state_model_idx, state_idxs, coefficients in raw_cls:
            if all_state_ids[target_state_model_idx] in eliminated_state_ids:
                # constants state, already eliminated
                continue
            # collect values for species engaged in the current CL
            state_ids = [all_state_ids[i_state] for i_state in state_idxs]
            compartment_sizes = [all_compartment_sizes[i] for i in state_idxs]

            target_state_id = all_state_ids[target_state_model_idx]
            total_abundance = symbol_with_assumptions(f"tcl_{target_state_id}")

            new_conservation_laws.append(
                {
                    "state": target_state_id,
                    "total_abundance": total_abundance,
                    "coefficients": {
                        state_id: coeff * compartment
                        for state_id, coeff, compartment in zip(
                            state_ids,
                            coefficients,
                            compartment_sizes,
                            strict=True,
                        )
                    },
                }
            )
            species_to_be_removed.add(target_state_model_idx)

        conservation_laws.extend(new_conservation_laws)

        # list of species that are not determined by conservation laws
        return [ix for ix in species_solver if ix not in species_to_be_removed]

    def _replace_compartments_with_volumes(self):
        """
        Replaces compartment symbols in expressions with their respective
        (possibly variable) volumes.
        """
        for comp, vol in self._compartments.items():
            if (
                comp in self._symbols[SymbolId.SPECIES]
                or comp in self._symbols[SymbolId.ALGEBRAIC_STATE]
            ):
                # for comps with rate rules volume is only initial
                for state in {
                    **self._symbols[SymbolId.SPECIES],
                    **self._symbols[SymbolId.ALGEBRAIC_STATE],
                }.values():
                    if isinstance(state["init"], sp.Expr):
                        state["init"] = smart_subs(state["init"], comp, vol)
                continue
            self._replace_in_all_expressions(comp, vol)

    def _replace_in_all_expressions(
        self, old: sp.Symbol, new: sp.Expr, replace_identifiers=False
    ) -> None:
        """
        Replace 'old' by 'new' in all _symbolic expressions.

        :param old:
            _symbolic variables to be replaced

        :param new:
            replacement _symbolic variables
        """
        fields = [
            "_stoichiometric_matrix",
            "_flux_vector",
        ]
        for field in fields:
            if field in dir(self):
                self.__setattr__(
                    field, smart_subs(self.__getattribute__(field), old, new)
                )

        dictfields = [
            "_compartment_assignment_rules",
            "_parameter_assignment_rules",
            "_initial_assignments",
        ]
        for dictfield in dictfields:
            d = getattr(self, dictfield)

            # replace identifiers
            if old in d and replace_identifiers:
                d[new] = d[old]
                del d[old]

            if dictfield == "_initial_assignments":
                tmp_new = self._make_initial(new)
            else:
                tmp_new = new

            # replace values
            for k in d:
                d[k] = smart_subs(d[k], old, tmp_new)

        # replace in identifiers
        if replace_identifiers:
            for symbol in [
                SymbolId.EXPRESSION,
                SymbolId.SPECIES,
                SymbolId.ALGEBRAIC_STATE,
            ]:
                # completely recreate the dict to keep ordering consistent
                if old not in self._symbols[symbol]:
                    continue
                self._symbols[symbol] = {
                    smart_subs(k, old, new): v
                    for k, v in self._symbols[symbol].items()
                }

            for symbol in [
                SymbolId.OBSERVABLE,
                SymbolId.LLHY,
                SymbolId.SIGMAY,
            ]:
                if old not in self._symbols[symbol]:
                    continue
                self._symbols[symbol][new] = self._symbols[symbol][old]
                del self._symbols[symbol][old]

        # replace in values
        for symbol in [
            SymbolId.OBSERVABLE,
            SymbolId.LLHY,
            SymbolId.LLHZ,
            SymbolId.SIGMAY,
            SymbolId.SIGMAZ,
            SymbolId.EXPRESSION,
            SymbolId.EVENT,
            SymbolId.EVENT_OBSERVABLE,
            SymbolId.ALGEBRAIC_EQUATION,
        ]:
            for element in self._symbols[symbol].values():
                element["value"] = smart_subs(element["value"], old, new)

        # replace in event assignments
        if self._symbols.get(SymbolId.EVENT, False):
            for event in self._symbols[SymbolId.EVENT].values():
                if event["assignments"] is not None:
                    event["assignments"] = {
                        smart_subs(target, old, new): smart_subs(
                            expr, old, new
                        )
                        for target, expr in event["assignments"].items()
                    }

        for state in {
            **self._symbols[SymbolId.SPECIES],
            **self._symbols[SymbolId.ALGEBRAIC_STATE],
        }.values():
            state["init"] = smart_subs(
                state["init"], old, self._make_initial(new)
            )

            if "dt" in state:
                state["dt"] = smart_subs(state["dt"], old, new)

        # Initial compartment volume may also be specified with an assignment
        # rule (at the end of the _process_species method), hence needs to be
        # processed here too.
        self._compartments = {
            smart_subs(c, old, new) if replace_identifiers else c: smart_subs(
                v, old, self._make_initial(new)
            )
            for c, v in self._compartments.items()
        }

        # Substitute inside spline definitions
        for spline in self._splines:
            spline._replace_in_all_expressions(old, new)

    def _clean_reserved_symbols(self) -> None:
        """
        Remove all reserved symbols from self.symbols
        """
        for sym in RESERVED_SYMBOLS:
            old_symbol = symbol_with_assumptions(sym)
            new_symbol = symbol_with_assumptions(f"amici_{sym}")
            self._replace_in_all_expressions(
                old_symbol, new_symbol, replace_identifiers=True
            )
            for symbols_ids, symbols in self._symbols.items():
                if old_symbol in symbols:
                    # reconstitute the whole dict in order to keep the ordering
                    self._symbols[symbols_ids] = {
                        new_symbol if k == old_symbol else k: v
                        for k, v in symbols.items()
                    }

    def _sympify(
        self,
        var_or_math: libsbml.SBase
        | libsbml.ASTNode
        | float
        | int
        | sp.Basic
        | str,
        piecewise_to_heaviside: bool = True,
        bool2num: bool = True,
    ) -> sp.Expr | None:
        """
        Sympify math expressions with all sanity checks and transformations.

        :param var_or_math:
            A math expression to sympify:

            * a string (see `sympy.sympify`)
            * a number
            * a libSBML ASTNode
            * a libSBML object that has a ``getMath()`` function
        :param piecewise_to_heaviside:
            If ``True``, piecewise expressions are transformed to Heaviside
            expressions. If ``False``, piecewise expressions are returned as
            is.
        :param bool2num:
            If ``True``, boolean expressions are transformed to numeric
            expressions. If ``False``, boolean expressions are returned as is.
        :raises SBMLException:
            In case of unsupported expressions.
        :returns:
            The sympy expression or None if the input was None or an empty string.
        """
        if var_or_math is None:
            return None

        # numeric values
        if isinstance(var_or_math, float | int):
            return (
                sp.Integer(var_or_math)
                if var_or_math.is_integer()
                else sp.Float(var_or_math)
            )

        def subs_locals(expr: sp.Basic) -> sp.Basic:
            """
            Substitute free symbols to match assumptions of other model
            entities where necessary, and replace parameters by values in case
            of hardcoded parameters.
            """
            return expr.subs(
                {
                    sym: local
                    for sym in expr.free_symbols
                    if (local := self._local_symbols.get(str(sym), sym)) != sym
                }
            )

        # already a sympy object
        if isinstance(var_or_math, sp.Basic):
            ele_name = "SymPy expression"
            expr = subs_locals(var_or_math)

        # an expression string
        elif isinstance(var_or_math, str):
            if not var_or_math:
                return None

            ele_name = "Expression string"
            try:
                expr = sp.sympify(
                    var_or_math,
                    locals={
                        "true": sp.Float(1.0),
                        "false": sp.Float(0.0),
                        **self._local_symbols,
                    },
                )
            except (sp.SympifyError, TypeError, ZeroDivisionError) as err:
                raise SBMLException(
                    f'{ele_name} "{var_or_math}" '
                    "contains an unsupported expression: "
                    f"{err}."
                )

        # actual libSBML objects
        elif isinstance(var_or_math, libsbml.ASTNode | libsbml.SBase):
            if isinstance(var_or_math, libsbml.ASTNode):
                ast_node = var_or_math
                sbml_obj = ast_node.getParentSBMLObject()
            elif isinstance(var_or_math, libsbml.SBase):
                sbml_obj = var_or_math
                ast_node = var_or_math.getMath()
                if ast_node is None:
                    return None
            ele_name = sbml_obj.getId()

            mathml = libsbml.writeMathMLWithNamespaceToString(
                ast_node,
                libsbml.SBMLNamespaces(
                    self.sbml_model.getLevel(), self.sbml_model.getVersion()
                ),
            )
            try:
                expr = self._mathml_parser.parse_str(mathml)
            except (ValueError, NotImplementedError) as err:
                raise SBMLException(
                    f'{ele_name} "{mathml}" '
                    "contains an unsupported expression: "
                    f"{err}."
                )
            # replace special time object by `sbml_time_symbol`
            #  which will later be replaced by `amici_time_symbol`
            expr = expr.replace(TimeSymbol, lambda *args: sbml_time_symbol)
            expr = expr.subs(avogadro, avogadro.evalf())
            expr = subs_locals(expr)
        else:
            raise ValueError(
                f"Unsupported input: {var_or_math}, type: {type(var_or_math)}"
            )

        try:
            _check_unsupported_functions_sbml(expr, expression_type=ele_name)
        except SBMLException:
            # try to (partially) evaluate expressions
            #  that would otherwise be unsupported
            expr = expr.simplify().evalf()
            _check_unsupported_functions_sbml(expr, expression_type=ele_name)

        # boolean to numeric piecewise
        #  (this may introduce piecewise functions and must happen before
        #   piecewise_to_heaviside)
        if bool2num and isinstance(expr, Boolean):
            from sbmlmath.mathml_parser import _bool2num

            expr = _bool2num(expr)

        # piecewise to heavisides
        if piecewise_to_heaviside:
            try:
                expr = expr.replace(sp.Piecewise, _piecewise_to_minmax)
                expr = expr.replace(
                    sp.Piecewise,
                    lambda *args: _parse_piecewise_to_heaviside(args),
                )
            except RuntimeError as err:
                raise SBMLException(str(err)) from err

        return expr

    def _get_element_initial_assignment(
        self, element_id: str
    ) -> sp.Expr | None:
        """
        Extract value of sbml variable according to its initial assignment

        :param element_id:
            sbml variable name
        :return:

        """
        assignment = self.sbml_model.getInitialAssignment(element_id)
        if assignment is None:
            return None
        sym = self._sympify(assignment, piecewise_to_heaviside=False)
        # this is an initial assignment so we need to use
        # initial conditions
        sym = self._make_initial(sym)
        return sym

    def _get_element_stoichiometry(self, ele: libsbml.SBase) -> sp.Expr:
        """
        Computes the stoichiometry of a reactant or product of a reaction

        :param ele:
            reactant or product
        :return:
            _symbolic variable that defines stoichiometry
        """
        if ele.isSetId():
            sym = self._get_element_initial_assignment(ele.getId())
            if sym is not None:
                return sym

            if self._is_assignment_rule_target(ele):
                return _get_identifier_symbol(ele)

        if hasattr(ele, "isSetStoichiometry") and ele.isSetStoichiometry():
            stoichiometry: float = ele.getStoichiometry()
            return (
                sp.Integer(stoichiometry)
                if stoichiometry.is_integer()
                else sp.Float(stoichiometry)
            )

        return sp.Integer(1)

    def _is_assignment_rule_target(self, element: libsbml.SBase) -> bool:
        """
        Checks if an element has a valid assignment rule in the specified
        model.

        :param element:
            SBML variable

        :return:
            boolean indicating truth of function name
        """
        a = self.sbml_model.getAssignmentRuleByVariable(element.getId())
        return a is not None and a.getMath() is not None

    def _is_rate_rule_target(self, element: libsbml.SBase) -> bool:
        """
        Checks if an element has a valid assignment rule in the specified
        model.

        :param element:
            SBML variable

        :return:
            boolean indicating truth of function name
        """
        a = self.sbml_model.getRateRuleByVariable(element.getId())
        return a is not None and a.getMath() is not None

    def _transform_dxdt_to_concentration(
        self, species_id: sp.Symbol, dxdt: sp.Expr
    ) -> sp.Expr:
        """
        Produces the appropriate expression for the first derivative of a
        species with respect to time, for species that reside in
        compartments with a constant volume, or a volume that is defined by
        an assignment or rate rule.

        :param species_id:
            The identifier of the species (generated in "__init__.py").

        :param dxdt:
            The element-wise product of the row in the stoichiometric
            matrix that corresponds to the species (row x_index) and the
            flux (kinetic laws) vector. Ignored in the case of rate rules.
        """
        # The derivation of the below return expressions can be found in
        # the documentation. They are found by rearranging
        # $\frac{d}{dt} (vx) = Sw$ for $\frac{dx}{dt}$, where $v$ is the
        # vector of species compartment volumes, $x$ is the vector of
        # species concentrations, $S$ is the stoichiometric matrix, and $w$
        # is the flux vector. The conditional below handles the cases of
        # species in (i) compartments with a rate rule, (ii) compartments
        # with an assignment rule, and (iii) compartments with a constant
        # volume, respectively.
        species = self._symbols[SymbolId.SPECIES][species_id]

        comp = species["compartment"]
        if comp in self._symbols[SymbolId.SPECIES]:
            dv_dt = self._symbols[SymbolId.SPECIES][comp]["dt"]
            xdot = (dxdt - dv_dt * species_id) / comp
            return xdot
        elif comp in self._compartment_assignment_rules:
            v = self._compartment_assignment_rules[comp]

            # we need to flatten out assignments in the compartment in
            # order to ensure that we catch all species dependencies
            v = smart_subs_dict(v, self._symbols[SymbolId.EXPRESSION], "value")
            dv_dt = v.diff(amici_time_symbol)
            # we may end up with a time derivative of the compartment
            # volume due to parameter rate rules
            comp_rate_vars = [
                p
                for p in v.free_symbols
                if p in self._symbols[SymbolId.SPECIES]
            ]
            for var in comp_rate_vars:
                dv_dt += (
                    v.diff(var) * self._symbols[SymbolId.SPECIES][var]["dt"]
                )
            dv_dx = v.diff(species_id)
            xdot = (dxdt - dv_dt * species_id) / (dv_dx * species_id + v)
            return xdot
        elif comp in self._symbols[SymbolId.ALGEBRAIC_STATE]:
            raise SBMLException(
                f"Species {species_id} is in a compartment {comp} that is"
                f" defined by an algebraic equation. This is not"
                f" supported."
            )
        else:
            v = self._compartments[comp]

            if v == 1.0:
                return dxdt

            return dxdt / v

    def _process_sbml_rate_of(self, de_model: DEModel) -> None:
        """Substitute any SBML-rateOf constructs in the model equations"""
        from sbmlmath import rate_of as rate_of_func

        # DAE models in general are supported, but not if they include rateOf()
        #  The is_ode check below will not trigger in case of a DAE model with
        #  rateOf() only in algebraic equations, so we need to check for that.
        for ae in de_model.algebraic_equations():
            if ae.get_val().find(rate_of_func):
                # Not implemented
                raise SBMLException(
                    "rateOf() is not supported in algebraic equations. "
                    f"Used in: {ae.get_id()}: {ae.get_val()}"
                )

        sym_to_state_var: dict[sp.Symbol, DifferentialState] = {
            state.get_sym(): state for state in de_model.differential_states()
        }

        is_ode = de_model.is_ode()

        def get_rate(symbol: sp.Symbol):
            """Get rate of change of the given symbol
            (symbol is argument of rateOf)"""
            if not is_ode:
                # Not implemented for DAE models
                #  sbml test case 01482
                raise SBMLException(
                    "rateOf() is only supported in ODE models."
                )

            if symbol.find(rate_of_func):
                raise SBMLException("Nesting rateOf() is not allowed.")

            # Replace all rateOf(some_state_var) by their respective
            #  xdot equation
            with contextlib.suppress(KeyError):
                return sym_to_state_var[symbol].get_dt()

            # For anything other than a state variable,
            #  rateOf(.) is 0 or invalid
            return sp.Integer(0)

        def do_subs(expr, rate_ofs) -> sp.Expr:
            """Substitute rateOf(...) in expr by appropriate expressions"""
            expr = expr.subs(
                {rate_of: get_rate(rate_of.args[0]) for rate_of in rate_ofs}
            )
            if rate_ofs := expr.find(rate_of_func):
                # recursively substitute until no rateOf remains
                expr = do_subs(expr, rate_ofs)
            return expr

        # replace rateOf-instances in xdot
        for state in de_model.differential_states():
            if rate_ofs := state.get_dt().find(rate_of_func):
                state.set_dt(do_subs(state.get_dt(), rate_ofs))

        # replace rateOf-instances in expressions which we will need for
        #  substitutions later
        expressions_changed = False
        for expr in de_model.expressions():
            if rate_ofs := expr.get_val().find(rate_of_func):
                expr.set_val(do_subs(expr.get_val(), rate_ofs))
                expressions_changed = True

        # get toposorted `w`, but only changed the order of expressions in the
        #  model if any expression changed.
        # (in principle, we could always reorder the expressions, but this
        #  will require regenerating all test oracle for
        #  `test_pregenerated_models.py`)
        w_toposorted = de_model.toposort_expressions(
            reorder=expressions_changed
        )

        # replace rateOf-instances in x0
        # indices of state variables whose x0 was modified
        changed_indices = []
        for i_state, state in enumerate(de_model.differential_states()):
            new, replacement = state.get_val().replace(
                rate_of_func, get_rate, map=True
            )
            if replacement:
                state.set_val(new)
                changed_indices.append(i_state)
        if changed_indices:
            # Replace any newly introduced state variables
            #  by their x0 expressions.
            subs = w_toposorted | toposort_symbols(
                {
                    state.get_x_rdata(): state.get_val()
                    for state in de_model.differential_states()
                }
            )
            states = de_model.differential_states()
            for i_state in changed_indices:
                states[i_state].set_val(
                    smart_subs_dict(states[i_state].get_val(), subs)
                )

        for component in chain(
            de_model.observables(),
            de_model.events(),
        ):
            if rate_ofs := component.get_val().find(rate_of_func):
                component.set_val(do_subs(component.get_val(), rate_ofs))

            if isinstance(component, Event):
                # for events, also substitute in state updates
                for target, assignment in component._assignments.items():
                    if rate_ofs := assignment.find(rate_of_func):
                        new_assignment = do_subs(assignment, rate_ofs)
                        # currently, deltax cannot depend on `w`.
                        #  this could be changed, or we could use xdot
                        #  symbols, but for now, we just flatten out w
                        #  expressions
                        new_assignment = smart_subs_dict(
                            new_assignment, w_toposorted
                        )
                        component._assignments[target] = new_assignment


def _check_lib_sbml_errors(
    sbml_doc: libsbml.SBMLDocument, show_warnings: bool = False
) -> None:
    """
    Checks the error log in the current self.sbml_doc.

    :param sbml_doc:
        SBML document

    :param show_warnings:
        display SBML warnings
    """
    num_warning = sbml_doc.getNumErrors(libsbml.LIBSBML_SEV_WARNING)
    num_error = sbml_doc.getNumErrors(libsbml.LIBSBML_SEV_ERROR)
    num_fatal = sbml_doc.getNumErrors(libsbml.LIBSBML_SEV_FATAL)

    if num_warning + num_error + num_fatal:
        for i_error in range(sbml_doc.getNumErrors()):
            error = sbml_doc.getError(i_error)
            # we ignore any info messages for now
            if error.getSeverity() >= libsbml.LIBSBML_SEV_ERROR or (
                show_warnings
                and error.getSeverity() >= libsbml.LIBSBML_SEV_WARNING
            ):
                logger.error(
                    f"libSBML {error.getCategoryAsString()} "
                    f"({error.getSeverityAsString()}):"
                    f" {error.getMessage()}"
                )

    if num_error + num_fatal:
        raise SBMLException(
            "SBML Document failed to load (see error messages above)"
        )


def _parse_event_trigger(trigger: sp.Expr) -> sp.Expr:
    """
    Recursively translates a boolean trigger function into a real valued
    root function

    :param trigger: The Boolean trigger expression.
        The event triggers when this expression changes from False to True.
    :return: real-valued root function expression
    """
    # Events can be defined without trigger, i.e., the event will never fire.
    # In this case, set a dummy trigger:
    if trigger is None or trigger is sp.false:
        return sp.Float(-1.0)

    if trigger is sp.true:
        return sp.Float(1.0)

    if trigger.is_Relational:
        root = trigger.args[0] - trigger.args[1]
        _check_unsupported_functions_sbml(root, "sympy.Expression")

        # convert relational expressions into trigger functions
        if isinstance(
            trigger,
            sp.core.relational.LessThan | sp.core.relational.StrictLessThan,
        ):
            # y < x or y <= x
            return -root
        if isinstance(
            trigger,
            sp.core.relational.GreaterThan
            | sp.core.relational.StrictGreaterThan,
        ):
            # y >= x or y > x
            return root

    # rewrite n-ary XOR to OR to be handled below:
    trigger = trigger.replace(sp.Xor, _xor_to_or)
    # rewrite ==, !=
    trigger = trigger.replace(sp.Eq, _eq_to_and)
    trigger = trigger.replace(sp.Ne, _ne_to_or)

    # or(x,y): any of {x,y} is > 0: sp.Max(x, y)
    if isinstance(trigger, sp.Or):
        return sp.Max(*[_parse_event_trigger(arg) for arg in trigger.args])
    # and(x,y): all out of {x,y} are > 0: sp.Min(x, y)
    if isinstance(trigger, sp.And):
        return sp.Min(*[_parse_event_trigger(arg) for arg in trigger.args])

    raise SBMLException(
        "AMICI can not parse piecewise/event trigger functions with argument "
        f"{trigger}."
    )


def assignment_rules_to_observables(
    sbml_model: libsbml.Model,
    filter_function: Callable = lambda *_: True,
    as_dict=False,
) -> list[MeasurementChannel] | dict[str, MeasurementChannel]:
    """
    Turn assignment rules into observables.

    The matching assignment rules and their targets are removed from the model.

    :param sbml_model:
        Model to operate on

    :param filter_function:
        Callback function taking assignment variable as input and returning
        ``True``/``False`` to indicate if the respective rule should be
        turned into an observable.

    :param as_dict:
        If ``True``, return a dict instead of a list, using the observable IDs
        as keys.

    :return:
        A list or dict of measurement channels.
    """
    observables = {}
    for rule in sbml_model.getListOfRules():
        if rule.getTypeCode() != libsbml.SBML_ASSIGNMENT_RULE:
            continue
        target_id = rule.getVariable()
        if (p := sbml_model.getParameter(target_id)) and filter_function(p):
            observables[target_id] = MeasurementChannel(
                id_=target_id,
                name=p.getName() if p.isSetName() else target_id,
                formula=sbml_model.getAssignmentRuleByVariable(
                    target_id
                ).getFormula(),
            )

    for target_id in observables:
        sbml_model.removeRuleByVariable(target_id)
        sbml_model.removeParameter(target_id)

    if as_dict:
        return observables

    return list(observables.values())


def _add_conservation_for_constant_species(
    ode_model: DEModel, conservation_laws: list[ConservationLaw]
) -> list[int]:
    """
    Adds constant species to conservations laws

    :param ode_model:
        DEModel object with basic definitions

    :param conservation_laws:
        List of already known conservation laws

    :returns species_solver:
        List of species indices which remain later in the DE solver
    """

    # decide which species to keep in stoichiometry
    species_solver = list(range(len(ode_model._differential_states)))

    # iterate over species, find constant ones
    for ix in reversed(range(len(ode_model._differential_states))):
        if ode_model.state_is_constant(ix):
            # dont use sym('x') here since conservation laws need to be
            # added before symbols are generated
            target_state = ode_model._differential_states[ix].get_sym()
            total_abundance = symbol_with_assumptions(f"tcl_{target_state}")
            conservation_laws.append(
                {
                    "state": target_state,
                    "total_abundance": total_abundance,
                    "coefficients": {target_state: 1.0},
                }
            )
            # mark species to delete from stoichiometric matrix
            species_solver.pop(ix)

    return species_solver


def _get_species_compartment_symbol(species: libsbml.Species) -> sp.Symbol:
    """
    Generate compartment symbol for the compartment of a specific species.
    This function will always return the same unique python object for a
    given species name.

    :param species:
        sbml species
    :return:
        compartment symbol
    """
    return symbol_with_assumptions(species.getCompartment())


def _get_identifier_symbol(var: libsbml.SBase) -> sp.Symbol:
    """
    Generate identifier symbol for a sbml variable.
    This function will always return the same unique python object for a
    given entity.

    :param var:
        sbml variable
    :return:
        identifier symbol
    """
    return symbol_with_assumptions(var.getId())


def get_species_initial(species: libsbml.Species) -> sp.Expr:
    """
    Extract the initial concentration from a given species

    :param species:
        species index

    :return:
        initial species concentration
    """
    if species.isSetInitialConcentration():
        conc = species.getInitialConcentration()
        if species.getHasOnlySubstanceUnits():
            return sp.Float(conc) * _get_species_compartment_symbol(species)
        else:
            return sp.Float(conc)

    if species.isSetInitialAmount():
        amt = species.getInitialAmount()
        if math.isnan(amt):
            return sp.Float(0.0)

        if species.getHasOnlySubstanceUnits():
            return sp.Float(amt)
        else:
            return sp.Float(amt) / _get_species_compartment_symbol(species)

    return sp.Float(0.0)


def _get_list_of_species_references(
    sbml_model: libsbml.Model,
) -> list[libsbml.SpeciesReference]:
    """
    Extracts list of species references as SBML doesn't provide a native
    function for this.

    :param sbml_model:
        SBML model instance

    :return:
        ListOfSpeciesReferences
    """
    return [
        reference
        for reaction in sbml_model.getListOfReactions()
        for reference in itt.chain(
            reaction.getListOfReactants(),
            reaction.getListOfProducts(),
            reaction.getListOfModifiers(),
        )
    ]


def _collect_event_assignment_parameter_targets(
    sbml_model: libsbml.Model,
) -> list[sp.Symbol]:
    targets = []
    sbml_parameters = sbml_model.getListOfParameters()
    sbml_parameter_ids = [p.getId() for p in sbml_parameters]
    for event in sbml_model.getListOfEvents():
        for event_assignment in event.getListOfEventAssignments():
            target_id = event_assignment.getVariable()
            if target_id in sbml_parameter_ids:
                targets.append(
                    _get_identifier_symbol(
                        sbml_parameters[sbml_parameter_ids.index(target_id)]
                    )
                )
    return targets


def _check_unsupported_functions_sbml(
    sym: sp.Expr, expression_type: str, full_sym: sp.Expr | None = None
):
    try:
        _check_unsupported_functions(sym, expression_type, full_sym)
    except RuntimeError as err:
        raise SBMLException(str(err))


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


def _non_const_conservation_laws_supported(sbml_model: libsbml.Model) -> bool:
    """Check whether non-constant conservation laws can be handled for the
    given model."""
    if any(
        rule.getTypeCode() == libsbml.SBML_RATE_RULE
        for rule in sbml_model.getListOfRules()
    ):
        # see SBML semantic test suite, case 33 for an example
        warnings.warn(
            "Conservation laws for non-constant species in "
            "models with RateRules are currently not supported "
            "and will be turned off.",
            stacklevel=1,
        )
        return False

    if any(
        rule.getTypeCode() == libsbml.SBML_ASSIGNMENT_RULE
        and sbml_model.getSpecies(rule.getVariable())
        for rule in sbml_model.getListOfRules()
    ):
        warnings.warn(
            "Conservation laws for non-constant species in "
            "models with Species-AssignmentRules are currently not "
            "supported and will be turned off.",
            stacklevel=1,
        )
        return False

    return True


def _rateof_to_dummy(
    sym_math: sp.Expr,
) -> tuple[sp.Expr, dict[sp.Function, sp.Dummy]]:
    """Replace rateOf(...) by dummy variable

    if `rateOf(some_species)` is used in an initial assignment, we don't want to substitute the species argument
    by its initial value.

    Usage:
            sym_math, rateof_to_dummy = _rateof_to_dummy(sym_math)
            [...substitute...]
            sym_math = _dummy_to_rateof(sym_math, rateof_to_dummy)
    """
    from sbmlmath import rate_of

    new_expr, rateof_to_dummies = sym_math.replace(
        # `x` is the rateOf argument
        rate_of,
        lambda x: sp.Dummy(f"Dummy_RateOf_{x.name}"),
        map=True,
    )
    return new_expr, rateof_to_dummies


def _dummy_to_rateof(
    sym_math: sp.Expr, rateof_dummies: dict[sp.Function, sp.Dummy]
) -> sp.Expr:
    """Back-substitution of dummies from `_rateof_to_dummy`"""
    if rateof_dummies:
        return sym_math.subs({v: k for k, v in rateof_dummies.items()})
    return sym_math


def _try_evalf(x: sp.Basic) -> sp.Basic:
    """Try to convert a sympy object to sympy.Float.

    :returns: The value of `x` as sympy.Float if there is a known conversion;
        `x` otherwise.
    """
    if hasattr(x, "evalf") and isinstance((x := x.evalf()), sp.Float):
        return x

    try:
        if (
            x := x.subs(
                {
                    BooleanTrue(): sp.Float(1.0),
                    BooleanFalse(): sp.Float(0.0),
                }
            )
        ) and isinstance(x, sp.Float):
            return x
    except TypeError:
        pass

    return x
