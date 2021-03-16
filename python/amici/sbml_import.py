"""
SBML Import
-----------
This module provides all necessary functionality to import a model specified
in the `Systems Biology Markup Language (SBML) <http://sbml.org/Main_Page>`_.
"""


import sympy as sp
import libsbml as sbml
import re
import math
import itertools as itt
import warnings
import logging
import copy
from typing import (
    Dict, List, Callable, Any, Iterable, Union, Optional, Tuple, Sequence
)

from .import_utils import (
    smart_subs, smart_subs_dict, toposort_symbols,
    _get_str_symbol_identifiers,
    noise_distribution_to_cost_function
)
from .ode_export import (
    ODEExporter, ODEModel, generate_measurement_symbol,
    symbol_with_assumptions
)
from .constants import SymbolId
from .logging import get_logger, log_execution_time, set_log_level
from . import has_clibs

from sympy.logic.boolalg import BooleanAtom


class SBMLException(Exception):
    pass


SymbolicFormula = Dict[sp.Symbol, sp.Expr]


default_symbols = {
    symbol: {} for symbol in SymbolId
}

ConservationLaw = Dict[str, Union[str, sp.Expr]]

logger = get_logger(__name__, logging.ERROR)


class SbmlImporter:
    """
    Class to generate AMICI C++ files for a model provided in the Systems
    Biology Markup Language (SBML).

    :ivar show_sbml_warnings:
        indicates whether libSBML warnings should be
        displayed

    :ivar symbols:
        dict carrying symbolic definitions

    :ivar sbml_reader:

        The libSBML sbml reader

        .. warning::
           Not storing this may result in a segfault.

    :ivar sbml_doc:
        document carrying the sbml definition

        .. warning::
           Not storing this may result in a segfault.

    :ivar sbml:
        SBML model to import

    :ivar compartments:
        dict of compartment ids and compartment volumes

    :ivar stoichiometric_matrix:
        stoichiometric matrix of the model

    :ivar flux_vector:
        reaction kinetic laws

    :ivar _local_symbols:
        model symbols for sympy to consider during sympification
        see `locals`argument in `sympy.sympify`

    :ivar species_assignment_rules:
        Assignment rules for species.
        Key is symbolic identifier and value is assignment value

    :ivar compartment_assignment_rules:
        Assignment rules for compartments.
        Key is symbolic identifier and value is assignment value

    :ivar parameter_assignment_rules:
        assignment rules for parameters, these parameters are not permissible
        for sensitivity analysis

    :ivar initial_assignments:
        initial assignments for parameters, these parameters are not
        permissible for sensitivity analysis

    :ivar sbml_parser_settings:
        sets behaviour of SBML Formula parsing

    """

    def __init__(self,
                 sbml_source: Union[str, sbml.Model],
                 show_sbml_warnings: bool = False,
                 from_file: bool = True) -> None:
        """
        Create a new Model instance.

        :param sbml_source:
            Either a path to SBML file where the model is specified,
            or a model string as created by sbml.sbmlWriter(
            ).writeSBMLToString() or an instance of `libsbml.Model`.

        :param show_sbml_warnings:
            Indicates whether libSBML warnings should be displayed.

        :param from_file:
            Whether `sbml_source` is a file name (True, default), or an SBML
            string
        """
        if isinstance(sbml_source, sbml.Model):
            self.sbml_doc: sbml.Document = sbml_source.getSBMLDocument()
        else:
            self.sbml_reader: sbml.SBMLReader = sbml.SBMLReader()
            if from_file:
                sbml_doc = self.sbml_reader.readSBMLFromFile(sbml_source)
            else:
                sbml_doc = self.sbml_reader.readSBMLFromString(sbml_source)
            self.sbml_doc = sbml_doc

        self.show_sbml_warnings: bool = show_sbml_warnings

        # process document
        self._process_document()

        self.sbml: sbml.Model = self.sbml_doc.getModel()

        # Long and short names for model components
        self.symbols: Dict[SymbolId, Dict[sp.Symbol, Dict[str, Any]]] = {}

        self._local_symbols: Dict[str, Union[sp.Expr, sp.Function]] = {}
        self.compartments: SymbolicFormula = {}
        self.compartment_assignment_rules: SymbolicFormula = {}
        self.species_assignment_rules: SymbolicFormula = {}
        self.parameter_assignment_rules: SymbolicFormula = {}
        self.initial_assignments: SymbolicFormula = {}

        self._reset_symbols()

        # http://sbml.org/Software/libSBML/5.18.0/docs/python-api/classlibsbml_1_1_l3_parser_settings.html#abcfedd34efd3cae2081ba8f42ea43f52
        # all defaults except disable unit parsing
        self.sbml_parser_settings = sbml.L3ParserSettings(
            self.sbml, sbml.L3P_PARSE_LOG_AS_LOG10,
            sbml.L3P_EXPAND_UNARY_MINUS, sbml.L3P_NO_UNITS,
            sbml.L3P_AVOGADRO_IS_CSYMBOL,
            sbml.L3P_COMPARE_BUILTINS_CASE_INSENSITIVE, None,
            sbml.L3P_MODULO_IS_PIECEWISE
        )

    def _process_document(self) -> None:
        """
        Validate and simplify document.
        """
        # Ensure we got a valid SBML model, otherwise further processing
        # might lead to undefined results
        self.sbml_doc.validateSBML()
        _check_lib_sbml_errors(self.sbml_doc, self.show_sbml_warnings)

        # apply several model simplifications that make our life substantially
        # easier
        if self.sbml_doc.getModel().getNumFunctionDefinitions():
            convert_config = sbml.SBMLFunctionDefinitionConverter()\
                .getDefaultProperties()
            self.sbml_doc.convert(convert_config)

        convert_config = sbml.SBMLLocalParameterConverter().\
            getDefaultProperties()
        self.sbml_doc.convert(convert_config)

        # If any of the above calls produces an error, this will be added to
        # the SBMLError log in the sbml document. Thus, it is sufficient to
        # check the error log just once after all conversion/validation calls.
        _check_lib_sbml_errors(self.sbml_doc, self.show_sbml_warnings)

    def _reset_symbols(self) -> None:
        """
        Reset the symbols attribute to default values
        """
        self.symbols = copy.deepcopy(default_symbols)
        self._local_symbols = dict()

    def sbml2amici(self,
                   model_name: str = None,
                   output_dir: str = None,
                   observables: Dict[str, Dict[str, str]] = None,
                   constant_parameters: Iterable[str] = None,
                   sigmas: Dict[str, Union[str, float]] = None,
                   noise_distributions: Dict[str, Union[str, Callable]] = None,
                   verbose: Union[int, bool] = logging.ERROR,
                   assume_pow_positivity: bool = False,
                   compiler: str = None,
                   allow_reinit_fixpar_initcond: bool = True,
                   compile: bool = True,
                   compute_conservation_laws: bool = True,
                   simplify: Callable = lambda x: sp.powsimp(x, deep=True),
                   log_as_log10: bool = True,
                   generate_sensitivity_code: bool = True,
                   **kwargs) -> None:
        """
        Generate and compile AMICI C++ files for the model provided to the
        constructor.

        The resulting model can be imported as a regular Python module (if
        `compile=True`), or used from Matlab or C++ as described in the
        documentation of the respective AMICI interface.

        Note that this generates model ODEs for changes in concentrations, not
        amounts unless the `hasOnlySubstanceUnits` attribute has been
        defined for a particular species.

        Sensitivity analysis for local parameters is enabled by creating
        global parameters _{reactionId}_{localParameterName}.

        :param model_name:
            name of the model/model directory

        :param output_dir:
            see :meth:`amici.ode_export.ODEExporter.set_paths`

        :param observables:
            dictionary( observableId:{'name':observableName
            (optional), 'formula':formulaString)}) to be added to the model

        :param constant_parameters:
            list of SBML Ids identifying constant parameters

        :param sigmas:
            dictionary(observableId: sigma value or (existing) parameter name)

        :param noise_distributions:
            dictionary(observableId: noise type).
            If nothing is passed for some observable id, a normal model is
            assumed as default. Either pass a noise type identifier, or a
            callable generating a custom noise string.

        :param verbose:
            verbosity level for logging, ``True``/``False`` default to
            ``logging.Error``/``logging.DEBUG``

        :param assume_pow_positivity:
            if set to ``True``, a special pow function is
            used to avoid problems with state variables that may become
            negative due to numerical errors

        :param compiler:
            distutils/setuptools compiler selection to build the
            python extension

        :param allow_reinit_fixpar_initcond:
            see :class:`amici.ode_export.ODEExporter`

        :param compile:
            If ``True``, compile the generated Python package,
            if ``False``, just generate code.

        :param compute_conservation_laws:
            if set to ``True``, conservation laws are automatically computed
            and applied such that the state-jacobian of the ODE
            right-hand-side has full rank. This option should be set to
            ``True`` when using the newton algorithm to compute steadystate
            sensitivities.

        :param simplify:
            see :attr:`ODEModel._simplify`

        :param log_as_log10:
            If ``True``, log in the SBML model will be parsed as ``log10``
            (default), if ``False``, log will be parsed as natural logarithm
            ``ln``

        :param generate_sensitivity_code:
            If ``False``, the code required for sensitivity computation will
            not be generated

        """
        set_log_level(logger, verbose)

        if 'constantParameters' in kwargs:
            logger.warning('Use of `constantParameters` as argument name '
                           'is deprecated and will be removed in a future '
                           'version. Please use `constant_parameters` as '
                           'argument name.')

            if constant_parameters is not None:
                raise ValueError('Cannot specify constant parameters using '
                                 'both `constantParameters` and '
                                 '`constant_parameters` as argument names.')

            constant_parameters = kwargs.pop('constantParameters', [])

        elif constant_parameters is None:
            constant_parameters = []
        constant_parameters = list(constant_parameters)

        if sigmas is None:
            sigmas = {}

        if noise_distributions is None:
            noise_distributions = {}

        if model_name is None:
            model_name = kwargs.pop('modelName', None)
            if model_name is None:
                raise ValueError('Missing argument: `model_name`')
            else:
                logger.warning('Use of `modelName` as argument name is '
                               'deprecated and will be removed in a future'
                               ' version. Please use `model_name` as '
                               'argument name.')
        else:
            if 'modelName' in kwargs:
                raise ValueError('Cannot specify model name using both '
                                 '`modelName` and `model_name` as argument '
                                 'names.')

        if len(kwargs):
            raise ValueError(f'Unknown arguments {kwargs.keys()}.')

        self._reset_symbols()
        self.sbml_parser_settings.setParseLog(
            sbml.L3P_PARSE_LOG_AS_LOG10 if log_as_log10 else
            sbml.L3P_PARSE_LOG_AS_LN
        )
        self._process_sbml(constant_parameters)
        if self.symbols.get(SymbolId.EVENT, False):
            if compute_conservation_laws:
                logger.warning(
                    'Conservation laws are currently not supported for models '
                    'with events, and will be turned off.'
                )
            compute_conservation_laws = False

        self._process_observables(observables, sigmas, noise_distributions)
        self._replace_compartments_with_volumes()

        self._clean_reserved_symbols()
        self._process_time()

        ode_model = ODEModel(verbose=verbose, simplify=simplify)
        ode_model.import_from_sbml_importer(
            self, compute_cls=compute_conservation_laws)
        exporter = ODEExporter(
            ode_model,
            outdir=output_dir,
            verbose=verbose,
            assume_pow_positivity=assume_pow_positivity,
            compiler=compiler,
            allow_reinit_fixpar_initcond=allow_reinit_fixpar_initcond,
            generate_sensitivity_code=generate_sensitivity_code
        )
        exporter.set_name(model_name)
        exporter.set_paths(output_dir)
        exporter.generate_model_code()

        if compile:
            if not has_clibs:
                warnings.warn('AMICI C++ extensions have not been built. '
                              'Generated model code, but unable to compile.')
            exporter.compile_model()

    @log_execution_time('importing SBML', logger)
    def _process_sbml(self, constant_parameters: List[str] = None) -> None:
        """
        Read parameters, species, reactions, and so on from SBML model

        :param constant_parameters:
            SBML Ids identifying constant parameters
        """
        self.check_support()
        self._gather_locals()
        self._process_parameters(constant_parameters)
        self._process_compartments()
        self._process_species()
        self._process_reactions()
        self._process_rules()
        self._process_initial_assignments()
        self._process_species_references()
        self._process_events()

    def check_support(self) -> None:
        """
        Check whether all required SBML features are supported.
        Also ensures that the SBML contains at least one reaction, or rate
        rule, or assignment rule, to produce change in the system over time.
        """
        if hasattr(self.sbml, 'all_elements_from_plugins') \
                and self.sbml.all_elements_from_plugins.getSize():
            raise SBMLException('SBML extensions are currently not supported!')

        if any([not rule.isAssignment() and not isinstance(
                    self.sbml.getElementBySId(rule.getVariable()),
                    (sbml.Compartment, sbml.Species, sbml.Parameter)
                ) for rule in self.sbml.getListOfRules()]):
            raise SBMLException('Algebraic rules are currently not supported, '
                                'and rate rules are only supported for '
                                'species, compartments, and parameters.')

        if any([not (rule.isAssignment() or rule.isRate())
                and isinstance(
                    self.sbml.getElementBySId(rule.getVariable()),
                    (sbml.Compartment, sbml.Species, sbml.Parameter)
                ) for rule in self.sbml.getListOfRules()]):
            raise SBMLException('Only assignment and rate rules are '
                                'currently supported for compartments, '
                                'species, and parameters!')

        if any([r.getFast() for r in self.sbml.getListOfReactions()]):
            raise SBMLException('Fast reactions are currently not supported!')

        # Check events for unsupported functionality
        self.check_event_support()

    def check_event_support(self) -> None:
        """
        Check possible events in the model, as AMICI does currently not support

        * delays in events
        * priorities of events
        * events fired at initial time

        Furthermore, event triggers are optional (e.g., if an event is fired at
        initial time, no trigger function is necessary).
        In this case, warn that this event will have no effect.
        """
        for event in self.sbml.getListOfEvents():
            event_id = event.getId()
            # Check for delays in events
            delay = event.getDelay()
            if delay is not None:
                try:
                    delay_time = float(self._sympy_from_sbml_math(delay))
                    if delay_time != 0:
                        raise ValueError
                # `TypeError` would be raised in the above `float(...)`
                # if the delay is not a fixed time
                except (TypeError, ValueError):
                    raise SBMLException('Events with execution delays are '
                                        'currently not supported in AMICI.')
            # Check for priorities
            if event.getPriority() is not None:
                raise SBMLException(f'Event {event_id} has a priority '
                                    'specified. This is currently not '
                                    'supported in AMICI.')

            # check trigger
            trigger_sbml = event.getTrigger()
            if trigger_sbml is None:
                logger.warning(f'Event {event_id} trigger has no trigger, '
                               'so will be skipped.')
                continue
            if trigger_sbml.getMath() is None:
                logger.warning(f'Event {event_id} trigger has no trigger '
                               'expression, so a dummy trigger will be set.')
            if not trigger_sbml.getInitialValue():
                # True: event not executed if triggered at time == 0
                # (corresponding to AMICI default). Raise if set to False.
                raise SBMLException(
                    f'Event {event_id} has a trigger that has an initial '
                    'value of False. This is currently not supported in AMICI.'
                )

            if not trigger_sbml.getPersistent():
                raise SBMLException(
                    f'Event {event_id} has a non-persistent trigger.'
                    'This is currently not supported in AMICI.'
                )

    @log_execution_time('gathering local SBML symbols', logger)
    def _gather_locals(self) -> None:
        """
        Populate self.local_symbols with all model entities.

        This is later used during sympifications to avoid sympy builtins
        shadowing model entities as well as to avoid possibly costly
        symbolic substitutions
        """
        self._gather_base_locals()
        self._gather_dependent_locals()

    def _gather_base_locals(self):
        """
        Populate self.local_symbols with pure symbol definitions that do not
        depend on any other symbol.
        """

        special_symbols_and_funs = {
            # oo is sympy infinity
            'INF': sp.oo,
            'NaN': sp.nan,
            'rem': sp.Mod,
            'time': symbol_with_assumptions('time'),
            # SBML L3 explicitly defines this value, which is not equal
            # to the most recent SI definition.
            'avogadro': sp.Float(6.02214179e23),
            'exponentiale': sp.E,
        }
        for s, v in special_symbols_and_funs.items():
            self.add_local_symbol(s, v)

        for c in itt.chain(self.sbml.getListOfSpecies(),
                           self.sbml.getListOfParameters(),
                           self.sbml.getListOfCompartments()):
            if not c.isSetId():
                continue

            self.add_local_symbol(c.getId(), _get_identifier_symbol(c))

        for x_ref in _get_list_of_species_references(self.sbml):
            if not x_ref.isSetId():
                continue
            if x_ref.isSetStoichiometry() and not \
                    self.is_assignment_rule_target(x_ref):
                value = sp.Float(x_ref.getStoichiometry())
            else:
                value = _get_identifier_symbol(x_ref)

            ia_sym = self._get_element_initial_assignment(x_ref.getId())
            if ia_sym is not None:
                value = ia_sym

            self.add_local_symbol(x_ref.getId(), value)

        for r in self.sbml.getListOfReactions():
            for e in itt.chain(r.getListOfReactants(), r.getListOfProducts()):
                if isinstance(e, sbml.SpeciesReference):
                    continue

                if not (e.isSetId() and e.isSetStoichiometry()) or \
                        self.is_assignment_rule_target(e):
                    continue

                self.add_local_symbol(e.getId(),
                                      sp.Float(e.getStoichiometry()))

    def _gather_dependent_locals(self):
        """
        Populate self.local_symbols with symbol definitions that may depend on
        other symbol definitions.
        """
        for r in self.sbml.getListOfReactions():
            if not r.isSetId():
                continue
            self.add_local_symbol(
                r.getId(),
                self._sympy_from_sbml_math(r.getKineticLaw())
            )

    def add_local_symbol(self, key: str, value: sp.Expr):
        """
        Add local symbols with some sanity checking for duplication which
        would indicate redefinition of internals, which SBML permits,
        but we don't.

        :param key:
            local symbol key

        :param value:
            local symbol value
        """
        if key in itt.chain(self._local_symbols.keys(),
                            ['True', 'False', 'pi']):
            raise SBMLException(
                f'AMICI tried to add a local symbol {key} with value {value}, '
                f'but {key} was already instantiated with '
                f'{self._local_symbols[key]}. This means either that there '
                f'are multiple SBML element with SId {key}, which is '
                f'invalid SBML, or that the employed SId {key} is a special '
                'reserved symbol in AMICI. This can be fixed by renaming '
                f'the element with SId {key}.'
            )
        self._local_symbols[key] = value

    @log_execution_time('processing SBML compartments', logger)
    def _process_compartments(self) -> None:
        """
        Get compartment information, stoichiometric matrix and fluxes from
        SBML model.
        """
        compartments = self.sbml.getListOfCompartments()
        self.compartments = {}
        for comp in compartments:
            init = sp.Float(1.0)

            if comp.isSetVolume():
                init = self._sympy_from_sbml_math(comp.getVolume())

            ia_sym = self._get_element_initial_assignment(comp.getId())
            if ia_sym is not None:
                init = ia_sym

            self.compartments[_get_identifier_symbol(comp)] = init

    @log_execution_time('processing SBML species', logger)
    def _process_species(self) -> None:
        """
        Get species information from SBML model.
        """
        if self.sbml.isSetConversionFactor():
            conversion_factor = symbol_with_assumptions(
                self.sbml.getConversionFactor()
            )
        else:
            conversion_factor = 1.0

        for s in self.sbml.getListOfSpecies():
            if self.is_assignment_rule_target(s):
                continue
            self.symbols[SymbolId.SPECIES][_get_identifier_symbol(s)] = {
                'name': s.getName() if s.isSetName() else s.getId(),
                'compartment': _get_species_compartment_symbol(s),
                'constant': s.getConstant() or s.getBoundaryCondition(),
                'amount': s.getHasOnlySubstanceUnits(),
                'conversion_factor': symbol_with_assumptions(
                    s.getConversionFactor()
                )
                if s.isSetConversionFactor()
                else conversion_factor,
                'index': len(self.symbols[SymbolId.SPECIES]),
            }

        self._convert_event_assignment_parameter_targets_to_species()
        self._process_species_initial()
        self._process_rate_rules()

    @log_execution_time('processing SBML species initials', logger)
    def _process_species_initial(self):
        """
        Extract initial values and initial assignments from species
        """
        for species_variable in self.sbml.getListOfSpecies():
            initial = get_species_initial(species_variable)

            species_id = _get_identifier_symbol(species_variable)
            # If species_id is a target of an AssignmentRule, species will be
            # None, but we don't have to account for the initial definition
            # of the species itself and SBML doesn't permit AssignmentRule
            # targets to have InitialAssignments.
            species = self.symbols[SymbolId.SPECIES].get(species_id, None)

            ia_initial = self._get_element_initial_assignment(
                species_variable.getId()
            )
            if ia_initial is not None:
                if species and species['amount'] \
                        and 'compartment' in species:
                    ia_initial *= self.compartments.get(
                        species['compartment'], species['compartment']
                    )
                initial = ia_initial
            if species:
                species['init'] = initial

        # don't assign this since they need to stay in order
        sorted_species = toposort_symbols(self.symbols[SymbolId.SPECIES],
                                          'init')
        for species in self.symbols[SymbolId.SPECIES].values():
            species['init'] = smart_subs_dict(species['init'],
                                              sorted_species,
                                              'init')

    @log_execution_time('processing SBML rate rules', logger)
    def _process_rate_rules(self):
        """
        Process rate rules for species, compartments and parameters.
        Compartments and parameters with rate rules are implemented as species.
        Note that, in the case of species, rate rules may describe the change
        in amount, not concentration, of a species.
        """
        rules = self.sbml.getListOfRules()
        # compartments with rules are replaced with constants in the relevant
        # equations during the _replace_in_all_expressions call inside
        # _process_rules
        for rule in rules:
            if rule.getTypeCode() != sbml.SBML_RATE_RULE:
                continue

            variable = symbol_with_assumptions(rule.getVariable())
            formula = self._sympy_from_sbml_math(rule)
            if formula is None:
                continue

            # Species rules are processed first, to avoid processing
            # compartments twice (as compartments with rate rules are
            # implemented as species).
            ia_init = self._get_element_initial_assignment(rule.getVariable())
            if variable in self.symbols[SymbolId.SPECIES]:
                init = self.symbols[SymbolId.SPECIES][variable]['init']
                name = None

            if variable in self.compartments:
                init = self.compartments[variable]
                name = str(variable)
                del self.compartments[variable]

            elif variable in self.symbols[SymbolId.PARAMETER]:
                init = self._sympy_from_sbml_math(
                    self.symbols[SymbolId.PARAMETER][variable]['value'],
                )
                name = self.symbols[SymbolId.PARAMETER][variable]['name']
                del self.symbols[SymbolId.PARAMETER][variable]

            # parameter with initial assignment, cannot use
            # self.initial_assignments as it is not filled at this
            # point
            elif ia_init is not None:
                init = ia_init
                par = self.sbml.getElementBySId(rule.getVariable())
                name = par.getName() if par.isSetName() else par.getId()

            self.add_d_dt(formula, variable, init, name)

    def add_d_dt(
            self,
            d_dt: sp.Expr,
            variable: sp.Symbol,
            variable0: Union[float, sp.Expr],
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
        if variable in self.symbols[SymbolId.SPECIES]:
            # only update dt if species was already generated
            self.symbols[SymbolId.SPECIES][variable]['dt'] = d_dt
        else:
            # update initial values
            for species_id, species in self.symbols[SymbolId.SPECIES].items():
                variable0 = smart_subs(variable0, species_id, species['init'])

            for species in self.symbols[SymbolId.SPECIES].values():
                species['init'] = smart_subs(species['init'],
                                             variable, variable0)

            # add compartment/parameter species
            self.symbols[SymbolId.SPECIES][variable] = {
                'name': name,
                'init': variable0,
                'amount': False,
                'conversion_factor': 1.0,
                'constant': False,
                'index': len(self.symbols[SymbolId.SPECIES]),
                'dt': d_dt,
            }

    @log_execution_time('processing SBML parameters', logger)
    def _process_parameters(self,
                            constant_parameters: List[str] = None) -> None:
        """
        Get parameter information from SBML model.

        :param constant_parameters:
            SBML Ids identifying constant parameters
        """

        if constant_parameters is None:
            constant_parameters = []

        # Ensure specified constant parameters exist in the model
        for parameter in constant_parameters:
            if not self.sbml.getParameter(parameter):
                raise KeyError('Cannot make %s a constant parameter: '
                               'Parameter does not exist.' % parameter)

        fixed_parameters = [
            parameter
            for parameter in self.sbml.getListOfParameters()
            if parameter.getId() in constant_parameters
        ]
        for parameter in fixed_parameters:
            if self._get_element_initial_assignment(parameter.getId()) is not \
                    None or self.is_assignment_rule_target(parameter) or \
                    self.is_rate_rule_target(parameter):
                raise SBMLException(
                    f'Cannot turn parameter {parameter.getId()} into a '
                    'constant/fixed parameter since it either has an '
                    'initial assignment or is the target of an assignment or '
                    'rate rule.'
                )

        parameters = [
            parameter for parameter
            in self.sbml.getListOfParameters()
            if parameter.getId() not in constant_parameters
            and self._get_element_initial_assignment(parameter.getId()) is None
            and not self.is_assignment_rule_target(parameter)
        ]

        loop_settings = {
            SymbolId.PARAMETER: {'var': parameters, 'name': 'parameter'},
            SymbolId.FIXED_PARAMETER: {'var': fixed_parameters,
                                       'name': 'fixed_parameter'}
        }

        for partype, settings in loop_settings.items():
            for par in settings['var']:
                self.symbols[partype][_get_identifier_symbol(par)] = {
                    'name': par.getName() if par.isSetName() else par.getId(),
                    'value': par.getValue()
                }

    @log_execution_time('processing SBML reactions', logger)
    def _process_reactions(self):
        """
        Get reactions from SBML model.
        """
        reactions = self.sbml.getListOfReactions()
        # nr (number of reactions) should have a minimum length of 1. This is
        # to ensure that, if there are no reactions, the stoichiometric matrix
        # and flux vector multiply to a zero vector with dimensions (nx, 1).
        nr = max(1, len(reactions))
        nx = len(self.symbols[SymbolId.SPECIES])
        # stoichiometric matrix
        self.stoichiometric_matrix = sp.SparseMatrix(sp.zeros(nx, nr))
        self.flux_vector = sp.zeros(nr, 1)

        reaction_ids = [
            reaction.getId() for reaction in reactions
            if reaction.isSetId()
        ]

        for reaction_index, reaction in enumerate(reactions):
            for element_list, sign in [(reaction.getListOfReactants(), -1.0),
                                       (reaction.getListOfProducts(), 1.0)]:
                for element in element_list:
                    stoichiometry = self._get_element_stoichiometry(
                        element
                    )
                    sbml_species = self.sbml.getSpecies(element.getSpecies())
                    if self.is_assignment_rule_target(sbml_species):
                        continue
                    species_id = _get_identifier_symbol(sbml_species)
                    species = self.symbols[SymbolId.SPECIES][species_id]

                    if species['constant']:
                        continue

                    # Division by species compartment size (to find the
                    # rate of change in species concentration) now occurs
                    # in the `dx_dt` method in "ode_export.py", which also
                    # accounts for possibly variable compartments.
                    self.stoichiometric_matrix[species['index'],
                                               reaction_index] += \
                        sign * stoichiometry * species['conversion_factor']

            if reaction.isSetId():
                sym_math = self._local_symbols[reaction.getId()]
            else:
                sym_math = self._sympy_from_sbml_math(reaction.getKineticLaw())

            self.flux_vector[reaction_index] = sym_math
            if any([
                str(symbol) in reaction_ids
                for symbol in self.flux_vector[reaction_index].free_symbols
            ]):
                raise SBMLException(
                    'Kinetic laws involving reaction ids are currently'
                    ' not supported!'
                )

    @log_execution_time('processing SBML rules', logger)
    def _process_rules(self) -> None:
        """
        Process Rules defined in the SBML model.
        """
        for rule in self.sbml.getListOfRules():
            # rate rules are processed in _process_species
            if rule.getTypeCode() == sbml.SBML_RATE_RULE:
                continue

            sbml_var = self.sbml.getElementBySId(rule.getVariable())
            sym_id = symbol_with_assumptions(rule.getVariable())
            formula = self._sympy_from_sbml_math(rule)
            if formula is None:
                continue

            if isinstance(sbml_var, sbml.Species):
                self.species_assignment_rules[sym_id] = formula

            elif isinstance(sbml_var, sbml.Compartment):
                self.compartment_assignment_rules[sym_id] = formula
                self.compartments[sym_id] = formula

            elif isinstance(sbml_var, sbml.Parameter):
                self.parameter_assignment_rules[sym_id] = formula

            self.symbols[SymbolId.EXPRESSION][sym_id] = {
                'name': str(sym_id),
                'value': formula
            }

        self.symbols[SymbolId.EXPRESSION] = toposort_symbols(
            self.symbols[SymbolId.EXPRESSION], 'value'
        )

        # expressions must not occur in definition of x0
        for species in self.symbols[SymbolId.SPECIES].values():
            species['init'] = self._make_initial(
                smart_subs_dict(species['init'],
                                self.symbols[SymbolId.EXPRESSION],
                                'value')
            )

    def _process_time(self) -> None:
        """
        Convert time_symbol into cpp variable.
        """
        sbml_time_symbol = symbol_with_assumptions('time')
        amici_time_symbol = symbol_with_assumptions('t')
        self.amici_time_symbol = amici_time_symbol

        self._replace_in_all_expressions(sbml_time_symbol, amici_time_symbol)

    def _convert_event_assignment_parameter_targets_to_species(self):
        """
        Convert parameters that are targets of event assignments to species.

        This is for the convenience of only implementing event assignments for
        "species".
        """
        parameter_targets = \
            _collect_event_assignment_parameter_targets(self.sbml)
        for parameter_target in parameter_targets:
            # Parameter rate rules already exist as species.
            if parameter_target in self.symbols[SymbolId.SPECIES]:
                continue
            if parameter_target in self.parameter_assignment_rules:
                raise SBMLException(
                    'AMICI does not currently support models with SBML events '
                    'that affect parameters that are also the target of '
                    'assignment rules.'
                )
            parameter_def = None
            for symbol_id in {SymbolId.PARAMETER, SymbolId.FIXED_PARAMETER}:
                if parameter_target in self.symbols[symbol_id]:
                    # `parameter_target` should only exist in one of the
                    # `symbol_id` dictionaries.
                    if parameter_def is not None:
                        raise AssertionError(
                            'Unexpected error. The parameter target of an '
                            'event assignment was processed twice.'
                        )
                    parameter_def = \
                        self.symbols[symbol_id].pop(parameter_target)
            if parameter_def is None:
                raise AssertionError(
                    'Unexpected error. The parameter target of an event '
                    'assignment could not be found.'
                )
            # Fixed parameters are added as species such that they can be
            # targets of events.
            self.symbols[SymbolId.SPECIES][parameter_target] = {
                'name': parameter_def['name'],
                'init': sp.Float(parameter_def['value']),
                # 'compartment': None,  # can ignore for amounts
                'constant': False,
                'amount': True,
                # 'conversion_factor': 1.0,  # can be ignored
                'index': len(self.symbols[SymbolId.SPECIES]),
                'dt': sp.Float(0),
            }

    @log_execution_time('processing SBML events', logger)
    def _process_events(self) -> None:
        """Process SBML events."""
        events = self.sbml.getListOfEvents()

        def get_empty_bolus_value() -> sp.Float:
            """
            Used in the event update vector for species that are not affected
            by the event.
            """
            return sp.Symbol('AMICI_EMTPY_BOLUS')

        # Used to update species concentrations when an event affects a
        # compartment.
        concentration_species_by_compartment = {
            symbol_with_assumptions(c.getId()): []
            for c in self.sbml.getListOfCompartments()
        }
        for species, species_def in self.symbols[SymbolId.SPECIES].items():
            if (
                    # Species is a concentration
                    not species_def.get('amount', True) and
                    # Species has a compartment
                    'compartment' in species_def
            ):
                concentration_species_by_compartment[
                    species_def['compartment']
                ].append(species)

        for ievent, event in enumerate(events):
            # get the event id (which is optional unfortunately)
            event_id = event.getId()
            if event_id is None or event_id == '':
                event_id = f'event_{ievent}'
            event_sym = sp.Symbol(event_id)

            # get and parse the trigger function
            trigger_sbml = event.getTrigger()
            trigger_sym = self._sympy_from_sbml_math(trigger_sbml)
            trigger = _parse_event_trigger(trigger_sym)

            # Currently, all event assignment targets must exist in
            # self.symbols[SymbolId.SPECIES]
            state_vector = list(self.symbols[SymbolId.SPECIES].keys())

            # parse the boluses / event assignments
            bolus = [get_empty_bolus_value() for _ in state_vector]
            event_assignments = event.getListOfEventAssignments()
            compartment_event_assignments = set()
            for event_assignment in event_assignments:
                variable_sym = \
                    symbol_with_assumptions(event_assignment.getVariable())
                if event_assignment.getMath() is None:
                    # Ignore event assignments with no change in value.
                    continue
                formula = self._sympy_from_sbml_math(event_assignment)
                try:
                    # Try to find the species in the state vector.
                    index = state_vector.index(variable_sym)
                    bolus[index] = formula
                except ValueError:
                    raise SBMLException(
                        'Could not process event assignment for '
                        f'{str(variable_sym)}. AMICI currently only allows '
                        'event assignments to species; parameters; or, '
                        'compartments with rate rules, at the moment.'
                    )
                try:
                    # Try working with the formula now to detect errors
                    # here instead of at multiple points downstream.
                    _ = formula - variable_sym
                except TypeError:
                    raise SBMLException(
                        'Could not process event assignment for '
                        f'{str(variable_sym)}. AMICI only allows symbolic '
                        'expressions as event assignments.'
                    )
                if variable_sym in concentration_species_by_compartment:
                    compartment_event_assignments.add(variable_sym)

                for comp, assignment in \
                        self.compartment_assignment_rules.items():
                    if variable_sym not in assignment.free_symbols:
                        continue
                    compartment_event_assignments.add(comp)

            # Update the concentration of species with concentration units
            # in compartments that were affected by the event assignments.
            for compartment_sym in compartment_event_assignments:
                for species_sym in concentration_species_by_compartment[
                        compartment_sym
                ]:
                    # If the species was not affected by an event assignment
                    # then the old value should be updated.
                    if (
                            bolus[state_vector.index(species_sym)]
                            == get_empty_bolus_value()
                    ):
                        species_value = species_sym
                    # else the species was affected by an event assignment,
                    # hence the updated value should be updated further.
                    else:
                        species_value = bolus[state_vector.index(species_sym)]
                    # New species value is old amount / new volume.
                    bolus[state_vector.index(species_sym)] = (
                        species_value * compartment_sym / formula
                    )

            # Subtract the current species value from each species with an
            # update, as the bolus will be added on to the current species
            # value during simulation.
            for index in range(len(bolus)):
                if bolus[index] != get_empty_bolus_value():
                    bolus[index] -= state_vector[index]
                bolus[index] = bolus[index].subs(get_empty_bolus_value(),
                                                 sp.Float(0.0))

            self.symbols[SymbolId.EVENT][event_sym] = {
                'name': event_id,
                'value': trigger,
                'state_update': sp.MutableDenseMatrix(bolus),
                'event_observable': None,
            }

    @log_execution_time('processing SBML observables', logger)
    def _process_observables(
            self,
            observables: Optional[Dict[str, Dict[str, str]]] = None,
            sigmas: Optional[Dict[str, Union[str, float]]] = None,
            noise_distributions: Optional[Dict[str, str]] = None
    ) -> None:
        """
        Perform symbolic computations required for observable and objective
        function evaluation.

        :param observables:
            dictionary(observableId: {'name':observableName
            (optional), 'formula':formulaString)})
            to be added to the model

        :param sigmas:
            dictionary(observableId: sigma value or (existing)
            parameter name)

        :param noise_distributions:
            dictionary(observableId: noise type)
            See :py:func:`sbml2amici`.
        """

        if sigmas is None or observables is None:
            sigmas = {}
        else:
            # Ensure no non-existing observableIds have been specified
            # (no problem here, but usually an upstream bug)
            unknown_ids = set(sigmas.keys()) - set(observables.keys())
            if unknown_ids:
                raise ValueError(
                    f"Sigma provided for unknown observableIds: "
                    f"{unknown_ids}.")

        if noise_distributions is None or observables is None:
            noise_distributions = {}
        else:
            # Ensure no non-existing observableIds have been specified
            # (no problem here, but usually an upstream bug)
            unknown_ids = set(noise_distributions.keys()) - \
                          set(observables.keys())
            if unknown_ids:
                raise ValueError(
                    f"Noise distribution provided for unknown observableIds: "
                    f"{unknown_ids}.")

        # add user-provided observables or make all species, and compartments
        # with assignment rules, observable
        if observables:
            self.symbols[SymbolId.OBSERVABLE] = {
                symbol_with_assumptions(obs): {
                    'name': definition.get('name', f'y{iobs}'),
                    # Replace logX(.) by log(., X) since sympy cannot parse the
                    # former.
                    'value': self._sympy_from_sbml_math(
                        definition['formula']
                    )
                }
                for iobs, (obs, definition) in enumerate(observables.items())
            }
        elif observables is None:
            self._generate_default_observables()

        self._process_log_likelihood(sigmas, noise_distributions)

    def _generate_default_observables(self):
        """
        Generate default observables from species, compartments and
        (initial) assignment rules.
        """
        self.symbols[SymbolId.OBSERVABLE] = {
            symbol_with_assumptions(f'y{species_id}'): {
                'name': specie['name'],
                'value': species_id
            }
            for ix, (species_id, specie)
            in enumerate(self.symbols[SymbolId.SPECIES].items())
        }

        for variable, formula in itt.chain(
                self.parameter_assignment_rules.items(),
                self.initial_assignments.items(),
                self.compartment_assignment_rules.items(),
                self.species_assignment_rules.items(),
                self.compartments.items()
        ):
            symbol = symbol_with_assumptions(f'y{variable}')
            # Assignment rules take precedence over compartment volume
            # definitions, so they need to be evaluated first.
            # Species assignment rules always overwrite.
            if symbol in self.symbols[SymbolId.OBSERVABLE] \
                    and variable not in self.species_assignment_rules:
                continue
            self.symbols[SymbolId.OBSERVABLE][symbol] = {
                'name': str(variable), 'value': formula
            }

    def _process_log_likelihood(self,
                                sigmas: Dict[str, Union[str, float]],
                                noise_distributions: Dict[str, str]):
        """
        Perform symbolic computations required for objective function
        evaluation.

        :param sigmas:
            See :py:func:`SBMLImporter._process_observables`

        :param noise_distributions:
            See :py:func:`SBMLImporter._process_observables`
        """

        for obs_id, obs in self.symbols[SymbolId.OBSERVABLE].items():
            obs['measurement_symbol'] = generate_measurement_symbol(obs_id)

        self.symbols[SymbolId.SIGMAY] = {
            symbol_with_assumptions(f'sigma_{obs_id}'): {
                'name': f'sigma_{obs["name"]}',
                'value': self._sympy_from_sbml_math(
                    sigmas.get(str(obs_id), '1.0')
                )
            }
            for obs_id, obs in self.symbols[SymbolId.OBSERVABLE].items()
        }

        self.symbols[SymbolId.LLHY] = {}
        for (obs_id, obs), (sigma_id, sigma) in zip(
                self.symbols[SymbolId.OBSERVABLE].items(),
                self.symbols[SymbolId.SIGMAY].items()
        ):
            symbol = symbol_with_assumptions(f'J{obs_id}')
            dist = noise_distributions.get(str(obs_id), 'normal')
            cost_fun = noise_distribution_to_cost_function(dist)(obs_id)
            value = sp.sympify(cost_fun, locals=dict(zip(
                _get_str_symbol_identifiers(obs_id),
                (obs_id, obs['measurement_symbol'], sigma_id)
            )))
            self.symbols[SymbolId.LLHY][symbol] = {
                    'name': f'J{obs["name"]}',
                    'value': value,
                    'dist': dist,
                }

    @log_execution_time('processing SBML initial assignments', logger)
    def _process_initial_assignments(self):
        """
        Accounts for initial assignments of parameters and species
        references. Initial assignments for species and compartments are
        processed in :py:func:`amici.SBMLImporter._process_initial_species` and
        :py:func:`amici.SBMLImporter._process_compartments` respectively.
        """
        for ia in self.sbml.getListOfInitialAssignments():
            identifier = _get_identifier_symbol(ia)
            if identifier in itt.chain(self.symbols[SymbolId.SPECIES],
                                       self.compartments):
                continue

            sym_math = self._get_element_initial_assignment(ia.getId())
            if sym_math is None:
                continue

            sym_math = self._make_initial(smart_subs_dict(
                sym_math, self.symbols[SymbolId.EXPRESSION], 'value'
            ))
            self.initial_assignments[_get_identifier_symbol(ia)] = sym_math

        # sort and flatten
        self.initial_assignments = toposort_symbols(self.initial_assignments)
        for ia_id, ia in self.initial_assignments.items():
            self.initial_assignments[ia_id] = smart_subs_dict(
                ia, self.initial_assignments
            )

        for identifier, sym_math in list(self.initial_assignments.items()):
            self._replace_in_all_expressions(identifier, sym_math)

    @log_execution_time('processing SBML species references', logger)
    def _process_species_references(self):
        """
        Replaces species references that define anything but stoichiometries.

        Species references for stoichiometries are processed in
        :py:func:`amici.SBMLImporter._process_reactions`.
        """
        # doesnt look like there is a better way to get hold of those lists:
        species_references = _get_list_of_species_references(self.sbml)
        for species_reference in species_references:
            if hasattr(species_reference, 'getStoichiometryMath') and \
                    species_reference.getStoichiometryMath() is not None:
                raise SBMLException('StoichiometryMath is currently not '
                                    'supported for species references.')
            if species_reference.getId() == '':
                continue

            stoich = self._get_element_stoichiometry(species_reference)
            self._replace_in_all_expressions(
                _get_identifier_symbol(species_reference),
                self._sympy_from_sbml_math(stoich)
            )

    def _make_initial(self, sym_math: Union[sp.Expr, None, float]
                      ) -> Union[sp.Expr, None, float]:
        """
        Transforms an expression to its value at the initial time point by
        replacing species by their initial values.

        :param sym_math:
            symbolic expression
        :return:
            transformed expression
        """

        if not isinstance(sym_math, sp.Expr):
            return sym_math

        for species_id, species in self.symbols[SymbolId.SPECIES].items():
            if 'init' in species:
                sym_math = smart_subs(sym_math, species_id, species['init'])

        sym_math = smart_subs(sym_math, self._local_symbols['time'],
                              sp.Float(0))

        return sym_math

    def process_conservation_laws(self, ode_model, volume_updates) -> List:
        """
        Find conservation laws in reactions and species.

        :param ode_model:
            ODEModel object with basic definitions

        :param volume_updates:
            List with updates for the stoichiometric matrix accounting for
            compartment volumes

        :returns volume_updates_solver:
            List (according to reduced stoichiometry) with updates for the
            stoichiometric matrix accounting for compartment volumes
        """
        conservation_laws = []

        # So far, only conservation laws for constant species are supported
        species_solver = _add_conservation_for_constant_species(
            ode_model, conservation_laws
        )

        # Check, whether species_solver is empty now. As currently, AMICI
        # cannot handle ODEs without species, CLs must switched in this case
        if len(species_solver) == 0:
            conservation_laws = []
            species_solver = list(range(ode_model.num_states_rdata()))

        # prune out species from stoichiometry and
        volume_updates_solver = self._reduce_stoichiometry(species_solver,
                                                           volume_updates)

        # add the found CLs to the ode_model
        for cl in conservation_laws:
            ode_model.add_conservation_law(**cl)

        return volume_updates_solver

    def _reduce_stoichiometry(self, species_solver, volume_updates) -> List:
        """
        Reduces the stoichiometry with respect to conserved quantities

        :param species_solver:
            List of species indices which remain later in the ODE solver

        :param volume_updates:
            List with updates for the stoichiometric matrix accounting for
            compartment volumes

        :returns volume_updates_solver:
            List (according to reduced stoichiometry) with updates for the
            stoichiometric matrix accounting for compartment volumes
        """

        # prune out constant species from stoichiometric matrix
        self.stoichiometric_matrix = \
            self.stoichiometric_matrix[species_solver, :]

        # updates of stoichiometry (later dxdotdw in ode_exporter) must be
        # corrected for conserved quantities:
        volume_updates_solver = [(species_solver.index(ix), iw, val)
                                 for (ix, iw, val) in volume_updates
                                 if ix in species_solver]

        return volume_updates_solver

    def _replace_compartments_with_volumes(self):
        """
        Replaces compartment symbols in expressions with their respective
        (possibly variable) volumes.
        """
        for comp, vol in self.compartments.items():
            if comp in self.symbols[SymbolId.SPECIES]:
                # for comps with rate rules volume is only initial
                for species in self.symbols[SymbolId.SPECIES].values():
                    if isinstance(species['init'], sp.Expr):
                        species['init'] = smart_subs(species['init'],
                                                     comp, vol)
                continue
            self._replace_in_all_expressions(comp, vol)

    def _replace_in_all_expressions(self,
                                    old: sp.Symbol,
                                    new: sp.Expr,
                                    replace_identifiers=False) -> None:
        """
        Replace 'old' by 'new' in all symbolic expressions.

        :param old:
            symbolic variables to be replaced

        :param new:
            replacement symbolic variables
        """
        fields = [
            'stoichiometric_matrix', 'flux_vector',
        ]
        for field in fields:
            if field in dir(self):
                self.__setattr__(field, smart_subs(
                    self.__getattribute__(field), old, new
                ))

        dictfields = [
            'compartment_assignment_rules', 'parameter_assignment_rules',
            'initial_assignments'
        ]
        for dictfield in dictfields:
            d = getattr(self, dictfield)

            # replace identifiers
            if old in d and replace_identifiers:
                d[new] = d[old]
                del d[old]

            if dictfield == 'initial_assignments':
                tmp_new = self._make_initial(new)
            else:
                tmp_new = new

            # replace values
            for k in d:
                d[k] = smart_subs(d[k], old, tmp_new)

        # replace in identifiers
        if replace_identifiers:
            for symbol in [SymbolId.EXPRESSION, SymbolId.SPECIES]:
                # completely recreate the dict to keep ordering consistent
                if old not in self.symbols[symbol]:
                    continue
                self.symbols[symbol] = {
                    smart_subs(k, old, new): v
                    for k, v in self.symbols[symbol].items()
                }

            for symbol in [SymbolId.OBSERVABLE, SymbolId.LLHY,
                           SymbolId.SIGMAY]:
                if old not in self.symbols[symbol]:
                    continue
                self.symbols[symbol][new] = self.symbols[symbol][old]
                del self.symbols[symbol][old]

        # replace in values
        for symbol in [SymbolId.OBSERVABLE, SymbolId.LLHY, SymbolId.SIGMAY,
                       SymbolId.EXPRESSION, SymbolId.EVENT]:
            if not self.symbols.get(symbol, None):
                continue
            for element in self.symbols[symbol].values():
                element['value'] = smart_subs(element['value'], old, new)

        # replace in event state updates (boluses)
        if self.symbols.get(SymbolId.EVENT, False):
            for event in self.symbols[SymbolId.EVENT].values():
                for index in range(len(event['state_update'])):
                    event['state_update'][index] = \
                        smart_subs(event['state_update'][index], old, new)

        if SymbolId.SPECIES in self.symbols:
            for species in self.symbols[SymbolId.SPECIES].values():
                species['init'] = smart_subs(species['init'],
                                             old, self._make_initial(new))

                fields = ['dt']
                if replace_identifiers:
                    fields.append('compartment')

                for field in ['dt']:
                    if field in species:
                        species[field] = smart_subs(species[field], old, new)

        # Initial compartment volume may also be specified with an assignment
        # rule (at the end of the _process_species method), hence needs to be
        # processed here too.
        self.compartments = {smart_subs(c, old, new) if replace_identifiers
                             else c:
                             smart_subs(v, old, self._make_initial(new))
                             for c, v in self.compartments.items()}

    def _clean_reserved_symbols(self) -> None:
        """
        Remove all reserved symbols from self.symbols
        """
        reserved_symbols = ['x', 'k', 'p', 'y', 'w', 'h', 't',
                            'AMICI_EMPTY_BOLUS']
        for sym in reserved_symbols:
            old_symbol = symbol_with_assumptions(sym)
            new_symbol = symbol_with_assumptions(f'amici_{sym}')
            self._replace_in_all_expressions(old_symbol, new_symbol,
                                             replace_identifiers=True)
            for symbols_ids, symbols in self.symbols.items():
                if old_symbol in symbols:
                    # reconstitute the whole dict in order to keep the ordering
                    self.symbols[symbols_ids] = {
                        new_symbol if k is old_symbol else k: v
                        for k, v in symbols.items()
                    }

    def _sympy_from_sbml_math(self, var_or_math: [sbml.SBase, str]
                              ) -> Union[sp.Expr, float, None]:
        """
        Sympify Math of SBML variables with all sanity checks and
        transformations

        :param var_or_math:
            SBML variable that has a getMath() function or math string
        :return:
            sympfified symbolic expression
        """
        if isinstance(var_or_math, sbml.SBase):
            math_string = sbml.formulaToL3StringWithSettings(
                var_or_math.getMath(),
                self.sbml_parser_settings
            )
            ele_name = var_or_math.element_name
        else:
            math_string = var_or_math
            ele_name = 'string'
        math_string = replace_logx(math_string)
        try:
            try:
                formula = sp.sympify(_parse_logical_operators(
                    math_string
                ), locals=self._local_symbols)
            except TypeError as err:
                if str(err) == 'BooleanAtom not allowed in this context.':
                    formula = sp.sympify(_parse_logical_operators(
                        math_string
                    ), locals={'true': sp.Float(1.0), 'false': sp.Float(0.0),
                               **self._local_symbols})
                else:
                    raise
        except (sp.SympifyError, TypeError, ZeroDivisionError) as err:
            raise SBMLException(f'{ele_name} "{math_string}" '
                                'contains an unsupported expression: '
                                f'{err}.')

        if isinstance(formula, sp.Expr):
            formula = _parse_special_functions(formula)
            _check_unsupported_functions(formula,
                                         expression_type=ele_name)
        return formula

    def _get_element_initial_assignment(self,
                                        element_id: str) -> Union[sp.Expr,
                                                                  None]:
        """
        Extract value of sbml variable according to its initial assignment

        :param element_id:
            sbml variable name
        :return:

        """
        assignment = self.sbml.getInitialAssignment(
            element_id
        )
        if assignment is None:
            return None
        sym = self._sympy_from_sbml_math(assignment)
        # this is an initial assignment so we need to use
        # initial conditions
        sym = self._make_initial(sym)
        return sym

    def _get_element_stoichiometry(self,
                                   ele: sbml.SBase) -> sp.Expr:
        """
        Computes the stoichiometry of a reactant or product of a reaction

        :param ele:
            reactant or product
        :return:
            symbolic variable that defines stoichiometry
        """
        if ele.isSetId():
            sym = self._get_element_initial_assignment(ele.getId())
            if sym is not None:
                return sym

            if self.is_assignment_rule_target(ele):
                return _get_identifier_symbol(ele)

        if ele.isSetStoichiometry():
            return sp.Float(ele.getStoichiometry())

        return sp.Float(1.0)

    def is_assignment_rule_target(self, element: sbml.SBase) -> bool:
        """
        Checks if an element has a valid assignment rule in the specified
        model.

        :param element:
            SBML variable

        :return:
            boolean indicating truth of function name
        """
        a = self.sbml.getAssignmentRuleByVariable(element.getId())
        if a is None or self._sympy_from_sbml_math(a) is None:
            return False

        return True

    def is_rate_rule_target(self, element: sbml.SBase) -> bool:
        """
        Checks if an element has a valid assignment rule in the specified
        model.

        :param element:
            SBML variable

        :return:
            boolean indicating truth of function name
        """
        a = self.sbml.getRateRuleByVariable(element.getId())
        if a is None or self._sympy_from_sbml_math(a) is None:
            return False

        return True


def _check_lib_sbml_errors(sbml_doc: sbml.SBMLDocument,
                           show_warnings: bool = False) -> None:
    """
    Checks the error log in the current self.sbml_doc.

    :param sbml_doc:
        SBML document

    :param show_warnings:
        display SBML warnings
    """
    num_warning = sbml_doc.getNumErrors(sbml.LIBSBML_SEV_WARNING)
    num_error = sbml_doc.getNumErrors(sbml.LIBSBML_SEV_ERROR)
    num_fatal = sbml_doc.getNumErrors(sbml.LIBSBML_SEV_FATAL)

    if num_warning + num_error + num_fatal:
        for i_error in range(0, sbml_doc.getNumErrors()):
            error = sbml_doc.getError(i_error)
            # we ignore any info messages for now
            if error.getSeverity() >= sbml.LIBSBML_SEV_ERROR \
                    or (show_warnings and
                        error.getSeverity() >= sbml.LIBSBML_SEV_WARNING):
                logger.error(f'libSBML {error.getCategoryAsString()} '
                             f'({error.getSeverityAsString()}):'
                             f' {error.getMessage()}')

    if num_error + num_fatal:
        raise SBMLException(
            'SBML Document failed to load (see error messages above)'
        )


def _check_unsupported_functions(sym: sp.Expr,
                                 expression_type: str,
                                 full_sym: Optional[sp.Expr] = None):
    """
    Recursively checks the symbolic expression for unsupported symbolic
    functions

    :param sym:
        symbolic expressions

    :param expression_type:
        type of expression, only used when throwing errors

    :param full sym:
        outermost symbolic expression in recursive checks, only used for errors
    """
    if full_sym is None:
        full_sym = sym

    # note that sp.functions.factorial, sp.functions.ceiling,
    # sp.functions.floor applied to numbers should be simplified out and
    # thus pass this test
    unsupported_functions = (
        sp.functions.factorial, sp.functions.ceiling, sp.functions.floor,
        sp.functions.sec, sp.functions.csc, sp.functions.cot,
        sp.functions.asec, sp.functions.acsc, sp.functions.acot,
        sp.functions.acsch, sp.functions.acoth,
        sp.Mod, sp.core.function.UndefinedFunction
    )

    if isinstance(sym.func, unsupported_functions) \
            or isinstance(sym, unsupported_functions):
        raise SBMLException(f'Encountered unsupported expression '
                            f'"{sym.func}" of type '
                            f'"{type(sym.func)}" as part of a '
                            f'{expression_type}: "{full_sym}"!')
    for arg in list(sym.args):
        _check_unsupported_functions(arg, expression_type)


def _parse_special_functions(sym: sp.Expr, toplevel: bool = True) -> sp.Expr:
    """
    Recursively checks the symbolic expression for functions which have be
    to parsed in a special way, such as piecewise functions

    :param sym:
        symbolic expressions

    :param toplevel:
        as this is called recursively, are we in the top level expression?
    """
    args = tuple(arg if arg.__class__.__name__ == 'piecewise'
                 and sym.__class__.__name__ == 'piecewise'
                 else _parse_special_functions(arg, False)
                 for arg in sym.args)

    fun_mappings = {
        'times': sp.Mul,
        'xor': sp.Xor,
        'abs': sp.Abs,
        'min': sp.Min,
        'max': sp.Max,
        'ceil': sp.functions.ceiling,
        'floor': sp.functions.floor,
        'factorial': sp.functions.factorial,
        'arcsin': sp.functions.asin,
        'arccos': sp.functions.acos,
        'arctan': sp.functions.atan,
        'arccot': sp.functions.acot,
        'arcsec': sp.functions.asec,
        'arccsc': sp.functions.acsc,
        'arcsinh': sp.functions.asinh,
        'arccosh': sp.functions.acosh,
        'arctanh': sp.functions.atanh,
        'arccoth': sp.functions.acoth,
        'arcsech': sp.functions.asech,
        'arccsch': sp.functions.acsch,
    }

    if sym.__class__.__name__ in fun_mappings:
        # c++ doesnt like mixing int and double for arguments of those
        # functions
        if sym.__class__.__name__ in ['min', 'max']:
            args = tuple([
                sp.Float(arg) if arg.is_number else arg
                for arg in sym.args
            ])
        else:
            args = sym.args
        return fun_mappings[sym.__class__.__name__](*args)

    elif sym.__class__.__name__ == 'piecewise':
        # We need to parse piecewise functions into Heavisides
        return _parse_piecewise_to_heaviside(
            _denest_piecewise(args)
        )

    if sym.__class__.__name__ == 'plus' and not sym.args:
        return sp.Float(0.0)

    if isinstance(sym, (sp.Function, sp.Mul, sp.Add)):
        sym._args = args

    elif toplevel and isinstance(sym, BooleanAtom):
        # Replace boolean constants by numbers so they can be differentiated
        #  must not replace in Piecewise function. Therefore, we only replace
        #  it the complete expression consists only of a Boolean value.
        sym = sp.Float(int(bool(sym)))

    return sym


def _denest_piecewise(
        args: Sequence[Union[sp.Expr, sp.logic.boolalg.Boolean, bool]]
) -> Tuple[Union[sp.Expr, sp.logic.boolalg.Boolean, bool]]:
    """
    Denest piecewise functions that contain piecewise as condition

    :param args:
        Arguments to the piecewise function

    :return:
        Arguments where conditions no longer contain piecewise functions and
        the conditional dependency is flattened out
    """
    args_out = []
    for coeff, cond in grouper(args, 2, True):
        # handling of this case is explicitely disabled in
        # _parse_special_functions as keeping track of coeff/cond
        # arguments is tricky. Simpler to just parse them out here
        if coeff.__class__.__name__ == 'piecewise':
            coeff = _parse_special_functions(coeff, False)

        # we can have conditions that are piecewise function
        # returning True or False
        if cond.__class__.__name__ == 'piecewise':
            # this keeps track of conditional that the previous
            # piece was picked
            previous_was_picked = sp.false
            # recursively denest those first
            for sub_coeff, sub_cond in grouper(
                    _denest_piecewise(cond.args), 2, True
            ):
                # flatten the individual pieces
                pick_this = sp.And(
                    sp.Not(previous_was_picked), sub_cond
                )
                if sub_coeff == sp.true:
                    args_out.extend([coeff, pick_this])
                previous_was_picked = pick_this

        else:
            args_out.extend([coeff, cond])
    # cut off last condition as that's the default
    return tuple(args_out[:-1])


def _parse_piecewise_to_heaviside(args: Iterable[sp.Expr]) -> sp.Expr:
    """
    Piecewise functions cannot be transformed into C++ right away, but AMICI
    has a special interface for Heaviside functions, so we transform them.

    :param args:
        symbolic expressions for arguments of the piecewise function
    """
    # how many condition-expression pairs will we have?
    formula = sp.Float(0.0)
    not_condition = sp.Float(1.0)

    for coeff, trigger in grouper(args, 2, True):
        if isinstance(coeff, BooleanAtom):
            coeff = sp.Float(int(bool(coeff)))

        if trigger == sp.true:
            return formula + coeff * not_condition

        if trigger == sp.false:
            continue

        tmp = _parse_heaviside_trigger(trigger)
        formula += coeff * sp.simplify(not_condition * tmp)
        not_condition *= (1-tmp)

    return formula


def _parse_heaviside_trigger(trigger: sp.Expr) -> sp.Expr:
    """
    Recursively translates a boolean trigger function into a real valued
    root function

    :param trigger:
    :return: real valued root function expression
    """
    if trigger.is_Relational:
        root = trigger.args[0] - trigger.args[1]
        _check_unsupported_functions(root, 'sympy.Expression')

        # normalize such that we always implement <,
        # this ensures that we can correctly evaluate the condition if
        # simulation starts at H(0). This is achieved by translating
        # conditionals into Heaviside functions H that is implemented as unit
        # step with H(0) = 1
        if isinstance(trigger, sp.core.relational.StrictLessThan):
            # x < y => x - y < 0 => r < 0
            return 1 - sp.Heaviside(root)
        if isinstance(trigger, sp.core.relational.LessThan):
            # x <= y => not(y < x) => not(y - x < 0) => not -r < 0
            return sp.Heaviside(-root)
        if isinstance(trigger, sp.core.relational.StrictGreaterThan):
            # y > x => y - x < 0 => -r < 0
            return 1 - sp.Heaviside(-root)
        if isinstance(trigger, sp.core.relational.GreaterThan):
            # y >= x => not(x < y) => not(x - y < 0) => not r < 0
            return sp.Heaviside(root)

    # or(x,y) = not(and(not(x),not(y))
    if isinstance(trigger, sp.Or):
        return 1-sp.Mul(*[1-_parse_heaviside_trigger(arg)
                        for arg in trigger.args])

    if isinstance(trigger, sp.And):
        return sp.Mul(*[_parse_heaviside_trigger(arg)
                        for arg in trigger.args])

    raise SBMLException(
        'AMICI can not parse piecewise/event trigger functions with argument '
        f'{trigger}.'
    )


def _parse_event_trigger(trigger: sp.Expr) -> sp.Expr:
    """
    Recursively translates a boolean trigger function into a real valued
    root function

    :param trigger:
    :return: real valued root function expression
    """
    # Events can be defined without trigger, i.e., the event will never fire.
    # In this case, set a dummy trigger:
    if trigger is None:
        return sp.Float(1.0)
    if trigger.is_Relational:
        root = trigger.args[0] - trigger.args[1]
        _check_unsupported_functions(root, 'sympy.Expression')

        # convert relational expressions into trigger functions
        if isinstance(trigger, sp.core.relational.LessThan) or \
                isinstance(trigger, sp.core.relational.StrictLessThan):
            # y < x or y <= x
            return -root
        if isinstance(trigger, sp.core.relational.GreaterThan) or \
                isinstance(trigger, sp.core.relational.StrictGreaterThan):
            # y >= x or y > x
            return root

    # or(x,y): any of {x,y} is > 0: sp.Max(x, y)
    if isinstance(trigger, sp.Or):
        return sp.Max(*[_parse_event_trigger(arg) for arg in trigger.args])
    # and(x,y): all out of {x,y} are > 0: sp.Min(x, y)
    if isinstance(trigger, sp.And):
        return sp.Min(*[_parse_event_trigger(arg) for arg in trigger.args])

    raise SBMLException(
        'AMICI can not parse piecewise/event trigger functions with argument '
        f'{trigger}.'
    )


def _parse_logical_operators(math_str: Union[str, float, None]
                             ) -> Union[str, float, None]:
    """
    Parses a math string in order to replace logical operators by a form
    parsable for sympy

    :param math_str:
        str with mathematical expression
    :param math_str:
        parsed math_str
    """
    if not isinstance(math_str, str):
        return math_str

    if ' xor(' in math_str or ' Xor(' in math_str:
        raise SBMLException('Xor is currently not supported as logical '
                            'operation.')

    return (math_str.replace('&&', '&')).replace('||', '|')


def grouper(iterable: Iterable, n: int,
            fillvalue: Any = None) -> Iterable[Tuple[Any]]:
    """
    Collect data into fixed-length chunks or blocks

    grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"

    :param iterable:
        any iterable

    :param n:
        chunk length

    :param fillvalue:
        padding for last chunk if length < n

    :return: itertools.zip_longest of requested chunks
    """
    args = [iter(iterable)] * n
    return itt.zip_longest(*args, fillvalue=fillvalue)


def assignmentRules2observables(sbml_model: sbml.Model,
                                filter_function: Callable = lambda *_: True):
    """
    Turn assignment rules into observables.

    :param sbml_model:
        Model to operate on

    :param filter_function:
        Callback function taking assignment variable as input and returning
        ``True``/``False`` to indicate if the respective rule should be
        turned into an observable.

    :return:
        A dictionary(observableId:{
        'name': observableName,
        'formula': formulaString
        })
    """
    observables = {}
    for p in sbml_model.getListOfParameters():
        parameter_id = p.getId()
        if filter_function(p):
            observables[parameter_id] = {
                'name': p.getName() if p.isSetName() else p.getId(),
                'formula': sbml_model.getAssignmentRuleByVariable(
                    parameter_id
                ).getFormula()
            }

    for parameter_id in observables:
        sbml_model.removeRuleByVariable(parameter_id)
        sbml_model.removeParameter(parameter_id)

    return observables


def _add_conservation_for_constant_species(
        ode_model: ODEModel,
        conservation_laws: List[ConservationLaw]
) -> List[int]:
    """
    Adds constant species to conservations laws

    :param ode_model:
        ODEModel object with basic definitions

    :param conservation_laws:
        List of already known conservation laws

    :returns species_solver:
        List of species indices which remain later in the ODE solver
    """

    # decide which species to keep in stoichiometry
    species_solver = list(range(ode_model.num_states_rdata()))

    # iterate over species, find constant ones
    for ix in reversed(range(ode_model.num_states_rdata())):
        if ode_model.state_is_constant(ix):
            # dont use sym('x') here since conservation laws need to be
            # added before symbols are generated
            target_state = ode_model._states[ix].get_id()
            total_abundance = symbol_with_assumptions(f'tcl_{target_state}')
            conservation_laws.append({
                'state': target_state,
                'total_abundance': total_abundance,
                'state_expr': total_abundance,
                'abundance_expr': target_state,
            })
            # mark species to delete from stoichiometric matrix
            species_solver.pop(ix)

    return species_solver


def _get_species_compartment_symbol(species: sbml.Species) -> sp.Symbol:
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


def _get_identifier_symbol(var: sbml.SBase) -> sp.Symbol:
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


def get_species_initial(species: sbml.Species) -> sp.Expr:
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


def _get_list_of_species_references(sbml_model: sbml.Model) \
        -> List[sbml.SpeciesReference]:
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
        for element in sbml_model.all_elements
        if isinstance(element, sbml.ListOfSpeciesReferences)
        for reference in element
    ]


def replace_logx(math_str: Union[str, float, None]) -> Union[str, float, None]:
    """
    Replace logX(.) by log(., X) since sympy cannot parse the former

    :param math_str:
        string for sympification

    :return:
        sympifiable string
    """
    if not isinstance(math_str, str):
        return math_str

    return re.sub(
        r'(^|\W)log(\d+)\(', r'\g<1>1/ln(\2)*ln(', math_str
    )


def _collect_event_assignment_parameter_targets(sbml_model: sbml.Model):
    targets = set()
    sbml_parameters = sbml_model.getListOfParameters()
    sbml_parameter_ids = [p.getId() for p in sbml_parameters]
    for event in sbml_model.getListOfEvents():
        for event_assignment in event.getListOfEventAssignments():
            target_id = event_assignment.getVariable()
            if target_id in sbml_parameter_ids:
                targets.add(_get_identifier_symbol(
                    sbml_parameters[sbml_parameter_ids.index(target_id)]
                ))
    return targets
