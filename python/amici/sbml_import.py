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
    Dict, List, Callable, Any, Iterable, Sequence, Union
)

from .ode_export import (
    ODEExporter, ODEModel, generate_measurement_symbol
)
from .constants import SymbolId
from .logging import get_logger, log_execution_time, set_log_level
from . import has_clibs

from sympy.logic.boolalg import BooleanTrue as spTrue
from sympy.logic.boolalg import BooleanFalse as spFalse
from sympy.printing.mathml import MathMLContentPrinter

# the following import can be removed when sympy 1.6.3 is released
from mpmath.libmp import repr_dps, to_str as mlib_to_str


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

    :ivar local_symbols:
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

    :ivar parameter_initial_assignments:
        initial assignments for parameters, these parameters are not
        permissible for sensitivity analysis

    :ivar reaction_ids:
        symbol definition as kinetic law of the respective reaction

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

        self.local_symbols: Dict[str, sp.Expr] = {}
        self.compartments: SymbolicFormula = {}
        self.compartment_assignment_rules: SymbolicFormula = {}
        self.species_assignment_rules: SymbolicFormula = {}
        self.parameter_assignment_rules: SymbolicFormula = {}
        self.parameter_initial_assignments: SymbolicFormula = {}

        self._reset_symbols()

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
        self.local_symbols = dict()

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
            verbosity level for logging, True/False default to
            logging.Error/logging.DEBUG

        :param assume_pow_positivity:
            if set to True, a special pow function is
            used to avoid problems with state variables that may become
            negative due to numerical errors

        :param compiler:
            distutils/setuptools compiler selection to build the
            python extension

        :param allow_reinit_fixpar_initcond:
            see :class:`amici.ode_export.ODEExporter`

        :param compile:
            If True, compile the generated Python package,
            if False, just generate code.

        :param compute_conservation_laws:
            if set to true, conservation laws are automatically computed and
            applied such that the state-jacobian of the ODE right-hand-side has
            full rank. This option should be set to True when using the newton
            algorithm to compute steadystate sensitivities.

        :param simplify:
            see :attr:`ODEModel._simplify`
        """
        set_log_level(logger, verbose)

        if observables is None:
            observables = {}

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
        self._process_sbml(constant_parameters)
        self._process_observables(observables, sigmas, noise_distributions)
        self._replace_compartments_with_volumes()

        self._process_time()
        self._clean_reserved_symbols()
        self._replace_special_constants()

        ode_model = ODEModel(verbose=verbose, simplify=simplify)
        ode_model.import_from_sbml_importer(
            self, compute_cls=compute_conservation_laws)
        exporter = ODEExporter(
            ode_model,
            outdir=output_dir,
            verbose=verbose,
            assume_pow_positivity=assume_pow_positivity,
            compiler=compiler,
            allow_reinit_fixpar_initcond=allow_reinit_fixpar_initcond
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

        if constant_parameters is None:
            constant_parameters = []

        self.check_support()
        self._gather_locals()
        self._process_parameters(constant_parameters)
        self._process_compartments()
        self._process_species()
        self._process_reactions()
        self._process_rules()
        self._process_initial_assignments()
        self._process_species_references()

    def check_support(self) -> None:
        """
        Check whether all required SBML features are supported.
        Also ensures that the SBML contains at least one reaction, or rate
        rule, or assignment rule, to produce change in the system over time.
        """
        if hasattr(self.sbml, 'all_elements_from_plugins') \
                and self.sbml.all_elements_from_plugins.getSize():
            raise SBMLException('SBML extensions are currently not supported!')

        if self.sbml.getNumEvents():
            raise SBMLException('Events are currently not supported!')

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
                                'species, and parameters! Error '
                                f'occurred with rule: {rule.getVariable()}')

        if any([r.getFast() for r in self.sbml.getListOfReactions()]):
            raise SBMLException('Fast reactions are currently not supported!')

        if any([any([not element.getStoichiometryMath() is None
                     for element in itt.chain(reaction.getListOfReactants(),
                                              reaction.getListOfProducts())])
                for reaction in self.sbml.getListOfReactions()]):
            raise SBMLException('Non-unity stoichiometry is'
                                ' currently not supported!')

    @log_execution_time('gathering local SBML symbols', logger)
    def _gather_locals(self) -> None:
        """
        Populate self.local_symbols with all model entities.

        This is later used during sympifications to avoid sympy builtins
        shadowing model entities as well as to avoid possibly costly
        symbolic substitutions/
        """
        self._gather_base_locals()
        self._gather_dependent_locals()

    def _gather_base_locals(self):
        """
        Populate self.local_symbols with pure symbol definitions that do not
        depend on any other symbol.
        """
        species_references = _get_list_of_species_references(self.sbml)
        for c in itt.chain(self.sbml.getListOfSpecies(),
                           self.sbml.getListOfParameters(),
                           self.sbml.getListOfCompartments(),
                           species_references):
            if c.getId():
                self.local_symbols[c.getId()] = _get_identifier_symbol(c)

        for r in self.sbml.getListOfRules():
            self.local_symbols[r.getVariable()] = symbol_with_assumptions(
                r.getVariable()
            )

        for r in self.sbml.getListOfReactions():
            if r.isSetId():
                self.local_symbols[r.getId()] = _get_identifier_symbol(r)

        # SBML time symbol + constants
        for symbol_name in ['time', 'avogadro']:
            if symbol_name in self.local_symbols:
                # Supporting this is probably kinda tricky and this sounds
                # like a stupid thing to do in the first place.
                raise SBMLException(
                    'AMICI does not support SBML models '
                    'containing variables with Id '
                    f'{symbol_name}.')
            self.local_symbols[symbol_name] = symbol_with_assumptions(
                symbol_name
            )

    def _gather_dependent_locals(self):
        """
        Populate self.local_symbols with symbol definitions that may depend on
        other symbol definitions.
        """
        for r in self.sbml.getListOfReactions():
            for e in itt.chain(r.getListOfReactants(), r.getListOfProducts()):
                if e.isSetId() and e.isSetStoichiometry() and \
                        e.getId() not in self.local_symbols:
                    self.local_symbols[e.getId()] = sp.sympify(
                        e.getStoichiometry(), locals=self.local_symbols
                    )

        for r in self.sbml.getListOfReactions():
            if r.isSetId():
                self.local_symbols[r.getId()] = self._sympy_from_sbml_math(
                    r.getKineticLaw()
                )

    @log_execution_time('processing SBML compartments', logger)
    def _process_compartments(self) -> None:
        """
        Get compartment information, stoichiometric matrix and fluxes from
        SBML model.
        """
        compartments = self.sbml.getListOfCompartments()
        self.compartments = {}
        for comp in compartments:
            init = sp.sympify(1.0)

            if comp.isSetVolume():
                init = sp.sympify(comp.getVolume(),
                                  locals=self.local_symbols)

            ia = self.sbml.getInitialAssignment(comp.getId())
            if ia is not None:
                init = self._sympy_from_sbml_math(ia)

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

        self.symbols[SymbolId.SPECIES] = {
            _get_identifier_symbol(s): {
                'name': s.getName() if s.isSetName() else s.getId(),
                'compartment': _get_species_compartment_symbol(s),
                'constant': s.getConstant() or s.getBoundaryCondition(),
                'amount': s.getHasOnlySubstanceUnits(),
                'conversion_factor': symbol_with_assumptions(
                    s.getConversionFactor()
                )
                if s.isSetConversionFactor()
                else conversion_factor,
                'index': ix,
            }
            for ix, s in enumerate(self.sbml.getListOfSpecies())
            if not is_assignment_rule_target(self.sbml, s)
        }

        self._process_species_initial()
        self._process_rate_rules()

    @log_execution_time('processing SBML species initials', logger)
    def _process_species_initial(self):
        """
        Extract initial values and initial assignments from species
        """
        for specie in self.sbml.getListOfSpecies():
            initial = _get_species_initial(specie)

            specie_id = _get_identifier_symbol(specie)
            # species is target of AssignmentRule targets, specie_def will be
            # None, but we don't have to account for the initial definition
            # of the species itself and SBML doesnt permit AssignmentRule
            # targets to have InitialAssignments.
            specie_def = self.symbols[SymbolId.SPECIES].get(specie_id, None)

            ia = self.sbml.getInitialAssignment(specie.getId())
            if ia is not None:
                ia_initial = self._sympy_from_sbml_math(ia)
                if ia_initial is not None:
                    if specie_def and specie_def['amount'] and 'compartment' \
                            in specie_def:
                        ia_initial *= self.compartments.get(
                            specie_def['compartment'],
                            specie_def['compartment']
                        )
                    initial = ia_initial
            if specie_def:
                specie_def['init'] = initial

        self._flatten_species_initial()

    @log_execution_time('flattening SBML species initials', logger)
    def _flatten_species_initial(self):
        """
        Flattens interdependency of species initial values
        """
        for specie in self.symbols[SymbolId.SPECIES].values():
            nested_species = True
            while nested_species:
                nested_species = False
                for symbol in specie['init'].free_symbols:
                    if is_assignment_rule_target(
                            self.sbml, self.sbml.getElementBySId(str(symbol))
                    ):
                        raise SBMLException('AMICI does not support '
                                            'species initial  values that '
                                            'are targets of assignment rules.'
                                            f'Initial value of {specie} '
                                            'depends on assignment of '
                                            f'{symbol}')

                    if symbol not in self.symbols[SymbolId.SPECIES]:
                        continue

                    nested_species = True
                    specie['init'] = specie['init'].subs(
                        symbol, self.symbols[SymbolId.SPECIES][symbol]['init']
                    )

    @log_execution_time('processing SBML rate rules', logger)
    def _process_rate_rules(self):
        """
        Process assignment and rate rules for species, compartments and
        parameters. Compartments and parameters with rate rules are
        implemented as species. Note that, in the case of species,
        rate rules may describe the change in amount, not concentration,
        of a species.
        """
        rules = self.sbml.getListOfRules()
        # compartments with rules are replaced with constants in the relevant
        # equations during the _replace_in_all_expressions call inside
        # _process_rules
        for rule in rules:
            if rule.getFormula() == '':
                continue
            if rule.getTypeCode() != sbml.SBML_RATE_RULE:
                continue

            variable = sp.sympify(rule.getVariable(),
                                  locals=self.local_symbols)
            formula = self._sympy_from_sbml_math(rule)

            # Species rules are processed first, to avoid processing
            # compartments twice (as compartments with rate rules are
            # implemented as species).
            if variable in self.symbols[SymbolId.SPECIES]:
                init = self.symbols[SymbolId.SPECIES][variable]['init']
                component_type = sbml.SBML_SPECIES
                name = None

            elif variable in self.compartments:
                init = self.compartments[variable]
                component_type = sbml.SBML_COMPARTMENT
                name = str(variable)
                del self.compartments[variable]

            elif variable in self.symbols[SymbolId.PARAMETER]:
                init = sp.sympify(
                    self.symbols[SymbolId.PARAMETER][variable]['value'],
                    locals=self.local_symbols
                )
                name = self.symbols[SymbolId.PARAMETER][variable]['name']
                del self.symbols[SymbolId.PARAMETER][variable]
                component_type = sbml.SBML_PARAMETER

            self.add_d_dt(formula, variable, init, component_type, name)

    def add_d_dt(
            self,
            d_dt: sp.Expr,
            variable: sp.Symbol,
            variable0: Union[float, sp.Expr],
            component_type: int,
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

        :param component_type:
            The type of SBML component. Currently, species and compartments
            are supported.

        :param name:
            Species name, only applicable if this function generates a new
            species
        """
        if component_type in [sbml.SBML_COMPARTMENT, sbml.SBML_PARAMETER]:
            # update initial values
            for specie_id, specie in self.symbols[SymbolId.SPECIES].items():
                variable0 = variable0.subs(specie_id, specie['init'])

            for specie in self.symbols[SymbolId.SPECIES].values():
                specie['init'] = specie['init'].subs(variable, variable0)

            # add compartment/parameter species
            self.symbols[SymbolId.SPECIES][variable] = {
                'name': name,
                'init': variable0,
                'amount': component_type == sbml.SBML_COMPARTMENT,
                'conversion_factor': 1.0,
                'constant': False,
                'index': len(self.symbols[SymbolId.SPECIES]),
                'dt': d_dt,
            }
        elif component_type == sbml.SBML_SPECIES:
            # SBML species are already in the species symbols
            self.symbols[SymbolId.SPECIES][variable]['dt'] = d_dt

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

        parameters = [
            parameter for parameter
            in self.sbml.getListOfParameters()
            if parameter.getId() not in constant_parameters
            and self.sbml.getInitialAssignment(parameter.getId()) is None
            and not is_assignment_rule_target(self.sbml, parameter)
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
                    species_id = _get_identifier_symbol(
                        self.sbml.getSpecies(element.getSpecies())
                    )
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
                sym_math = self.local_symbols[reaction.getId()]
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
            if rule.getFormula() == '':
                continue

            sbml_var = self.sbml.getElementBySId(rule.getVariable())
            sym_id = sp.sympify(rule.getVariable(),
                                locals=self.local_symbols)
            formula = self._sympy_from_sbml_math(rule)

            if isinstance(sbml_var, sbml.Species):
                self.species_assignment_rules[sym_id] = formula
                self.symbols[SymbolId.EXPRESSION][sym_id] = {
                    'name': str(sym_id),
                    'value': formula
                }

            elif isinstance(sbml_var, sbml.Compartment):
                self.compartment_assignment_rules[sym_id] = formula
                self.symbols[SymbolId.EXPRESSION][sym_id] = {
                    'name': str(sym_id),
                    'value': formula
                }

            elif isinstance(sbml_var, sbml.Parameter):
                self.parameter_assignment_rules[sym_id] = formula
                self.symbols[SymbolId.EXPRESSION][sym_id] = {
                    'name': str(sym_id),
                    'value': formula
                }

    def _process_time(self) -> None:
        """
        Convert time_symbol into cpp variable.
        """
        sbml_time_symbol = symbol_with_assumptions('time')
        amici_time_symbol = symbol_with_assumptions('t')
        self.amici_time_symbol = amici_time_symbol

        self._replace_in_all_expressions(sbml_time_symbol, amici_time_symbol)

    @log_execution_time('processing SBML observables', logger)
    def _process_observables(self,
                             observables: Dict[str, Dict[str, str]],
                             sigmas: Dict[str, Union[str, float]],
                             noise_distributions: Dict[str, str]) -> None:
        """
        Perform symbolic computations required for objective function
        evaluation.

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

        if observables is None:
            observables = {}

        if sigmas is None:
            sigmas = {}
        else:
            # Ensure no non-existing observableIds have been specified
            # (no problem here, but usually an upstream bug)
            unknown_ids = set(sigmas.keys()) - set(observables.keys())
            if unknown_ids:
                raise ValueError(
                    f"Sigma provided for unknown observableIds: "
                    f"{unknown_ids}.")

        if noise_distributions is None:
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
                    'value': re.sub(
                        r'(^|\W)log(\d+)\(', r'\g<1>1/ln(\2)*ln(',
                        definition['formula']
                    )
                }
                for iobs, (obs, definition) in enumerate(observables.items())
            }
        else:
            self.symbols[SymbolId.OBSERVABLE] = {
                symbol_with_assumptions(f'y{specie["name"]}'): {
                    'name': specie['name'],
                    'value': specie_id
                }
                for ix, (specie_id, specie)
                in enumerate(self.symbols[SymbolId.SPECIES].items())
            }

            # Assignment rules take precedence over compartment volume
            # definitions, so they need to be evaluated first.
            # Species assignment rules always overwrite
            for variable, formula in (
                *self.parameter_assignment_rules.items(),
                *self.parameter_initial_assignments.items(),
                *self.compartment_assignment_rules.items(),
                *self.species_assignment_rules.items(),
                *self.compartments.items()
            ):
                symbol = symbol_with_assumptions(f'y{variable}')
                if symbol in self.symbols[SymbolId.OBSERVABLE] \
                        and variable not in self.species_assignment_rules:
                    continue
                self.symbols[SymbolId.OBSERVABLE][symbol] = {
                    'name': str(variable), 'value': formula
                }

        for obs_id, obs in self.symbols[SymbolId.OBSERVABLE].items():
            obs['measurement_symbol'] = generate_measurement_symbol(obs_id)

        self.symbols[SymbolId.SIGMAY] = {
            symbol_with_assumptions(f'sigma_{obs_id}'): {
                'name': f'sigma_{obs["name"]}',
                'value': sigmas.get(str(obs_id), '1.0')
            }
            for obs_id, obs in self.symbols[SymbolId.OBSERVABLE].items()
        }

        self.symbols[SymbolId.LLHY] = {
            symbol_with_assumptions(f'J{obs_id}'): {
                'name': f'J{obs["name"]}',
                'value': sp.sympify(noise_distribution_to_cost_function(
                    noise_distributions.get(str(obs_id), 'normal')
                )(obs_id), locals=dict(zip(
                    _get_str_symbol_identifiers(obs_id),
                    (obs_id, obs['measurement_symbol'], sigma_id)
                )))
            }
            for (obs_id, obs), (sigma_id, sigma) in zip(
                    self.symbols[SymbolId.OBSERVABLE].items(),
                    self.symbols[SymbolId.SIGMAY].items()
            )
        }

    @log_execution_time('processing SBML initial assignments', logger)
    def _process_initial_assignments(self):
        """
        Accounts for initial assignments of parameters and species
        references. Initial assignments for species and compartments are
        processed in :py:func:`amici.SBMLImporter._process_initial_species` and
        :py:func:`amici.SBMLImporter._process_compartments` respectively.
        """
        parameter_ids = [p.getId() for p in self.sbml.getListOfParameters()]
        species_ids = [s.getId() for s in self.sbml.getListOfSpecies()]
        comp_ids = [c.getId() for c in self.sbml.getListOfCompartments()]
        for ia in self.sbml.getListOfInitialAssignments():
            if ia.getId() in species_ids + comp_ids:
                continue

            sym_math = self._sympy_from_sbml_math(ia)
            sym_math = self._make_initial(sym_math)

            identifier = _get_identifier_symbol(ia)
            if ia.getId() in parameter_ids:
                self.parameter_initial_assignments[identifier] = sym_math
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
                sp.sympify(stoich, locals=self.local_symbols)
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
            sym_math = sym_math.subs(species_id, species['init'])

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
                for specie in self.symbols[SymbolId.SPECIES].values():
                    if isinstance(specie['init'], sp.Expr):
                        specie['init'] = specie['init'].subs(comp, vol)
                continue
            self._replace_in_all_expressions(comp, vol)

    def _replace_in_all_expressions(self,
                                    old: sp.Symbol,
                                    new: sp.Expr) -> None:
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
                self.__setattr__(field, self.__getattribute__(field).subs(
                    old, new
                ))

        dictfields = [
            'compartment_assignment_rules', 'parameter_assignment_rules',
            'parameter_initial_assignments'
        ]
        for dictfield in dictfields:
            d = getattr(self, dictfield)

            if dictfield == 'parameter_initial_assignments':
                new = self._make_initial(new)

            for k in d:
                d[k] = d[k].subs(old, new)

        for symbol in [SymbolId.OBSERVABLE, SymbolId.LLHY, SymbolId.SIGMAY]:
            if not self.symbols.get(symbol, None):
                continue
            for element in self.symbols[symbol].values():
                element['value'] = element['value'].subs(old, new)

        if SymbolId.SPECIES in self.symbols:
            for specie in self.symbols[SymbolId.SPECIES].values():
                specie['init'] = specie['init'].subs(old,
                                                     self._make_initial(new))
                if 'dt' in specie:
                    specie['dt'] = specie['dt'].subs(old, new)


        # Initial compartment volume may also be specified with an assignment
        # rule (at the end of the _process_species method), hence needs to be
        # processed here too.
        subs = 0 if getattr(self, 'amici_time_symbol', sp.nan) == new else new
        self.compartments = {c: v.subs(old, subs)
                             for c, v in self.compartments.items()}

    def _clean_reserved_symbols(self) -> None:
        """
        Remove all reserved symbols from self.symbols
        """
        reserved_symbols = ['x', 'k', 'p', 'y', 'w']
        for sym in reserved_symbols:
            old_symbol = symbol_with_assumptions(sym)
            new_symbol = symbol_with_assumptions(f'amici_{sym}')
            self._replace_in_all_expressions(old_symbol, new_symbol)
            for symbols in self.symbols.values():
                if old_symbol in symbols:
                    symbols[new_symbol] = symbols[old_symbol]
                    del symbols[old_symbol]

    def _replace_special_constants(self) -> None:
        """
        Replace all special constants by their respective SBML
        csymbol definition
        """
        constants = [
            (symbol_with_assumptions('avogadro'), sp.sympify(6.02214179e23)),
        ]
        for constant, value in constants:
            # do not replace if any symbol is shadowing default definition
            if not any([constant in self.symbols[symbol]['identifier']
                        for symbol in self.symbols.keys()
                        if 'identifier' in self.symbols[symbol].keys()]):
                self._replace_in_all_expressions(constant, value)
            else:
                # yes sbml supports this but we wont. Are you really expecting
                # to be saved if you are trying to shoot yourself in the foot?
                raise SBMLException(
                    f'Encountered currently unsupported element id {constant}!'
                )

    def _sympy_from_sbml_math(self, var: sbml.SBase) -> sp.Expr:
        """
        Sympify Math of SBML variables with all sanity checks and
        transformations

        :param var:
            SBML variable that has a getMath() function
        :return:
            sympfified symbolic expression
        """

        math_string = sbml.formulaToL3String(var.getMath())
        try:
            formula = sp.sympify(_parse_logical_operators(
                math_string
            ), locals=self.local_symbols)
        except sp.SympifyError:
            raise SBMLException(f'{var.element_name} "{math_string}" '
                                f'contains an unsupported expression!')

        if formula is not None:
            formula = _parse_special_functions(formula)
            _check_unsupported_functions(formula,
                                         expression_type=var.element_name)
        return formula

    def _get_element_from_assignment(self, element_id: str) -> sp.Expr:
        """
        Extract value of sbml variable according to its initial assignment

        :param element_id:
            sbml variable name
        :return:

        """
        assignment = self.sbml.getInitialAssignment(
            element_id
        )
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
            if self.sbml.getInitialAssignment(ele.getId()) is not None:
                sym = self._get_element_from_assignment(ele.getId())
                if sym is None:
                    sym = sp.sympify(ele.getStoichiometry(),
                                     locals=self.local_symbols)
            elif is_assignment_rule_target(self.sbml, ele):
                sym = _get_identifier_symbol(ele)
            else:
                # dont put the symbol if it wont get replaced by a
                # rule
                sym = sp.sympify(ele.getStoichiometry(),
                                 locals=self.local_symbols)
        elif ele.isSetStoichiometry():
            sym = sp.sympify(ele.getStoichiometry(),
                             locals=self.local_symbols)
        else:
            sym = sp.sympify(1.0)
        sym = _parse_special_functions(sym)
        _check_unsupported_functions(sym, 'Stoichiometry')
        return sym


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
                                 full_sym: sp.Expr = None):
    """
    Recursively checks the symbolic expression for unsupported symbolic
    functions

    :param sym:
        symbolic expressions

    :param expression_type:
        type of expression, only used when throwing errors
    """
    if full_sym is None:
        full_sym = sym

    unsupported_functions = [
        sp.functions.factorial, sp.functions.ceiling, sp.functions.floor,
        sp.function.UndefinedFunction
    ]

    unsupp_fun_type = next(
        (
            fun_type
            for fun_type in unsupported_functions
            if isinstance(sym.func, fun_type)
        ),
        None
    )
    if unsupp_fun_type:
        raise SBMLException(f'Encountered unsupported expression '
                            f'"{sym.func}" of type '
                            f'"{unsupp_fun_type}" as part of a '
                            f'{expression_type}: "{full_sym}"!')
    for fun in list(sym._args) + [sym]:
        unsupp_fun_type = next(
            (
                fun_type
                for fun_type in unsupported_functions
                if isinstance(fun, fun_type)
            ),
            None
        )
        if unsupp_fun_type:
            raise SBMLException(f'Encountered unsupported expression '
                                f'"{fun}" of type '
                                f'"{unsupp_fun_type}" as part of a '
                                f'{expression_type}: "{full_sym}"!')
        if fun is not sym:
            _check_unsupported_functions(fun, expression_type)


def _parse_special_functions(sym: sp.Expr, toplevel: bool = True) -> sp.Expr:
    """
    Recursively checks the symbolic expression for functions which have be
    to parsed in a special way, such as piecewise functions

    :param sym:
        symbolic expressions

    :param toplevel:
        as this is called recursively, are we in the top level expression?
    """
    args = tuple(_parse_special_functions(arg, False) for arg in sym._args)

    if sym.__class__.__name__ == 'abs':
        return sp.Abs(sym._args[0])
    elif sym.__class__.__name__ == 'xor':
        return sp.Xor(*sym.args)
    elif sym.__class__.__name__ == 'piecewise':
        # how many condition-expression pairs will we have?
        return sp.Piecewise(*grouper(args, 2, True))
    elif isinstance(sym, (sp.Function, sp.Mul, sp.Add)):
        sym._args = args
    elif toplevel:
        # Replace boolean constants by numbers so they can be differentiated
        #  must not replace in Piecewise function. Therefore, we only replace
        #  it the complete expression consists only of a Boolean value.
        if isinstance(sym, spTrue):
            sym = sp.Float(1.0)
        elif isinstance(sym, spFalse):
            sym = sp.Float(0.0)

    return sym


def _parse_logical_operators(math_str: str) -> Union[str, None]:
    """
    Parses a math string in order to replace logical operators by a form
    parsable for sympy

    :param math_str:
        str with mathematical expression
    :param math_str:
        parsed math_str
    """
    if math_str is None:
        return None

    if ' xor(' in math_str or ' Xor(' in math_str:
        raise SBMLException('Xor is currently not supported as logical '
                            'operation.')

    return (math_str.replace('&&', '&')).replace('||', '|')


def grouper(iterable: Iterable, n: int,
            fillvalue: Any = None) -> Iterable[Iterable]:
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


def noise_distribution_to_cost_function(
        noise_distribution: str
) -> Callable[[str], str]:
    """
    Parse noise distribution string to a cost function definition amici can
    work with.

    The noise distributions listed in the following are supported. :math:`m`
    denotes the measurement, :math:`y` the simulation, and :math:`\\sigma` a
    distribution scale parameter
    (currently, AMICI only supports a single distribution parameter).

    - `'normal'`, `'lin-normal'`: A normal distribution:

      .. math::
         \\pi(m|y,\\sigma) = \\frac{1}{\\sqrt{2\\pi}\\sigma}\\
         exp\\left(-\\frac{(m-y)^2}{2\\sigma^2}\\right)

    - `'log-normal'`: A log-normal distribution (i.e. log(m) is
      normally distributed):

      .. math::
         \\pi(m|y,\\sigma) = \\frac{1}{\\sqrt{2\\pi}\\sigma m}\\
         exp\\left(-\\frac{(\\log m - \\log y)^2}{2\\sigma^2}\\right)

    - `'log10-normal'`: A log10-normal distribution (i.e. log10(m) is
      normally distributed):

      .. math::
         \\pi(m|y,\\sigma) = \\frac{1}{\\sqrt{2\\pi}\\sigma m \\log(10)}\\
         exp\\left(-\\frac{(\\log_{10} m - \\log_{10} y)^2}{2\\sigma^2}\\right)

    - `'laplace'`, `'lin-laplace'`: A laplace distribution:

      .. math::
         \\pi(m|y,\\sigma) = \\frac{1}{2\\sigma}
         \\exp\\left(-\\frac{|m-y|}{\\sigma}\\right)

    - `'log-laplace'`: A log-Laplace distribution (i.e. log(m) is Laplace
      distributed):

      .. math::
         \\pi(m|y,\\sigma) = \\frac{1}{2\\sigma m}
         \\exp\\left(-\\frac{|\\log m - \\log y|}{\\sigma}\\right)

    - `'log10-laplace'`: A log10-Laplace distribution (i.e. log10(m) is
      Laplace distributed):

      .. math::
         \\pi(m|y,\\sigma) = \\frac{1}{2\\sigma m \\log(10)}
         \\exp\\left(-\\frac{|\\log_{10} m - \\log_{10} y|}{\\sigma}\\right)

    - `'binomial'`, `'lin-binomial'`: A (continuation of a discrete) binomial
      distribution, parameterized via the success probability
      :math:`p=\\sigma`:

      .. math::
         \\pi(m|y,\\sigma) = \\operatorname{Heaviside}(y-m) \\cdot
                \\frac{\\Gamma(y+1)}{\\Gamma(m+1) \\Gamma(y-m+1)}
                \\sigma^m (1-\\sigma)^{(y-m)}

    - `'negative-binomial'`, `'lin-negative-binomial'`: A (continuation of a
      discrete) negative binomial distribution, with with `mean = y`,
      parameterized via success probability `p`:

      .. math::

         \\pi(m|y,\\sigma) = \\frac{\\Gamma(m+r)}{\\Gamma(m+1) \\Gamma(r)}
            (1-\\sigma)^m \\sigma^r

      where

      .. math::
         r = \\frac{1-\\sigma}{\\sigma} y

    The distributions above are for a single data point.
    For a collection :math:`D=\\{m_i\\}_i` of data points and corresponding
    simulations :math:`Y=\\{y_i\\}_i` and noise parameters
    :math:`\\Sigma=\\{\\sigma_i\\}_i`, AMICI assumes independence,
    i.e. the full distributions is

    .. math::
       \\pi(D|Y,\\Sigma) = \\prod_i\\pi(m_i|y_i,\\sigma_i)

    AMICI uses the logarithm :math:`\\log(\\pi(m|y,\\sigma)`.

    In addition to the above mentioned distributions, it is also possible to
    pass a function taking a symbol string and returning a log-distribution
    string with variables '{str_symbol}', 'm{str_symbol}', 'sigma{str_symbol}'
    for y, m, sigma, respectively.

    :param noise_distribution: An identifier specifying a noise model.
        Possible values are

        {`'normal'`, `'lin-normal'`, `'log-normal'`, `'log10-normal'`,
        `'laplace'`, `'lin-laplace'`, `'log-laplace'`, `'log10-laplace'`,
        `'binomial'`, `'lin-binomial'`, `'negative-binomial'`,
        `'lin-negative-binomial'`, `<Callable>`}

        For the meaning of the values see above.

    :return: A function that takes a strSymbol and then creates a cost
        function string (negative log-likelihood) from it, which can be
        sympified.
    """

    if isinstance(noise_distribution, Callable):
        return noise_distribution

    if noise_distribution in ['normal', 'lin-normal']:
        y_string = '0.5*log(2*pi*{sigma}**2) + 0.5*(({y} - {m}) / {sigma})**2'
    elif noise_distribution == 'log-normal':
        y_string = '0.5*log(2*pi*{sigma}**2*{m}**2) ' \
                   '+ 0.5*((log({y}) - log({m})) / {sigma})**2'
    elif noise_distribution == 'log10-normal':
        y_string = '0.5*log(2*pi*{sigma}**2*{m}**2*log(10)**2) ' \
                   '+ 0.5*((log({y}, 10) - log({m}, 10)) / {sigma})**2'
    elif noise_distribution in ['laplace', 'lin-laplace']:
        y_string = 'log(2*{sigma}) + Abs({y} - {m}) / {sigma}'
    elif noise_distribution == 'log-laplace':
        y_string = 'log(2*{sigma}*{m}) + Abs(log({y}) - log({m})) / {sigma}'
    elif noise_distribution == 'log10-laplace':
        y_string = 'log(2*{sigma}*{m}*log(10)) ' \
                   '+ Abs(log({y}, 10) - log({m}, 10)) / {sigma}'
    elif noise_distribution in ['binomial', 'lin-binomial']:
        # Binomial noise model parameterized via success probability p
        y_string = '- log(Heaviside({y} - {m})) - loggamma({y}+1) ' \
                   '+ loggamma({m}+1) + loggamma({y}-{m}+1) ' \
                   '- {m} * log({sigma}) - ({y} - {m}) * log(1-{sigma})'
    elif noise_distribution in ['negative-binomial', 'lin-negative-binomial']:
        # Negative binomial noise model of the number of successes m
        # (data) before r=(1-sigma)/sigma * y failures occur,
        # with mean number of successes y (simulation),
        # parameterized via success probability p = sigma.
        r = '{y} * (1-{sigma}) / {sigma}'
        y_string = f'- loggamma({{m}}+{r}) + loggamma({{m}}+1) ' \
                   f'+ loggamma({r}) - {r} * log(1-{{sigma}}) ' \
                   f'- {{m}} * log({{sigma}})'
    else:
        raise ValueError(
            f"Cost identifier {noise_distribution} not recognized.")

    def nllh_y_string(str_symbol):
        y, m, sigma = _get_str_symbol_identifiers(str_symbol)
        return y_string.format(y=y, m=m, sigma=sigma)

    return nllh_y_string


def _get_str_symbol_identifiers(str_symbol: str) -> tuple:
    """Get identifiers for simulation, measurement, and sigma."""
    y, m, sigma = f"{str_symbol}", f"m{str_symbol}", f"sigma{str_symbol}"
    return y, m, sigma


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


def _get_species_initial(species: sbml.Species) -> sp.Expr:
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
            return sp.sympify(conc) * _get_species_compartment_symbol(species)
        else:
            return sp.sympify(conc)

    if species.isSetInitialAmount():
        amt = species.getInitialAmount()
        if math.isnan(amt):
            return sp.sympify(0.0)

        if species.getHasOnlySubstanceUnits():
            return sp.sympify(amt)
        else:
            return sp.sympify(amt) / _get_species_compartment_symbol(species)

    return sp.sympify(0.0)


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


def symbol_with_assumptions(name: str):
    """
    Central function to create symbols with consistent, canonical assumptions

    :param name:
        name of the symbol

    :return:
        symbol with canonical assumptions
    """
    return sp.Symbol(name, real=True)


class MathMLSbmlPrinter(MathMLContentPrinter):
    """
    Prints a SymPy expression to a MathML expression parsable by libSBML.

    Differences from :class:`sympy.MathMLContentPrinter`:

    1. underscores in symbol names are not converted to subscripts
    2. symbols with name 'time' are converted to the SBML time symbol
    """
    def _print_Symbol(self, sym):
        ci = self.dom.createElement(self.mathml_tag(sym))
        ci.appendChild(self.dom.createTextNode(sym.name))
        return ci

    # _print_Float can be removed when sympy 1.6.3 is released
    def _print_Float(self, expr):
        x = self.dom.createElement(self.mathml_tag(expr))
        repr_expr = mlib_to_str(expr._mpf_, repr_dps(expr._prec))
        x.appendChild(self.dom.createTextNode(repr_expr))
        return x

    def doprint(self, expr):
        mathml_str = super().doprint(expr)
        mathml_str = '<math xmlns="http://www.w3.org/1998/Math/MathML">' + \
                     mathml_str + '</math>'
        mathml_str = mathml_str.replace(
            '<ci>time</ci>',
            '<csymbol encoding="text" definitionURL='
            '"http://www.sbml.org/sbml/symbols/time"> time </csymbol>'
        )
        return mathml_str


def mathml(expr, **settings):
    return MathMLSbmlPrinter(settings).doprint(expr)


def is_assignment_rule_target(model: sbml.Model, element: sbml.SBase) -> bool:
    """
    Checks if an element has a valid assignment rule in the specified model

    :param model:
        SBML model

    :param element:
        SBML variable

    :return:
        boolean indicating truth of function name
    """
    a = model.getAssignmentRuleByVariable(element.getId())
    if a is None:
        return False

    if not a.getFormula():
        return False

    return True
