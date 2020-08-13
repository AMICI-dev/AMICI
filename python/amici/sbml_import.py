"""
SBML Import
-----------
This module provides all necessary functionality to import a model specified
in the System Biology Markup Language (SBML)
"""


import sympy as sp
import libsbml as sbml
import re
import math
import itertools as itt
import warnings
import logging
from typing import Dict, Union, List, Callable, Any, Iterable

from .ode_export import ODEExporter, ODEModel
from .logging import get_logger, log_execution_time, set_log_level
from . import has_clibs

from sympy.logic.boolalg import BooleanTrue as spTrue
from sympy.logic.boolalg import BooleanFalse as spFalse
from sympy.printing.mathml import MathMLContentPrinter

# the following import can be removed if sympy PR #19958 is merged
from mpmath.libmp import repr_dps, to_str as mlib_to_str


class SBMLException(Exception):
    pass


default_symbols = {
    'species': {},
    'parameter': {},
    'fixed_parameter': {},
    'observable': {},
    'expression': {},
    'sigmay': {},
    'my': {},
    'llhy': {},
}

ConservationLaw = Dict[str, Union[str, sp.Basic]]

logger = get_logger(__name__, logging.ERROR)


class SbmlImporter:
    """
    Class to generate AMICI C++ files for a model provided in the Systems
    Biology Markup Language (SBML).

    :ivar show_sbml_warnings: bool
        indicates whether libSBML warnings should be
        displayed

    :ivar symbols:
        dict carrying symbolic definitions

    :ivar sbml_reader:
        the libSBML sbml reader [!not storing this will result
        in a segfault!]

    :ivar sbml_doc:
        document carrying the sbml definition [!not storing this
        will result in a segfault!]

    :ivar sbml:
        sbml definition [!not storing this will result in a segfault!]

    :ivar species_index:
        maps species names to indices

    :ivar species_compartment: sympy.Matrix
        compartment for each species

    :ivar constant_species:
        ids of species that are marked as constant

    :ivar boundary_condition_species:
        ids of species that are marked as boundary
        condition

    :ivar species_has_only_substance_units:
        flags indicating whether a species has only substance units

    :ivar species_conversion_factor:
        conversion factors for every species

    :ivar compartment_symbols
        compartment ids

    :ivar compartment_volume:
        numeric/symbolic compartment volumes

    :ivar stoichiometric_matrix:
        stoichiometric matrix of the model

    :ivar flux_vector:
        reaction kinetic laws

    :ivar local_symbols:
        model symbols for sympy to consider during sympification
        see `locals`argument in `sympy.sympify`

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
            self.sbml_doc = sbml_source.getSBMLDocument()
        else:
            self.sbml_reader = sbml.SBMLReader()
            if from_file:
                sbml_doc = self.sbml_reader.readSBMLFromFile(sbml_source)
            else:
                sbml_doc = self.sbml_reader.readSBMLFromString(sbml_source)
            self.sbml_doc = sbml_doc

        self.show_sbml_warnings : bool = show_sbml_warnings

        # process document
        self._process_document()

        self.sbml = self.sbml_doc.getModel()

        # Long and short names for model components
        self.symbols = dict()
        self._reset_symbols()

        self.local_symbols : dict = {}

        self.compartment_rate_rules : dict = {}
        self.species_rate_rules : dict = {}
        self.compartment_assignment_rules : dict = {}
        self.species_assignment_rules : dict = {}

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
        if len(self.sbml_doc.getModel().getListOfFunctionDefinitions()) > 0:
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
        self.symbols = default_symbols

    def sbml2amici(self,
                   model_name: str = None,
                   output_dir: str = None,
                   observables: Dict[str, Dict[str, str]] = None,
                   constant_parameters: List[str] = None,
                   sigmas: Dict[str, Union[str, float]] = None,
                   noise_distributions: Dict[str, str] = None,
                   verbose: Union[int, bool] = logging.ERROR,
                   assume_pow_positivity: bool = False,
                   compiler: str = None,
                   allow_reinit_fixpar_initcond: bool = True,
                   compile: bool = True,
                   compute_conservation_laws: bool = True,
                   simplify: Callable = lambda x: sp.powsimp(x, deep=True),
                   **kwargs) -> None:
        """
        Generate AMICI C++ files for the model provided to the constructor.

        The resulting model can be imported as a regular Python module (if
        `compile=True`), or used from Matlab or C++ as described in the
        documentation of the respective AMICI interface.

        Note that this generates model ODEs for changes in concentrations, not
        amounts. The simulation results obtained from the model will be
        concentrations, independently of the SBML `hasOnlySubstanceUnits`
        attribute.

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
            assumed as default.

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

        self._process_time()
        self._clean_reserved_symbols()
        self._replace_special_constants()

        ode_model = ODEModel(simplify=simplify)
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
        self._process_volume_conversion()

    def check_support(self) -> None:
        """
        Check whether all required SBML features are supported.
        Also ensures that the SBML contains at least one reaction, or rate
        rule, or assignment rule, to produce change in the system over time.
        """
        if len(self.sbml.getListOfSpecies()) == 0:
            raise SBMLException('Models without species '
                                'are currently not supported!')

        if hasattr(self.sbml, 'all_elements_from_plugins') \
                and self.sbml.all_elements_from_plugins.getSize() > 0:
            raise SBMLException('SBML extensions are currently not supported!')

        if len(self.sbml.getListOfEvents()) > 0:
            raise SBMLException('Events are currently not supported!')

        # Contains condition to allow compartment rate rules
        compartment_ids = list(map(lambda x: x.getId(),
                                   self.sbml.getListOfCompartments()))
        species_ids = list(map(lambda x: x.getId(),
                               self.sbml.getListOfSpecies()))
        if any([not(rule.isAssignment()) and
                not(rule.getVariable() in compartment_ids) and
                not(rule.getVariable() in species_ids)
                for rule in self.sbml.getListOfRules()]):
            raise SBMLException('Algebraic rules are currently not supported, '
                                'and rate rules are only supported for '
                                'species and compartments.')

        if any([not(rule.isAssignment() or rule.isRate()) and
                (rule.getVariable() in compartment_ids)
                for rule in self.sbml.getListOfRules()]):
            raise SBMLException('Only assignment and rate rules are currently '
                                'supported for compartments!')

        if any([not(rule.isAssignment() or rule.isRate()) and
                (rule.getVariable() in species_ids)
                for rule in self.sbml.getListOfRules()]):
            raise SBMLException('Only assignment and rate rules are currently '
                                'supported for species!')

        if any([reaction.getFast()
                for reaction in self.sbml.getListOfReactions()]):
            raise SBMLException('Fast reactions are currently not supported!')

        if any([any([not element.getStoichiometryMath() is None
                for element in list(reaction.getListOfReactants())
                + list(reaction.getListOfProducts())])
                for reaction in self.sbml.getListOfReactions()]):
            raise SBMLException('Non-unity stoichiometry is'
                                ' currently not supported!')

    def _gather_locals(self) -> None:
        """
        Populate self.local_symbols with all model entities.

        This is later used during sympifications to avoid sympy builtins
        shadowing model entities.
        """
        for s in self.sbml.getListOfSpecies():
            self.local_symbols[s.getId()] = sp.Symbol(s.getId(), real=True)

        for p in self.sbml.getListOfParameters():
            self.local_symbols[p.getId()] = sp.Symbol(p.getId(), real=True)

        for c in self.sbml.getListOfCompartments():
            self.local_symbols[c.getId()] = sp.Symbol(c.getId(), real=True)

        for r in self.sbml.getListOfRules():
            self.local_symbols[r.getVariable()] = sp.Symbol(r.getVariable(),
                                                            real=True)

        # SBML time symbol + constants
        self.local_symbols['time'] = sp.Symbol('time', real=True)
        self.local_symbols['avogadro'] = sp.Symbol('avogadro', real=True)

    @log_execution_time('processing SBML compartments', logger)
    def _process_compartments(self) -> None:
        """
        Get compartment information, stoichiometric matrix and fluxes from
        SBML model.
        """
        compartments = self.sbml.getListOfCompartments()
        self.compartment_symbols = sp.Matrix(
            [sp.Symbol(comp.getId(), real=True) for comp in compartments]
        )

        # Initial volumes may be overridden at the end of _process_species,
        # where compartment assignment rules are processed.
        self.compartment_volume = sp.Matrix([
            sp.sympify(comp.getVolume()) if comp.isSetVolume()
            else sp.sympify(1.0) for comp in compartments
        ])

        compartment_ids = [comp.getId() for comp in compartments]
        for initial_assignment in self.sbml.getListOfInitialAssignments():
            if initial_assignment.getId() in compartment_ids:
                index = compartment_ids.index(
                        initial_assignment.getId()
                    )
                self.compartment_volume[index] = sp.sympify(
                    sbml.formulaToL3String(initial_assignment.getMath()),
                    locals=self.local_symbols
                )

    @log_execution_time('processing SBML species', logger)
    def _process_species(self) -> None:
        """
        Get species information from SBML model.
        """
        species = self.sbml.getListOfSpecies()

        self.species_index = {
            species_element.getId(): species_index
            for species_index, species_element in enumerate(species)
        }

        self.symbols['species']['identifier'] = sp.Matrix(
            [sp.Symbol(spec.getId(), real=True) for spec in species]
        )
        self.symbols['species']['name'] = [spec.getName() for spec in species]

        self.species_compartment = sp.Matrix(
            [sp.Symbol(spec.getCompartment(), real=True) for spec in species]
        )

        self.constant_species = [species_element.getId()
                                 for species_element in species
                                 if species_element.getConstant()]

        self.boundary_condition_species = [
            species_element.getId()
            for species_element in species
            if species_element.getBoundaryCondition()
        ]
        self.species_has_only_substance_units = [
            specie.getHasOnlySubstanceUnits() for specie in species
        ]

        concentrations = [spec.getInitialConcentration() for spec in species]
        amounts = [spec.getInitialAmount() for spec in species]

        def get_species_initial(index, conc):
            # We always simulate concentrations!
            if self.species_has_only_substance_units[index]:
                if species[index].isSetInitialAmount() \
                        and not math.isnan(amounts[index]):
                    return sp.sympify(amounts[index]) \
                           / self.species_compartment[index]
                if species[index].isSetInitialConcentration():
                    return sp.sympify(conc)
            else:
                if species[index].isSetInitialConcentration():
                    return sp.sympify(conc)

                if species[index].isSetInitialAmount() \
                        and not math.isnan(amounts[index]):
                    return sp.sympify(amounts[index]) \
                           / self.species_compartment[index]

            return self.symbols['species']['identifier'][index]

        species_initial = sp.Matrix(
            [get_species_initial(index, conc)
             for index, conc in enumerate(concentrations)]
        )

        species_ids = [spec.getId() for spec in self.sbml.getListOfSpecies()]
        for initial_assignment in self.sbml.getListOfInitialAssignments():
            if initial_assignment.getId() in species_ids:
                index = species_ids.index(
                        initial_assignment.getId()
                    )
                sym_math = sp.sympify(_parse_logical_operators(
                    sbml.formulaToL3String(initial_assignment.getMath())),
                    locals=self.local_symbols
                )
                if sym_math is not None:
                    sym_math = _parse_special_functions(sym_math)
                    _check_unsupported_functions(sym_math, 'InitialAssignment')
                    species_initial[index] = sym_math

        for ix, (symbol, init) in enumerate(zip(
            self.symbols['species']['identifier'], species_initial
        )):
            if symbol == init:
                species_initial[ix] = sp.sympify(0.0)

        # flatten initSpecies
        while any([species in species_initial.free_symbols
                   for species in self.symbols['species']['identifier']]):
            species_initial = species_initial.subs([
                (symbol, init)
                for symbol, init in zip(
                    self.symbols['species']['identifier'], species_initial
                )
            ])

        self.symbols['species']['value'] = species_initial

        if self.sbml.isSetConversionFactor():
            conversion_factor = sp.Symbol(self.sbml.getConversionFactor(),
                                          real=True)
        else:
            conversion_factor = 1.0

        self.species_conversion_factor = sp.Matrix([
             sp.sympify(specie.getConversionFactor())
             if specie.isSetConversionFactor()
             else conversion_factor
             for specie in species
        ])

        # Process assignment and rate rules for species and compartments.
        # Compartments with rate rules are implemented as species. Species and
        # compartments with assignments are implemented as observables (and
        # replaced with their assignment in all expressions). Note that, in the
        # case of species, rate rules may describe the change in amount, not
        # concentration, of a species.
        rules = self.sbml.getListOfRules()
        compartmentvars = self.compartment_symbols.free_symbols
        # compartments with rules are replaced with constants in the relevant
        # equations during the _replace_in_all_expressions call inside
        # _process_rules
        for rule in rules:
            ### Taken from _process_rules
            if rule.getFormula() == '':
                continue
            variable = sp.sympify(rule.getVariable(),
                                  locals=self.local_symbols)
            formula = sp.sympify(_parse_logical_operators(
                sbml.formulaToL3String(rule.getMath())),
                locals=self.local_symbols)
            formula = _parse_special_functions(formula)
            _check_unsupported_functions(formula, 'Rule')
            formula = self._replace_reactions_in_rule_formula(rule, formula)
            ###

            # Species rules are processed first, to avoid processing
            # compartments twice (as compartments with rate rules are
            # implemented as species). Could also be avoided with a
            # `not in self.compartment_rate_rules` condition.
            if variable in self.symbols['species']['identifier']:
                if rule.getTypeCode() == sbml.SBML_ASSIGNMENT_RULE:
                    # Handled in _process_rules and _process_observables.
                    pass
                elif rule.getTypeCode() == sbml.SBML_RATE_RULE:
                    self.add_d_dt(
                        formula,
                        variable,
                        self.symbols['species']['value'],
                        sbml.SBML_SPECIES)
                else:
                    raise SBMLException('The only rules currently supported '
                                        'for species are assignment and rate '
                                        'rules!')

            if variable in compartmentvars:
                if rule.getTypeCode() == sbml.SBML_ASSIGNMENT_RULE:
                    # Handled in _process_rules and _process_observables
                    # SBML Assignment Rules can be used to specify initial
                    # values (see SBML L3V2 manual, Section 3.4.8).
                    # Priority appears to be above InitialAssignment.
                    self.compartment_volume[list(
                        self.compartment_symbols
                    ).index(variable)] = formula
                elif rule.getTypeCode() == sbml.SBML_RATE_RULE:
                    self.add_d_dt(
                        formula,
                        variable,
                        self.compartment_volume[list(
                            self.compartment_symbols
                        ).index(variable)],
                        sbml.SBML_COMPARTMENT)
                else:
                    raise SBMLException('The only rules currently supported '
                                        'for compartments are assignment and '
                                        'rate rules!')

    def add_d_dt(
            self,
            d_dt: sp.Expr,
            variable: sp.Symbol,
            variable0: Union[float, sp.Expr],
            component: int,
            name: str = None
    ) -> None:
        '''
        Creates or modifies species, to implement rate rules for compartments
        and species, respectively.

        :param d_dt:
            The rate rule (or, right-hand side of an ODE).

        :param variable:
            The subject of the rate rule.

        :param variable0:
            The initial value of the variable.

        :param component:
            The type of SBML component. Currently, only libsbml.SBML_SPECIES
            and libsbml.SBML_COMPARTMENT are supported.
        '''
        if name is None:
            name = ''

        ### d_dt may contain speciesReference symbols, that may be defined in
        # an initial assignment (e.g. see SBML test suite case 1498, which
        # uses a speciesReference Id in a species rate rule).
        # Here, such speciesReference symbols are replaced with the initial
        # assignment expression, if the expression is a constant (time-
        # dependent expression symbols are not evaluated at zero, rather raise
        # an error).
        # One method to implement expressions with time-dependent symbols
        # may be to produce a dictionary of speciesReference symbols and
        # their initial assignment expressions here, then add this dictionary
        # to the _replace_in_all_expressions method. After _process_sbml,
        # substitute initial values in for any remaining symbols, evaluate the
        # the expressions at $t=0$ (self.amici_time_symbol), then substitute
        # them into d_dt.

        # Initial assignment symbols may be compartments, species, parameters,
        # speciesReferences, or an (extension?) package element. Here, it is
        # assumed that a symbol is a speciesReference if it is not a
        # compartment, species, or parameter, and is the symbol of an initial
        # assignment.
        alternative_components = [s.getId() for s in
            list(self.sbml.getListOfCompartments()) +\
            list(self.sbml.getListOfSpecies()) +\
            list(self.sbml.getListOfParameters())]
        initial_assignments = {ia.getId(): ia for ia in
                                       self.sbml.getListOfInitialAssignments()}
        for symbol in d_dt.free_symbols:
            if str(symbol) not in alternative_components and\
               str(symbol) in initial_assignments:
                # Taken from _process_species
                sym_math = sp.sympify(_parse_logical_operators(
                    sbml.formulaToL3String(initial_assignments[
                        str(symbol)
                    ].getMath())),
                    locals=self.local_symbols
                )
                if sym_math is not None:
                    sym_math = _parse_special_functions(sym_math)
                    _check_unsupported_functions(sym_math, 'InitialAssignment')

                if not isinstance(sym_math, sp.Float):
                    raise SBMLException('Rate rules that contain '
                                        'speciesReferences, defined in '
                                        'initial assignments that contain '
                                        'symbols, are currently not '
                                        'supported! Rate rule symbol: '
                                        f'{variable}, species reference '
                                        f'symbol: {symbol}, initial '
                                        f'assignment: {sym_math}, type: '
                                        f'{type(sym_math)}.')
                else:
                    d_dt = d_dt.subs(symbol, sym_math)
        ###

        if component == sbml.SBML_COMPARTMENT:
            self.symbols['species']['identifier'] = \
                self.symbols['species']['identifier'].col_join(
                    sp.Matrix([variable]))
            self.symbols['species']['name'].append(name)
            self.symbols['species']['value'] = \
                self.symbols['species']['value'].col_join(
                    sp.Matrix([variable0]))
            self.species_index[str(variable)] = len(self.species_index)
            self.compartment_rate_rules[variable] = d_dt
        elif component == sbml.SBML_SPECIES:
            if str(variable) in self.species_index:
                # SBML species are already in the species symbols
                x_index = self.species_index[str(variable)]
                if self.species_has_only_substance_units[x_index]:
                    self.symbols['species']['value'][x_index] *= \
                        self.species_compartment[x_index]
                self.species_rate_rules[variable] = d_dt
            else:
                # For Jakob to specify rules independently of SBML (untested)
                # Remove assert statement to use
                # Currently untested...
                assert False, "Error while processing species rate rules."
                self.symbols['species']['identifier'] = \
                    self.symbols['species']['identifier'].col_join(
                        sp.Matrix([variable]))
                self.symbols['species']['name'].append(name)
                self.symbols['species']['value'] = \
                    self.symbols['species']['value'].col_join(
                        sp.Matrix([variable0]))
                self.species_index[str(variable)] = len(self.species_index)
                self.species_rate_rules[variable] = d_dt
                # May need to modify this function (`add_d_dt`) to allow
                # specification of whether the argument `d_dt` describes change
                # in amount or concentration.
                #if self.species_has_only_substance_units[x_index]:
                #    self.symbols['species']['value'][x_index] *= \
                #    self.species_compartment[x_index]
        else:
            raise TypeError(f'Rate rules are currently only supported for '
                            'libsbml.SBML_COMPARTMENT and '
                            'libsbml.SBML_SPECIES components.')

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

        parameter_ids = [par.getId() for par
                         in self.sbml.getListOfParameters()]
        for initial_assignment in self.sbml.getListOfInitialAssignments():
            if initial_assignment.getId() in parameter_ids:
                raise SBMLException('Initial assignments for parameters are'
                                    ' currently not supported')

        fixed_parameters = [
            parameter
            for parameter in self.sbml.getListOfParameters()
            if parameter.getId() in constant_parameters
        ]

        rulevars = [rule.getVariable() for rule in self.sbml.getListOfRules()]

        parameters = [parameter for parameter
                      in self.sbml.getListOfParameters()
                      if parameter.getId() not in constant_parameters
                      and parameter.getId() not in rulevars]

        loop_settings = {
            'parameter': {
                'var': parameters,
                'name': 'parameter',

            },
            'fixed_parameter': {
                'var': fixed_parameters,
                'name': 'fixed_parameter'
            }

        }

        for partype, settings in loop_settings.items():
            self.symbols[partype]['identifier'] = sp.Matrix(
                [sp.Symbol(par.getId(), real=True) for par in settings['var']]
            )
            self.symbols[partype]['name'] = [
                par.getName() for par in settings['var']
            ]
            self.symbols[partype]['value'] = [
                par.getValue() for par in settings['var']
            ]
            setattr(
                self,
                f'{settings["name"]}_index',
                {
                    parameter_element.getId(): parameter_index
                    for parameter_index, parameter_element
                    in enumerate(settings['var'])
                }
            )

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
        nx = len(self.symbols['species']['name'])
        # stoichiometric matrix
        self.stoichiometric_matrix = sp.SparseMatrix(sp.zeros(nx, nr))
        self.flux_vector = sp.zeros(nr, 1)

        assignment_ids = [ass.getId()
                          for ass in self.sbml.getListOfInitialAssignments()]
        rulevars = [rule.getVariable()
                    for rule in self.sbml.getListOfRules()
                    if rule.getFormula() != '']

        reaction_ids = [
            reaction.getId() for reaction in reactions
            if reaction.isSetId()
        ]

        def get_element_from_assignment(element_id):
            assignment = self.sbml.getInitialAssignment(
                element_id
            )
            sym = sp.sympify(sbml.formulaToL3String(assignment.getMath()),
                             locals=self.local_symbols)
            # this is an initial assignment so we need to use
            # initial conditions
            if sym is not None:
                sym = sym.subs(
                    self.symbols['species']['identifier'],
                    self.symbols['species']['value']
                )
            return sym

        def get_element_stoichiometry(ele):
            if ele.isSetId():
                if ele.getId() in assignment_ids:
                    sym = get_element_from_assignment(ele.getId())
                    if sym is None:
                        sym = sp.sympify(ele.getStoichiometry())
                elif ele.getId() in rulevars:
                    return sp.Symbol(ele.getId(), real=True)
                else:
                    # dont put the symbol if it wont get replaced by a
                    # rule
                    sym = sp.sympify(ele.getStoichiometry())
            elif ele.isSetStoichiometry():
                sym = sp.sympify(ele.getStoichiometry())
            else:
                return sp.sympify(1.0)
            sym = _parse_special_functions(sym)
            _check_unsupported_functions(sym, 'Stoichiometry')
            return sym

        def is_constant(specie):
            return specie in self.constant_species or \
                   specie in self.boundary_condition_species

        for reaction_index, reaction in enumerate(reactions):
            for elementList, sign in [(reaction.getListOfReactants(), -1.0),
                                      (reaction.getListOfProducts(), 1.0)]:
                elements = {}
                for index, element in enumerate(elementList):
                    # we need the index here as we might have multiple elements
                    # for the same species
                    elements[index] = {'species': element.getSpecies()}
                    elements[index]['stoichiometry'] = \
                        get_element_stoichiometry(element)

                for index in elements.keys():
                    if not is_constant(elements[index]['species']):
                        specie_index = self.species_index[
                            elements[index]['species']
                        ]
                        # Division by species compartment size (to find the
                        # rate of change in species concentration) now occurs
                        # in the `dx_dt` method in "ode_export.py", which also
                        # accounts for possibly variable compartments.
                        self.stoichiometric_matrix[specie_index,
                                                   reaction_index] += \
                            sign \
                            * elements[index]['stoichiometry'] \
                            * self.species_conversion_factor[specie_index]

            # usage of formulaToL3String ensures that we get "time" as time
            # symbol
            kmath = sbml.formulaToL3String(reaction.getKineticLaw().getMath())
            try:
                sym_math = sp.sympify(_parse_logical_operators(kmath),
                                      locals=self.local_symbols)
            except SBMLException as Ex:
                raise Ex
            except sp.SympifyError:
                raise SBMLException(f'Kinetic law "{kmath}" contains an '
                                    'unsupported expression!')
            sym_math = _parse_special_functions(sym_math)
            _check_unsupported_functions(sym_math, 'KineticLaw')
            for r in reactions:
                elements = list(r.getListOfReactants()) \
                           + list(r.getListOfProducts())
                for element in elements:
                    if element.isSetId() & element.isSetStoichiometry():
                        sym_math = sym_math.subs(
                            sp.sympify(element.getId(),
                                       locals=self.local_symbols),
                            sp.sympify(element.getStoichiometry())
                        )

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
        rules = self.sbml.getListOfRules()

        rulevars = get_rule_vars(rules, local_symbols=self.local_symbols)
        fluxvars = self.flux_vector.free_symbols
        specvars = self.symbols['species']['identifier'].free_symbols
        volumevars = self.compartment_volume.free_symbols
        compartmentvars = self.compartment_symbols.free_symbols
        parametervars = sp.Matrix([
            sp.Symbol(par.getId(), real=True)
            for par in self.sbml.getListOfParameters()
        ])
        stoichvars = self.stoichiometric_matrix.free_symbols

        assignments = {}

        for rule in rules:
            # Rate rules should not be substituted for the target of the rate
            # rule.
            if rule.getTypeCode() == sbml.SBML_RATE_RULE:
                continue
            if rule.getFormula() == '':
                continue
            variable = sp.sympify(rule.getVariable(),
                                  locals=self.local_symbols)
            # avoid incorrect parsing of pow(x, -1) in symengine
            formula = sp.sympify(_parse_logical_operators(
                sbml.formulaToL3String(rule.getMath())),
                locals=self.local_symbols)
            formula = _parse_special_functions(formula)
            _check_unsupported_functions(formula, 'Rule')
            formula = self._replace_reactions_in_rule_formula(rule, formula)

            if variable in stoichvars:
                self.stoichiometric_matrix = \
                    self.stoichiometric_matrix.subs(variable, formula)

            if variable in specvars:
                if rule.getTypeCode() == sbml.SBML_ASSIGNMENT_RULE:
                    self.species_assignment_rules[variable] = formula
                    assignments[str(variable)] = formula
                else:
                    # Rate rules are handled in _process_species, and are
                    # skipped in this loop
                    raise KeyError('Only assignment and rate rules are '
                                   'currently supported for species!')

            if variable in compartmentvars:
                if rule.getTypeCode() == sbml.SBML_ASSIGNMENT_RULE:
                    self.compartment_assignment_rules[variable] = formula
                    assignments[str(variable)] = formula
                else:
                    # Rate rules are handled in _process_species, and are
                    # skipped in this loop
                    raise KeyError('Only assignment and rate rules are '
                                   'currently supported for compartments!')

            if variable in parametervars:
                if str(variable) in self.parameter_index:
                    idx = self.parameter_index[str(variable)]
                    self.symbols['parameter']['value'][idx] \
                        = float(formula)
                else:
                    self.sbml.removeParameter(str(variable))
                    assignments[str(variable)] = formula

            if variable in fluxvars:
                self.flux_vector = self.flux_vector.subs(variable, formula)

            if variable in volumevars:
                self.compartment_volume = \
                    self.compartment_volume.subs(variable, formula)

            if variable in rulevars:
                for nested_rule in rules:

                    nested_formula = sp.sympify(
                        _parse_logical_operators(sbml.formulaToL3String(nested_rule.getMath())),
                        locals=self.local_symbols).subs(variable, formula)
                    nested_formula = _parse_special_functions(nested_formula)
                    _check_unsupported_functions(nested_formula, 'Rule')

                    nested_rule_math_ml = mathml(nested_formula)
                    nested_rule_math_ml_ast_node = sbml.readMathMLFromString(nested_rule_math_ml)

                    if nested_rule_math_ml_ast_node is None:
                        raise SBMLException(f'Formula {sbml.formulaToL3String(nested_rule.getMath())}'
                                            f' cannot be correctly read by SymPy'
                                            f' or cannot be converted to valid MathML by SymPy!')

                    elif nested_rule.setMath(nested_rule_math_ml_ast_node) != sbml.LIBSBML_OPERATION_SUCCESS:
                        raise SBMLException(f'Formula {sbml.formulaToL3String(nested_rule.getMath())}'
                                            f' cannot be parsed by libSBML!')

                for assignment in assignments:
                    assignments[assignment] = assignments[assignment].subs(variable, formula)

        # do this at the very end to ensure we have flattened all recursive
        # rules
        for variable in assignments.keys():
            self._replace_in_all_expressions(
                sp.Symbol(variable, real=True),
                assignments[variable]
            )
        for comp, vol in zip(self.compartment_symbols,
                             self.compartment_volume):
            self._replace_in_all_expressions(
                comp, vol
            )

    def _process_volume_conversion(self) -> None:
        """
        Convert equations from amount to volume.
        """
        compartments = self.species_compartment
        for comp, vol in zip(self.compartment_symbols,
                             self.compartment_volume):
            if comp not in self.compartment_rate_rules:
                compartments = compartments.subs(comp, vol)
        for index, sunits in enumerate(self.species_has_only_substance_units):
            if sunits:
                self.flux_vector = \
                    self.flux_vector.subs(
                        self.symbols['species']['identifier'][index],
                        self.symbols['species']['identifier'][index]
                        * compartments[index]
                    )

    def _process_time(self) -> None:
        """
        Convert time_symbol into cpp variable.
        """
        sbml_time_symbol = sp.Symbol('time', real=True)
        amici_time_symbol = sp.Symbol('t', real=True)
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
            See :func:`sbml2amici`.
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

        species_syms = self.symbols['species']['identifier']
        assignments = {str(c): str(r)
                       for c, r in self.compartment_assignment_rules.items()}
        assignments.update({str(s): str(r)
                            for s, r in self.species_assignment_rules.items()})

        def replace_assignments(formula: str) -> sp.Basic:
            """
            Replace assignment rules in observables

            :param formula:
                algebraic formula of the observable

            :return:
                observable formula with assignment rules replaced
            """
            formula = sp.sympify(formula, locals=self.local_symbols)
            for s in formula.free_symbols:
                r = self.sbml.getAssignmentRuleByVariable(str(s))
                if r is not None:
                    rule_formula = _parse_logical_operators(
                        sbml.formulaToL3String(r.getMath()))
                    rule_formula = sp.sympify(
                        rule_formula, locals=self.local_symbols)
                    rule_formula = _parse_special_functions(rule_formula)
                    _check_unsupported_functions(rule_formula, 'Rule')
                    formula = formula.replace(s, rule_formula)
            return formula

        # add user-provided observables or make all species, and compartments
        # with assignment rules, observable
        if observables:
            # Replace logX(.) by log(., X) since symengine cannot parse the
            # former. Also replace symengine-incompatible sbml log(basis, x)
            for observable in observables:
                observables[observable]['formula'] = re.sub(
                    r'(^|\W)log(\d+)\(', r'\g<1>1/ln(\2)*ln(',
                    observables[observable]['formula']
                )
                repl = replaceLogAB(observables[observable]['formula'])
                if repl != observables[observable]['formula']:
                    warnings.warn(
                        f'Replaced "{observables[observable]["formula"]}" by '
                        f'"{repl}", assuming first argument to log() was the '
                        f'basis.'
                    )
                    observables[observable]['formula'] = repl

            observable_values = sp.Matrix([
                replace_assignments(observables[observable]['formula'])
                for observable in observables
            ])
            observable_names = [
                observables[observable]['name'] if 'name' in observables[
                    observable].keys()
                else f'y{index}'
                for index, observable in enumerate(observables)
            ]
            observable_syms = sp.Matrix([
                sp.symbols(obs, real=True) for obs in observables.keys()
            ])
            observable_ids = observables.keys()
        else:
            observable_values = species_syms.copy() # prefer sympy's copy over deepcopy, see sympy issue #7672
            observable_ids = [
                f'x{index}' for index in range(len(species_syms))
            ]
            observable_names = observable_ids[:]
            observable_syms = sp.Matrix(
                [sp.symbols(f'y{index}', real=True)
                 for index in range(len(species_syms))]
            )
            # Add compartment and species assignment rules as observables
            # Useful for passing the SBML Test Suite (compartment volumes are
            # used to calculate species amounts).
            # The id's and names below may conflict with the automatically
            # generated id's and names above.
            for compartment in self.compartment_assignment_rules:
                observable_values = observable_values.col_join(sp.Matrix(
                    [self.compartment_assignment_rules[compartment]]))
                observable_ids.append(str(compartment))
                observable_names.append(str(compartment))
                observable_syms = observable_syms.col_join(sp.Matrix(
                    [compartment]))
            for species in self.species_assignment_rules:
                x_index = self.species_index[str(species)]
                observable_values[x_index] = sp.Matrix(
                    [replace_assignments(str(species))])
                observable_ids[x_index] = str(species)
                observable_names[x_index] = str(species)
                observable_syms[x_index] = sp.Matrix([species])

        sigma_y_syms = sp.Matrix(
            [sp.symbols(f'sigma{symbol}', real=True)
             for symbol in observable_syms]
        )
        sigma_y_values = sp.Matrix(
            [1.0] * len(observable_syms)
        )

        # set user-provided sigmas
        for iy, obs_name in enumerate(observables):
            if obs_name in sigmas:
                sigma_y_values[iy] = replace_assignments(sigmas[obs_name])

        measurement_y_syms = sp.Matrix(
            [sp.symbols(f'm{symbol}', real=True) for symbol in observable_syms]
        )
        measurement_y_values = sp.Matrix(
            [0.0] * len(observable_syms)
        )

        # set cost functions
        llh_y_strings = []
        for y_name in observable_ids:
            llh_y_strings.append(noise_distribution_to_cost_function(
                noise_distributions.get(y_name, 'normal')))

        llh_y_values = []
        for llhYString, o_sym, m_sym, s_sym \
                in zip(llh_y_strings, observable_syms,
                       measurement_y_syms, sigma_y_syms):
            f = sp.sympify(llhYString(o_sym), locals={str(o_sym): o_sym,
                                                      str(m_sym): m_sym,
                                                      str(s_sym): s_sym})
            llh_y_values.append(f)
        llh_y_values = sp.Matrix(llh_y_values)

        llh_y_syms = sp.Matrix(
            [sp.Symbol(f'J{symbol}', real=True) for symbol in observable_syms]
        )

        # set symbols
        self.symbols['observable']['identifier'] = observable_syms
        self.symbols['observable']['name'] = l2s(observable_names)
        self.symbols['observable']['value'] = observable_values
        self.symbols['sigmay']['identifier'] = sigma_y_syms
        self.symbols['sigmay']['name'] = l2s(sigma_y_syms)
        self.symbols['sigmay']['value'] = sigma_y_values
        self.symbols['my']['identifier'] = measurement_y_syms
        self.symbols['my']['name'] = l2s(measurement_y_syms)
        self.symbols['my']['value'] = measurement_y_values
        self.symbols['llhy']['value'] = llh_y_values
        self.symbols['llhy']['name'] = l2s(llh_y_syms)
        self.symbols['llhy']['identifier'] = llh_y_syms

    def process_conservation_laws(self, ode_model, volume_updates) -> List:
        """
        Find conservation laws in reactions and species.

        :param ode_model:
            ODEModel object with basic definitions

        :param volume_updates:
            List with updates for the stoichiometrix matrix accounting for
            compartment volumes

        :returns volume_updates_solver:
            List (according to reduced stoichiometry) with updates for the
            stoichiometrix matrix accounting for compartment volumes
        """
        conservation_laws = []

        # So far, only conservation laws for constant species are supported
        species_solver = self._add_conservation_for_constant_species(
            ode_model, conservation_laws)

        # Check, whether species_solver is empty now. As currently, AMICI
        # cannot handle ODEs without species, CLs must switched in this case
        if len(species_solver) == 0:
            conservation_laws = []
            species_solver = list(range(ode_model.nx_rdata()))

        # prune out species from stoichiometry and
        volume_updates_solver = self._reduce_stoichiometry(species_solver,
                                                           volume_updates)

        # add the found CLs to the ode_model
        for cl in conservation_laws:
            ode_model.add_conservation_law(**cl)

        return volume_updates_solver

    def _add_conservation_for_constant_species(
            self,
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
        species_solver = list(range(ode_model.nx_rdata()))

        # iterate over species, find constant ones
        for ix in reversed(range(ode_model.nx_rdata())):
            if ode_model.state_is_constant(ix):
                # dont use sym('x') here since conservation laws need to be
                # added before symbols are generated
                target_state = ode_model._states[ix].get_id()
                total_abundance = sp.Symbol(f'tcl_{target_state}')
                conservation_laws.append({
                    'state': target_state,
                    'total_abundance': total_abundance,
                    'state_expr': total_abundance,
                    'abundance_expr': target_state,
                })
                # mark species to delete from stoichiometrix matrix
                species_solver.pop(ix)

        return species_solver

    def _reduce_stoichiometry(self, species_solver, volume_updates) -> List:
        """
        Reduces the stoichiometry with respect to conserved quantities

        :param species_solver:
            List of species indices which remain later in the ODE solver

        :param volume_updates:
            List with updates for the stoichiometrix matrix accounting for
            compartment volumes

        :returns volume_updates_solver:
            List (according to reduced stoichiometry) with updates for the
            stoichiometrix matrix accounting for compartment volumes
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


    def _replace_in_all_expressions(self,
                                    old: sp.Symbol,
                                    new: sp.Symbol) -> None:
        """
        Replace 'old' by 'new' in all symbolic expressions.

        :param old:
            symbolic variables to be replaced

        :param new:
            replacement symbolic variables
        """
        # Avoid replacing variables with rates
        if old not in set((*self.compartment_rate_rules.keys(),
                           *self.species_rate_rules.keys())):
            fields = [
                'stoichiometric_matrix', 'flux_vector',
            ]
            for field in fields:
                if field in dir(self):
                    self.__setattr__(field, self.__getattribute__(field).subs(
                        old, new
                    ))
            for compartment in self.compartment_rate_rules:
                self.compartment_rate_rules[compartment] = \
                    self.compartment_rate_rules[compartment].subs(old, new)
            for species in self.species_rate_rules:
                self.species_rate_rules[species] = \
                    self.species_rate_rules[species].subs(old, new)

        symbols = [
            'species', 'observable',
        ]
        for symbol in symbols:
            if symbol in self.symbols:
                # Initial species values that are specified as amounts need to
                # be divided by their compartment volume to obtain
                # concentration. The condition below ensures that species
                # initial amount is divided by the initial compartment size,
                # and not the expression for a compartment assignment rule.
                if symbol == 'species' and (
                        old in self.compartment_assignment_rules):
                    comp_v0 = self.compartment_volume[
                        list(self.compartment_symbols).index(old)]
                    self.symbols[symbol]['value'] = \
                        self.symbols[symbol]['value'].subs(old, comp_v0)
                else:
                    # self.symbols['observable'] may not yet be defined.
                    if not self.symbols[symbol]:
                        continue
                    self.symbols[symbol]['value'] = \
                        self.symbols[symbol]['value'].subs(old, new)

        for compartment in self.compartment_assignment_rules:
            self.compartment_assignment_rules[compartment] = \
                self.compartment_assignment_rules[compartment].subs(old, new)
        # Initial compartment volume may also be specified with an assignment
        # rule (at the end of the _process_species method), hence needs to be
        # processed here too.
        for index in range(len(self.compartment_volume)):
            if 'amici_time_symbol' in dir(self) and (
                    new == self.amici_time_symbol):
                self.compartment_volume[index] = \
                    self.compartment_volume[index].subs(old, 0)
            else:
                self.compartment_volume[index] = \
                    self.compartment_volume[index].subs(old, new)

    def _clean_reserved_symbols(self) -> None:
        """
        Remove all reserved symbols from self.symbols
        """
        reserved_symbols = ['k', 'p', 'y', 'w']
        for sym in reserved_symbols:
            old_symbol = sp.Symbol(sym, real=True)
            new_symbol = sp.Symbol('amici_' + sym, real=True)
            self._replace_in_all_expressions(old_symbol, new_symbol)
            for symbol in self.symbols.keys():
                if 'identifier' in self.symbols[symbol].keys():
                    self.symbols[symbol]['identifier'] = \
                        self.symbols[symbol]['identifier'].subs(old_symbol,
                                                                new_symbol)

    def _replace_special_constants(self) -> None:
        """
        Replace all special constants by their respective SBML
        csymbol definition
        """
        constants = [
            (sp.Symbol('avogadro', real=True), sp.Symbol('6.02214179e23')),
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

    def _replace_reactions_in_rule_formula(self,
                                           rule: sbml.Rule,
                                           formula: sp.Expr):
        """
        SBML allows reaction IDs in rules, which should be interpreted as the
        reaction rate.

        An assignment or rate "...rule cannot be defined for a species that is
        created or destroyed in a reaction, unless that species is defined as
        a boundary condition in the model." Here, valid SBML is assumed, so
        this restriction is not checked.

        :param rule:
            The SBML rule.

        :param formula:
            The `rule` formula that has already been parsed.
            TODO create a function to parse rule formulae, as this logic is
                 repeated a few times.

        :return:
            The formula, but reaction IDs are replaced with respective
            reaction rate symbols.
        """
        reaction_ids = [r.getId() for r in self.sbml.getListOfReactions()]
        reactions_in_rule_formula = {s
                                     for s in formula.free_symbols
                                     if str(s) in reaction_ids}
        if reactions_in_rule_formula:
            if rule.getTypeCode() not in (sbml.SBML_ASSIGNMENT_RULE,
                                          sbml.SBML_RATE_RULE):
                raise SBMLException('Currently, only assignment and rate'
                                    ' rules have reaction replacement'
                                    ' implemented.')

        # Reactions are assigned indices in
        # `sbml_import.py:_process_reactions()`, and these indices are used to
        # generate flux variables in
        # `ode_export.py:import_from_sbml_importer()`.
        # These flux variables are anticipated here, as the symbols that
        # represent the rates of reactions in the model.
        subs = {r_sym: sp.Symbol(f'flux_r{reaction_ids.index(str(r_sym))}',
                                 real=True)
                for r_sym in reactions_in_rule_formula}
        return formula.subs(subs)


def get_rule_vars(rules: List[sbml.Rule],
                  local_symbols: Dict[str, sp.Symbol] = None) -> sp.Matrix:
    """
    Extract free symbols in SBML rule formulas.

    :param rules:
        sbml definitions of rules

    :param local_symbols:
        locals to pass to sympy.sympify

    :return:
        Tuple of free symbolic variables in the formulas all provided rules
    """
    return sp.Matrix(
        [sp.sympify(_parse_logical_operators(
                    sbml.formulaToL3String(rule.getMath())),
                    locals=local_symbols)
         for rule in rules if rule.getFormula() != '']
    ).free_symbols


def replaceLogAB(x: str) -> str:
    """
    Replace log(a, b) in the given string by ln(b)/ln(a).

    Works for nested parentheses and nested 'log's. This can be used to
    circumvent the incompatible argument order between symengine (log(x,
    basis)) and libsbml (log(basis, x)).

    :param x:
        string to replace

    :return:
        string with replaced 'log's
    """

    match = re.search(r'(^|\W)log\(', x)
    if not match:
        return x

    # index of 'l' of 'log'
    log_start = match.start() if match.end() - match.start() == 4 \
        else match.start() + 1
    level = 0  # parenthesis level
    pos_comma = -1  # position of comma in log(a,b)
    for i in range(log_start + 4, len(x)):
        if x[i] == '(':
            level += 1
        elif x[i] == ')':
            level -= 1
            if level == -1:
                break
        elif x[i] == ',' and level == 0:
            pos_comma = i

    if pos_comma < 0:
        # was log(a), not log(a,b), so nothing to replace
        return x

    prefix = x[:log_start]
    suffix = x[i+1:]
    basis = x[log_start+4: pos_comma]
    a = x[pos_comma+1: i]

    replacement = f'{prefix}ln({a})/ln({basis}){suffix}'

    return replaceLogAB(replacement)


def l2s(inputs: List) -> List[str]:
    """
    Transforms a list into list of strings.

    :param inputs:
        objects

    :return: list of str(object)
    """
    return [str(inp) for inp in inputs]


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
        for iError in range(0, sbml_doc.getNumErrors()):
            error = sbml_doc.getError(iError)
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


def _check_unsupported_functions(sym: sp.Basic,
                                 expression_type: str,
                                 full_sym: sp.Basic = None):
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


def _parse_special_functions(sym: sp.Basic, toplevel: bool = True) -> sp.Basic:
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


def assignmentRules2observables(sbml_model,
                                filter_function=lambda *_: True):
    """
    Turn assignment rules into observables.

    :param sbml_model:
        an sbml Model instance

    :param filter_function:
        callback function taking assignment variable as input and returning
        True/False to indicate if the respective rule should be
        turned into an observable

    :return:
        A dictionary(observableId:{
        'name': observableName,
        'formula': formulaString
        })
    """
    warnings.warn("This function will be removed in future releases.",
                  DeprecationWarning)
    observables = {}
    for p in sbml_model.getListOfParameters():
        parameter_id = p.getId()
        if filter_function(p):
            observables[parameter_id] = {
                'name': p.getName(),
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

    :param noise_distribution: An identifier specifying a noise model.
        Possible values are {'normal', 'log-normal', 'log10-normal', 'laplace',
        'log-laplace', 'log10-laplace'}

    :return: A function that takes a strSymbol and then creates a cost
        function string (negative log-likelihood) from it, which can be
        sympified.
    """
    if noise_distribution in ['normal', 'lin-normal']:
        def nllh_y_string(str_symbol):
            return f'0.5*log(2*pi*sigma{str_symbol}**2) ' \
                f'+ 0.5*(({str_symbol} - m{str_symbol}) ' \
                f'/ sigma{str_symbol})**2'
    elif noise_distribution == 'log-normal':
        def nllh_y_string(str_symbol):
            return f'0.5*log(2*pi*sigma{str_symbol}**2*m{str_symbol}**2) ' \
                f'+ 0.5*((log({str_symbol}) - log(m{str_symbol})) ' \
                f'/ sigma{str_symbol})**2'
    elif noise_distribution == 'log10-normal':
        def nllh_y_string(str_symbol):
            return f'0.5*log(2*pi*sigma{str_symbol}**2' \
                f'*m{str_symbol}**2*log(10)**2) ' \
                f'+ 0.5*((log({str_symbol}, 10) - log(m{str_symbol}, 10)) ' \
                f'/ sigma{str_symbol})**2'
    elif noise_distribution in ['laplace', 'lin-laplace']:
        def nllh_y_string(str_symbol):
            return f'log(2*sigma{str_symbol}) ' \
                f'+ Abs({str_symbol} - m{str_symbol}) ' \
                f'/ sigma{str_symbol}'
    elif noise_distribution == 'log-laplace':
        def nllh_y_string(str_symbol):
            return f'log(2*sigma{str_symbol}*m{str_symbol}) ' \
                f'+ Abs(log({str_symbol}) - log(m{str_symbol})) ' \
                f'/ sigma{str_symbol}'
    elif noise_distribution == 'log10-laplace':
        def nllh_y_string(str_symbol):
            return f'log(2*sigma{str_symbol}*m{str_symbol}*log(10)) ' \
                f'+ Abs(log({str_symbol}, 10) - log(m{str_symbol}, 10)) ' \
                f'/ sigma{str_symbol}'
    else:
        raise ValueError(
            f"Cost identifier {noise_distribution} not recognized.")

    return nllh_y_string

class MathMLSbmlPrinter(MathMLContentPrinter):
    """Prints a SymPy expression to a MathML expression parsable by libSBML.
    Differences from `sympy.MathMLContentPrinter`:
    1. underscores in symbol names are not converted to subscripts
    2. symbols with name 'time' are converted to the SBML time symbol
    """
    def _print_Symbol(self, sym):
        ci = self.dom.createElement(self.mathml_tag(sym))
        ci.appendChild(self.dom.createTextNode(sym.name))
        return ci
    # _print_Float can be removed if sympy PR #19958 is merged 
    def _print_Float(self, expr):
        x = self.dom.createElement(self.mathml_tag(expr))
        repr_expr = mlib_to_str(expr._mpf_, repr_dps(expr._prec))
        x.appendChild(self.dom.createTextNode(repr_expr))
        return x
    def doprint(self, expr):
        mathml = super().doprint(expr)
        mathml = '<math xmlns="http://www.w3.org/1998/Math/MathML">' + mathml + '</math>'
        mathml = mathml.replace(f'<ci>time</ci>', '<csymbol encoding="text" definitionURL="http://www.sbml.org/sbml/symbols/time"> time </csymbol>')
        return mathml

def mathml(expr, **settings):
    return MathMLSbmlPrinter(settings).doprint(expr)
