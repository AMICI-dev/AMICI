"""@package amici.sbml_import The SBML import module for Python."""

import sympy as sp
import libsbml as sbml
import re
import math
import itertools as itt
import warnings
import logging
from typing import Dict, Union, List, Callable, Any, Iterable

from .ode_export import ODEExporter, ODEModel
from .logging import get_logger, log_execution_time
from . import has_clibs

from sympy.logic.boolalg import BooleanTrue as spTrue
from sympy.logic.boolalg import BooleanFalse as spFalse

class SBMLException(Exception):
    pass


## default dict for symbols
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


## python log manager
logger = get_logger(__name__, logging.ERROR)


class SbmlImporter:
    """The SbmlImporter class generates AMICI C++ files for a model provided in
    the Systems Biology Markup Language (SBML).

    Attributes:

        show_sbml_warnings: indicates whether libSBML warnings should be
        displayed @type bool

        symbols: dict carrying symbolic definitions @type dict

        sbml_reader: the libSBML sbml reader [!not storing this will result
        in a segfault!]

        sbml_doc: document carrying the sbml definition [!not storing this
        will result in a segfault!]

        sbml: sbml definition [!not storing this will result in a segfault!]

        speciesIndex: maps species names to indices @type dict

        speciesCompartment: compartment for each species @type
        sympy.Matrix

        constantSpecies: ids of species that are marked as constant @type list

        boundaryConditionSpecies: ids of species that are marked as boundary
        condition @type list

        speciesHasOnlySubstanceUnits: flags indicating whether a species has
        only substance units @type list

        speciesConversionFactor: conversion factors for every species @type
        sympy.Matrix

        compartmentSymbols: compartment ids @type sympy.Matrix

        compartmentVolume: numeric/symbolic compartment volumes @type
        sympy.Matrix

        stoichiometricMatrix: stoichiometric matrix of the model @type
        sympy.Matrix

        fluxVector: reaction kinetic laws @type sympy.Matrix

        local_symbols: model symbols for sympy to consider during sympification
        see `locals`argument in `sympy.sympify` @type dict

    """

    def __init__(
            self,
            sbml_source: Union[str, sbml.Model],
            show_sbml_warnings: bool = False,
            from_file: bool = True):
        """Create a new Model instance.

        Arguments:

            sbml_source: Either a path to SBML file where the model is
                specified, or a model string as created by
                sbml.sbmlWriter().writeSBMLToString() or an instance of
                `libsbml.Model`.

            show_sbml_warnings: Indicates whether libSBML warnings should be displayed.

            from_file: Whether `sbml_source` is a file name (True, default), or
                an SBML string

        Raises:

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

        self.show_sbml_warnings = show_sbml_warnings

        # process document
        self.process_document()

        self.sbml = self.sbml_doc.getModel()

        # Long and short names for model components
        self.symbols = dict()
        self.reset_symbols()

        self.local_symbols = {}

    def process_document(self):
        """Validate and simplify document.

        Arguments:

        Returns:

        Raises:

        """
        # Ensure we got a valid SBML model, otherwise further processing
        # might lead to undefined results
        self.sbml_doc.validateSBML()
        checkLibSBMLErrors(self.sbml_doc, self.show_sbml_warnings)

        # apply several model simplifications that make our life substantially
        # easier
        if len(self.sbml_doc.getModel().getListOfFunctionDefinitions()) > 0:
            convertConfig = sbml.SBMLFunctionDefinitionConverter()\
                .getDefaultProperties()
            self.sbml_doc.convert(convertConfig)

        convertConfig = sbml.SBMLLocalParameterConverter().\
            getDefaultProperties()
        self.sbml_doc.convert(convertConfig)

        # If any of the above calls produces an error, this will be added to
        # the SBMLError log in the sbml document. Thus, it is sufficient to
        # check the error log just once after all conversion/validation calls.
        checkLibSBMLErrors(self.sbml_doc, self.show_sbml_warnings)

    def reset_symbols(self):
        """Reset the symbols attribute to default values

        Arguments:

        Returns:

        Raises:

        """
        self.symbols = default_symbols

    def sbml2amici(self,
                   modelName: str,
                   output_dir: str = None,
                   observables: Dict[str, Dict[str, str]] = None,
                   constantParameters: List[str] = None,
                   sigmas: Dict[str, Union[str, float]] = None,
                   noise_distributions: Dict[str, str] = None,
                   verbose: Union[int, bool] = logging.ERROR,
                   assume_pow_positivity: bool = False,
                   compiler: str = None,
                   allow_reinit_fixpar_initcond: bool = True,
                   compile: bool = True
                   ) -> None:
        """Generate AMICI C++ files for the model provided to the constructor.

        The resulting model can be imported as a regular Python module (if
        `compile=True`), or used from Matlab or C++ as described in the
        documentation of the respective AMICI interface.

        Note that this generates model ODEs for changes in concentrations, not
        amounts. The simulation results obtained from the model will be
        concentrations, independently of the SBML `hasOnlySubstanceUnits`
        attribute.

        Arguments:
            modelName: name of the model/model directory

            output_dir: see sbml_import.setPaths()

            observables: dictionary( observableId:{'name':observableName
                (optional), 'formula':formulaString)}) to be added to the model

            constantParameters: list of SBML Ids identifying constant parameters

            sigmas: dictionary(observableId:
                    sigma value or (existing) parameter name)

            noise_distributions: dictionary(observableId: noise type).
                If nothing is passed
                for some observable id, a normal model is assumed as default.

            verbose: verbosity level for logging, True/False default to
                logging.Error/logging.DEBUG

            assume_pow_positivity: if set to True, a special pow function is
                used to avoid problems with state variables that may become
                negative due to numerical errors

            compiler: distutils/setuptools compiler selection to build the
                python extension

            allow_reinit_fixpar_initcond: see ode_export.ODEExporter

            compile: If True, compile the generated Python package,
                if False, just generate code.

        Returns:

        Raises:

        """
        if observables is None:
            observables = {}

        if constantParameters is None:
            constantParameters = []

        if sigmas is None:
            sigmas = {}

        if noise_distributions is None:
            noise_distributions = {}

        logger.setLevel(verbose)

        self.reset_symbols()
        self.processSBML(constantParameters)
        self.processObservables(observables, sigmas, noise_distributions)
        ode_model = ODEModel(simplify=sp.powsimp)
        ode_model.import_from_sbml_importer(self)
        exporter = ODEExporter(
            ode_model,
            outdir=output_dir,
            verbose=verbose,
            assume_pow_positivity=assume_pow_positivity,
            compiler=compiler,
            allow_reinit_fixpar_initcond=allow_reinit_fixpar_initcond
        )
        exporter.setName(modelName)
        exporter.setPaths(output_dir)
        exporter.generateModelCode()

        if compile:
            if not has_clibs:
                warnings.warn('AMICI C++ extensions have not been built. '
                              'Generated model code, but unable to compile.')
            exporter.compileModel()

    def processSBML(self, constantParameters: List[str] = None):
        """Read parameters, species, reactions, and so on from SBML model

        Arguments:
            constantParameters: SBML Ids identifying constant parameters

        Returns:

        Raises:

        """

        if constantParameters is None:
            constantParameters = []

        self.checkSupport()
        self._gather_locals()
        self.processParameters(constantParameters)
        self.processSpecies()
        self.processReactions()
        self.processCompartments()
        self.processRules()
        self.processVolumeConversion()
        self.processTime()
        self.cleanReservedSymbols()
        self.replaceSpecialConstants()

    def checkSupport(self):
        """Check whether all required SBML features are supported.

        Arguments:

        Returns:

        Raises:

        """
        if len(self.sbml.getListOfSpecies()) == 0:
            raise SBMLException('Models without species '
                                'are currently not supported!')

        if hasattr(self.sbml, 'all_elements_from_plugins') \
                and self.sbml.all_elements_from_plugins.getSize() > 0:
            raise SBMLException('SBML extensions are currently not supported!')

        if len(self.sbml.getListOfEvents()) > 0:
            raise SBMLException('Events are currently not supported!')

        if any([not(rule.isAssignment())
                for rule in self.sbml.getListOfRules()]):
            raise SBMLException('Algebraic and rate '
                                'rules are currently not supported!')

        if any([reaction.getFast()
                for reaction in self.sbml.getListOfReactions()]):
            raise SBMLException('Fast reactions are currently not supported!')

        if any([any([not element.getStoichiometryMath() is None
                for element in list(reaction.getListOfReactants())
                + list(reaction.getListOfProducts())])
                for reaction in self.sbml.getListOfReactions()]):
            raise SBMLException('Non-unity stoichiometry is'
                                ' currently not supported!')

    def _gather_locals(self):
        """Populate self.local_symbols with all model entities.

        This is later used during sympifications to avoid sympy builtins
        shadowing model entities.

        Arguments:

        Returns:

        Raises:

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

    @log_execution_time('processing SBML species', logger)
    def processSpecies(self):
        """Get species information from SBML model.

        Arguments:

        Returns:

        Raises:

        """
        species = self.sbml.getListOfSpecies()

        self.speciesIndex = {
            species_element.getId(): species_index
            for species_index, species_element in enumerate(species)
        }

        self.symbols['species']['identifier'] = sp.Matrix(
            [sp.Symbol(spec.getId(), real=True) for spec in species]
        )
        self.symbols['species']['name'] = [spec.getName() for spec in species]

        self.speciesCompartment = sp.Matrix(
            [sp.Symbol(spec.getCompartment(), real=True) for spec in species]
        )

        self.constantSpecies = [species_element.getId()
                                for species_element in species
                                if species_element.getConstant()]

        self.boundaryConditionSpecies = [
            species_element.getId()
            for species_element in species
            if species_element.getBoundaryCondition()
        ]
        self.speciesHasOnlySubstanceUnits = [
            specie.getHasOnlySubstanceUnits() for specie in species
        ]

        concentrations = [spec.getInitialConcentration() for spec in species]
        amounts = [spec.getInitialAmount() for spec in species]

        def get_species_initial(index, conc):
            # We always simulate concentrations!
            if self.speciesHasOnlySubstanceUnits[index]:
                if species[index].isSetInitialAmount() \
                        and not math.isnan(amounts[index]):
                    return sp.sympify(amounts[index]) \
                           / self.speciesCompartment[index]
                if species[index].isSetInitialConcentration():
                    return sp.sympify(conc)
            else:
                if species[index].isSetInitialConcentration():
                    return sp.sympify(conc)

                if species[index].isSetInitialAmount() \
                        and not math.isnan(amounts[index]):
                    return sp.sympify(amounts[index]) \
                           / self.speciesCompartment[index]

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
                symMath = sp.sympify(_parse_logical_operators(
                    sbml.formulaToL3String(initial_assignment.getMath())),
                    locals=self.local_symbols
                )
                if symMath is not None:
                    symMath = _parse_special_functions(symMath)
                    _check_unsupported_functions(symMath, 'InitialAssignment')
                    species_initial[index] = symMath

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

        self.speciesConversionFactor = sp.Matrix([
             sp.sympify(specie.getConversionFactor())
             if specie.isSetConversionFactor()
             else conversion_factor
             for specie in species
        ])

    @log_execution_time('processing SBML parameters', logger)
    def processParameters(self, constantParameters: List[str] = None):
        """Get parameter information from SBML model.

        Arguments:
            constantParameters: SBML Ids identifying constant parameters

        Returns:

        Raises:

        """

        if constantParameters is None:
            constantParameters = []

        # Ensure specified constant parameters exist in the model
        for parameter in constantParameters:
            if not self.sbml.getParameter(parameter):
                raise KeyError('Cannot make %s a constant parameter: '
                               'Parameter does not exist.' % parameter)

        parameter_ids = [par.getId() for par
                         in self.sbml.getListOfParameters()]
        for initial_assignment in self.sbml.getListOfInitialAssignments():
            if initial_assignment.getId() in parameter_ids:
                raise SBMLException('Initial assignments for parameters are'
                                    ' currently not supported')

        fixedParameters = [parameter for parameter
                           in self.sbml.getListOfParameters()
                           if parameter.getId() in constantParameters
                           ]

        rulevars = [rule.getVariable() for rule in self.sbml.getListOfRules()]

        parameters = [ parameter for parameter
                       in self.sbml.getListOfParameters()
                       if parameter.getId() not in constantParameters
                       and parameter.getId() not in rulevars]

        loop_settings = {
            'parameter': {
                'var': parameters,
                'name': 'parameter',

            },
            'fixed_parameter': {
                'var': fixedParameters,
                'name': 'fixedParameter'
            }

        }

        for partype in loop_settings:
            settings = loop_settings[partype]

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
                f'{settings["name"]}Index',
                {
                    parameter_element.getId(): parameter_index
                    for parameter_index, parameter_element
                    in enumerate(settings['var'])
                }
            )

    @log_execution_time('processing SBML compartments', logger)
    def processCompartments(self):
        """Get compartment information, stoichiometric matrix and fluxes from
        SBML model.

        Arguments:

        Returns:

        Raises:

        """
        compartments = self.sbml.getListOfCompartments()
        self.compartmentSymbols = sp.Matrix(
            [sp.Symbol(comp.getId(), real=True) for comp in compartments]
        )
        self.compartmentVolume = sp.Matrix(
            [sp.sympify(comp.getVolume()) if comp.isSetVolume()
            else sp.sympify(1.0) for comp in compartments]
        )

        compartment_ids = [comp.getId() for comp in compartments]
        for initial_assignment in self.sbml.getListOfInitialAssignments():
            if initial_assignment.getId() in compartment_ids:
                index = compartment_ids.index(
                        initial_assignment.getId()
                    )
                self.compartmentVolume[index] = sp.sympify(
                    sbml.formulaToL3String(initial_assignment.getMath()),
                    locals=self.local_symbols
                )

    @log_execution_time('processing SBML reactions', logger)
    def processReactions(self):
        """Get reactions from SBML model.

        Arguments:

        Returns:

        Raises:

        """
        reactions = self.sbml.getListOfReactions()
        nr = len(reactions)
        nx = len(self.symbols['species']['name'])
        # stoichiometric matrix
        self.stoichiometricMatrix = sp.SparseMatrix(sp.zeros(nx, nr))
        self.fluxVector = sp.zeros(nr, 1)

        assignment_ids = [ass.getId()
                          for ass in self.sbml.getListOfInitialAssignments()]
        rulevars = [rule.getVariable()
                                for rule in self.sbml.getListOfRules()
                                if rule.getFormula() != '']

        reaction_ids = [
            reaction.getId() for reaction in reactions
            if reaction.isSetId()
        ]

        def getElementFromAssignment(element_id):
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

        def getElementStoichiometry(element):
            if element.isSetId():
                if element.getId() in assignment_ids:
                    symMath = getElementFromAssignment(element.getId())
                    if symMath is None:
                        symMath = sp.sympify(element.getStoichiometry())
                elif element.getId() in rulevars:
                    return sp.Symbol(element.getId(), real=True)
                else:
                    # dont put the symbol if it wont get replaced by a
                    # rule
                    symMath = sp.sympify(element.getStoichiometry())
            elif element.isSetStoichiometry():
                symMath = sp.sympify(element.getStoichiometry())
            else:
                return sp.sympify(1.0)
            symMath = _parse_special_functions(symMath)
            _check_unsupported_functions(symMath, 'Stoichiometry')
            return symMath

        def isConstant(specie):
            return specie in self.constantSpecies or \
                specie in self.boundaryConditionSpecies

        for reactionIndex, reaction in enumerate(reactions):
            for elementList, sign in [(reaction.getListOfReactants(), -1.0),
                                       (reaction.getListOfProducts(), 1.0)]:
                elements = {}
                for index, element in enumerate(elementList):
                    # we need the index here as we might have multiple elements
                    # for the same species
                    elements[index] = {'species': element.getSpecies()}
                    elements[index]['stoichiometry'] = getElementStoichiometry(
                        element
                    )

                for index in elements.keys():
                    if not isConstant(elements[index]['species']):
                        specieIndex = self.speciesIndex[
                            elements[index]['species']
                        ]
                        self.stoichiometricMatrix[specieIndex, reactionIndex] \
                            += sign \
                            * elements[index]['stoichiometry'] \
                            * self.speciesConversionFactor[specieIndex] \
                            / self.speciesCompartment[specieIndex]

            # usage of formulaToL3String ensures that we get "time" as time
            # symbol
            math = sbml.formulaToL3String(reaction.getKineticLaw().getMath())
            try:
                symMath = sp.sympify(_parse_logical_operators(math),
                                     locals=self.local_symbols)
            except SBMLException as Ex:
                raise Ex
            except:
                raise SBMLException(f'Kinetic law "{math}" contains an '
                                    'unsupported expression!')
            symMath = _parse_special_functions(symMath)
            _check_unsupported_functions(symMath, 'KineticLaw')
            for r in reactions:
                elements = list(r.getListOfReactants()) \
                           + list(r.getListOfProducts())
                for element in elements:
                    if element.isSetId() & element.isSetStoichiometry():
                        symMath = symMath.subs(
                            sp.sympify(element.getId(),
                                       locals=self.local_symbols),
                            sp.sympify(element.getStoichiometry())
                        )

            self.fluxVector[reactionIndex] = symMath
            if any([
                str(symbol) in reaction_ids
                for symbol in self.fluxVector[reactionIndex].free_symbols
            ]):
                raise SBMLException(
                    'Kinetic laws involving reaction ids are currently'
                    ' not supported!'
                )

    @log_execution_time('processing SBML rules', logger)
    def processRules(self):
        """Process Rules defined in the SBML model.

        Arguments:

        Returns:

        Raises:

        """
        rules = self.sbml.getListOfRules()

        rulevars = getRuleVars(rules, local_symbols=self.local_symbols)
        fluxvars = self.fluxVector.free_symbols
        specvars = self.symbols['species']['identifier'].free_symbols
        volumevars = self.compartmentVolume.free_symbols
        compartmentvars = self.compartmentSymbols.free_symbols
        parametervars = sp.Matrix([
            sp.Symbol(par.getId(), real=True)
            for par in self.sbml.getListOfParameters()
        ])
        stoichvars = self.stoichiometricMatrix.free_symbols

        assignments = {}

        for rule in rules:
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

            if variable in stoichvars:
                self.stoichiometricMatrix = \
                    self.stoichiometricMatrix.subs(variable, formula)

            if variable in specvars:
                raise SBMLException('Species assignment rules are currently'
                                    ' not supported!')

            if variable in compartmentvars:
                raise SBMLException('Compartment assignment rules are'
                                    ' currently not supported!')

            if variable in parametervars:
                try:
                    idx = self.parameterIndex[str(variable)]
                    self.symbols['parameter']['value'][idx] \
                        = float(formula)
                except:
                    self.sbml.removeParameter(str(variable))
                    assignments[str(variable)] = formula

            if variable in fluxvars:
                self.fluxVector = self.fluxVector.subs(variable, formula)

            if variable in volumevars:
                self.compartmentVolume = \
                    self.compartmentVolume.subs(variable, formula)

            if variable in rulevars:
                for nested_rule in rules:
                    nested_formula = sp.sympify(
                        sbml.formulaToL3String(nested_rule.getMath()),
                        locals=self.local_symbols)
                    nested_formula = \
                        nested_formula.subs(variable, formula)
                    nested_rule.setFormula(str(nested_formula))

                for variable in assignments:
                    assignments[variable].subs(variable, formula)

        # do this at the very end to ensure we have flattened all recursive
        # rules
        for variable in assignments.keys():
            self.replaceInAllExpressions(
                sp.Symbol(variable, real=True),
                assignments[variable]
            )
        for comp, vol in zip(self.compartmentSymbols, self.compartmentVolume):
            self.replaceInAllExpressions(
               comp, vol
            )

    def processVolumeConversion(self):
        """Convert equations from amount to volume.

        Arguments:

        Returns:

        Raises:

        """
        compartments = self.speciesCompartment
        for comp, vol in zip(self.compartmentSymbols, self.compartmentVolume):
            compartments = compartments.subs(comp, vol)
        for index, bool in enumerate(self.speciesHasOnlySubstanceUnits):
            if bool:
                self.fluxVector = \
                    self.fluxVector.subs(
                        self.symbols['species']['identifier'][index],
                        self.symbols['species']['identifier'][index]
                        * compartments[index]
                    )

    def processTime(self):
        """Convert time_symbol into cpp variable.

        Arguments:

        Returns:

        Raises:

        """
        sbmlTimeSymbol = sp.Symbol('time', real=True)
        amiciTimeSymbol = sp.Symbol('t', real=True)

        self.replaceInAllExpressions(sbmlTimeSymbol, amiciTimeSymbol)

    @log_execution_time('processing SBML observables', logger)
    def processObservables(self, observables: Dict[str, Dict[str, str]],
                           sigmas: Dict[str, Union[str, float]],
                           noise_distributions: Dict[str, str]):
        """Perform symbolic computations required for objective function
        evaluation.

        Arguments:
            observables: dictionary(observableId: {'name':observableName
                (optional), 'formula':formulaString)})
                to be added to the model

            sigmas: dictionary(observableId: sigma value or (existing)
                parameter name)

            noise_distributions: dictionary(observableId: noise type)
                See `sbml2amici`.

        Returns:

        Raises:

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
            unknown_ids = set(noise_distributions.keys()) - set(observables.keys())
            if unknown_ids:
                raise ValueError(
                    f"Noise distribution provided for unknown observableIds: "
                    f"{unknown_ids}.")

        speciesSyms = self.symbols['species']['identifier']

        # add user-provided observables or make all species observable
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

            def replace_assignments(formula):
                """Replace assignment rules in observables"""
                formula = sp.sympify(formula, locals=self.local_symbols)
                for s in formula.free_symbols:
                    r = self.sbml.getAssignmentRuleByVariable(str(s))
                    if r is not None:
                        formula = formula.replace(s, sp.sympify(
                            sbml.formulaToL3String(r.getMath()),
                             locals=self.local_symbols))
                return formula

            observableValues = sp.Matrix([
                replace_assignments(observables[observable]['formula'])
                for observable in observables
            ])
            observableNames = [
                observables[observable]['name'] if 'name' in observables[
                    observable].keys()
                else f'y{index}'
                for index, observable in enumerate(observables)
            ]
            observableSyms = sp.Matrix([
                sp.symbols(obs, real=True) for obs in observables.keys()
            ])
            observable_ids = observables.keys()
        else:
            observableValues = speciesSyms
            observable_ids = [
                f'x{index}' for index in range(len(speciesSyms))
            ]
            observableNames = observable_ids[:]
            observableSyms = sp.Matrix(
                [sp.symbols(f'y{index}', real=True)
                 for index in range(len(speciesSyms))]
            )

        sigmaYSyms = sp.Matrix(
            [sp.symbols(f'sigma{symbol}', real=True)
             for symbol in observableSyms]
        )
        sigmaYValues = sp.Matrix(
            [1.0] * len(observableSyms)
        )

        # set user-provided sigmas
        for iy, obsName in enumerate(observables):
            if obsName in sigmas:
                sigmaYValues[iy] = sp.Symbol(sigmas[obsName], real=True)

        measurementYSyms = sp.Matrix(
            [sp.symbols(f'm{symbol}', real=True) for symbol in observableSyms]
        )
        measurementYValues = sp.Matrix(
            [0.0] * len(observableSyms)
        )

        # set cost functions
        llhYStrings = []
        for y_name in observable_ids:
            llhYStrings.append(noise_distribution_to_cost_function(
                noise_distributions.get(y_name, 'normal')))

        llhYValues = []
        for llhYString, o_sym, m_sym, s_sym \
                in zip(llhYStrings, observableSyms,
                       measurementYSyms, sigmaYSyms):
            f = sp.sympify(llhYString(o_sym), locals={str(o_sym): o_sym,
                                                      str(m_sym): m_sym,
                                                      str(s_sym): s_sym})
            llhYValues.append(f)
        llhYValues = sp.Matrix(llhYValues)

        llhYSyms = sp.Matrix(
            [sp.Symbol(f'J{symbol}', real=True) for symbol in observableSyms]
        )

        # set symbols
        self.symbols['observable']['identifier'] = observableSyms
        self.symbols['observable']['name'] = l2s(observableNames)
        self.symbols['observable']['value'] = observableValues
        self.symbols['sigmay']['identifier'] = sigmaYSyms
        self.symbols['sigmay']['name'] = l2s(sigmaYSyms)
        self.symbols['sigmay']['value'] = sigmaYValues
        self.symbols['my']['identifier'] = measurementYSyms
        self.symbols['my']['name'] = l2s(measurementYSyms)
        self.symbols['my']['value'] = measurementYValues
        self.symbols['llhy']['value'] = llhYValues
        self.symbols['llhy']['name'] = l2s(llhYSyms)
        self.symbols['llhy']['identifier'] = llhYSyms

    def replaceInAllExpressions(self, old: sp.Symbol, new: sp.Symbol):
        """Replace 'old' by 'new' in all symbolic expressions.

        Arguments:
            old: symbolic variables to be replaced

            new: replacement symbolic variables


        Returns:

        Raises:

        """
        fields = [
            'stoichiometricMatrix', 'fluxVector',
        ]
        for field in fields:
            if field in dir(self):
                self.__setattr__(field, self.__getattribute__(field).subs(
                    old, new
                ))
        symbols = [
            'species', 'observables',
        ]
        for symbol in symbols:
            if symbol in self.symbols:
                self.symbols[symbol]['value'] = \
                    self.symbols[symbol]['value'].subs(old, new)

    def cleanReservedSymbols(self):
        """Remove all reserved symbols from self.symbols

        Arguments:

        Returns:

        Raises:

        """
        reservedSymbols = ['k','p','y','w']
        for str in reservedSymbols:
            old_symbol = sp.Symbol(str, real=True)
            new_symbol = sp.Symbol('amici_' + str, real=True)
            self.replaceInAllExpressions(old_symbol, new_symbol)
            for symbol in self.symbols.keys():
                if 'identifier' in self.symbols[symbol].keys():
                    self.symbols[symbol]['identifier'] = \
                        self.symbols[symbol]['identifier'].subs(old_symbol,new_symbol)

    def replaceSpecialConstants(self):
        """Replace all special constants by their respective SBML
        csymbol definition

        Arguments:

        Returns:

        Raises:

        """
        constants = [
            (sp.Symbol('avogadro', real=True), sp.Symbol('6.02214179e23')),
        ]
        for constant, value in constants:
            # do not replace if any symbol is shadowing default definition
            if not any([constant in self.symbols[symbol]['identifier']
                        for symbol in self.symbols.keys()
                        if 'identifier' in self.symbols[symbol].keys()]):
                self.replaceInAllExpressions(constant, value)
            else:
                # yes sbml supports this but we wont. Are you really expecting
                # to be saved if you are trying to shoot yourself in the foot?
                raise SBMLException(
                    f'Encountered currently unsupported element id {constant}!'
                )


def getRuleVars(rules, local_symbols=None):
    """Extract free symbols in SBML rule formulas.

    Arguments:
        rules: sbml definitions of rules @type list
        local_symbols: locals to pass to sympy.sympify

    Returns:
    Tuple of free symbolic variables in the formulas all provided rules

    Raises:

    """
    return sp.Matrix(
        [sp.sympify(sbml.formulaToL3String(rule.getMath()),
                    locals=local_symbols)
         for rule in rules if rule.getFormula() != '']
    ).free_symbols


def replaceLogAB(x):
    """Replace log(a, b) in the given string by ln(b)/ln(a)

    Works for nested parentheses and nested 'log's. This can be used to
    circumvent the incompatible argument order between symengine (log(x,
    basis)) and libsbml (log(basis, x)).

    Arguments:
        x: string to replace @type str

    Returns:
    string with replaced 'log's

    Raises:

    """

    match = re.search(r'(^|\W)log\(', x)
    if not match:
        return x

    # index of 'l' of 'log'
    logStart = match.start() \
            if match.end() - match.start() == 4 \
            else match.start() + 1
    level = 0 # parenthesis level
    posComma = -1 # position of comma in log(a,b)
    for i in range(logStart + 4, len(x)):
        if x[i] == '(':
            level += 1
        elif x[i] == ')':
            level -= 1
            if level == -1: break
        elif x[i] == ',' and level == 0:
            posComma = i

    if posComma < 0:
        # was log(a), not log(a,b), so nothing to replace
        return x

    prefix = x[:logStart]
    suffix = x[i+1:]
    basis = x[logStart+4: posComma]
    a = x[posComma+1: i]

    replacement = f'{prefix}ln({a})/ln({basis}){suffix}'

    return replaceLogAB(replacement)


def l2s(inputs):
    """Transforms a list into list of strings.

    Arguments:
        inputs: objects @type list

    Returns:
    list of str(object)

    Raises:

    """
    return [str(inp) for inp in inputs]


def checkLibSBMLErrors(sbml_doc, show_warnings=False):
    """Checks the error log in the current self.sbml_doc.

    Arguments:
        sbml_doc: SBML document @type libsbml.SBMLDocument
        show_warnings: display SBML warnings @type bool

    Returns:

    Raises:
        raises SBMLException if errors with severity ERROR or FATAL have
        occurred
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


def _check_unsupported_functions(sym, expression_type, full_sym=None):
    """Recursively checks the symbolic expression for unsupported symbolic
    functions

        Arguments:
            sym: symbolic expressions @type sympy.Basic
            expression_type: type of expression

        Returns:

        Raises:
            raises SBMLException if an unsupported function is encountered
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


def _parse_special_functions(sym, toplevel=True):
    """Recursively checks the symbolic expression for functions which have be
    to parsed in a special way, such as piecewise functions

        Arguments:
            sym: symbolic expressions @type sympy.Basic
            toplevel: as this is called recursively,
                are we in the top level expression?
        Returns:

        Raises:
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


def _parse_logical_operators(math_str: str) -> str:
    """Parses a math string in order to replace logical operators by a form
    parsable for sympy

        Arguments:
            math_str: str with mathematical expression

        Returns:
            math_str: parsed math_str

        Raises:
    """
    if math_str is None:
        return None

    if ' xor(' in math_str or ' Xor(' in math_str:
        raise SBMLException('Xor is currently not supported as logical '
                            'operation.')

    return (math_str.replace('&&', '&')).replace('||', '|')


def grouper(iterable: Iterable, n: int, fillvalue: Any = None):
    """Collect data into fixed-length chunks or blocks

    E.g. grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"

    Arguments:
        iterable: any iterable
        n: chunk length
        fillvalue: padding for last chunk if length < n

    Returns:
        itertools.zip_longest of requested chunks
    """
    args = [iter(iterable)] * n
    return itt.zip_longest(*args, fillvalue=fillvalue)


def assignmentRules2observables(sbml_model,
                                filter_function=lambda *_: True):
    """Turn assignment rules into observables.
    Arguments:
    sbml_model: an sbml Model instance
    filter_function: callback function taking assignment variable as input
    and returning True/False to indicate if the respective rule should be
    turned into an observable
    Returns:
    A dictionary(observableId:{
        'name': observableName,
        'formula': formulaString
    })
    Raises:
    """
    warnings.warn("This function will be removed in future releases. "
                  "This functionality is now included in "
                  "https://github.com/ICB-DCM/PEtab .", DeprecationWarning)
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
        noise_distribution: str) -> Callable[[str], str]:
    """Parse noise distribution string to a cost function definition amici can
    work with.

    Arguments:

    noise_distribution: A code specifying a noise model. Can be any of
    [normal, log-normal, log10-normal, laplace, log-laplace, log10-laplace].

    Returns:

    A function that takes a strSymbol and then creates a cost function string
    (negative log-likelihood) from it, which can be sympified.

    Raises:
        ValueError: in case of invalid ``noise_distribution``
    """
    if noise_distribution in ['normal', 'lin-normal']:
        nllh_y_string = lambda str_symbol: \
            f'0.5*log(2*pi*sigma{str_symbol}**2) ' \
            f'+ 0.5*(({str_symbol} - m{str_symbol}) ' \
            f'/ sigma{str_symbol})**2'
    elif noise_distribution == 'log-normal':
        nllh_y_string = lambda str_symbol: \
            f'0.5*log(2*pi*sigma{str_symbol}**2*m{str_symbol}**2) ' \
            f'+ 0.5*((log({str_symbol}) - log(m{str_symbol})) ' \
            f'/ sigma{str_symbol})**2'
    elif noise_distribution == 'log10-normal':
        nllh_y_string = lambda str_symbol: \
            f'0.5*log(2*pi*sigma{str_symbol}**2*m{str_symbol}**2) ' \
            f'+ 0.5*((log({str_symbol}, 10) - log(m{str_symbol}, 10)) ' \
            f'/ sigma{str_symbol})**2'
    elif noise_distribution in ['laplace', 'lin-laplace']:
        nllh_y_string = lambda str_symbol: \
            f'log(2*sigma{str_symbol}) ' \
            f'+ Abs({str_symbol} - m{str_symbol}) ' \
            f'/ sigma{str_symbol}'
    elif noise_distribution == 'log-laplace':
        nllh_y_string = lambda str_symbol: \
            f'log(2*sigma{str_symbol}*m{str_symbol}) ' \
            f'+ Abs(log({str_symbol}) - log(m{str_symbol})) ' \
            f'/ sigma{str_symbol}'
    elif noise_distribution == 'log10-laplace':
        nllh_y_string = lambda str_symbol: \
            f'log(2*sigma{str_symbol}*m{str_symbol}) ' \
            f'+ Abs(log({str_symbol}, 10) - log(m{str_symbol}, 10)) ' \
            f'/ sigma{str_symbol}'
    else:
        raise ValueError(
            f"Cost type {noise_distribution} not recognized.")

    return nllh_y_string
