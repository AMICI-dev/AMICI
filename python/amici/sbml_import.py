#!/usr/bin/env python3
""" @package amici.sbml_import The python sbml import module for python
"""

import sympy as sp
import libsbml as sbml
import re
import math
import warnings
from sympy.logic.boolalg import BooleanTrue as spTrue
from sympy.logic.boolalg import BooleanFalse as spFalse

from .ode_export import ODEExporter, ODEModel
from . import has_clibs

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


class SbmlImporter:
    """The SbmlImporter class generates AMICI C++ files for a model provided in
    the Systems Biology Markup Language (SBML).

    Attributes:

        show_sbml_warnings: indicates whether libSBML warnings should be
        displayed @type bool

        symbols: dict carrying symbolic definitions @type dict

        SBMLreader: the libSBML sbml reader [!not storing this will result
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

    """

    def __init__(self, SBMLFile, show_sbml_warnings=False):
        """Create a new Model instance.

        Arguments:

        SBMLFile: Path to SBML file where the model is specified @type string

        show_sbml_warnings: indicates whether libSBML warnings should be
        displayed @type bool

        Returns:
        SbmlImporter instance with attached SBML document

        Raises:

        """

        self.show_sbml_warnings = show_sbml_warnings

        self.loadSBMLFile(SBMLFile)

        """Long and short names for model components"""
        self.symbols = dict()
        self.reset_symbols()

    def reset_symbols(self):
        """Reset the symbols attribute to default values

        Arguments:

        Returns:

        Raises:

        """
        self.symbols = default_symbols

    def loadSBMLFile(self, SBMLFile):
        """Parse the provided SBML file.

        Arguments:
            SBMLFile: path to SBML file @type str

        Returns:

        Raises:

        """

        self.SBMLreader = sbml.SBMLReader()
        self.sbml_doc = self.SBMLreader.readSBML(SBMLFile)

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

        self.sbml = self.sbml_doc.getModel()

    def sbml2amici(self,
                   modelName,
                   output_dir = None,
                   observables = None,
                   constantParameters = None,
                   sigmas = None,
                   verbose = False,
                   assume_pow_positivity = False,
                   compiler = None,
                   allow_reinit_fixpar_initcond = True,
                   compile = True
                   ):
        """Generate AMICI C++ files for the model provided to the constructor.

        Arguments:
            modelName: name of the model/model directory @type str

            output_dir: see sbml_import.setPaths() @type str

            observables: dictionary( observableId:{'name':observableName
                (optional), 'formula':formulaString)}) to be added to the model
                @type dict

            sigmas: dictionary(observableId: sigma value or (existing) parameter name)
                @type dict

            constantParameters: list of SBML Ids identifying constant parameters
                @type list

            verbose: more verbose output if True @type bool

            assume_pow_positivity: if set to true, a special pow function is
                used to avoid problems with state variables that may become
                negative due to numerical errors @type bool

            compiler: distutils/setuptools compiler selection to build the
                python extension @type str

            allow_reinit_fixpar_initcond: see ode_export.ODEExporter @type bool

            compile: If True, compile the generated Python package, if False, just
                generate code. @type bool

        Returns:

        Raises:

        """
        if observables is None:
            observables = {}

        if constantParameters is None:
            constantParameters = []

        if sigmas is None:
            sigmas = {}

        self.reset_symbols()
        self.processSBML(constantParameters)
        self.processObservables(observables, sigmas)
        ode_model = ODEModel()
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

    def processSBML(self, constantParameters=None):
        """Read parameters, species, reactions, and so on from SBML model

        Arguments:
            constantParameters: SBML Ids identifying constant parameters
            @type list

        Returns:

        Raises:

        """

        if constantParameters is None:
            constantParameters = []

        self.checkSupport()
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
                and len(self.sbml.all_elements_from_plugins) > 0:
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
            [sp.Symbol(spec.getId()) for spec in species]
        )
        self.symbols['species']['name'] = [spec.getName() for spec in species]

        self.speciesCompartment = sp.Matrix(
            [sp.Symbol(spec.getCompartment()) for spec in species]
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

        def getSpeciesInitial(index, conc):
            if not math.isnan(conc):
                return sp.sympify(conc)
            if not math.isnan(amounts[index]):
                return \
                    sp.sympify(amounts[index]) / self.speciesCompartment[index]
            return self.symbols['species']['identifier'][index]

        speciesInitial = sp.Matrix(
            [getSpeciesInitial(index, conc)
            for index, conc in enumerate(concentrations)]
        )

        species_ids = [spec.getId() for spec in self.sbml.getListOfSpecies()]
        for initial_assignment in self.sbml.getListOfInitialAssignments():
            if initial_assignment.getId() in species_ids:
                index = species_ids.index(
                        initial_assignment.getId()
                    )
                symMath = sp.sympify(
                    sbml.formulaToL3String(initial_assignment.getMath())
                )
                if symMath is not None:
                    _check_unsupported_functions(symMath, 'InitialAssignment')
                    speciesInitial[index] = symMath

        for ix, (symbol, init) in enumerate(zip(
                    self.symbols['species']['identifier'], speciesInitial
        )):
            if symbol == init:
                speciesInitial[ix] = sp.sympify(0.0)

        # flatten initSpecies
        while any([species in speciesInitial.free_symbols
                   for species in self.symbols['species']['identifier']]):
            speciesInitial = speciesInitial.subs([
                (symbol, init)
                for symbol, init in zip(
                    self.symbols['species']['identifier'], speciesInitial
                )
            ])

        self.symbols['species']['value'] = speciesInitial

        if self.sbml.isSetConversionFactor():
            conversion_factor = sp.Symbol(self.sbml.getConversionFactor())
        else:
            conversion_factor = 1.0

        self.speciesConversionFactor = sp.Matrix([
             sp.sympify(specie.getConversionFactor())
             if specie.isSetConversionFactor()
             else conversion_factor
             for specie in species
        ])




    def processParameters(self, constantParameters=None):
        """Get parameter information from SBML model.

        Arguments:
            constantParameters: SBML Ids identifying constant parameters
            @type list

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
                [sp.Symbol(par.getId()) for par in settings['var']]
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

    def processCompartments(self):
        """Get compartment information, stoichiometric matrix and fluxes from
        SBML model.

        Arguments:

        Returns:

        Raises:

        """
        compartments = self.sbml.getListOfCompartments()
        self.compartmentSymbols = sp.Matrix(
            [sp.Symbol(comp.getId()) for comp in compartments]
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
                    sbml.formulaToL3String(initial_assignment.getMath())
                )

        for comp, vol in zip(self.compartmentSymbols, self.compartmentVolume):
            self.replaceInAllExpressions(
               comp, vol
            )


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
            sym = sp.sympify(sbml.formulaToL3String(assignment.getMath()))
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
                    return sp.Symbol(element.getId())
                else:
                    # dont put the symbol if it wont get replaced by a
                    # rule
                    symMath = sp.sympify(element.getStoichiometry())
            elif element.isSetStoichiometry():
                symMath = sp.sympify(element.getStoichiometry())
            else:
                return sp.sympify(1.0)
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
                symMath = sp.sympify(math)
            except:
                raise SBMLException(f'Kinetic law "{math}" contains an '
                                    'unsupported expression!')
            _check_unsupported_functions(symMath, 'KineticLaw')
            for r in reactions:
                elements = list(r.getListOfReactants()) \
                           + list(r.getListOfProducts())
                for element in elements:
                    if element.isSetId() & element.isSetStoichiometry():
                        symMath = symMath.subs(
                            sp.sympify(element.getId()),
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

    def processRules(self):
        """Process Rules defined in the SBML model.

        Arguments:

        Returns:

        Raises:

        """
        rules = self.sbml.getListOfRules()

        rulevars = getRuleVars(rules)
        fluxvars = self.fluxVector.free_symbols
        specvars = self.symbols['species']['identifier'].free_symbols
        volumevars = self.compartmentVolume.free_symbols
        compartmentvars = self.compartmentSymbols.free_symbols
        parametervars = sp.Matrix([
            sp.Symbol(par.getId()) for par in self.sbml.getListOfParameters()
        ])
        stoichvars = self.stoichiometricMatrix.free_symbols

        assignments = {}

        for rule in rules:
            if rule.getFormula() == '':
                continue
            variable = sp.sympify(rule.getVariable())
            # avoid incorrect parsing of pow(x, -1) in symengine
            formula = sp.sympify(sbml.formulaToL3String(rule.getMath()))
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
                        sbml.formulaToL3String(nested_rule.getMath()))
                    nested_formula = \
                        nested_formula.subs(variable, formula)
                    nested_rule.setFormula(str(nested_formula))

                for variable in assignments:
                    assignments[variable].subs(variable, formula)

        # do this at the very end to ensure we have flattened all recursive
        # rules
        for variable in assignments.keys():
            self.replaceInAllExpressions(
                sp.sympify(variable),
                assignments[variable]
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
        sbmlTimeSymbol = sp.Symbol('time')
        amiciTimeSymbol = sp.Symbol('t')

        self.replaceInAllExpressions(sbmlTimeSymbol, amiciTimeSymbol)

    def processObservables(self, observables, sigmas):
        """Perform symbolic computations required for objective function
        evaluation.

        Arguments:
            observables: dictionary( observableId:{'name':observableName
            (optional), 'formula':formulaString)}) to be added to the model
            @type dict

            sigmas: dictionary(observableId: sigma value or (existing)
            parameter name) @type dict

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
            unknown_observables = set(sigmas.keys()) - set(observables.keys())
            if unknown_observables:
                raise ValueError('Sigma provided for an unknown observableId: '
                                 + str(unknown_observables))

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
                    print(
                        f'Replaced "{observables[observable]["formula"]}" by '
                        f'"{repl}", assuming first argument to log() was the '
                        f'basis.'
                    )
                    observables[observable]['formula'] = repl

            observableValues = sp.Matrix([
                sp.sympify(observables[observable]['formula'])
                for observable in observables
            ])
            observableNames = [
                observables[observable]['name'] if 'name' in observables[
                    observable].keys()
                else f'y{index}'
                for index, observable in enumerate(observables)
            ]
            observableSyms = sp.Matrix([
                sp.Symbol(obs) for obs in observables.keys()
            ])
        else:
            observableValues = speciesSyms
            observableNames = [
                f'x{index}' for index in range(len(speciesSyms))
            ]
            observableSyms = sp.Matrix(
                [sp.Symbol(f'y{index}') for index in range(len(speciesSyms))]
            )

        sigmaYSyms = sp.Matrix(
            [sp.Symbol(f'sigma{symbol}') for symbol in observableSyms]
        )
        sigmaYValues = sp.Matrix(
            [1.0] * len(observableSyms)
        )

        # set user-provided sigmas
        for iy, obsName in enumerate(observables):
            if obsName in sigmas:
                sigmaYValues[iy] = sigmas[obsName]

        measurementYSyms = sp.Matrix(
            [sp.Symbol(f'm{symbol}') for symbol in observableSyms]
        )
        measurementYValues = sp.Matrix(
            [0.0] * len(observableSyms)
        )

        llhYString = lambda \
            strSymbol: f'0.5*log(2*pi*sigma{strSymbol}**2)' \
                       f'+ 0.5*(({strSymbol} - m{strSymbol})' \
                       f'/ sigma{strSymbol})**2'
        llhYValues = sp.Matrix(
            [sp.sympify(llhYString(symbol))
             for symbol in observableSyms]
        )
        llhYSyms = sp.Matrix(
            [sp.Symbol(f'J{symbol}') for symbol in observableSyms]
        )

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

    def replaceInAllExpressions(self, old, new):
        """Replace 'old' by 'new' in all symbolic expressions.

        Arguments:
            old: symbolic variables to be replaced @type symengine.Symbol

            new: replacement symbolic variables @type symengine.Symbol


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
            old_symbol = sp.Symbol(str)
            new_symbol = sp.Symbol('amici_' + str)
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
            (sp.Symbol('avogadro'), sp.Symbol('6.02214179e23')),
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


def getRuleVars(rules):
    """Extract free symbols in SBML rule formulas.

    Arguments:
        rules: sbml definitions of rules @type list

    Returns:
    Tuple of free symbolic variables in the formulas all provided rules

    Raises:

    """
    return sp.Matrix(
        [sp.sympify(sbml.formulaToL3String(rule.getMath())) for rule in rules
         if rule.getFormula() != '']
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
    """transforms an list into list of strings

    Arguments:
        inputs: objects @type list

    Returns:
    list of str(object)

    Raises:

    """
    return [str(inp) for inp in inputs]


def checkLibSBMLErrors(sbml_doc, show_warnings=False):
    """Checks the error log in the current self.sbml_doc

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
                category = error.getCategoryAsString()
                severity = error.getSeverityAsString()
                error_message = error.getMessage()
                print(f'libSBML {severity} ({category}): {error_message}')

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
        sp.functions.Piecewise, spTrue, spFalse, sp.function.UndefinedFunction
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
