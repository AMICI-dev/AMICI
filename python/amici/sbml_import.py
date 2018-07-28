""" @package amici.sbml_import The python sbml import module for python 

"""
#!/usr/bin/env python3

import symengine as sp
from symengine.printing import CCodePrinter
from symengine import symbols
import libsbml as sbml
import os
import re
import math
import shutil
import subprocess
import sys
from string import Template

from . import amici_path, amiciSwigPath, amiciSrcPath, amiciModulePath

class SBMLException(Exception):
    pass

class SbmlImporter:
    """The SbmlImporter class generates AMICI C++ files for a model provided in the Systems Biology Markup Language (SBML).
    
    Attributes:
        check_validity: flag indicating whether the validity of the SBML document should be checked
        codeprinter: codeprinter that allows export of symbolic variables as C++ code
        functions: dict carrying function specific definitions
        symbols: symbolic definitions 
        SBMLreader: the libSBML sbml reader [!not storing this will result in a segfault!]
        sbml_doc: document carrying the sbml defintion [!not storing this will result in a segfault!]
        sbml: sbml definition [!not storing this will result in a segfault!]
        modelName: name of the model that will be used for compilation
        modelPath: path to the generated model specific files
        modelSwigPath: path to the generated swig files
        n_species: number of species
        speciesIndex: dict that maps species names to indices
        speciesCompartment: array of compartment for each species
        constantSpecies: ids of species that are marked as constant
        boundaryConditionSpecies: ids of species that are marked as boundary condition
        speciesHasOnlySubstanceUnits: array of flags indicating whether a species has only substance units
        speciesInitial: array of initial concentrations
        speciesConversionFactor: array of conversion factors for every
        parameterValues: array of parameter values
        n_parameters: number of parameters
        parameterIndex: dict that maps parameter names to indices
        compartmentSymbols: array of compartment ids
        compartmentVolume: array of compartment volumes
        n_reactions: number of reactions
        stoichiometricMatrix: stoichiometrix matrix of the model
        fluxVector: vector of reaction kinetic laws 
        observables: array of observable definitions
        n_observables: number of observables
        fixedParameterValues: array of fixed parameter values
        fixedParameterIndex: dict that maps fixed parameter names to indices
        n_fixed_parameters: number of fixed parameters

    """

    def __init__(self, SBMLFile, check_validity=True):
        """Create a new Model instance.
        
        Arguments:
            SBMLFile: Path to SBML file where the model is specified
            check_validity: Flag indicating whether the validity of the SBML document should be checked

        Returns:
            SbmlIMporter instance with attached SBML document

        Raises:

        """

        self.check_validity = check_validity

        self.loadSBMLFile(SBMLFile)

        self.codeprinter = CCodePrinter()

        """Signatures and properties of generated model functions (see include/amici/model.h for details)."""
        self.functions = {
            'J': {
                'signature': '(realtype *J, const realtype t, const realtype *x, const double *p,'
                             ' const double *k, const realtype *h, const realtype *w, const realtype *dwdx)'},
            'JB': {
                'signature': '(realtype *JB, const realtype t, const realtype *x, const realtype *p,'
                             ' const realtype *k, const realtype *h, const realtype *xB, const realtype *w,'
                             ' const realtype *dwdx)'},
            'JDiag': {
                'signature': '(realtype *JDiag, const realtype t, const realtype *x, const realtype *p,'
                             ' const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx)'},
            'JSparse': {
                'signature': '(SlsMat JSparse, const realtype t, const realtype *x, const realtype *p,'
                             ' const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx)',
                'symbol': 'sparseList'},
            'JSparseB': {
                'signature': '(SlsMat JSparseB, const realtype t, const realtype *x, const realtype *p,'
                             ' const realtype *k, const realtype *h, const realtype *xB, const realtype *w,'
                             ' const realtype *dwdx)',
                'symbol': 'sparseList'},
            'Jv': {
                'signature': '(realtype *Jv, const realtype t, const realtype *x, const realtype *p,'
                             ' const realtype *k, const realtype *h, const realtype *v, const realtype *w,'
                             ' const realtype *dwdx)'},
            'JvB': {
                'signature': '(realtype *JvB, const realtype t, const realtype *x, const realtype *p,'
                             ' const realtype *k, const realtype *h, const realtype *xB, const realtype *vB,'
                             ' const realtype *w, const realtype *dwdx)'},
            'Jy': {
                'signature': '(double *nllh, const int iy, const realtype *p, const realtype *k, const double *y,'
                             ' const double *sigmay, const double *my)',
                'variable': 'nllh',
                'multiobs': True},
            'dJydsigma': {
                'signature': '(double *dJydsigma, const int iy, const realtype *p, const realtype *k, const double *y,'
                             ' const double *sigmay, const double *my)',
                'multiobs': True},
            'dJydy': {
                'signature': '(double *dJydy, const int iy, const realtype *p, const realtype *k, const double *y,'
                             ' const double *sigmay, const double *my)',
                'multiobs': True},
            'dwdp': {
                'signature': '(realtype *dwdp, const realtype t, const realtype *x, const realtype *p,'
                             ' const realtype *k, const realtype *h, const realtype *w)',
                'symbol': 'sparseList'},
            'dwdx': {
                'signature': '(realtype *dwdx, const realtype t, const realtype *x, const realtype *p,'
                             ' const realtype *k, const realtype *h, const realtype *w)',
                'symbol': 'sparseList'},
            'dxdotdp': {
                'signature': '(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p,'
                             ' const realtype *k, const realtype *h, const int ip, const realtype *w,'
                             ' const realtype *dwdp)',
                'sensitivity': True},
            'dydx': {
                'signature': '(double *dydx, const realtype t, const realtype *x, const realtype *p,'
                             ' const realtype *k, const realtype *h)'},
            'dydp': {
                'signature': '(double *dydp, const realtype t, const realtype *x, const realtype *p,'
                             ' const realtype *k, const realtype *h, const int ip)',
                'sensitivity': True},
            'dsigmaydp': {
                'signature': '(double *dsigmaydp, const realtype t, const realtype *p,'
                             ' const realtype *k, const int ip)',
                'sensitivity': True},
            'qBdot': {
                'signature': '(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p,'
                             ' const realtype *k, const realtype *h, const realtype *xB, const realtype *w,'
                             ' const realtype *dwdp)',
                'sensitivity': True},
            'sigmay': {'signature': '(double *sigmay, const realtype t, const realtype *p, const realtype *k)',
                        'variable': 'sigmay'},
            'sxdot': {
                'signature': '(realtype *sxdot, const realtype t, const realtype *x, const realtype *p,'
                             ' const realtype *k, const realtype *h, const int ip, const realtype *sx,'
                             ' const realtype *w, const realtype *dwdx, const realtype *J, const realtype *dxdotdp)'},
            'w': {
                'signature': '(realtype *w, const realtype t, const realtype *x, const realtype *p,'
                             ' const realtype *k, const realtype *h)'},
            'x0': {'signature': '(realtype *x0, const realtype t, const realtype *p, const realtype *k)'},
            'sx0': {
                'signature': '(realtype *sx0, const realtype t,const realtype *x0, const realtype *p,'
                             ' const realtype *k, const int ip)'},
            'xBdot': {
                'signature': '(realtype *xBdot, const realtype t, const realtype *x, const realtype *p,'
                             ' const realtype *k, const realtype *h, const realtype *xB, const realtype *w,'
                             ' const realtype *dwdx)'},
            'xdot': {
                'signature': '(realtype *xdot, const realtype t, const realtype *x, const realtype *p,'
                             ' const realtype *k, const realtype *h, const realtype *w)'},
            'y': {
                'signature': '(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k,'
                             ' const realtype *h)'}
        }

        """Long and short names for model components"""
        self.symbols = {
            'species': {'shortName': 'x'},
            'sensitivity': {'shortName': 'sx'},
            'vector': {'shortName': 'v'},
            'vectorB': {'shortName': 'vB'},
            'parameter': {'shortName': 'p'},
            'fixed_parameter': {'shortName': 'k'},
            'observable': {'shortName': 'y'},
            'adjoint': {'shortName': 'xB'},
            'flux': {'shortName': 'w'},
            'dxdotdp': {'shortName': 'dxdotdp'},
            'dwdx': {'shortName': 'dwdx'},
            'dwdp': {'shortName': 'dwdp'},
            'JSparse': {'shortName': 'J'},
            'JSparseB': {'shortName': 'JB'},
            'my': {'shortName': 'my'},
            'sigmay': {'shortName': 'sigmay'}
        }

    def loadSBMLFile(self, SBMLFile):
        """Parse the provided SBML file.
        
        Arguments:
            SBMLFile: path to SBML file

        Returns:

        Raises:
        
        """

        self.SBMLreader = sbml.SBMLReader()
        self.sbml_doc = self.SBMLreader.readSBML(SBMLFile)
        self.checkLibSBMLErrors()
        # If any of the above calls produces an error, this will be added to the SBMLError log in the sbml document.
        # Thus, it is sufficient to check the error log just once after all conversion/validation calls.

        # apply several model simplifications that make our life substantially easier
        if len(self.sbml_doc.getModel().getListOfFunctionDefinitions()) > 0:
            convertConfig = sbml.SBMLFunctionDefinitionConverter().getDefaultProperties()
            self.sbml_doc.convert(convertConfig)

        convertConfig = sbml.SBMLInitialAssignmentConverter().getDefaultProperties()
        self.sbml_doc.convert(convertConfig)

        convertConfig = sbml.SBMLLocalParameterConverter().getDefaultProperties()
        self.sbml_doc.convert(convertConfig)

        if self.check_validity:
            self.sbml_doc.validateSBML()

        # If any of the above calls produces an error, this will be added to the SBMLError log in the sbml document.
        # Thus, it is sufficient to check the error log just once after all conversion/validation calls.
        self.checkLibSBMLErrors()

        self.sbml = self.sbml_doc.getModel()


    def checkLibSBMLErrors(self):
        """Generate AMICI C++ files for the model provided to the constructor.

        Arguments:

        Returns:

        Raises:
            Checks the error log in the current self.sbml_doc and raises an exception if errors with severity ERROR
            or FATAL have occured

        """
        num_warning = self.sbml_doc.getNumErrors(sbml.LIBSBML_SEV_WARNING)
        num_error = self.sbml_doc.getNumErrors(sbml.LIBSBML_SEV_ERROR)
        num_fatal = self.sbml_doc.getNumErrors(sbml.LIBSBML_SEV_FATAL)
        if num_warning + num_error + num_fatal:
            for iError in range(0, self.sbml_doc.getNumErrors()):
                error = self.sbml_doc.getError(iError)
                if error.getSeverity() >= sbml.LIBSBML_SEV_WARNING: # we ignore any info messages for now
                    category = error.getCategoryAsString()
                    severity = error.getSeverityAsString()
                    error_message = error.getMessage()
                    print('libSBML ' + severity + ' (' + category + '): ' + error_message)
        if num_error + num_fatal:
            raise SBMLException('SBML Document failed to load (see error messages above)')


    def sbml2amici(self, modelName, output_dir=None, observables={}, constantParameters=[], sigmas={}, verbose=False):
        """Generate AMICI C++ files for the model provided to the constructor.
        
        Arguments:
            modelName: name of the model/model directory
            output_dir: see sbml_import.setPaths()
            observables: dictionary(observableId:{'name':observableName (optional) ,'formula':formulaString)
                         to be added to the model
            sigmas: dictionary(observableId: sigma value or (existing) parameter name)
            constantParameters: list of SBML Ids identifying constant parameters
            verbose: more verbose output if True
        Returns:

        Raises:

        """

        self.setName(modelName)
        self.setPaths(output_dir)
        self.processSBML(constantParameters)
        self.computeModelEquations(observables, sigmas)
        self.prepareModelFolder()
        self.generateCCode()
        self.compileCCode(verbose)


    def setName(self, modelName):
        """Sets the model name

        Arguments:
            modelName: name of the model (must only contain valid filename characters)

        Returns:

        Raises:

        """
        self.modelName = modelName


    def setPaths(self, output_dir=None):
        """Set output paths for the model and create if necessary

        Arguments:
            output_dir: relative or absolute path where the generated model code is to be placed. will be created if does not exists. defaults to `pwd`/amici-$modelname.

        Returns:

        Raises:

        """

        if not output_dir:
            output_dir = '%s/amici-%s' % (os.getcwd(), self.modelName)

        self.modelPath = os.path.abspath(output_dir)
        self.modelSwigPath = os.path.join(self.modelPath, 'swig')

        for dir in [self.modelPath, self.modelSwigPath]:
            if not os.path.exists(dir):
                os.makedirs(dir)


    def processSBML(self, constantParameters=[]):
        """Read parameters, species, reactions, and so on from SBML model

        Arguments:
            constantParameters: list of SBML Ids identifying constant parameters

        Returns:

        Raises:

        """
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
            raise SBMLException('Models without species are currently not supported!')

        if len(self.sbml.all_elements_from_plugins) > 0:
            raise SBMLException('SBML extensions are currently not supported!')

        if len(self.sbml.getListOfEvents()) > 0:
            raise SBMLException('Events are currently not supported!')

        if any([not(rule.isAssignment()) for rule in self.sbml.getListOfRules()]):
            raise SBMLException('Algebraic and rate rules are currently not supported!')

        if any([reaction.getFast() for reaction in self.sbml.getListOfReactions()]):
            raise SBMLException('Fast reactions are currently not supported!')

        if any([any([not element.getStoichiometryMath() is None
                for element in list(reaction.getListOfReactants()) + list(reaction.getListOfProducts())])
                for reaction in self.sbml.getListOfReactions()]):
            raise SBMLException('Non-unity stoichiometry is currently not supported!')


    def processSpecies(self):
        """Get species information from SBML model.

        Arguments:

        Returns:

        Raises:

        """
        species = self.sbml.getListOfSpecies()
        self.n_species = len(species)
        self.speciesIndex = {species_element.getId(): species_index
                             for species_index, species_element in enumerate(species)}
        self.symbols['species']['sym'] = sp.DenseMatrix([symbols(spec.getId()) for spec in species])
        self.symbols['species']['name'] = [spec.getName() for spec in species]
        self.speciesCompartment = sp.DenseMatrix([symbols(spec.getCompartment()) for spec in species])
        self.constantSpecies = [species_element.getId() if species_element.getConstant() else None
                                for species_element in species]
        self.boundaryConditionSpecies = [species_element.getId() if species_element.getBoundaryCondition() else None
                                         for species_element in species]
        self.speciesHasOnlySubstanceUnits = [specie.getHasOnlySubstanceUnits() for specie in species]

        concentrations = [spec.getInitialConcentration() for spec in species]
        amounts = [spec.getInitialAmount() for spec in species]

        self.speciesInitial = sp.DenseMatrix([sp.sympify(conc) if not math.isnan(conc) else
                                              sp.sympify(amounts[index])/self.speciesCompartment[index] if not
                                              math.isnan(amounts[index]) else
                                              self.symbols['species']['sym'][index]
                                              for index, conc in enumerate(concentrations)])

        if self.sbml.isSetConversionFactor():
            conversionFactor = self.sbml.getConversionFactor()
        else:
            conversionFactor = 1.0
        self.speciesConversionFactor = sp.DenseMatrix([sp.sympify(specie.getConversionFactor()) if
                                                       specie.isSetConversionFactor() else conversionFactor
                                                       for specie in species])


    def processParameters(self, constantParameters=[]):
        """Get parameter information from SBML model.

        Arguments:
            constantParameters: list of SBML Ids identifying constant parameters
        Returns:

        Raises:

        """

        fixedParameters = [ parameter for parameter in self.sbml.getListOfParameters() if parameter.getId() in constantParameters ]
        parameters = [ parameter for parameter in self.sbml.getListOfParameters() if parameter.getId() not in constantParameters ]

        self.symbols['parameter']['sym'] = sp.DenseMatrix([symbols(par.getId()) for par in parameters])
        self.symbols['parameter']['name'] = [par.getName() for par in parameters]
        self.parameterValues = [par.getValue() for par in parameters]
        self.n_parameters = len(self.symbols['parameter']['sym'])
        self.parameterIndex = {parameter_element.getId(): parameter_index
                             for parameter_index, parameter_element in enumerate(parameters)}

        self.symbols['fixed_parameter']['sym'] = sp.DenseMatrix([symbols(par.getId()) for par in fixedParameters])
        self.symbols['fixed_parameter']['name'] = [par.getName() for par in fixedParameters]
        self.fixedParameterValues = [par.getValue() for par in fixedParameters]
        self.n_fixed_parameters = len(self.symbols['fixed_parameter']['sym'])
        self.fixedParameterIndex = {parameter_element.getId(): parameter_index
                             for parameter_index, parameter_element in enumerate(fixedParameters)}


    def processCompartments(self):
        """Get compartment information, stoichiometric matrix and fluxes from SBML model.

        Arguments:

        Returns:

        Raises:

        """
        compartments = self.sbml.getListOfCompartments()
        self.compartmentSymbols = sp.DenseMatrix([symbols(comp.getId()) for comp in compartments])
        self.compartmentVolume = sp.DenseMatrix([sp.sympify(comp.getVolume()) if comp.isSetVolume()
                                                 else sp.sympify(1.0) for comp in compartments])

        self.replaceInAllExpressions(self.compartmentSymbols,self.compartmentVolume)


    def processReactions(self):
        """Get reactions from SBML model.

        Arguments:

        Returns:

        Raises:

        """
        reactions = self.sbml.getListOfReactions()
        self.n_reactions = len(reactions)

        # stoichiometric matrix
        self.stoichiometricMatrix = sp.zeros(self.n_species, self.n_reactions)
        self.fluxVector = sp.zeros(self.n_reactions, 1)
        self.symbols['flux']['sym'] = sp.zeros(self.n_reactions, 1)

        for reactionIndex, reaction in enumerate(reactions):

            for elementList, sign  in [(reaction.getListOfReactants(), -1.0),
                                       (reaction.getListOfProducts(), 1.0)]:
                elements = {}
                for index, element in enumerate(elementList):
                    # we need the index here as we might have multiple elements for the same species
                    elements[index] = {'species': element.getSpecies()}
                    if element.isSetId():
                        if element.getId() in [rule.getVariable()
                                               for rule in self.sbml.getListOfRules() if rule.getFormula() != '']:
                            elements[index]['stoichiometry'] = sp.sympify(element.getId())
                        else:
                            # dont put the symbol if it wont get replaced by a rule
                            elements[index]['stoichiometry'] = sp.sympify(element.getStoichiometry())
                    else:
                        elements[index]['stoichiometry'] = sp.sympify(element.getStoichiometry())

                for index in elements.keys():
                    if not (elements[index]['species'] in self.constantSpecies
                            or elements[index]['species'] in self.boundaryConditionSpecies):
                        specieIndex = self.speciesIndex[elements[index]['species']]
                        self.stoichiometricMatrix[specieIndex, reactionIndex] += sign * \
                            elements[index]['stoichiometry'] * self.speciesConversionFactor[specieIndex]/ \
                            self.speciesCompartment[specieIndex]

            # usage of formulaToL3String ensures that we get "time" as time symbol
            math = sbml.formulaToL3String(reaction.getKineticLaw().getMath())
            try:
                symMath = sp.sympify(math)
            except:
                raise SBMLException('Kinetic law "' + math + '" contains an unsupported expression!')
            for r in reactions:
                for element in list(r.getListOfReactants()) + list(r.getListOfProducts()):
                    if element.isSetId() & element.isSetStoichiometry():
                        symMath = symMath.subs(sp.sympify(element.getId()),sp.sympify(element.getStoichiometry()))
            self.fluxVector[reactionIndex] = symMath
            self.symbols['flux']['sym'][reactionIndex] = sp.sympify('w' + str(reactionIndex))
            if any([str(symbol) in [reaction.getId() for reaction in reactions if reaction.isSetId()]
                         for symbol in self.fluxVector[reactionIndex].free_symbols]):
                raise SBMLException('Kinetic laws involving reaction ids are currently not supported!')


    def processRules(self):
        """Process *Rules defined in the SBML model.

        Arguments:

        Returns:

        Raises:

        """
        rules = self.sbml.getListOfRules()

        rulevars = getRuleVars(rules)
        fluxvars = self.fluxVector.free_symbols
        specvars = self.symbols['species']['sym'].free_symbols
        volumevars = self.compartmentVolume.free_symbols
        compartmentvars = self.compartmentSymbols.free_symbols
        parametervars = self.symbols['parameter']['sym'].free_symbols
        stoichvars = self.stoichiometricMatrix.free_symbols

        for rule in rules:
            if rule.getFormula() == '':
                continue
            variable = sp.sympify(rule.getVariable())
            formula = sp.sympify(rule.getFormula())

            if variable in stoichvars:
                self.stoichiometricMatrix = self.stoichiometricMatrix.subs(variable, formula)

            if variable in specvars:
                raise SBMLException('Species assignment rules are currently not supported!')

            if variable in compartmentvars:
                raise SBMLException('Compartment assignment rules are currently not supported!')

            if variable in parametervars:
                try:
                    self.parameterValues[self.parameterIndex[str(variable)]] = float(formula)
                except:
                    raise SBMLException('Non-float parameter assignment rules are currently not supported!')

            if variable in fluxvars:
                self.fluxVector = self.fluxVector.subs(variable, formula)

            if variable in volumevars:
                self.compartmentVolume = self.compartmentVolume.subs(variable, formula)

            if variable in rulevars:
                for nested_rule in rules:
                    nested_formula = sp.sympify(nested_rule.getFormula())
                    nested_formula.subs(variable, formula)
                    nested_rule.setFormula(str(nested_formula))

    def processVolumeConversion(self):
        """Convert equations from amount to volume.

        Arguments:

        Returns:

        Raises:

        """
        for index, bool in enumerate(self.speciesHasOnlySubstanceUnits):
            if bool:
                self.fluxVector = self.fluxVector.subs(self.symbols['species']['sym'][index],
                                                self.symbols['species']['sym'][index] * self.speciesCompartment[index]
                                                            .subs(self.compartmentSymbols,
                                                                  self.compartmentVolume))


    def processTime(self):
        """Convert time_symbol into cpp variable.

        Arguments:

        Returns:

        Raises:

        """
        sbmlTimeSymbol = sp.sympify('time')
        amiciTimeSymbol = sp.sympify('t')

        self.replaceInAllExpressions(sbmlTimeSymbol, amiciTimeSymbol)



    def replaceInAllExpressions(self,old,new):
        """Replace 'old' by 'new' in all symbolic expressions.

        Arguments:
            old: symbolic variables to be replaced
            new: replacement symbolic variables

        Returns:

        Raises:

        """
        fields = ['observables', 'stoichiometricMatrix', 'fluxVector', 'speciesInitial']
        for field in fields:
            if field in dir(self):
                self.__setattr__(field, self.__getattribute__(field).subs(old, new))


    def cleanReservedSymbols(self):
        """Remove all reserved symbols from self.symbols

        Arguments:

        Returns:

        Raises:

        """
        reservedSymbols = ['k','p','y','w']
        for str in reservedSymbols:
            old_symbol = sp.sympify(str)
            new_symbol = sp.sympify('amici_' + str)
            self.replaceInAllExpressions(old_symbol, new_symbol)
            for symbol in self.symbols.keys():
                if 'sym' in self.symbols[symbol].keys():
                    self.symbols[symbol]['sym'] = self.symbols[symbol]['sym'].subs(old_symbol,new_symbol)

    def replaceSpecialConstants(self):
        """Replace all special constants by their respective SBML csymbol definition

        Arguments:

        Returns:

        Raises:

        """
        constants = [(sp.sympify('avogadro'),sp.sympify('6.02214179*1e23')),]
        for constant, value in constants:
            # do not replace if any symbol is shadowing default definition
            if not any([constant in self.symbols[symbol]['sym']
                        for symbol in self.symbols.keys() if 'sym' in self.symbols[symbol].keys()]):
                self.replaceInAllExpressions(constant, value)
            else:
                # yes sbml supports this but we wont, are you really expecting to be saved if you are trying to shoot
                # yourself in the foot?
                raise SBMLException('Ecountered currently unsupported element id ' + str(constant) + '!')


    def getSparseSymbols(self,symbolName):
        """Create sparse symbolic matrix.

        Arguments:
            symbolName: name of the function

        Returns:
            sparseMatrix: sparse matrix containing symbolic entries
            symbolList: symbolic vector containing the list of symbol names
            sparseList: symbolic vector containing the list of symbol formulas
            symbolColPtrs: Column pointer as specified in the SlsMat definition in CVODES
            symbolRowVals: Row Values as specified in the SlsMat definition in CVODES

        Raises:

        """
        matrix = self.functions[symbolName]['sym']
        symbolIndex = 0
        sparseMatrix = sp.zeros(matrix.rows,matrix.cols)
        symbolList = []
        sparseList = []
        symbolColPtrs = []
        symbolRowVals = []
        for col in range(0,matrix.cols):
            symbolColPtrs.append(symbolIndex)
            for row in range(0, matrix.rows):
                if(not(matrix[row, col]==0)):
                    name = symbolName + str(symbolIndex)
                    sparseMatrix[row, col] = sp.sympify(name)
                    symbolList.append(name)
                    sparseList.append(matrix[row, col])
                    symbolRowVals.append(row)
                    symbolIndex += 1
        symbolColPtrs.append(symbolIndex)
        sparseList = sp.DenseMatrix(sparseList)
        return sparseMatrix, symbolList, sparseList, symbolColPtrs, symbolRowVals

    def computeModelEquations(self, observables={}, sigmas={}):
        """Perform symbolic computations required to populate functions in `self.functions`.

        Arguments:
            observables: dictionary(observableId:{'name':observableName (optional) ,'formula':formulaString)
                         to be added to the model
            sigmas: dictionary(observableId: sigma value or (existing) parameter name)

        Returns:

        Raises:

        """

        # core
        self.functions['xdot']['sym'] = self.stoichiometricMatrix * self.symbols['flux']['sym']
        self.computeModelEquationsLinearSolver()
        self.computeModelEquationsObjectiveFunction(observables, sigmas)

        # sensitivity
        self.computeModelEquationsSensitivitesCore()
        self.computeModelEquationsForwardSensitivites()
        self.computeModelEquationsAdjointSensitivites()


    def computeModelEquationsLinearSolver(self):
        """Perform symbolic computations required for the use of various linear solvers.

        Arguments:

        Returns:

        Raises:

        """
        self.functions['w']['sym'] = self.fluxVector
        self.functions['dwdx']['sym'] = self.fluxVector.jacobian(self.symbols['species']['sym'])

        self.functions['dwdx']['sparseSym'],\
        self.symbols['dwdx']['sym'],\
        self.functions['dwdx']['sparseList']  = self.getSparseSymbols('dwdx')[0:3]

        self.functions['J']['sym'] = self.functions['xdot']['sym'].jacobian(self.symbols['species']['sym']) \
                                        + self.stoichiometricMatrix * self.functions['dwdx']['sparseSym']
        self.symbols['vector']['sym'] = getSymbols('v',self.n_species)
        self.functions['Jv']['sym'] = self.functions['J']['sym']*self.symbols['vector']['sym']

        self.functions['JSparse']['sym'],\
        self.symbols['JSparse']['sym'],\
        self.functions['JSparse']['sparseList'],\
        self.functions['JSparse']['colPtrs'],\
        self.functions['JSparse']['rowVals'] = self.getSparseSymbols('J')

        self.functions['x0']['sym'] = self.speciesInitial
        self.functions['JDiag']['sym'] = getSymbolicDiagonal(self.functions['J']['sym'])


    def computeModelEquationsObjectiveFunction(self, observables={}, sigmas={}):
        """Perform symbolic computations required for objective function evaluation.

        Arguments:
            observables: dictionary(observableId:{'name':observableName (optional) ,'formula':formulaString)
                         to be added to the model
            sigmas: dictionary(observableId: sigma value or (existing) parameter name)

        Returns:

        Raises:

        """
        speciesSyms = self.symbols['species']['sym']

        # add user-provided observables or make all species observable
        if(observables):
            self.observables = sp.DenseMatrix([observables[observable]['formula'] for observable in observables])
            observableNames = [observables[observable]['name'] if 'name' in observables[observable].keys()
                               else 'y' + str(index)
                               for index, observable in enumerate(observables)]
            observableSyms = sp.DenseMatrix(observables.keys())
            self.n_observables = len(observables)
        else:
            self.observables = speciesSyms
            observableNames = ['x' + str(index) for index in range(0,len(speciesSyms))]
            observableSyms = sp.DenseMatrix([sp.sympify('y' + str(index)) for index in range(0,len(speciesSyms))])
            self.n_observables = len(speciesSyms)
        self.functions['y']['sym'] = self.observables
        sigmaYSyms = sp.DenseMatrix([sp.sympify('sigma' + str(symbol)) for symbol in observableSyms])
        self.functions['sigmay']['sym'] = sp.DenseMatrix([1.0] * len(observableSyms))

        # set user-provided sigmas
        for iy, obsName in enumerate(observables):
            if obsName in sigmas:
                self.functions['sigmay']['sym'][iy] = sigmas[obsName]

        self.symbols['my']['sym'] = sp.DenseMatrix([sp.sympify('m' + str(symbol)) for symbol in observableSyms])

        loglikelihoodString = lambda strSymbol: '0.5*log(2*pi*sigma{symbol}**2) + 0.5*(({symbol}-m{symbol})/sigma{symbol})**2'.format(symbol=strSymbol)
        self.functions['Jy']['sym'] = sp.DenseMatrix([sp.sympify(loglikelihoodString(str(symbol))) for symbol in observableSyms])
        self.functions['dJydy']['sym']     = self.functions['Jy']['sym'].jacobian(observableSyms).transpose()
        self.functions['dJydsigma']['sym'] = self.functions['Jy']['sym'].jacobian(sigmaYSyms).transpose()
        self.functions['Jy']['sym']        = self.functions['Jy']['sym'].transpose()

        self.symbols['observable']['sym'] = observableSyms
        self.symbols['observable']['name'] = observableNames
        self.symbols['sigmay']['sym'] = sigmaYSyms


    def computeModelEquationsSensitivitesCore(self):
        """Perform symbolic computations required for any sensitivity analysis.

        Arguments:

        Returns:

        Raises:

        """
        self.functions['dydp']['sym'] = self.functions['y']['sym']\
                                                .jacobian(self.symbols['parameter']['sym'])
        self.functions['dsigmaydp']['sym'] = self.functions['sigmay']['sym']\
                                                .jacobian(self.symbols['parameter']['sym'])
        self.functions['dydx']['sym'] = self.functions['y']['sym']\
                                                .jacobian(self.symbols['species']['sym']).transpose()

        self.functions['dwdp']['sym'] = self.fluxVector.jacobian(self.symbols['parameter']['sym'])

        self.functions['dwdp']['sparseSym'],\
        self.symbols['dwdp']['sym'],\
        self.functions['dwdp']['sparseList'] = self.getSparseSymbols('dwdp')[0:3]

        self.functions['dxdotdp']['sym'] = self.functions['xdot']['sym']\
                                                      .jacobian(self.symbols['parameter']['sym'])\
                                                  + self.stoichiometricMatrix\
                                                    * self.functions['dwdp']['sparseSym']
        self.symbols['dxdotdp']['sym'] = getSymbols('dxdotdp',self.n_species)
        self.functions['sx0']['sym'] = self.speciesInitial.jacobian(self.symbols['parameter']['sym'])


    def computeModelEquationsForwardSensitivites(self):
        """Perform symbolic computations required for forward sensitivity analysis.

        Arguments:

        Returns:

        Raises:

        """
        self.symbols['sensitivity']['sym'] = getSymbols('sx',self.n_species)
        self.functions['sxdot']['sym'] = self.functions['JSparse']['sym']\
                                                * self.symbols['sensitivity']['sym'] \
                                                + self.symbols['dxdotdp']['sym']

    def computeModelEquationsAdjointSensitivites(self):
        """Perform symbolic computations required for adjoint sensitivity analysis.

        Arguments:

        Returns:

        Raises:

        """
        self.functions['JB']['sym'] = self.functions['J']['sym'].transpose()
        self.symbols['vectorB']['sym'] = getSymbols('vB', self.n_species)
        self.functions['JvB']['sym'] = self.functions['JB']['sym'] \
                                              * self.symbols['vectorB']['sym']
        self.functions['JSparseB']['sym'], self.symbols['JSparseB']['sym'],\
            self.functions['JSparseB']['sparseList'], self.functions['JSparseB']['colPtrs'],\
            self.functions['JSparseB']['rowVals'] = self.getSparseSymbols('JB')

        self.symbols['adjoint']['sym'] = getSymbols('xB',self.n_species)
        self.functions['xBdot']['sym'] = - self.functions['JB']['sym']\
                                                  * self.symbols['adjoint']['sym']
        self.functions['qBdot']['sym'] = - self.symbols['adjoint']['sym'].transpose()\
                                                  * self.functions['dxdotdp']['sym']


    def prepareModelFolder(self):
        """Remove all files from the model folder.

        Arguments:

        Returns:

        Raises:

        """
        for file in os.listdir(self.modelPath):
            file_path = os.path.join(self.modelPath, file)
            if os.path.isfile(file_path):
                os.remove(file_path)



    def generateCCode(self):
        """Create C++ code files for the model based on.

        Arguments:

        Returns:

        Raises:

        """
        for name in self.symbols.keys():
            self.writeIndexFiles(name)

        for function in self.functions.keys():
            self.writeFunctionFile(function)

        self.writeWrapfunctionsCPP()
        self.writeWrapfunctionsHeader()
        self.writeModelHeader()
        self.writeCMakeFile()
        self.writeSwigFiles()
        self.writeModuleSetup()

        shutil.copy(os.path.join(amiciSrcPath, 'main.template.cpp'),
                    os.path.join(self.modelPath, 'main.cpp'))


    def compileCCode(self, verbose=False):
        """Compile the generated model code

        Arguments:
            verbose: Make model compilation verbose
        Returns:

        Raises:

        """

        # setup.py assumes it is run from within the model directory
        moduleDir = self.modelPath

        script_args = [sys.executable, '%s%ssetup.py' % (moduleDir, os.sep)]

        if verbose:
            script_args.append('--verbose')
        else:
            script_args.append('--quiet')

        script_args.extend(['build_ext', '--build-lib=%s' % moduleDir])

        # distutils.core.run_setup looks nicer, but does not let us check the
        # result easily
        try:
            result = subprocess.run(script_args,
                                    cwd=moduleDir,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT,
                                    check=True)
        except subprocess.CalledProcessError as e:
            print(e.output.decode('utf-8'))
            raise

        if verbose:
            print(result.stdout.decode('utf-8'))

    def writeIndexFiles(self,name):
        """Write index file for a symbolic array.

        Arguments:
            name: key in self.symbols for which the respective file should be written

        Returns:

        Raises:


        """
        lines = []
        [lines.append('#define ' + str(symbol) + ' ' + str(self.symbols[name]['shortName']) + '[' + str(index) + ']')
            for index, symbol in enumerate(self.symbols[name]['sym'])]

        with open(os.path.join(self.modelPath,name + '.h'), 'w') as fileout:
            fileout.write('\n'.join(lines))


    def writeFunctionFile(self,function):
        """Write the function `function`.

        Arguments:
            function: name of the function to be written (see self.functions)

        Returns:

        Raises:


        """

        # function header
        lines = ['#include "amici/symbolic_functions.h"',
            '#include "amici/defines.h" //realtype definition',
            'using amici::realtype;',
            '#include <cmath> ',
            '']

        # function signature
        signature = self.functions[function]['signature']

        if(not signature.find('SlsMat') == -1):
            lines.append('#include <sundials/sundials_sparse.h>')

        lines.append('')

        for symbol in self.symbols.keys():
            # added |double for data
            # added '[0]*' for initial conditions
            if not re.search('const (realtype|double) \*' + self.symbols[symbol]['shortName'] + '[0]*[,)]+',
                             signature) is None :
                lines.append('#include "' + symbol + '.h"')

        lines.append('')

        lines.append('void ' + function + '_' + self.modelName + signature + '{')

        # function body
        body = self.getFunctionBody(function)
        self.functions[function]['body'] = body
        lines += body
        lines.append('}')
        #if not body is None:
        with open(os.path.join(self.modelPath, self.modelName + '_' + function + '.cpp'), 'w') as fileout:
            fileout.write('\n'.join(lines))


    def getFunctionBody(self,function):
        """Generate C++ code for body of function `function`.

        Arguments:
            function: name of the function to be written (see self.functions)

        Returns:

        Raises:

        """

        if('variable' in self.functions[function].keys()):
            variableName = self.functions[function]['variable']
        else:
            variableName = function

        if ('symbol' in self.functions[function].keys()):
            symbol = self.functions[function][self.functions[function]['symbol']]
        else:
            symbol = self.functions[function]['sym']
        lines = []


        if('sensitivity' in self.functions[function].keys()):
            lines.append(' '*4 + 'switch(ip) {')
            for ipar in range(0,self.n_parameters):
                lines.append(' ' * 8 + 'case ' + str(ipar) + ':')
                lines += self.getSymLines(symbol[:, ipar], variableName, 12)
                lines.append(' ' * 12 + 'break;')
            lines.append('}')
        elif('multiobs' in self.functions[function].keys()):
            lines.append(' '*4 + 'switch(iy) {')
            for iobs in range(0,self.n_observables):
                lines.append(' ' * 8 + 'case ' + str(iobs) + ':')
                lines += self.getSymLines(symbol[:, iobs], variableName, 12)
                lines.append(' ' * 12 + 'break;')
            lines.append('}')
        else:
            if('colPtrs' in self.functions[function].keys()):
                rowVals = self.functions[function]['rowVals']
                colPtrs = self.functions[function]['colPtrs']
                lines += self.getSparseSymLines(symbol, rowVals, colPtrs, variableName, 4)
            else:
                lines += self.getSymLines(symbol, variableName, 4)

        return [line for line in lines if line]


    def writeWrapfunctionsCPP(self):
        """Write model-specific 'wrapper' file (wrapfunctions.cpp).

        Arguments:

        Returns:

        Raises:

        """
        templateData = {'MODELNAME': self.modelName}
        applyTemplate(os.path.join(amiciSrcPath, 'wrapfunctions.template.cpp'),
                      os.path.join(self.modelPath, 'wrapfunctions.cpp'), templateData)


    def writeWrapfunctionsHeader(self):
        """Write model-specific header file (wrapfunctions.h).

        Arguments:

        Returns:

        Raises:

        """
        templateData = {'MODELNAME': str(self.modelName)}
        applyTemplate(os.path.join(amiciSrcPath, 'wrapfunctions.ODE_template.h'),
                      os.path.join(self.modelPath, 'wrapfunctions.h'), templateData)

    def writeModelHeader(self):
        """Write model-specific header file (MODELNAME.h).

        Arguments:

        Returns:

        Raises:

        """
        templateData = {'MODELNAME': str(self.modelName),
                        'NX': str(self.n_species),
                        'NXTRUE': str(self.n_species),
                        'NY': str(self.n_observables),
                        'NYTRUE': str(self.n_observables),
                        'NZ': '0',
                        'NZTRUE': '0',
                        'NEVENT': '0',
                        'NOBJECTIVE': '1',
                        'NW': str(len(self.symbols['flux']['sym'])),
                        'NDWDDX': str(len(self.functions['dwdx']['sparseList'])),
                        'NDWDP': str(len(self.functions['dwdp']['sparseList'])),
                        'NNZ': str(len(self.functions['JSparse']['sparseList'])),
                        'UBW': str(self.n_species),
                        'LBW': str(self.n_species),
                        'NP': str(self.n_parameters),
                        'NK': str(self.n_fixed_parameters),
                        'O2MODE': 'amici::AMICI_O2MODE_NONE',
                        'PARAMETERS': str(self.parameterValues)[1:-1],
                        'FIXED_PARAMETERS': str(self.fixedParameterValues)[1:-1],
                        'PARAMETER_NAMES_INITIALIZER_LIST': self.getSymbolNameInitializerList('parameter'),
                        'STATE_NAMES_INITIALIZER_LIST': self.getSymbolNameInitializerList('species'),
                        'FIXED_PARAMETER_NAMES_INITIALIZER_LIST': self.getSymbolNameInitializerList('fixed_parameter'),
                        'OBSERVABLE_NAMES_INITIALIZER_LIST': self.getSymbolNameInitializerList('observable'),
                        'PARAMETER_IDS_INITIALIZER_LIST': self.getSymbolIDInitializerList('parameter'),
                        'STATE_IDS_INITIALIZER_LIST': self.getSymbolIDInitializerList('species'),
                        'FIXED_PARAMETER_IDS_INITIALIZER_LIST': self.getSymbolIDInitializerList('fixed_parameter'),
                        'OBSERVABLE_IDS_INITIALIZER_LIST': self.getSymbolIDInitializerList('observable'),
                        }
        applyTemplate(os.path.join(amiciSrcPath, 'model_header.ODE_template.h'),
                      os.path.join(self.modelPath, str(self.modelName) + '.h'), templateData)

    def getSymbolNameInitializerList(self, name):
        """Get SBML name initializer list for vector of names for the given model entity

        Arguments:
        name: string, any key present in self.symbols

        Returns:
            Template initializer list of names

        Raises:

        """
        return '\n'.join([ '"%s",' % str(symbol) for symbol in self.symbols[name]['name'] ])

    def getSymbolIDInitializerList(self, name):
        """Get C++ initializer list for vector of names for the given model entity

        Arguments:
        name: string, any key present in self.symbols

        Returns:
            Template initializer list of ids

        Raises:

        """
        return '\n'.join([ '"%s",' % str(symbol) for symbol in self.symbols[name]['sym'] ])

    def writeCMakeFile(self):
        """Write CMake CMakeLists.txt file for this model.

        Arguments:

        Returns:

        Raises:

        """
        sources = [self.modelName + '_' + function + '.cpp ' if self.functions[function]['body'] is not None else ''
                    for function in self.functions.keys() ]
        try:
            sources.remove('')
        except:
            pass
        templateData = {'MODELNAME': self.modelName,
                        'SOURCES': '\n'.join(sources)}
        applyTemplate(os.path.join(amiciSrcPath, 'CMakeLists.template.txt'),
                      os.path.join(self.modelPath, 'CMakeLists.txt'), templateData)

    def writeSwigFiles(self):
        """Write SWIG interface files for this model.

        Arguments:

        Returns:

        Raises:

        """
        if not os.path.exists(self.modelSwigPath):
            os.makedirs(self.modelSwigPath)
        templateData = {'MODELNAME': self.modelName}
        applyTemplate(os.path.join(amiciSwigPath, 'modelname.template.i'),
                      os.path.join(self.modelSwigPath, self.modelName + '.i'), templateData)
        shutil.copy(os.path.join(amiciSwigPath, 'CMakeLists_model.txt'),
                    os.path.join(self.modelSwigPath, 'CMakeLists.txt'))

    def writeModuleSetup(self):
        """Create a distutils setup.py file for compile the model module.

        Arguments:

        Returns:

        Raises:

        """

        templateData = {'MODELNAME': self.modelName,
                        'VERSION': '0.1.0'}
        applyTemplate(os.path.join(amiciModulePath, 'setup.template.py'),
                      os.path.join(self.modelPath, 'setup.py'), templateData)

        # write __init__.py for the model module
        if not os.path.exists(os.path.join(self.modelPath, self.modelName)):
            os.makedirs(os.path.join(self.modelPath, self.modelName))

        applyTemplate(os.path.join(amiciModulePath, '__init__.template.py'),
                      os.path.join(self.modelPath, self.modelName, '__init__.py'), templateData)

    def getSymLines(self, symbols, variable, indentLevel):
        """Generate C++ code for assigning symbolic terms in symbols to C++ array `variable`.

        Arguments:
            symbols: vectors of symbolic terms
            variable: name of the C++ array to assign to
            indentLevel: indentation level (number of leading blanks)

        Returns:
            C++ code as list of lines

        Raises:

        """

        lines = [' ' * indentLevel + variable + '[' + str(index) + '] = ' + self.printWithException(math) + ';'
                if not math == 0 else '' for index, math in enumerate(symbols)]

        try:
            lines.remove('')
        except:
            pass

        return lines

    def getSparseSymLines(self, symbolList, RowVals, ColPtrs, variable, indentLevel):
        """Generate C++ code for assigning sparse symbolic matrix to a C++ array `variable`.

        Arguments:
            symbolList: vectors of symbolic terms
            RowVals: list of row indices of each nonzero entry (see CVODES SlsMat documentation for details)
            ColPtrs: list of indices of the first column entries (see CVODES SlsMat documentation for details)
            variable: name of the C++ array to assign to
            indentLevel: indentation level (number of leading blanks)

        Returns:
            C++ code as list of lines

        Raises:

        """
        lines = [
            ' ' * indentLevel + variable + '->indexvals[' + str(index) + '] = ' + self.printWithException(math) + ';'
            for index, math in enumerate(RowVals)]
        lines.extend([' ' * indentLevel + variable + '->indexptrs[' + str(index) + '] = ' + self.printWithException(math) + ';'
                      for index, math in enumerate(ColPtrs)])
        lines.extend([' ' * indentLevel + variable + '->data[' + str(index) + '] = ' + self.printWithException(math) + ';'
                      for index, math in enumerate(symbolList)])

        return lines

    def printWithException(self,math):
        """Generate C++ code for a symbolic expression

        Arguments:
            math: symbolic expression

        Returns:
            C++ code for the specified expression

        Raises:
            SBMLException: The specified expression containted an unsupported function

        """
        try:
            return self.codeprinter.doprint(math)
        except:
            raise SBMLException('Encountered unsupported function in expression "' + str(math) + '"!')

def applyTemplate(sourceFile,targetFile,templateData):
    """Load source file, apply template substitution as provided in templateData and save as targetFile.

    Arguments:
        sourceFile: relative or absolute path to template file
        targetFile: relative or absolute path to output file
        templateData: dictionary with template keywords to substitute (key is template veriable without TemplateAmici.delimiter)

    Returns:

    Raises:

    """
    with open(sourceFile) as filein:
        src = TemplateAmici(filein.read())
    result = src.safe_substitute(templateData)
    with open(targetFile, 'w') as fileout:
        fileout.write(result)

def getSymbols(prefix,length):
    """Get symbolic matrix with symbols prefix0..prefix(length-1).

    Arguments:
        prefix: variable name
        length: number of symbolic variables + 1

    Returns:
        A symbolic matrix with symbols prefix0..prefix(length-1)

    Raises:

    """
    return sp.DenseMatrix([sp.sympify(prefix + str(i)) for i in range(0, length)])


def getSymbolicDiagonal(matrix):
    """Get symbolic matrix with diagonal of matrix `matrix`.

    Arguments:
        matrix: Matrix from which to return the diagonal

    Returns:
        A Symbolic matrix with the diagonal of `matrix`.

    Raises:
        Exception: The provided matrix was not square
    """
    if not matrix.cols == matrix.rows:
        raise Exception('Provided matrix is not square!')

    diagonal = [matrix[index,index] for index in range(matrix.cols)]

    return sp.DenseMatrix(diagonal)

def getRuleVars(rules):
    """Extract free symbols in SBML rule formulas.

    Argumentss:
        rules: list of rules

    Returns:
        Vector of free symbolic variables in the formulas all provided rules

    Raises:

    """
    return sp.DenseMatrix([sp.sympify(rule.getFormula()) for rule in rules if rule.getFormula() != '']).free_symbols

class TemplateAmici(Template):
    """Template format used in AMICI (see string.template for more details).

    Attributes:
        delimiter: delimiter that identifies template variables

    Returns:

    Raises:

    """
    delimiter = 'TPL_'


def assignmentRules2observables(sbml, filter = lambda *_: True):
    """Turn assignment rules into observables.

    Arguments:
        sbml: an sbml Model instance
        filter: callback function taking assignment variable as input and returning True/False to indicate if the respective rule should be turned into an observable

    Returns:
        A dictionary(observableId:{'name':observableNamem,'formula':formulaString})

    Raises:

    """
    observables = {}
    for p in sbml.getListOfParameters():
        parameterId = p.getId()
        if filter(parameterId):
            observables[parameterId] = {'name': p.getName(), 'formula': sbml.getAssignmentRuleByVariable(parameterId).getFormula()}

    for parameterId in observables:
        sbml.removeRuleByVariable(parameterId)
        sbml.removeParameter(parameterId)

    return observables