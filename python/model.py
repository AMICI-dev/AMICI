#!/usr/bin/env python3

import symengine as sp
from symengine.printing import CCodePrinter
import libsbml as sbml
import os
import re
import math
import shutil
import subprocess
from symengine import symbols
from string import Template


class Model:
    """The Model class generates AMICI C++ files for a model provided in the Systems Biology Markup Language (SBML).
    
    Attributes:
    -----------
    SBMLreader = sbml.SBMLReader()
        sbml_doc: private?
        sbml: private?
        modelname:
        amiciPath:
        amiciSwigPath:
        amiciSrcPath:
        modelPath:
        modelSwigPath:
        self.functionBodies:
    """

    # TODO: are units respected on sbml import? if not convert; at least throw if differ?
    # TODO: camelCase?


    # TODO: move arguments to wrap model? 
    def __init__(self, SBMLFile):
        """Create a new Model instance.
        
        Arguments:
        ----------
        SBMLFile: Path to SBML file where the model is specified
        modelname: Name of Model to be generated
        """
        self.loadSBMLFile(SBMLFile)

        # TODO: move amici*path to amici module?
        dirname, filename = os.path.split(os.path.abspath(__file__))
        self.amiciPath = os.path.split(dirname)[0]
        self.amiciSwigPath = os.path.join(self.amiciPath, 'swig')
        self.amiciSrcPath = os.path.join(self.amiciPath, 'src')

        self.functionBodies = {} # TODO: "private" ?
        self.Codeprinter = CCodePrinter()

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
                             ' const realtype *k, const realtype *h)'},
            'qBdot': {
                'signature': '(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p,'
                             ' const realtype *k, const realtype *h, const realtype *xB, const realtype *w,'
                             ' const realtype *dwdp)',
                'sensitivity': True},
            'sigma_y': {'signature': '(double *sigmay, const realtype t, const realtype *p, const realtype *k)',
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
            'observable': {'shortName': 'y'},
            'adjoint': {'shortName': 'xB'},
            'flux': {'shortName': 'w'},
            'dxdotdp': {'shortName': 'dxdotdp'},
            'dwdx': {'shortName': 'dwdx'},
            'dwdp': {'shortName': 'dwdp'},
            'JSparse': {'shortName': 'J'},
            'JSparseB': {'shortName': 'JB'},
            'my': {'shortName': 'my'},
            'sigma_y': {'shortName': 'sigmay'}
        }

    def loadSBMLFile(self, SBMLFile):
        """Parse the provided SBML file."""
        self.SBMLreader = sbml.SBMLReader()
        self.sbml_doc = self.SBMLreader.readSBML(SBMLFile)

        if (self.sbml_doc.getNumErrors() > 0):
            raise Exception('Provided SBML file does not exists or is invalid!')

        # apply several model simplifications that make our life substantially easier
        if len(self.sbml_doc.getModel().getListOfFunctionDefinitions())>0:
            convertConfig = sbml.SBMLFunctionDefinitionConverter().getDefaultProperties()
            status = self.sbml_doc.convert(convertConfig)
            if status != sbml.LIBSBML_OPERATION_SUCCESS:
                raise Exception('Could not flatten function definitions!')

        convertConfig = sbml.SBMLInitialAssignmentConverter().getDefaultProperties()
        status = self.sbml_doc.convert(convertConfig)
        if status != sbml.LIBSBML_OPERATION_SUCCESS:
            raise Exception('Could not flatten initial assignments!')

        convertConfig = sbml.SBMLLocalParameterConverter().getDefaultProperties()
        status = self.sbml_doc.convert(convertConfig)
        if status != sbml.LIBSBML_OPERATION_SUCCESS:
            raise Exception('Could not flatten local parameters!')

        self.sbml = self.sbml_doc.getModel()

   

    def setPaths(self):
        """Deduce paths input and output paths from script file name."""
        self.modelPath = os.path.join(self.amiciPath,'models', self.modelName)
        self.modelSwigPath = os.path.join(self.modelPath, 'swig')

    def wrapModel(self, modelName):
        """Generate AMICI C++ files for the model provided to the constructor."""
        self.setName(modelName)
        self.processSBML()
        self.computeModelEquations()
        self.prepareModelFolder()
        self.generateCCode()
        self.compileCCode()
        
    def setName(self, modelName):
        """Sets the model name and adapts paths accordingly.

        Args:
        modelName: name of the model/model directory
        """
        self.modelName = modelName
        self.setPaths()
        if not os.path.exists(self.modelPath):
            os.makedirs(self.modelPath)

    def processSBML(self):
        """Read parameters, species, reactions, and so on from SBML model"""
        self.checkSupport()
        self.processParameters()
        self.processSpecies()
        self.processReactions()
        self.processCompartments()
        self.processRules()
        self.processVolumeConversion()
        self.processTime()
        self.cleanReservedSymbols()


    def checkSupport(self):
        """Check whether all required SBML features are supported."""
        if len(self.sbml.getListOfSpecies()) == 0:
            raise Exception('Models without species are currently not supported!')

        if len(self.sbml.getListOfEvents()) > 0:
            raise Exception('Events are currently not supported!')

        if any([not(rule.isAssignment()) for rule in self.sbml.getListOfRules()]):
            raise Exception('Algebraic and rate rules are currently not supported!')

        if any([reaction.isSetFast() for reaction in self.sbml.getListOfReactions()]):
            raise Exception('Fast reactions are currently not supported!')

        if any([any([not element.getStoichiometryMath() is None
                for element in list(reaction.getListOfReactants()) + list(reaction.getListOfProducts())])
                for reaction in self.sbml.getListOfReactions()]):
            raise Exception('Non-unity stoichiometry is currently not supported!')


    def processSpecies(self):
        """Get species information from SBML model."""
        species = self.sbml.getListOfSpecies()
        self.n_species = len(species)
        self.speciesIndex = {species_element.getId(): species_index
                             for species_index, species_element in enumerate(species)}
        self.symbols['species']['sym'] = sp.DenseMatrix([symbols(spec.getId()) for spec in species])
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


    def processParameters(self):
        """Get parameter information from SBML model."""
        parameters = self.sbml.getListOfParameters()
        self.symbols['parameter']['sym'] = sp.DenseMatrix([symbols(par.getId()) for par in parameters])
        self.parameterValues = [par.getValue() for par in parameters]
        self.n_parameters = len(self.symbols['parameter']['sym'])
        self.parameterIndex = {parameter_element.getId(): parameter_index
                             for parameter_index, parameter_element in enumerate(parameters)}

    def processCompartments(self):
        """Get compartment information, stoichiometric matrix and fluxes from SBML model."""
        compartments = self.sbml.getListOfCompartments()
        self.compartmentSymbols = sp.DenseMatrix([symbols(comp.getId()) for comp in compartments])
        self.compartmentVolume = sp.DenseMatrix([sp.sympify(comp.getVolume()) if comp.isSetVolume()
                                                 else sp.sympify(1.0) for comp in compartments])

        self.replaceInAllExpressions(self.compartmentSymbols,self.compartmentVolume)


    def processReactions(self):
        """Get reactions from SBML model."""
        reactions = self.sbml.getListOfReactions()
        self.n_reactions = len(reactions)

        # stoichiometric matrix
        self.stoichiometricMatrix = sp.zeros(self.n_species, self.n_reactions)
        self.fluxVector = sp.zeros(self.n_reactions, 1)
        self.symbols['flux']['sym'] = sp.zeros(self.n_reactions, 1)

        for reaction_index, reaction in enumerate(reactions):

            for type in ['reactant','product']:
                if type == 'reactant':
                    elementList = reaction.getListOfReactants()
                    sign = -1.0
                else:
                    elementList = reaction.getListOfProducts()
                    sign = +1.0

                elements = {}

                for e in elementList:
                    if e.element_name == 'species':
                        elements[e.getSpecies()] = e.getStoichiometry()
                    elif e.element_name == 'speciesReference':
                        if e.isSetStoichiometry():
                            elements[e.getSpecies()] = e.getStoichiometry()
                        else:
                            elements[e.getSpecies()] = e.getId()

                for element in elements.keys():
                    if not (element in self.constantSpecies or element in self.boundaryConditionSpecies):
                        self.stoichiometricMatrix[self.speciesIndex[element], reaction_index] += sign * \
                            sp.sympify(elements[element])*self.speciesConversionFactor[self.speciesIndex[element]]/ \
                            self.speciesCompartment[self.speciesIndex[element]]

            # usage of formulaToL3String ensures that we get "time" as time symbol
            self.fluxVector[reaction_index] = sp.sympify(sbml.formulaToL3String(reaction.getKineticLaw().getMath()))
            self.symbols['flux']['sym'][reaction_index] = sp.sympify('w' + str(reaction_index))


    def processRules(self):
        """Process *Rules defined in the SBML model."""
        rules = self.sbml.getListOfRules()

        rulevars = sp.DenseMatrix([sp.sympify(rule.getFormula()) for rule in rules]).free_symbols
        fluxvars = self.fluxVector.free_symbols
        specvars = self.symbols['species']['sym'].free_symbols
        volumevars = self.compartmentVolume.free_symbols
        compartmentvars = self.compartmentSymbols.free_symbols
        parametervars = self.symbols['parameter']['sym'].free_symbols
        stoichvars = self.stoichiometricMatrix.free_symbols

        observables = []
        observableSymbols = []
        self.observableNames = []
        for rule in rules:
            variable = sp.sympify(rule.getVariable())
            formula = sp.sympify(rule.getFormula())

            isObservable = True

            if variable in stoichvars:
                self.stoichiometricMatrix = self.stoichiometricMatrix.subs(variable, formula)
                isObservable = False

            if variable in specvars:
                raise Exception('Species assignment rules are currently not supported!')

            if variable in compartmentvars:
                raise Exception('Compartment assignment rules are currently not supported!')

            if variable in parametervars:
                try:
                    self.parameterValues[self.parameterIndex[str(variable)]] = float(formula)
                except:
                    raise Exception('Non-float parameter assignment rules are currently not supported!')
            if variable in fluxvars:
                self.fluxVector = self.fluxVector.subs(variable, formula)
                isObservable = False

            if variable in volumevars:
                self.compartmentVolume = self.compartmentVolume.subs(variable, formula)
                isObservable = False

            if variable in rulevars:
                for nested_rule in rules:
                    nested_formula = sp.sympify(nested_rule.getFormula())
                    nested_formula.subs(variable, formula)
                    nested_rule.setFormula(str(nested_formula))
                isObservable = False

            if isObservable:
                observables.append(formula)
                observableSymbols.append(sp.sympify('y' + str(len(observableSymbols))))
                self.observableNames.append(str(variable))

        if(len(observables)>0):
            self.observables = sp.DenseMatrix(observables)
            self.symbols['observable']['sym'] = sp.DenseMatrix(observableSymbols)
            self.n_observables = len(observables)
        else:
            self.observables = self.symbols['species']['sym']
            self.symbols['observable']['sym'] = sp.DenseMatrix([sp.sympify('y' + str(index))
                                                     for index in range(0,len(self.symbols['species']['sym']))])
            self.n_observables = len(self.symbols['species']['sym'])

    def processVolumeConversion(self):
        """Convert equations from amount to volume."""
        for index, bool in enumerate(self.speciesHasOnlySubstanceUnits):
            if bool:
                self.fluxVector = self.fluxVector.subs(self.symbols['species']['sym'][index],
                                                self.symbols['species']['sym'][index] * self.speciesCompartment[index]
                                                            .subs(self.compartmentSymbols,
                                                                  self.compartmentVolume))


    def processTime(self):
        """Converts time_symbol into cpp variable."""
        sbmlTimeSymbol = sp.sympify('time')
        amiciTimeSymbol = sp.sympify('t')

        self.replaceInAllExpressions(sbmlTimeSymbol, amiciTimeSymbol)



    def replaceInAllExpressions(self,old,new):
        """Replaces 'old' by 'new' in all symbolic expressions.

        Args:
        old: symbolic variables to be replaced
        new: replacement symbolic variables

        """
        self.stoichiometricMatrix = self.stoichiometricMatrix.subs(old,new)
        self.speciesInitial = self.speciesInitial.subs(old,new)
        self.fluxVector = self.fluxVector.subs(old,new)


    def cleanReservedSymbols(self):
        reservedSymbols = ['k','p']
        for str in reservedSymbols:
            old_symbol = sp.sympify(str)
            new_symbol = sp.sympify('amici_' + str)
            fields = ['observables','stoichiometricMatrix','compartmentSymbols','fluxVector']
            for field in fields:
                self.__setattr__(field,self.__getattribute__(field).subs(old_symbol,new_symbol))
            for symbol in self.symbols.keys():
                if 'sym' in self.symbols[symbol].keys():
                    self.symbols[symbol]['sym'] = self.symbols[symbol]['sym'].subs(old_symbol,new_symbol)


    def getSparseSymbols(self,symbolName):
        """Create sparse symbolic matrix.
        
        Args:
        symbolName:
        
        Returns:
        sparseMatrix:
        symbolList:
        sparseList:
        symbolColPtrs:
        symbolRowVals:
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

    def computeModelEquations(self):
        """Perform symbolic computations required to populate functions in `self.functions`."""
        
        # core
        self.functions['xdot']['sym'] = self.stoichiometricMatrix * self.symbols['flux']['sym']
        self.computeModelEquationsLinearSolver()
        self.computeModelEquationsObjectiveFunction()
        
        # sensitivity
        self.computeModelEquationsSensitivitesCore()
        self.computeModelEquationsForwardSensitivites()
        self.computeModelEquationsAdjointSensitivites()
        
        
    def computeModelEquationsLinearSolver(self):
        """Perform symbolic computations required for the use of various linear solvers."""
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


    def computeModelEquationsObjectiveFunction(self):
        """Perform symbolic computations required for objective function evaluation."""
        self.functions['y']['sym'] = self.observables

        self.symbols['sigma_y']['sym'] = sp.DenseMatrix([sp.sympify('sigma' + str(symbol))
                                                                for symbol in self.symbols['observable']['sym']])
        self.functions['sigma_y']['sym'] = sp.zeros(self.symbols['sigma_y']['sym'].cols,
                                                           self.symbols['sigma_y']['sym'].rows)
        self.symbols['my']['sym'] = sp.DenseMatrix([sp.sympify('m' + str(symbol))
                                                           for symbol in self.symbols['observable']['sym']])
        self.functions['Jy']['sym'] = sp.DenseMatrix([sp.sympify('0.5*sqrt(2*pi*sigma' + str(symbol) + '**2) ' +
                                    '+ 0.5*((' + str(symbol) + '-m' + str(symbol) + ')/sigma' + str(symbol) + ')**2')
                                    for iy, symbol in enumerate(self.symbols['observable']['sym'])])
        self.functions['dJydy']['sym'] = self.functions['Jy']['sym'].jacobian(self.symbols['observable']['sym'])
        self.functions['dJydsigma']['sym'] = self.functions['Jy']['sym'].jacobian(self.symbols['sigma_y']['sym'])
        self.functions['Jy']['sym'] = self.functions['Jy']['sym'].transpose()
        self.functions['dJydy']['sym'] = self.functions['dJydy']['sym'].transpose()
        self.functions['dJydsigma']['sym'] = self.functions['dJydsigma']['sym'].transpose()


    def computeModelEquationsSensitivitesCore(self):
        """Perform symbolic computations required for any sensitivity analysis."""
        self.functions['dydp']['sym'] = self.functions['y']['sym']\
                                                .jacobian(self.symbols['parameter']['sym'])
        self.functions['dydx']['sym'] = self.functions['y']['sym']\
                                                .jacobian(self.symbols['species']['sym'])
        
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
        """Perform symbolic computations required for forward sensitivity analysis."""
        self.symbols['sensitivity']['sym'] = getSymbols('sx',self.n_species)
        self.functions['sxdot']['sym'] = self.functions['JSparse']['sym']\
                                                * self.symbols['sensitivity']['sym'] \
                                                + self.symbols['dxdotdp']['sym']

    def computeModelEquationsAdjointSensitivites(self):
        """Perform symbolic computations required for adjoint sensitivity analysis."""
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

        for file in os.listdir(self.modelPath):
            file_path = os.path.join(self.modelPath, file)
            if os.path.isfile(file_path):
                os.remove(file_path)



    def generateCCode(self):
        """Create C++ code files for the model based on."""
        for name in self.symbols.keys():
            self.writeIndexFiles(name)

        for function in self.functions.keys():
            self.writeFunctionFile(function)

        self.writeWrapfunctionsCPP()
        self.writeWrapfunctionsHeader()
        self.writeCMakeFile()
        self.writeSwigFiles()

        shutil.copy(os.path.join(self.amiciSrcPath, 'main.template.cpp'),
                    os.path.join(self.modelPath, 'main.cpp'))


    def compileCCode(self):
        """Compile the generated model code"""
        subprocess.call([os.path.join(self.amiciPath, 'scripts', 'buildModel.sh'), self.modelName])



    def writeIndexFiles(self,name):
        """Write index file for a symbolic array.
        
        Args:
        Symbols: symbolic variables
        CVariableName: Name of the C++ array the symbols  
        fileName: filename without path
        """
        lines = []
        [lines.append('#define ' + str(symbol) + ' ' + str(self.symbols[name]['shortName']) + '[' + str(index) + ']')
            for index, symbol in enumerate(self.symbols[name]['sym'])]
        open(os.path.join(self.modelPath,name + '.h'), 'w').write('\n'.join(lines))


    def writeFunctionFile(self,function):
        """Write the function `function`.
        
        Args:
        function: name of the function to be written (see self.functions)"""
        
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
        open(os.path.join(self.modelPath, self.modelName + '_' + function + '.cpp'), 'w').write('\n'.join(lines))


    def getFunctionBody(self,function):
        """Generate C++ code for body of function `function`.
        
        Args:
        function: name of the function to be written (see self.functions)
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

        return lines


    def writeWrapfunctionsCPP(self):
        """Write model-specific 'wrapper' file (wrapfunctions.cpp)."""
        templateData = {'MODELNAME': self.modelName}
        applyTemplate(os.path.join(self.amiciSrcPath, 'wrapfunctions.template.cpp'),
                      os.path.join(self.modelPath, 'wrapfunctions.cpp'), templateData)


    def writeWrapfunctionsHeader(self):
        """Write model-specific header file (wrapfunctions.h)."""
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
                        'NK': '0',
                        'O2MODE': 'amici::AMICI_O2MODE_NONE',
                        'PARAMETERS': str(self.parameterValues)[1:-1]}
        applyTemplate(os.path.join(self.amiciSrcPath, 'wrapfunctions.ODE_template.h'),
                      os.path.join(self.modelPath, 'wrapfunctions.h'), templateData)

    def writeCMakeFile(self):
        """Write CMake CMakeLists.txt file for this model."""
        sources = [self.modelName + '_' + function + '.cpp ' if self.functions[function]['body'] is not None else ''
                    for function in self.functions.keys() ]
        try:
            sources.remove('')
        except:
            pass
        templateData = {'MODELNAME': self.modelName,
                        'SOURCES': '\n'.join(sources)}
        applyTemplate(os.path.join(self.amiciSrcPath, 'CMakeLists.template.txt'),
                      os.path.join(self.modelPath, 'CMakeLists.txt'), templateData)

    def writeSwigFiles(self):
        """Write SWIG interface files for this model."""
        if not os.path.exists(self.modelSwigPath):
            os.makedirs(self.modelSwigPath)
        templateData = {'MODELNAME': self.modelName}
        applyTemplate(os.path.join(self.amiciSwigPath, 'modelname.template.i'),
                      os.path.join(self.modelSwigPath, self.modelName + '.i'), templateData)
        shutil.copy(os.path.join(self.amiciSwigPath, 'CMakeLists_model.txt'),
                    os.path.join(self.modelSwigPath, 'CMakeLists.txt'))

    def getSymLines(self, symbols, variable, indentLevel):
        """Generate C++ code for assigning symbolic terms in symbols to C++ array `variable`.
        
        Args: 
        symbols: vectors of symbolic terms
        variable: name of the C++ array to assign to
        indentLevel: indentation level (number of leading blanks)
        
        Returns: 
        C++ code as list of lines"""

        lines = [' ' * indentLevel + variable + '[' + str(index) + '] = ' + self.printWithException(math) + ';'
                if not math == 0 else '' for index, math in enumerate(symbols)]

        try:
            lines.remove('')
        except:
            pass

        return lines

    def getSparseSymLines(self, symbolList, RowVals, ColPtrs, variable, indentLevel):
        """Generate C++ code for assigning sparse symbolic matrix to a C++ array `variable`.
     
        Args: 
        symbolList: vectors of symbolic terms
        RowVals: list of row indices of each nonzero entry (see CVODES SlsMat documentation for details)
        ColPtrs: list of indices of the first column entries (see CVODES SlsMat documentation for details)
        variable: name of the C++ array to assign to
        indentLevel: indentation level (number of leading blanks)
        
        Returns: 
        C++ code as list of lines"""
        lines = [
            ' ' * indentLevel + variable + '->indexvals[' + str(index) + '] = ' + self.printWithException(math) + ';'
            for index, math in enumerate(RowVals)]
        lines.extend([' ' * indentLevel + variable + '->indexptrs[' + str(index) + '] = ' + self.printWithException(math) + ';'
                      for index, math in enumerate(ColPtrs)])
        lines.extend([' ' * indentLevel + variable + '->data[' + str(index) + '] = ' + self.printWithException(math) + ';'
                      for index, math in enumerate(symbolList)])

        return lines

    def printWithException(self,math):
        try:
            return self.Codeprinter.doprint(math)
        except:
            raise Exception('Encountered unsupported function in expression "' + str(math) + '"!')

def applyTemplate(sourceFile,targetFile,templateData):
    """Load source file, apply template substitution as provided in templateData and save as targetFile.
    
    Args:
    sourceFile: relative or absolute path to template file
    targetFile: relative or absolute path to output file
    templateData: dictionary with template keywords to substitute (key is template veriable without TemplateAmici.delimiter)
    """
    with open(sourceFile) as filein:
        src = TemplateAmici(filein.read())
    result = src.safe_substitute(templateData)
    with open(targetFile, 'w') as fileout:
        fileout.write(result)

def getSymbols(prefix,length):
    """Get symbolic matrix with symbols prefix0..prefix(length-1).
    
    Args:
    prefix: variable name
    length: number of symbolic variables + 1
    """ 
    return sp.DenseMatrix([sp.sympify(prefix + str(i)) for i in range(0, length)])


def getSymbolicDiagonal(matrix):
    """Get symbolic matrix with diagonal of matrix `matrix`.
    
    Args:
    matrix: Matrix from which to return the diagonal
    Returns:
    Symbolic matrix with the diagonal of `matrix`.
    """
    assert matrix.cols == matrix.rows
    diagonal = [matrix[index,index] for index in range(matrix.cols)]
    
    return sp.DenseMatrix(diagonal)

class TemplateAmici(Template):
    """Template format used in AMICI (see string.template for more details)."""
    delimiter = 'TPL_'
