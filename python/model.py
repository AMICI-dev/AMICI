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
        amici_path:
        amici_swig_path:
        amici_src_path:
        model_path:
        model_swig_path:
        self.functionBodies:
    """

    # TODO: are units respected on sbml import? if not convert; at least throw if differ?
    # TODO: camelCase?
    # TODO: dupe:  ('JSparse', 'J'); to dict?
    """Long and short names for model components"""
    shortenedVariables = [('species','x'), ('sensitivity', 'sx'), ('vector','v'), ('vectorB','vB'), ('parameter','p'),
                          ('observable','y'), ('adjoint','xB'), ('flux', 'w'), ('dxdotdp','dxdotdp'),
                          ('dwdxSparse', 'dwdx'), ('dwdpSparse', 'dwdp'), ('JSparse', 'J'), ('JSparse', 'J'),
                          ('my', 'my'), ('sigmay', 'sigmay')]

    """Signatures and properties of generated model functions (see include/amici/model.h for details)."""
    functions = {
        'J': {'signature': '(realtype *J, const realtype t, const realtype *x, const double *p, const double *k, const realtype *h, const realtype *w, const realtype *dwdx)'},
        'JB': {'signature': '(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx)'},
        'JDiag': {'signature': '(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx)',
                  'symbol': 'J',
                  'diagonality': True},
        'JSparse': {'signature': '(SlsMat JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx)',
                    'symbol': 'JSparseList',
                    'sparsity': True},
        'JSparseB': {'signature': '(SlsMat JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx)',
                     'symbol': 'JSparseBList',
                     'sparsity': True},
        'Jv': {'signature': '(realtype *Jv, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *v, const realtype *w, const realtype *dwdx)'},
        'JvB': {'signature': '(realtype *JvB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *vB, const realtype *w, const realtype *dwdx)'},
        'Jy': {'signature': '(double *nllh, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my)',
               'variable': 'nllh',
               'multiobs': True},
        'dJydsigma': {'signature': '(double *dJydsigma, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my)',
                      'multiobs': True},
        'dJydy': {'signature': '(double *dJydy, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my)',
                  'multiobs': True},
        'dwdp': {'signature': '(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w)',
                 'symbol': 'dwdxSparseList'},
        'dwdx': {'signature': '(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w)',
                 'symbol': 'dwdpSparseList'},
        'dxdotdp': {'signature': '(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dwdp)',
                    'sensitivity': True},
        'dydx': {'signature': '(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h)'},
        'qBdot': {'signature': '(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdp)',
                  'sensitivity': True},
        'sigma_y': {'signature': '(double *sigmay, const realtype t, const realtype *p, const realtype *k)',
                    'symbol': 'sigmay',
                    'variable': 'sigmay'},
        'sxdot': {'signature': '(realtype *sxdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *sx, const realtype *w, const realtype *dwdx, const realtype *J, const realtype *dxdotdp)'},
        'w': {'signature': '(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h)'},
        'x0': {'signature': '(realtype *x0, const realtype t, const realtype *p, const realtype *k)'},
        'sx0': {'signature': '(realtype *sx0, const realtype t,const realtype *x0, const realtype *p, const realtype *k, const int ip)'},
        'xBdot': {'signature': '(realtype *xBdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx)'},
        'xdot': {'signature': '(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w)'},
        'y': {'signature': '(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h)'}
    }


    # TODO: move arguments to wrap model? 
    def __init__(self, SBMLFile, modelname):
        """Create a new Model instance.
        
        Arguments:
        ----------
        SBMLFile: Path to SBML file where the model is specified
        modelname: Name of Model to be generated
        """
        self.loadSBMLFile(SBMLFile)
        self.modelname = modelname
        
        self.setPaths()
        if not os.path.exists(self.model_path):
            os.makedirs(self.model_path)

        self.functionBodies = {} # TODO: "private" ?
        self.Codeprinter = CCodePrinter()

    def loadSBMLFile(self, SBMLFile):
        """Parse the provided SBML file."""
        self.SBMLreader = sbml.SBMLReader()
        self.sbml_doc = self.SBMLreader.readSBML(SBMLFile)

        if (self.sbml_doc is None):
            raise Exception('Provided SBML file does not exit or is invalid!')

        # apply several model simplifications that make our life substantially easier
        convertConfig = sbml.SBMLFunctionDefinitionConverter().getDefaultProperties()
        status = self.sbml_doc.convert(convertConfig)
        if status != sbml.LIBSBML_OPERATION_SUCCESS:
            raise Exception('Could not flatten function defintions!')

        convertConfig = sbml.SBMLInitialAssignmentConverter().getDefaultProperties()
        status = self.sbml_doc.convert(convertConfig)
        if status != sbml.LIBSBML_OPERATION_SUCCESS:
            raise Exception('Could not flatten initial assignments!')

        convertConfig = sbml.SBMLLocalParameterConverter().getDefaultProperties()
        status = self.sbml_doc.convert(convertConfig)
        if status != sbml.LIBSBML_OPERATION_SUCCESS:
            raise Exception('Could not flatten local parameters!')

        self.sbml = self.sbml_doc.getModel()
        if (self.sbml is None):
            raise Exception('Provided SBML file is invalid!')

        self.SBMLWriter = sbml.SBMLWriter()
        self.SBMLWriter.writeSBML(self.sbml_doc,SBMLFile.replace('.xml','_converted.xml'))

   

    def setPaths(self):
        """Deduce paths input and output paths from script file name."""
        # TODO: move amici_*_path to amici module?
        dirname, filename = os.path.split(os.path.abspath(__file__))
        self.amici_path = os.path.split(dirname)[0]
        self.amici_swig_path = os.path.join(self.amici_path, 'swig')
        self.amici_src_path = os.path.join(self.amici_path, 'src')
        self.model_path = os.path.join(self.amici_path,'models', self.modelname)
        self.model_swig_path = os.path.join(self.model_path, 'swig')

    def wrapModel(self):
        """Generate AMICI C++ files for the model provided to the constructor."""
        self.processSBML()
        self.checkModelDimensions()
        self.getModelEquations()
        self.generateCCode()
        self.compileCCode()


    def checkModelDimensions(self):
        """Ensure the model has species to simulate"""
        assert self.n_species > 0


    def processSBML(self):
        """Read parameters, species, reactions, and so on from SBML model"""
        self.checkSupport()
        self.processParameters()
        self.processSpecies()
        self.processReactions()
        self.processCompartments()
        self.processRules()
        self.processVolumeConversion()


    def checkSupport(self):
        if len(self.sbml.getListOfEvents()):
            raise Exception('Events are currently not supported.')


    def processSpecies(self):
        """Get species information from SBML model."""
        species = self.sbml.getListOfSpecies()
        self.n_species = len(species)
        self.speciesIndex = {species_element.getId(): species_index for species_index, species_element in enumerate(species)}
        self.speciesSymbols = sp.DenseMatrix([symbols(spec.getId()) for spec in species])
        self.speciesCompartment = sp.DenseMatrix([symbols(spec.getCompartment()) for spec in species])
        self.constantSpecies = [species_element.getId() if species_element.getConstant() else None
                                for species_element in species]
        self.boundaryConditionSpecies = [species_element.getId() if species_element.getBoundaryCondition() else None
                                         for species_element in species]
        concentrations = [spec.getInitialConcentration() for spec in species]
        amounts = [spec.getInitialAmount() for spec in species]

        self.speciesInitial = sp.DenseMatrix([sp.sympify(conc) if not math.isnan(conc) else
                                              sp.sympify(amounts[index])/self.speciesCompartment[index] if not
                                              math.isnan(amounts[index]) else
                                              self.speciesSymbols[index]
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
        self.parameterSymbols = sp.DenseMatrix([symbols(par.getId()) for par in parameters])
        self.parameterValues = [par.getValue() for par in parameters]
        self.n_parameters = len(self.parameterSymbols)

    def processCompartments(self):
        """Get compartment information, stoichiometric matrix and fluxes from SBML model."""
        compartments = self.sbml.getListOfCompartments()
        self.compartmentSymbols = sp.DenseMatrix([symbols(comp.getId()) for comp in compartments])
        self.compartmentVolume = sp.DenseMatrix([sp.sympify(comp.getVolume()) for comp in compartments])

        #TODO: ?? overwritten in processReactions
        self.stoichiometricMatrix = self.stoichiometricMatrix.subs(self.compartmentSymbols, self.compartmentVolume)
        self.speciesInitial = self.speciesInitial.subs(self.compartmentSymbols, self.compartmentVolume)
        self.fluxVector = self.fluxVector.subs(self.compartmentSymbols, self.compartmentVolume)


    def processReactions(self):
        """Get reactions from SBML model."""
        reactions = self.sbml.getListOfReactions()
        self.n_reactions = len(reactions)

        # stoichiometric matrix
        self.stoichiometricMatrix = sp.zeros(self.n_species, self.n_reactions)
        self.fluxVector = sp.zeros(self.n_reactions, 1)
        self.fluxSymbols = sp.zeros(self.n_reactions, 1)

        for reaction_index, reaction in enumerate(reactions):
            reactants = {r.getSpecies(): r.getStoichiometry() for r in reaction.getListOfReactants()}
            products = {p.getSpecies(): p.getStoichiometry() for p in reaction.getListOfProducts()}

            for reactant in reactants.keys():
                if not (reactant in self.constantSpecies or reactant in self.boundaryConditionSpecies):
                    self.stoichiometricMatrix[self.speciesIndex[reactant], reaction_index] -= \
                        sp.sympify(reactants[reactant])*self.speciesConversionFactor[self.speciesIndex[reactant]]/ \
                        self.speciesCompartment[self.speciesIndex[reactant]]


            for product in products.keys():
                if not (product in self.constantSpecies or product in self.boundaryConditionSpecies):
                    self.stoichiometricMatrix[self.speciesIndex[product], reaction_index] += \
                        sp.sympify(products[product])*self.speciesConversionFactor[self.speciesIndex[product]]/ \
                        self.speciesCompartment[self.speciesIndex[product]]

            self.fluxVector[reaction_index] = sp.sympify(reaction.getKineticLaw().getFormula())
            self.fluxSymbols[reaction_index] = sp.sympify('w' + str(reaction_index))


    def processRules(self):
        """Process *Rules defined in the SBML model."""
        rules = self.sbml.getListOfRules()

        rulevars = sp.DenseMatrix([sp.sympify(rule.getFormula()) for rule in rules]).free_symbols
        fluxvars = self.fluxVector.free_symbols
        initvars = self.speciesInitial.free_symbols
        volumevars = self.compartmentVolume.free_symbols
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

            if variable in initvars:
                self.speciesInitial = self.speciesInitial.subs(variable, formula)
                isObservable = False

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
            self.observableSymbols = sp.DenseMatrix(observableSymbols)
            self.n_observables = len(observables)
        else:
            self.observables = self.speciesSymbols
            self.observableSymbols = sp.DenseMatrix([sp.sympify('y' + str(index))
                                                     for index in range(0,len(self.speciesSymbols))])
            self.n_observables = len(self.speciesSymbols)

    def processVolumeConversion(self):
        if __name__ == '__main__':
            self.fluxVector = self.fluxVector.subs(self.speciesSymbols,
                                                   self.speciesSymbols.mul_matrix(
                                                       self.speciesCompartment.applyfunc(lambda x: 1/x).
                                                           subs(self.compartmentSymbols,self.compartmentVolume)))


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
        matrix = self.__getattribute__(symbolName)
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


    # TODO: misleading name, does not return
    # TODO: write results to self.functions instead of to instance? saves the need for getattr
    def getModelEquations(self):
        """Perform symbolic computations required to populate functions in `self.functions`."""
        # core
        self.xdot = self.stoichiometricMatrix * self.fluxSymbols

        self.w = self.fluxVector

        self.dwdx = self.fluxVector.jacobian(self.speciesSymbols)
        self.dwdxSparse, self.dwdxSparseSymbols, self.dwdxSparseList  = self.getSparseSymbols('dwdx')[0:3]

        self.J = self.xdot.jacobian(self.speciesSymbols) + self.stoichiometricMatrix * self.dwdxSparse
        self.vectorSymbols = getSymbols('v',self.n_species)
        self.Jv = self.J*self.vectorSymbols
        self.JSparse, self.JSparseSymbols, self.JSparseList, self.JSparseColPtrs, self.JSparseRowVals \
            = self.getSparseSymbols('J')

        self.x0 = self.speciesInitial

        # sensitivity
        self.dwdp = self.fluxVector.jacobian(self.parameterSymbols)
        self.dwdpSparse, self.dwdpSparseSymbols, self.dwdpSparseList = self.getSparseSymbols('dwdp')[0:3]

        self.dxdotdp = self.xdot.jacobian(self.parameterSymbols) + self.stoichiometricMatrix * self.dwdpSparse
        self.dxdotdpSymbols = getSymbols('dxdotdp',self.n_species)
        self.sx0 = self.speciesInitial.jacobian(self.parameterSymbols)

        # forward
        self.sensitivitySymbols = getSymbols('sx',self.n_species)
        self.sxdot = self.JSparse*self.sensitivitySymbols + self.dxdotdpSymbols

        # adjoint
        self.JB = self.J.transpose()
        self.vectorBSymbols = getSymbols('vB', self.n_species)
        self.JvB = self.JB * self.vectorBSymbols
        self.JSparseB, self.JSparseBSymbols, self.JSparseBList, self.JSparseBColPtrs, self.JSparseBRowVals \
            = self.getSparseSymbols('JB')

        self.adjointSymbols = getSymbols('xB',self.n_species)
        self.xBdot = - self.JB*self.adjointSymbols
        self.qBdot = - self.adjointSymbols.transpose()*self.dxdotdp

        # observables
        self.y = self.observables
        self.dydp = self.y.jacobian(self.parameterSymbols)
        self.dydx = self.y.jacobian(self.speciesSymbols)

        # objective function
        self.sigmaySymbols = sp.DenseMatrix([sp.sympify('sigma' + str(symbol)) for symbol in self.observableSymbols])
        self.sigmay = sp.zeros(self.sigmaySymbols.cols, self.sigmaySymbols.rows)
        self.mySymbols = sp.DenseMatrix([sp.sympify('m' + str(symbol)) for symbol in self.observableSymbols])
        self.Jy = sp.DenseMatrix([sp.sympify('0.5*sqrt(2*pi*sigma' + str(symbol) + '**2) ' +
                                        '+ 0.5*((' + str(symbol) + '-m' + str(symbol) + ')/sigma' + str(symbol) + ')**2')
                             for iy, symbol in enumerate(self.observableSymbols)])
        self.dJydy = self.Jy.jacobian(self.observableSymbols)
        self.dJydsigma = self.Jy.jacobian(self.sigmaySymbols)
        self.Jy = self.Jy.transpose()
        self.dJydy = self.dJydy.transpose()
        self.dJydsigma = self.dJydsigma.transpose()


    def generateCCode(self):
        """Create C++ code files for the model based on."""
        for field,name in self.shortenedVariables:
            self.writeIndexFiles(self.__getattribute__(field + 'Symbols'), name, field + '.h')

        for function in self.functions.keys():
            self.writeFunctionFile(function)

        self.writeWrapfunctionsCPP()
        self.writeWrapfunctionsHeader()
        self.writeCMakeFile()
        self.writeSwigFiles()

        shutil.copy(os.path.join(self.amici_src_path, 'main.template.cpp'),
                    os.path.join(self.model_path, 'main.cpp'))


    def compileCCode(self):
        """Compile the generated model code"""
        subprocess.call([os.path.join(self.amici_path, 'scripts', 'buildModel.sh'), self.modelname])



    def writeIndexFiles(self,Symbols,CVariableName,fileName):
        """Write index file for a symbolic array.
        
        Args:
        Symbols: symbolic variables
        CVariableName: Name of the C++ array the symbols  
        fileName: filename without path
        """
        lines = []
        [lines.append('#define ' + str(symbol) + ' ' + CVariableName + '[' + str(index) + ']')
         for index, symbol in enumerate(Symbols)]
        open(os.path.join(self.model_path,fileName), 'w').write('\n'.join(lines))


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

        for variable, shortenedVariable in self.shortenedVariables:
            if(not re.search('const (realtype|double) \*' + shortenedVariable + '[,)]+',signature) is None):
                lines.append('#include "' + variable + '.h"')

        lines.append('')

        lines.append('void ' + function + '_' + self.modelname + signature + '{')
        
        # function body
        body = self.getFunctionBody(function)
        self.functionBodies[function] = body
        lines += body
        lines.append('}')
        #if not body is None:
        open(os.path.join(self.model_path, self.modelname + '_' + function + '.cpp'), 'w').write('\n'.join(lines))


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
            symbol = self.__getattribute__(self.functions[function]['symbol'])
        else:
            symbol = self.__getattribute__(function)
        lines = []

        if ('diagonality' in self.functions[function].keys()):
            symbol = getSymbolicDiagonal(symbol)

        if('sensitivity' in self.functions[function].keys()):
            lines.append(' '*4 + 'switch(ip) {')
            for ipar in range(0,self.n_parameters):
                lines.append(' ' * 8 + 'case ' + str(ipar) + ':')
                lines += self.writeSym(symbol[:,ipar], variableName, 12)
                lines.append(' ' * 12 + 'break;')
            lines.append('}')
        elif('multiobs' in self.functions[function].keys()):
            lines.append(' '*4 + 'switch(iy) {')
            for iobs in range(0,self.n_observables):
                lines.append(' ' * 8 + 'case ' + str(iobs) + ':')
                lines += self.writeSym(symbol[:,iobs], variableName, 12)
                lines.append(' ' * 12 + 'break;')
            lines.append('}')
        else:
            if('sparsity' in self.functions[function].keys()):
                rowVals = self.__getattribute__(function + 'RowVals')
                colPtrs = self.__getattribute__(function + 'ColPtrs')
                lines += self.writeSparseSym(symbol, rowVals, colPtrs, variableName, 4)
            else:
                lines += self.writeSym(symbol, variableName, 4)

        return lines


    def writeWrapfunctionsCPP(self):
        """Write model-specific 'wrapper' file (wrapfunctions.cpp)."""
        templateData = {'MODELNAME': self.modelname}
        applyTemplate(os.path.join(self.amici_src_path, 'wrapfunctions.template.cpp'),
                      os.path.join(self.model_path, 'wrapfunctions.cpp'), templateData)


    def writeWrapfunctionsHeader(self):
        """Write model-specific header file (wrapfunctions.h)."""
        templateData = {'MODELNAME': str(self.modelname),
                        'NX': str(self.n_species),
                        'NXTRUE': str(self.n_species),
                        'NY': str(self.n_observables),
                        'NYTRUE': str(self.n_observables),
                        'NZ': '0',
                        'NZTRUE': '0',
                        'NEVENT': '0',
                        'NOBJECTIVE': '1',
                        'NW': str(len(self.fluxSymbols)),
                        'NDWDDX': str(len(self.dwdxSparseSymbols)),
                        'NDWDP': str(len(self.dwdpSparseSymbols)),
                        'NNZ': str(len(self.JSparseSymbols)),
                        'UBW': str(self.n_species),
                        'LBW': str(self.n_species),
                        'NP': str(self.n_parameters),
                        'NK': '0',
                        'O2MODE': 'amici::AMICI_O2MODE_NONE',
                        'PARAMETERS': str(self.parameterValues)[1:-1]}
        applyTemplate(os.path.join(self.amici_src_path, 'wrapfunctions.ODE_template.h'),
                      os.path.join(self.model_path, 'wrapfunctions.h'), templateData)

    def writeCMakeFile(self):
        """Write CMake CMakeLists.txt file for this model."""
        sources = [ self.modelname + '_' + function + '.cpp ' if self.functionBodies[function] is not None else ''
                    for function in self.functionBodies.keys() ]
        try:
            sources.remove('')
        except:
            pass
        templateData = {'MODELNAME': self.modelname,
                        'SOURCES': '\n'.join(sources)}
        applyTemplate(os.path.join(self.amici_src_path, 'CMakeLists.template.txt'),
                      os.path.join(self.model_path, 'CMakeLists.txt'), templateData)

    def writeSwigFiles(self):
        """Write SWIG interface files for this model."""
        if not os.path.exists(self.model_swig_path):
            os.makedirs(self.model_swig_path)
        templateData = {'MODELNAME': self.modelname}
        applyTemplate(os.path.join(self.amici_swig_path, 'modelname.template.i'),
                      os.path.join(self.model_swig_path, self.modelname + '.i'), templateData)
        shutil.copy(os.path.join(self.amici_swig_path, 'CMakeLists_model.txt'),
                    os.path.join(self.model_swig_path, 'CMakeLists.txt'))

    def writeSym(self,symbols,variable,indentLevel):
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
    
    # TODO: rename, does not *write*
    def writeSparseSym(self, symbolList, RowVals, ColPtrs, variable, indentLevel):
        """Generate C++ code for assigning sparse symbolic matrix to a C++ array `variable`.
     
        Args: 
        symbolList: vectors of symbolic terms
        RowVals:
        ColPtrs:
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
            raise Exception('Encountered unsupported function in expression "' + str(math) + '"')

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
