""" @package amici.ode_export The python ode export module for python
"""
#!/usr/bin/env python3


import symengine as sp

import re
import shutil
import subprocess
import sys
import os
import copy


from symengine.printing import CCodePrinter
from symengine import symbols
from string import Template

from . import amiciSwigPath, amiciSrcPath, amiciModulePath


functions = {
    'J': {
        'signature':
            '(realtype *J, const realtype t, const realtype *x,'
            ' const double *p, const double *k, const realtype *h,'
            ' const realtype *w, const realtype *dwdx)',
        'assume_pow_positivity':
            True,
    },
    'JB': {
        'signature':
            '(realtype *JB, const realtype t, const realtype *x,'
            ' const realtype *p, const realtype *k, const realtype *h,'
            ' const realtype *xB, const realtype *w,'
            ' const realtype *dwdx)',
        'assume_pow_positivity':
            True,
    },
    'JDiag': {
        'signature':
            '(realtype *JDiag, const realtype t, const realtype *x,'
            ' const realtype *p, const realtype *k, const realtype *h,'
            ' const realtype *w, const realtype *dwdx)',
        'assume_pow_positivity':
            True,
    },
    'JSparse': {
        'signature':
            '(SlsMat JSparse, const realtype t, const realtype *x,'
            ' const realtype *p, const realtype *k, const realtype *h,'
            ' const realtype *w, const realtype *dwdx)',
        'sparse':
            True,
        'assume_pow_positivity':
            True,
    },
    'JSparseB': {
        'signature':
            '(SlsMat JSparseB, const realtype t, const realtype *x,'
            ' const realtype *p, const realtype *k, const realtype *h,'
            ' const realtype *xB, const realtype *w,'
            ' const realtype *dwdx)',
        'sparse':
            True,
        'assume_pow_positivity':
            True,
    },
    'Jv': {
        'signature':
            '(realtype *Jv, const realtype t, const realtype *x,'
            ' const realtype *p, const realtype *k, const realtype *h,'
            ' const realtype *v, const realtype *w,'
            ' const realtype *dwdx)',
        'assume_pow_positivity':
            True,
    },
    'JvB': {
        'signature':
            '(realtype *JvB, const realtype t, const realtype *x,'
            ' const realtype *p, const realtype *k, const realtype *h,'
            ' const realtype *xB, const realtype *vB,'
            ' const realtype *w, const realtype *dwdx)',
        'assume_pow_positivity':
            True,
    },
    'Jy': {
        'signature':
            '(double *nllh, const int iy, const realtype *p,'
            ' const realtype *k, const double *y,'
            ' const double *sigmay, const double *my)',
        'variable':
            'nllh',
        'multiobs':
            True,
    },
    'dJydsigmay': {
        'signature':
            '(double *dJydsigma, const int iy, const realtype *p,'
            ' const realtype *k, const double *y,'
            ' const double *sigmay, const double *my)',
        'variable':
            'dJydsigma',
        'multiobs':
            True,
    },
    'dJydy': {
        'signature':
            '(double *dJydy, const int iy, const realtype *p,'
            ' const realtype *k, const double *y,'
            ' const double *sigmay, const double *my)',
        'multiobs':
            True,
    },
    'dwdp': {
        'signature':
            '(realtype *dwdp, const realtype t, const realtype *x,'
            ' const realtype *p, const realtype *k, const realtype *h,'
            ' const realtype *w)',
        'sparse':
            True,
        'assume_pow_positivity':
            True,
    },
    'dwdx': {
        'signature':
            '(realtype *dwdx, const realtype t, const realtype *x,'
            ' const realtype *p, const realtype *k, const realtype *h,'
            ' const realtype *w)',
        'sparse':
            True,
        'assume_pow_positivity':
            True,
    },
    'dxdotdp': {
        'signature':
            '(realtype *dxdotdp, const realtype t, const realtype *x,'
            ' const realtype *p, const realtype *k, const realtype *h,'
            ' const int ip, const realtype *w, const realtype *dwdp)',
        'sensitivity':
            True,
        'assume_pow_positivity':
            True,
    },
    'dydx': {
        'signature':
            '(double *dydx, const realtype t, const realtype *x,'
            ' const realtype *p, const realtype *k,'
            ' const realtype *h)',
    },
    'dydp': {
        'signature':
            '(double *dydp, const realtype t, const realtype *x,'
            ' const realtype *p, const realtype *k, const realtype *h,'
            ' const int ip)',
        'sensitivity':
            True,
    },
    'dsigmaydp': {
        'signature':
            '(double *dsigmaydp, const realtype t, const realtype *p,'
            ' const realtype *k, const int ip)',
        'sensitivity':
            True,
    },
    'qBdot': {
        'signature':
            '(realtype *qBdot, const int ip, const realtype t,'
            ' const realtype *x, const realtype *p, const realtype *k,'
            ' const realtype *h, const realtype *xB,'
            ' const realtype *w, const realtype *dwdp)',
        'sensitivity':
            True,
        'assume_pow_positivity':
            True,
    },
    'sigmay': {
        'signature':
            '(double *sigmay, const realtype t, const realtype *p,'
            ' const realtype *k)',
    },
    'sxdot': {
        'signature':
            '(realtype *sxdot, const realtype t, const realtype *x,'
            ' const realtype *p, const realtype *k, const realtype *h,'
            ' const int ip, const realtype *sx, const realtype *w,'
            ' const realtype *dwdx, const realtype *J,'
            ' const realtype *dxdotdp)',
        'assume_pow_positivity':
            True,
    },
    'w': {
        'signature':
            '(realtype *w, const realtype t, const realtype *x,'
            ' const realtype *p, const realtype *k,'
            ' const realtype *h)',
        'assume_pow_positivity':
            True,
    },
    'x0': {
        'signature':
            '(realtype *x0, const realtype t, const realtype *p,'
            ' const realtype *k)',
    },
    'x0_fixedParameters': {
        'signature':
            '(realtype *x0, const realtype t, const realtype *p,'
            ' const realtype *k)',
        'variable': 'x0',
    },
    'sx0': {
        'signature':
            '(realtype *sx0, const realtype t,const realtype *x,'
            ' const realtype *p, const realtype *k, const int ip)',
        'sensitivity':
            True,
    },
    'sx0_fixedParameters': {
        'signature':
            '(realtype *sx0, const realtype t,const realtype *x0,'
            ' const realtype *p, const realtype *k, const int ip)',
        'variable':
            'sx0',
        'sensitivity':
            True,
    },
    'xBdot': {
        'signature':
            '(realtype *xBdot, const realtype t, const realtype *x,'
            ' const realtype *p, const realtype *k, const realtype *h,'
            ' const realtype *xB, const realtype *w,'
            ' const realtype *dwdx)',
        'assume_pow_positivity':
            True,
    },
    'xdot': {
        'signature':
            '(realtype *xdot, const realtype t, const realtype *x,'
            ' const realtype *p, const realtype *k, const realtype *h,'
            ' const realtype *w)',
        'assume_pow_positivity':
            True,
    },
    'y': {
        'signature':
            '(double *y, const realtype t, const realtype *x,'
            ' const realtype *p, const realtype *k,'
            ' const realtype *h)',
    }
}

sparse_functions = [
    function for function in functions
    if 'sparse' in functions[function]
       and functions[function]['sparse']
]

sensi_functions = [
    function for function in functions
    if 'sensitivity' in functions[function]
       and functions[function]['sensitivity']
]

class ModelQuantity:
    def __init__(self, identifier, name, value):
        self.identifier = identifier
        self.name = name
        self.value = value

    def __repr__(self):
        return str(self.identifier)


class State(ModelQuantity):
    def __init__(self, identifier, name, value, dt):
        super(State, self).__init__(identifier, name, value)
        self.dt = dt


class Observable(ModelQuantity):
    pass


class SigmaY(ModelQuantity):
    pass


class Expression(ModelQuantity):
    pass


class Parameter(ModelQuantity):
    pass


class Constant(ModelQuantity):
    pass


class LogLikelihood(ModelQuantity):
    pass

symbol_to_type = {
    'species': State,
    'parameter': Parameter,
    'fixed_parameter': Constant,
    'observable': Observable,
    'sigmay': SigmaY,
    'llhy': LogLikelihood,
    'expression': Expression,
}


class ODEModel:
    def __init__(self, flux_as_expressions=True):
        self._states = []
        self._observables = []
        self._sigmays = []
        self._parameters = []
        self._constants = []
        self._loglikelihoods = []
        self._expressions = []
        self.imported_from = ''
        self.sdim_funs = {
            'sx': self.nx,
            'v': self.nx,
            'vB': self.nx,
            'xB': self.nx,
            'sigmay': self.ny,
        }
        self.flux_as_expressions = flux_as_expressions
        self._eqs = dict()
        self._sparseeqs = dict()
        self._vals = dict()
        self._syms = dict()
        self._names = dict()
        self._sparsesyms = dict()
        self._colptrs = dict()
        self._rowvals = dict()

        self.equation_prototype = {
            'x0': '_states',
            'y': '_observables',
            'Jy': '_loglikelihoods',
            'w': '_expressions',
            'sigmay': '_sigmays',
        }
        self.variable_prototype = {
            'x': '_states',
            'y': '_observables',
            'p': '_parameters',
            'k': '_constants',
            'sigmay': '_sigmays'
        }
        self.value_prototype = {
            'p': '_parameters',
            'k': '_constants',
        }
        self.total_derivative_prototypes = {
            'J': {
                'eq': 'xdot',
                'chainvar': 'w',
                'var': 'x',
            },
            'sxdot': {
                'eq': 'xdot',
                'chainvar': 'x',
                'var': 'p',
                'dydx_name': 'JSparse',
                'dxdz_name': 'sx',
            },
        }
        self.multiplication_prototypes = {
            'Jv': {
                'x': 'J',
                'y': 'v',
            },
            'JvB': {
                'x': 'JB',
                'y': 'vB',
            },
            'qBdot': {
                'x': 'xB',
                'transpose_x': True,
                'y': 'dxdotdp',
            },
            'xBdot': {
                'x': 'JB',
                'y': 'xB',
            },
        }

    def import_from_sbml_importer(self, si, flux_as_expressions=True):
        self.imported_from = 'amici_sbml'

        if flux_as_expressions:
            self._eqs['dxdotdw'] = si.stoichiometricMatrix
            self._eqs['w'] = si.fluxVector

        self.symbols = copy.copy(si.symbols)
        self.symbols['species']['dt'] = si.stoichiometricMatrix * self.sym('w')

        for symbol in [s for s in self.symbols if s != 'my']:
            # transform dict of lists into a list of dicts
            protos = [dict(zip(self.symbols[symbol], t))
                      for t in zip(*self.symbols[symbol].values())]
            for proto in protos:
                self.add_component(symbol_to_type[symbol](**proto))

        self.generateBasicVariables()

    def add_component(self, component):
        for comp_type in [Observable, Expression, Parameter, Constant, State,
                          LogLikelihood, SigmaY]:
            if isinstance(component, comp_type):
                getattr(self, f'_{type(component).__name__.lower()}s').append(
                    component
                )
                return
        Exception(f'Invalid component type {type(component)}')

    def nx(self):
        return len(self.sym('x'))

    def ny(self):
        return len(self.sym('y'))

    def nk(self):
        return len(self.sym('k'))

    def np(self):
        return len(self.sym('p'))

    def sym(self, name):
        if name not in self._syms:
            self.generateSymbol(name)
        return self._syms[name]

    def sparsesym(self, name):
        if name not in sparse_functions:
            raise Exception(f'{name} is not marked as sparse')
        if name not in self._sparsesyms:
            self.generateSparseSymbol(name)
        return self._sparsesyms[name]

    def eq(self, name):
        if name not in self._eqs:
            self.computeEquation(name)
        return self._eqs[name]

    def sparseeq(self, name):
        if name not in sparse_functions:
            raise Exception(f'{name} is not marked as sparse')
        if name not in self._sparseeqs:
            self.generateSparseSymbol(name)
        return self._sparseeqs[name]

    def colptr(self, name):
        if name not in sparse_functions:
            raise Exception(f'{name} is not marked as sparse')
        if name not in self._sparseeqs:
            self.generateSparseSymbol(name)
        return self._colptrs[name]

    def rowval(self, name):
        if name not in sparse_functions:
            raise Exception(f'{name} is not marked as sparse')
        if name not in self._sparseeqs:
            self.generateSparseSymbol(name)
        return self._rowvals[name]

    def val(self, name):
        if name not in self._vals:
            self.generateValue(name)
        return self._vals[name]

    def name(self, name):
        if name not in self._names:
            self.generateName(name)
        return self._names[name]

    def generateSymbol(self, name):
        if name in self.variable_prototype:
            component = self.variable_prototype[name]
            self._syms[name] = sp.DenseMatrix(
                [comp.identifier for comp in getattr(self, component)]
            )
            if name == 'y':
                self._syms['my'] = sp.DenseMatrix(
                    [sp.Symbol(f'm{comp.identifier}')
                    for comp in getattr(self, component)]
                )
            return
        elif name in sparse_functions:
            self.generateSparseSymbol(name)
            return
        elif name in self.sdim_funs:
            length = self.sdim_funs[name]()
        elif name in sensi_functions:
            length = self.eq(name).shape[0]
        else:
            length = len(self.eq(name))

        self._syms[name] = sp.DenseMatrix([
            sp.Symbol(f'{name}{i}') for i in range(length)
        ])

    def generateBasicVariables(self):
        for var in self.variable_prototype:
            self.generateSymbol(var)

    def generateSparseSymbol(self, name):
        """
        Create sparse symbolic matrix.

        sparseMatrix: sparse matrix containing symbolic entries
        symbolList: symbolic vector containing the list of symbol names
        sparseList: symbolic vector containing the list of symbol formulas
        symbolColPtrs: Column pointer as specified in the SlsMat definition in CVODES
        symbolRowVals: Row Values as specified in the SlsMat definition in CVODES

        Arguments:
        name: name of the equation

        Returns:

        Raises:
        """
        matrix = self.eq(name)
        idx = 0
        sparseMatrix = sp.zeros(matrix.rows, matrix.cols)
        symbolList = []
        sparseList = []
        symbolColPtrs = []
        symbolRowVals = []
        for col in range(0, matrix.cols):
            symbolColPtrs.append(idx)
            for row in range(0, matrix.rows):
                if not (matrix[row, col] == 0):
                    symbolName = f'{name}{idx}'
                    sparseMatrix[row, col] = sp.Symbol(symbolName)
                    symbolList.append(symbolName)
                    sparseList.append(matrix[row, col])
                    symbolRowVals.append(row)
                    idx += 1

        symbolColPtrs.append(idx)
        sparseList = sp.DenseMatrix(sparseList)

        self._colptrs[name] = symbolColPtrs
        self._rowvals[name] = symbolRowVals
        self._sparseeqs[name] = sparseList
        self._sparsesyms[name] = symbolList
        self._syms[name] = sparseMatrix

    def computeEquation(self, name):
        match_deriv = re.match(r'd([\w]+)d([a-z]+)', name)

        if match_deriv:
            self.partialDerivative(match_deriv.group(1), match_deriv.group(2))

        elif name in self.equation_prototype:
            self.equationFromComponent(name, self.equation_prototype[name])

        elif name in self.total_derivative_prototypes:
            args = self.total_derivative_prototypes[name]
            args['name'] = name
            self.totalDerivative(**args)

        elif name in self.multiplication_prototypes:
            args = self.multiplication_prototypes[name]
            args['name'] = name
            self.multiplication(**args)

        elif name == 'xdot':
            self._eqs['xdot'] = sp.DenseMatrix(
                [comp.dt for comp in self._states]
            )

        elif name in ['sx0', 'sx0_fixedParameters']:
            self.partialDerivative(name[1:], 'p', name=name)

        elif name == 'JB':
            self._eqs[name] = self.eq('J').transpose()

        elif name == 'JDiag':
            self._eqs[name] = getSymbolicDiagonal(self.eq('J'))

        elif name == 'x0_fixedParameters':
            k = self.sym('k')
            self._eqs[name] = sp.DenseMatrix([
                eq
                if any([sym in eq.free_symbols
                        for sym in k])
                else 0.0
                for eq in self.eq('x0')
            ])

        elif name in ['JSparse', 'JSparseB']:
            self._eqs[name] = self.eq(name.replace('Sparse', ''))

        else:
            raise Exception(f'Unknown equation {name}')

        if name in ['Jy', 'dydx']:
            self._eqs[name] = self._eqs[name].transpose()

    def symNames(self):
        return list(self._syms.keys())

    def partialDerivative(self, eq, var, name=None):
        if not name:
            name = f'd{eq}d{var}'
        self._eqs[name] = self.eq(eq).jacobian(self.sym(var))

    def multiplication(self, name, x, y,
                       transpose_x=False, sign=1):
        vars = dict()
        for varname in [x, y]:
            if var_in_function_signature(name, varname):
                vars[varname] = self.sym(varname)
            else:
                vars[varname] = self.eq(varname)

        if transpose_x:
            xx = vars[x].transpose()
        else:
            xx = vars[x]

        self._eqs[name] = sign * xx * vars[y]

    def equationFromComponent(self, name, component):
        self._eqs[name] = sp.DenseMatrix(
            [comp.value for comp in getattr(self, component)]
        )

    def generateValue(self, name):
        if name in self.value_prototype:
            component = self.value_prototype[name]
        else:
            raise Exception(f'No values for {name}')

        self._vals[name] = [comp.value for comp in getattr(self, component)]

    def generateName(self, name):
        if name in self.variable_prototype:
            component = self.variable_prototype[name]
        elif name in self.equation_prototype:
            component = self.equation_prototype[name]
        else:
            raise Exception(f'No names for {name}')

        self._names[name] = [comp.name for comp in getattr(self, component)]

    def totalDerivative(self, name, eq, chainvar, var,
                        dydx_name=None, dxdz_name=None):

        # compute total derivative according to chainrule
        # Dydz = dydx*dxdz + dydz
        vars = dict()
        vars['dydx'] = dict()
        vars['dxdz'] = dict()
        vars['dydz'] = dict()

        vars['dydx']['name'] = f'd{eq}d{chainvar}' \
            if dydx_name is None else dydx_name
        vars['dxdz']['name'] = f'd{chainvar}d{var}' \
            if dxdz_name is None else dxdz_name
        vars['dydz']['name'] = f'd{eq}d{var}'

        # if the symbol appears in the function signature,
        # we can use the respective symbol instead of the equation
        for var in vars:
            varname = vars[var]["name"]
            if var_in_function_signature(name, varname):
                vars[var]['sym'] = self.sym(varname)
            else:
                vars[var]['sym'] = self.eq(varname)

        self._eqs[name] = \
            vars['dydx']['sym'] * vars['dxdz']['sym'] + vars['dydz']['sym']

def var_in_function_signature(name, varname):
    varname = varname.replace('sparse', '')
    return name in functions \
           and re.search(
                    f'const (realtype|double) \*{varname}[0]*[,)]+',
                    functions[name]['signature']
                )


class ODEExporter:
    """
    The ODEExporter class generates AMICI C++ files for ODE model as
    defined in symbolic expressions.

    Attributes:
    ----------
    codeprinter: symengine.printing.CCodePrinter
        codeprinter that allows export of symbolic variables as C++ code

    functions: dict
        carries C++ function signatures and other specifications

    modelName: string
        name of the model that will be used for compilation

    modelPath: string
        path to the generated model specific files

    modelSwigPath: string
        path to the generated swig files

    allow_reinit_fixpar_initcond: bool
        indicates whether reinitialization of initial states depending on
        fixedParmeters is allowed for this model

    """

    def __init__(
            self,
            ode_model,
            outdir=None,
            verbose=False,
            assume_pow_positivity=False,
            compiler=None
    ):
        """
        Generate AMICI C++ files for the ODE provided to the constructor.

        Arguments:
        ----------
        ode_model: ODEModel
            ODE definition

        output_dir: string
            see sbml_import.setPaths()

        verbose: bool
            more verbose output if True

        assume_pow_positivity: bool
            if set to true, a special pow function is used to avoid
            problems with state variables that may become negative due
            to numerical errors

        compiler: string
            distutils/setuptools compiler selection to build the python
            extension

        Returns:
        ----------

        Raises:
        ----------

        """

        self.outdir = outdir
        self.verbose = verbose
        self.assume_pow_positivity = assume_pow_positivity
        self.compiler = compiler

        self.codeprinter = CCodePrinter()

        self.modelName = 'model'
        output_dir = os.path.join(os.getcwd(),
                                  f'amici-{self.modelName}')
        self.modelPath = os.path.abspath(output_dir)
        self.modelSwigPath = os.path.join(self.modelPath, 'swig')

        # Signatures and properties of generated model functions (see
        # include/amici/model.h for details)

        self.model = ode_model

        self.functions = functions

        self.allow_reinit_fixpar_initcond = False

    def compileODE(self):
        self.prepareModelFolder()
        self.generateCCode()
        self.compileCCode(compiler=self.compiler, verbose=self.verbose)


    def prepareModelFolder(self):
        """Remove all files from the model folder.

        Arguments:
        ----------

        Returns:

        Raises:

        """
        for file in os.listdir(self.modelPath):
            file_path = os.path.join(self.modelPath, file)
            if os.path.isfile(file_path):
                os.remove(file_path)



    def generateCCode(self):
        """
        Create C++ code files for the model based on.

        Arguments:
        ----------

        Returns:
        ----------

        Raises:
        ----------

        """
        for function in self.functions.keys():
            self.writeFunctionFile(function)

        for name in self.model.symNames():
            self.writeIndexFiles(name)

        self.writeWrapfunctionsCPP()
        self.writeWrapfunctionsHeader()
        self.writeModelHeader()
        self.writeCMakeFile()
        self.writeSwigFiles()
        self.writeModuleSetup()

        shutil.copy(os.path.join(amiciSrcPath, 'main.template.cpp'),
                    os.path.join(self.modelPath, 'main.cpp'))


    def compileCCode(self, verbose=False, compiler=None):
        """
        Compile the generated model code

        Arguments:
        ----------
        verbose: bool
            Make model compilation verbose
        compiler: string
            distutils/setuptools compiler selection to build the python
            extension

        Returns:
        ----------

        Raises:
        ----------

        """

        # setup.py assumes it is run from within the model directory
        moduleDir = self.modelPath
        script_args = [sys.executable, os.path.join(moduleDir, 'setup.py')]

        if verbose:
            script_args.append('--verbose')
        else:
            script_args.append('--quiet')

        script_args.extend(['build_ext', f'--build-lib={moduleDir}'])

        if compiler is not None:
            script_args.extend([f'--compiler={compiler}'])

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

    def writeIndexFiles(self, name):
        """
        Write index file for a symbolic array.

        Arguments:
        ----------
        name: string
            key in self.symbols for which the respective file should be written

        Returns:
        ----------

        Raises:
        ----------


        """
        lines = []
        if name in self.model.symNames():
            if name in sparse_functions:
                symbols = self.model.sparsesym(name)
            else:
                symbols = self.model.sym(name)
        else:
            raise Exception('Unknown symbolic array')

        for index, symbol in enumerate(symbols):
            lines.append(
                f'#define {symbol} {name}[{index}]'
            )

        with open(os.path.join(self.modelPath,f'{name}.h'), 'w') as fileout:
            fileout.write('\n'.join(lines))

    def writeFunctionFile(self, function):
        """
        Write the function `function`.

        Arguments:
        ----------
        function: str
            name of the function to be written (see self.functions)

        Returns:
        ----------

        Raises:
        ----------

        """

        # first generate the equations to make sure we have everything we
        # need in subsequent steps
        if 'sparse' in self.functions[function] and \
                self.functions[function]['sparse']:
            symbol = self.model.sparseeq(function)
        else:
            symbol = self.model.eq(function)

        # function header
        lines = [
            '#include "amici/symbolic_functions.h"',
            '#include "amici/defines.h" //realtype definition',
            'using amici::realtype;',
            '#include <cmath> ',
            ''
        ]

        # function signature
        signature = self.functions[function]['signature']

        if not signature.find('SlsMat') == -1:
            lines.append('#include <sundials/sundials_sparse.h>')

        lines.append('')

        for sym in self.model.symNames():
            # added |double for data
            # added '[0]*' for initial conditions
            if re.search(
                    f'const (realtype|double) \*{sym}[0]*[,)]+', signature
            ):
                lines.append(f'#include "{sym}.h"')

        lines.append('')

        lines.append(f'void {function}_{self.modelName}{signature}{{')

        # function body
        body = self.getFunctionBody(function, symbol)
        if self.assume_pow_positivity \
                and 'assume_pow_positivity' in self.functions[function].keys()\
                and self.functions[function]['assume_pow_positivity']:
            body = [re.sub(r'(^|\W)pow\(', r'\1amici::pos_pow(', line)
                    for line in body]
            # execute this twice to catch cases where the ending ( would be the
            # starting (^|\W) for the following match
            body = [re.sub(r'(^|\W)pow\(', r'\1amici::pos_pow(', line)
                    for line in body]
        self.functions[function]['body'] = body
        lines += body
        lines.append('}')
        #if not body is None:
        with open(os.path.join(
                self.modelPath, f'{self.modelName}_{function}.cpp'), 'w'
        ) as fileout:
            fileout.write('\n'.join(lines))


    def getFunctionBody(self, function, symbol):
        """
        Generate C++ code for body of function `function`.

        Arguments:
        ----------
        function: str
            name of the function to be written (see self.functions)

        Returns:
        ----------

        Raises:
        ----------

        """

        if 'variable' in self.functions[function]:
            variableName = self.functions[function]['variable']
        else:
            variableName = function

        lines = []

        if function == 'sx0_fixedParameters':
            # here we specifically want to overwrite some of the values with 0
            lines.append(' ' * 4 + 'switch(ip) {')
            for ipar in range(self.model.np()):
                lines.append(' ' * 8 + f'case {ipar}:')
                for index, formula in enumerate(
                        self.model.eq('x0_fixedParameters')
                ):
                    if formula != 0 and formula != 0.0:
                        lines.append(' ' * 12 + f'sx0[{index}] = 0.0;')
                lines.append(' ' * 12 + 'break;')
            lines.append('}')
        elif 'sensitivity' in self.functions[function]:
            lines.append(' '*4 + 'switch(ip) {')
            for ipar in range(self.model.np()):
                lines.append(' ' * 8 + f'case {ipar}:')
                lines += self.getSymLines(symbol[:, ipar], variableName, 12)
                lines.append(' ' * 12 + 'break;')
            lines.append('}')
        elif 'multiobs' in self.functions[function]:
            lines.append(' '*4 + 'switch(iy) {')
            for iobs in range(self.model.ny()):
                lines.append(' ' * 8 + f'case {iobs}:')
                lines += self.getSymLines(symbol[:, iobs], variableName, 12)
                lines.append(' ' * 12 + 'break;')
            lines.append('}')
        else:
            if function in ['JSparse', 'JSparseB']:
                rowVals = self.model.rowval(function)
                colPtrs = self.model.colptr(function)
                lines += self.getSparseSymLines(
                    symbol, rowVals, colPtrs, variableName, 4
                )
            else:
                lines += self.getSymLines(symbol, variableName, 4)

        return [line for line in lines if line]


    def writeWrapfunctionsCPP(self):
        """
        Write model-specific 'wrapper' file (wrapfunctions.cpp).

        Arguments:
        ----------

        Returns:
        ----------

        Raises:
        ----------

        """
        templateData = {'MODELNAME': self.modelName}
        applyTemplate(
            os.path.join(amiciSrcPath, 'wrapfunctions.template.cpp'),
            os.path.join(self.modelPath, 'wrapfunctions.cpp'),
            templateData
        )


    def writeWrapfunctionsHeader(self):
        """
        Write model-specific header file (wrapfunctions.h).

        Arguments:
        ----------

        Returns:
        ----------

        Raises:
        ----------

        """
        templateData = {'MODELNAME': str(self.modelName)}
        applyTemplate(
            os.path.join(amiciSrcPath, 'wrapfunctions.ODE_template.h'),
            os.path.join(self.modelPath, 'wrapfunctions.h'),
            templateData
        )

    def writeModelHeader(self):
        """
        Write model-specific header file (MODELNAME.h).

        Arguments:
        ----------

        Returns:
        ----------

        Raises:
        ----------

        """
        if any([math != 0 and math != 0.0 for math in
                self.model.eq('sx0_fixedParameters')]):
            self.allow_reinit_fixpar_initcond = False
        else:
            self.allow_reinit_fixpar_initcond = True

        templateData = {
            'MODELNAME': str(self.modelName),
            'NX': str(self.model.nx()),
            'NXTRUE': str(self.model.nx()),
            'NY': str(self.model.ny()),
            'NYTRUE': str(self.model.ny()),
            'NZ': '0',
            'NZTRUE': '0',
            'NEVENT': '0',
            'NOBJECTIVE': '1',
            'NW': str(len(self.model.sym('w'))),
            'NDWDDX': str(len(self.model.sparsesym('dwdx'))),
            'NDWDP': str(len(self.model.sparsesym('dwdp'))),
            'NNZ': str(len(self.model.sparsesym('JSparse'))),
            'UBW': str(self.model.nx()),
            'LBW': str(self.model.nx()),
            'NP': str(self.model.np()),
            'NK': str(self.model.nk()),
            'O2MODE': 'amici::SecondOrderMode::none',
            'PARAMETERS': str(self.model.val('p'))[1:-1],
            'FIXED_PARAMETERS': str(self.model.val('k'))[1:-1],
            'PARAMETER_NAMES_INITIALIZER_LIST':
                self.getSymbolNameInitializerList('p'),
            'STATE_NAMES_INITIALIZER_LIST':
                self.getSymbolNameInitializerList('x'),
            'FIXED_PARAMETER_NAMES_INITIALIZER_LIST':
                self.getSymbolNameInitializerList('k'),
            'OBSERVABLE_NAMES_INITIALIZER_LIST':
                self.getSymbolNameInitializerList('y'),
            'PARAMETER_IDS_INITIALIZER_LIST':
                self.getSymbolIDInitializerList('p'),
            'STATE_IDS_INITIALIZER_LIST':
                self.getSymbolIDInitializerList('x'),
            'FIXED_PARAMETER_IDS_INITIALIZER_LIST':
                self.getSymbolIDInitializerList('k'),
            'OBSERVABLE_IDS_INITIALIZER_LIST':
                self.getSymbolIDInitializerList('y'),
            'REINIT_FIXPAR_INITCOND':
                'true' if self.allow_reinit_fixpar_initcond else
                'false',
        }
        applyTemplate(
            os.path.join(amiciSrcPath,'model_header.ODE_template.h'),
            os.path.join(self.modelPath,f'{self.modelName}.h'),
            templateData
        )

    def getSymbolNameInitializerList(self, name):
        """
        Get SBML name initializer list for vector of names for the given model
        entity

        Arguments:
        ----------
        name: string
            any key present in self.symbols

        Returns:
        ----------
        Template initializer list of names

        Raises:
        ----------

        """
        return '\n'.join(
            [f'"{symbol}",' for symbol in self.model.name(name)]
        )

    def getSymbolIDInitializerList(self, name):
        """
        Get C++ initializer list for vector of names for the given model entity

        Arguments:
        ----------
        name: string
            any key present in self.symbols

        Returns:
        ----------
        Template initializer list of ids

        Raises:
        ----------

        """
        return '\n'.join(
            [f'"{symbol}",' for symbol in self.model.sym(name)]
        )

    def writeCMakeFile(self):
        """
        Write CMake CMakeLists.txt file for this model.

        Arguments:
        ----------

        Returns:
        ----------

        Raises:
        ----------

        """
        sources = [self.modelName + '_' + function + '.cpp '
                   for function in self.functions.keys()
                   if self.functions[function]['body'] is not None]
        templateData = {'MODELNAME': self.modelName,
                        'SOURCES': '\n'.join(sources)}
        applyTemplate(
            os.path.join(amiciSrcPath, 'CMakeLists.template.txt'),
            os.path.join(self.modelPath, 'CMakeLists.txt'),
            templateData
        )

    def writeSwigFiles(self):
        """
        Write SWIG interface files for this model.

        Arguments:
        ----------

        Returns:
        ----------

        Raises:
        ----------

        """
        if not os.path.exists(self.modelSwigPath):
            os.makedirs(self.modelSwigPath)
        templateData = {'MODELNAME': self.modelName}
        applyTemplate(
            os.path.join(amiciSwigPath, 'modelname.template.i'),
            os.path.join(self.modelSwigPath, self.modelName + '.i'),
            templateData
        )
        shutil.copy(os.path.join(amiciSwigPath, 'CMakeLists_model.txt'),
                    os.path.join(self.modelSwigPath, 'CMakeLists.txt'))

    def writeModuleSetup(self):
        """
        Create a distutils setup.py file for compile the model module.

        Arguments:
        ----------

        Returns:
        ----------

        Raises:
        ----------

        """

        templateData = {'MODELNAME': self.modelName,
                        'VERSION': '0.1.0'}
        applyTemplate(os.path.join(amiciModulePath, 'setup.template.py'),
                      os.path.join(self.modelPath, 'setup.py'), templateData)

        # write __init__.py for the model module
        if not os.path.exists(os.path.join(self.modelPath, self.modelName)):
            os.makedirs(os.path.join(self.modelPath, self.modelName))

        applyTemplate(
            os.path.join(amiciModulePath, '__init__.template.py'),
            os.path.join(self.modelPath, self.modelName, '__init__.py'),
            templateData
        )

    def getSymLines(self, symbols, variable, indentLevel):
        """
        Generate C++ code for assigning symbolic terms in symbols to C++ array
        `variable`.

        Arguments:
        ----------
        symbols: list
            vectors of symbolic terms

        variable: string
            name of the C++ array to assign to

        indentLevel: int
            indentation level (number of leading blanks)


        Returns:
        ----------
            C++ code as list of lines

        Raises:
        ----------

        """

        lines = [' ' * indentLevel + f'{variable}[{index}] = '
                 f'{self.printWithException(math)};'
                 if not (math == 0 or math == 0.0)
                 else ''
                 for index, math in enumerate(symbols)]

        try:
            lines.remove('')
        except:
            pass

        return lines

    def getSparseSymLines(
            self, symbolList, RowVals, ColPtrs, variable, indentLevel
    ):
        """
        Generate C++ code for assigning sparse symbolic matrix to a C++ array
        `variable`.

        Arguments:
        ----------
        symbolList: list
            symbolic terms

        RowVals: list
            row indices of each nonzero entry
            (see CVODES SlsMat documentation for details)

        ColPtrs: list
            indices of the first column entries
            (see CVODES SlsMat documentation for details)

        variable: string
            name of the C++ array to assign to

        indentLevel: int
            indentation level (number of leading blanks)

        Returns:
        ----------
        C++ code as list of lines

        Raises:
        ----------

        """
        lines = [
            ' ' * indentLevel + f'{variable}->indexvals[{index}] = ' 
            f'{self.printWithException(math)};'
            for index, math in enumerate(RowVals)
        ]

        lines.extend(
            [' ' * indentLevel + f'{variable}->indexptrs[{index}] = '
            f'{self.printWithException(math)};'
            for index, math in enumerate(ColPtrs)]
        )
        lines.extend(
            [' ' * indentLevel + f'{variable}->data[{index}] = '
            f'{self.printWithException(math)};'
            for index, math in enumerate(symbolList)]
        )

        return lines

    def printWithException(self, math):
        """Generate C++ code for a symbolic expression

        Arguments:
        ----------
        math:
            symbolic expression

        Returns:
        ----------
        C++ code for the specified expression

        Raises:
        ----------
        Exception:
            The specified expression contained an unsupported function

        """
        try:
            ret = self.codeprinter.doprint(math)
            ret = re.sub(r'(^|\W)M_PI(\W|$)', r'\1amici::pi\2', ret)
            return ret
        except:
            raise Exception(
                f'Encountered unsupported function in expression "{math}"!'
            )

    def setPaths(self, output_dir):
        """
        Set output paths for the model and create if necessary

        Arguments:
        ----------
        output_dir: string
            relative or absolute path where the generated model code is to
            be placed. will be created if does not exists.

        Returns:
        ----------

        Raises:
        ----------

        """
        self.modelPath = os.path.abspath(output_dir)
        self.modelSwigPath = os.path.join(self.modelPath, 'swig')

        for directory in [self.modelPath, self.modelSwigPath]:
            if not os.path.exists(directory):
                os.makedirs(directory)

    def setName(self, modelName):
        """
        Sets the model name

        Arguments:
        ----------
        modelName: string
            name of the model (must only contain valid filename characters)

        Returns:
        ----------

        Raises:
        ----------

        """
        self.modelName = modelName

def getSymbols(prefix, length):
    """
    Get symbolic matrix with symbols prefix0..prefix(length-1).

    Arguments:
    ----------
    prefix: string
        variable name

    length: int
        number of symbolic variables + 1

    Returns:
    ----------
    A symbolic matrix with symbols prefix0..prefix(length-1)

    Raises:
    ----------

    """
    return sp.DenseMatrix([
        sp.Symbol(prefix + str(i)) for i in range(0, length)
    ])

def getSymbolicDiagonal(matrix):
    """
    Get symbolic matrix with diagonal of matrix `matrix`.

    Arguments:
    ----------
    matrix: symengine.DenseMatrix
        Matrix from which to return the diagonal

    Returns:
    ----------
    A Symbolic matrix with the diagonal of `matrix`.

    Raises:
    ----------
    Exception: The provided matrix was not square
    """
    if not matrix.cols == matrix.rows:
        raise Exception('Provided matrix is not square!')

    diagonal = [matrix[index,index] for index in range(matrix.cols)]

    return sp.DenseMatrix(diagonal)

class TemplateAmici(Template):
    """
    Template format used in AMICI (see string.template for more details).

    Attributes:
    ----------
    delimiter: string
        delimiter that identifies template variables

    """
    delimiter = 'TPL_'


def applyTemplate(sourceFile,targetFile,templateData):
    """Load source file, apply template substitution as provided in templateData and save as targetFile.

    Arguments:
    ----------
    sourceFile: string
        relative or absolute path to template file

    targetFile: string
        relative or absolute path to output file

    templateData: dict
        template keywords to substitute (key is template variable without
        TemplateAmici.delimiter)

    Returns:
    ----------

    Raises:
    ----------

    """
    with open(sourceFile) as filein:
        src = TemplateAmici(filein.read())
    result = src.safe_substitute(templateData)
    with open(targetFile, 'w') as fileout:
        fileout.write(result)

def getSparseSymbols(functions, symbolName):
    """
    Create sparse symbolic matrix.

    Arguments:
    ----------
    functions: dict
        carries function definition

    symbolName: str
        name of the function

    Returns:
    ----------
    sparseMatrix:
        sparse matrix containing symbolic entries

    symbolList:
        symbolic vector containing the list of symbol names

    sparseList:
        symbolic vector containing the list of symbol formulas

    symbolColPtrs:
        Column pointer as specified in the SlsMat definition in CVODES

    symbolRowVals:
        Row Values as specified in the SlsMat definition in CVODES

    Raises:
    ----------

    """
    matrix = functions[symbolName]['sym']
    symbolIndex = 0
    sparseMatrix = sp.zeros(matrix.rows,matrix.cols)
    symbolList = []
    sparseList = []
    symbolColPtrs = []
    symbolRowVals = []
    for col in range(0,matrix.cols):
        symbolColPtrs.append(symbolIndex)
        for row in range(0, matrix.rows):
            if not matrix[row, col] == 0:
                name = f'{symbolName}{symbolIndex}'
                sparseMatrix[row, col] = sp.sympify(name)
                symbolList.append(name)
                sparseList.append(matrix[row, col])
                symbolRowVals.append(row)
                symbolIndex += 1
    symbolColPtrs.append(symbolIndex)
    sparseList = sp.DenseMatrix(sparseList)

    return (
        sparseMatrix, symbolList, sparseList, symbolColPtrs, symbolRowVals
    )
