""" @package amici.ode_export The C++ ODE export module for python
"""
#!/usr/bin/env python3

import symengine as sp
import re
import shutil
import subprocess
import sys
import os
import copy
import numbers

from symengine.printing import CCodePrinter
from string import Template

from . import amiciSwigPath, amiciSrcPath, amiciModulePath

## prototype for generated C++ functions, keys are the names of functions
#
# signature: defines the argument part of the function signature,
# input variables
# should have a const flag
#
# assume_pow_positivity: identifies the functions on which
# assume_pow_positivity will have an effect when specified during model
# generation. generally these are functions that are used for solving the
# ODE, where negative values may negatively affect convergence of the
# integration algorithm
#
# sparse: specifies whether the result of this function will be stored in
# sparse format. sparse format means that the function will only return an
# array of nonzero values and not a full matrix.
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
            '(double *Jy, const int iy, const realtype *p,'
            ' const realtype *k, const double *y,'
            ' const double *sigmay, const double *my)',
    },
    'dJydsigmay': {
        'signature':
            '(double *dJydsigmay, const int iy, const realtype *p,'
            ' const realtype *k, const double *y,'
            ' const double *sigmay, const double *my)',
    },
    'dJydy': {
        'signature':
            '(double *dJydy, const int iy, const realtype *p,'
            ' const realtype *k, const double *y,'
            ' const double *sigmay, const double *my)',
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
        'assume_pow_positivity':
            True,
    },
    'dydx': {
        'signature':
            '(double *dydx, const realtype t, const realtype *x,'
            ' const realtype *p, const realtype *k,'
            ' const realtype *h, const realtype *w, const realtype *dwdx)',
    },
    'dydp': {
        'signature':
            '(double *dydp, const realtype t, const realtype *x,'
            ' const realtype *p, const realtype *k, const realtype *h,'
            ' const int ip, const realtype *w, const realtype *dwdp)',
    },
    'dsigmaydp': {
        'signature':
            '(double *dsigmaydp, const realtype t, const realtype *p,'
            ' const realtype *k, const int ip)',
    },
    'qBdot': {
        'signature':
            '(realtype *qBdot, const int ip, const realtype t,'
            ' const realtype *x, const realtype *p, const realtype *k,'
            ' const realtype *h, const realtype *xB,'
            ' const realtype *w, const realtype *dwdp)',
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
            '(realtype *x0_fixedParameters, const realtype t, const realtype *p,'
            ' const realtype *k)',
    },
    'sx0': {
        'signature':
            '(realtype *sx0, const realtype t,const realtype *x,'
            ' const realtype *p, const realtype *k, const int ip)',
    },
    'sx0_fixedParameters': {
        'signature':
            '(realtype *sx0_fixedParameters, const realtype t,const realtype *x0,'
            ' const realtype *p, const realtype *k, const int ip)',
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
            ' const realtype *h, const realtype *w)',
    }
}

## list of sparse functions
sparse_functions = [
    function for function in functions
    if 'sparse' in functions[function]
    and functions[function]['sparse']
]
## list of sensitivity functions
sensi_functions = [
    function for function in functions
    if 'const int ip' in functions[function]['signature']
    and function is not 'sxdot'
]
## list of multiobs functions
multiobs_functions = [
    function for function in functions
    if 'const int iy' in functions[function]['signature']
]

def var_in_function_signature(name, varname):
    """Checks if the values for a symbolic variable is passed in the signature
    of a function

    Arguments:

        name: name of the function @type str

        varname: name of the symbolic variable @type str

    Returns:

    Raises:

    """
    varname = varname.replace('sparse', '')
    return name in functions \
           and re.search(
                    f'const (realtype|double) \*{varname}[0]*[,)]+',
                    functions[name]['signature']
                )

class ModelQuantity:
    """Base class for model components

    Attributes:

    """
    def __init__(self, identifier,  name, value):
        """Create a new ModelQuantity instance.

        Arguments:
            identifier: unique identifier of the quantity @type sympy.Symbol

            name: individual name of the quantity (does not need to be unique)
            @type str

            value: either formula, numeric value or initial value

        Returns:
        ModelQuantity instance

        Raises:
        TypeError:
            is thrown if input types do not match documented types
        """
        if not isinstance(identifier, sp.Symbol):
            raise TypeError(f'identifier must be sympy.Symbol, was '
                            f'{type(identifier)}')
        self._identifier = identifier

        if not isinstance(name, str):
            raise TypeError(f'name must be str, was {type(name)}')
        self._name = name

        if isinstance(value, sp.RealNumber) or isinstance(value,
                                                          numbers.Number):
            value = float(value)
        if not isinstance(value, sp.Basic) and not isinstance(value, float):
            raise TypeError(f'value must be sympy.Symbol or float, was '
                            f'{type(value)}')
        self._value = value

    def __repr__(self):
        """Representation of the ModelQuantity object

        Arguments:

        Returns:
        string representation of the ModelQuantity

        Raises:

        """
        return str(self._identifier)


class State(ModelQuantity):
    """A State variable defines an entity that evolves with time according to the
    provided time derivative, abbreviated by `x`

    """

    def __init__(self, identifier, name, value, dt):
        """Create a new State instance. Extends ModelQuantity.__init__ by dt

        Arguments:
            identifier: unique identifier of the state @type sympy.Symbol

            name: individual name of the state (does not need to be unique)
            @type str

            value: initial value @type symengine.Basic

            dt: time derivative @type symengine.Basic

        Returns:
        ModelQuantity instance

        Raises:
        TypeError:
            is thrown if input types do not match documented types
        """
        super(State, self).__init__(identifier, name, value)
        if not isinstance(dt, sp.Basic):
            raise TypeError(f'dt must be sympy.Symbol, was '
                            f'{type(dt)}')
        self._dt = dt

class Observable(ModelQuantity):
    """An Observable links model simulations to experimental measurements,
    abbreviated by `y`

    """
    def __init__(self, identifier, name, value):
        """Create a new Observable instance.

        Arguments:
            identifier: unique identifier of the Observable @type sympy.Symbol

            name: individual name of the Observable (does not need to be
            unique) @type str

            value: formula @type symengine.Basic

        Returns:
        ModelQuantity instance

        Raises:
        TypeError:
            is thrown if input types do not match documented types
        """
        super(Observable, self).__init__(identifier, name, value)


class SigmaY(ModelQuantity):
    """A Standard Deviation SigmaY rescales the distance between simulations and
    measurements when computing residuals, abbreviated by `sigmay`

    """
    def __init__(self, identifier, name, value):
        """Create a new Standard Deviation instance.

        Arguments:
            identifier: unique identifier of the Standard Deviation
            @type sympy.Symbol

            name: individual name of the Standard Deviation (does not need to
            be unique) @type str

            value: formula @type symengine.Basic

        Returns:
        ModelQuantity instance

        Raises:
        TypeError:
            is thrown if input types do not match documented types
        """
        super(SigmaY, self).__init__(identifier, name, value)

class Expression(ModelQuantity):
    """An Expressions is a recurring elements in symbolic formulas. Specifying
    this may yield more compact expression which may lead to substantially
    shorter model compilation times, but may also reduce model simulation time,
    abbreviated by `w`

    """
    def __init__(self, identifier, name, value):
        """Create a new Expression instance.

        Arguments:
            identifier: unique identifier of the Expression @type sympy.Symbol

            name: individual name of the Expression (does not need to be
            unique) @type str

            value: formula @type symengine.Basic

        Returns:
        ModelQuantity instance

        Raises:
        TypeError:
            is thrown if input types do not match documented types
        """
        super(Expression, self).__init__(identifier, name, value)


class Parameter(ModelQuantity):
    """A Parameter is a free variable in the model with respect to which
    sensitivites may be computed, abbreviated by `p`

    """

    def __init__(self, identifier, name, value):
        """Create a new Expression instance.

        Arguments:
            identifier: unique identifier of the Parameter @type sympy.Symbol

            name: individual name of the Parameter (does not need to be
            unique) @type str

            value: numeric value @type float

        Returns:
        ModelQuantity instance

        Raises:
        TypeError:
            is thrown if input types do not match documented types
        """
        super(Parameter, self).__init__(identifier, name, value)


class Constant(ModelQuantity):
    """A Constant is a fixed variable in the model with respect to which
    sensitivites cannot be computed, abbreviated by `k`

    """

    def __init__(self, identifier, name, value):
        """Create a new Expression instance.

        Arguments:
            identifier: unique identifier of the Constant @type sympy.Symbol

            name: individual name of the Constant (does not need to be unique)
             @type str

            value: numeric value @type float

        Returns:
        ModelQuantity instance

        Raises:
        TypeError:
            is thrown if input types do not match documented types
        """
        super(Constant, self).__init__(identifier, name, value)


class LogLikelihood(ModelQuantity):
    """A LogLikelihood defines the distance between measurements and
    experiments for a particular observable. The final LogLikelihood value
    in the simulation will be the sum of all specified LogLikelihood
    instances evaluated at all timepoints, abbreviated by `Jy`

    """

    def __init__(self, identifier, name, value):
        """Create a new Expression instance.

        Arguments:
            identifier: unique identifier of the LogLikelihood
            @type sympy.Symbol

            name: individual name of the LogLikelihood (does not need to be
             unique) @type str

            value: formula @type symengine.Basic

        Returns:
        ModelQuantity instance

        Raises:
        TypeError:
            is thrown if input types do not match documented types
        """
        super(LogLikelihood, self).__init__(identifier, name, value)

## indicates which type some attributes in ODEModel should have
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
    """An ODEModel defines an Ordinay Differential Equation as set of
    ModelQuantities. This class provides general purpose interfaces to
    compute arbitrary symbolic derivatives that are necessary for model
    simulation or sensitivity computation

    Attributes:
        _states: list of State instances  @type list

        _observables: list of Observable instances  @type list

        _sigmays: list of SigmaY instances  @type list

        _parameters: list of Parameter instances  @type list

        _loglikelihoods: list of LogLikelihood instances  @type list

        _expressions: list of Expression instances  @type list

        _symboldim_funs: define functions that compute model dimensions, these
        are functions as the underlying symbolic expressions have not been
        populated at compile time @type dict

        _eqs: carries symbolic formulas of the symbolic variables of the model
        @type dict

        _sparseeqs: carries linear list of all symbolic formulas for sparsified
        variables  @type dict

        _vals: carries numeric values of symbolic identifiers of the symbolic
        variables of the model  @type dict

        _names: carries names of symbolic identifiers of the symbolic variables
        of the model  @type dict

        _syms: carries symbolic identifiers of the symbolic variables of the
        model  @type dict

        _sparsesyms: carries linear list of all symbolic identifiers for
        sparsified variables  @type dict

        _colptrs: carries column pointers for sparsified variables. See
        SlsMat definition in CVODES for more details about ColPtrs @type dict

        _rowvals: carries row values for sparsified variables. See SlsMat
        definition in CVODES for more details about RowVals @type dict

        _equation_prototype: defines the attribute from which an equation
        should be generated via list comprehension (see
        ODEModel.generateEquation()) @type dict

        _variable_prototype: defines the attribute from which a variable should
        be generated via list comprehension (see ODEModel.generateSymbol())
        @type dict

        _value_prototype: defines the attribute from which a value should be
        generated via list comprehension (see ODEModel.generateValue())
        @type dict

        _total_derivative_prototypes: defines how a total derivative equation
        is computed for an equation, key defines the name and values should
        be arguments for ODEModel.totalDerivative() @type dict

        _multiplication_prototypes: defines how a multiplication equation is
        computed for an equation, key defines the name and values should be
        arguments for ODEModel.multiplication() @type dict

        _lock_total_derivative: set this to true when computing a total
        derivative from a partial derivative call to enforce a partial
        derivative in the next recursion. prevents infinite recursion

    """

    def __init__(self):
        """Create a new ODEModel instance.

        Arguments:

        Returns:
        New ODEModel instance

        Raises:

        """
        self._states = []
        self._observables = []
        self._sigmays = []
        self._parameters = []
        self._constants = []
        self._loglikelihoods = []
        self._expressions = []
        self._symboldim_funs = {
            'sx': self.nx,
            'v': self.nx,
            'vB': self.nx,
            'xB': self.nx,
            'sigmay': self.ny,
        }
        self._eqs = dict()
        self._sparseeqs = dict()
        self._vals = dict()
        self._names = dict()
        self._syms = dict()
        self._sparsesyms = dict()
        self._colptrs = dict()
        self._rowvals = dict()

        self._equation_prototype = {
            'x0': '_states',
            'y': '_observables',
            'Jy': '_loglikelihoods',
            'w': '_expressions',
            'sigmay': '_sigmays',
        }
        self._variable_prototype = {
            'x': '_states',
            'y': '_observables',
            'p': '_parameters',
            'k': '_constants',
            'w': '_expressions',
            'sigmay': '_sigmays'
        }
        self._value_prototype = {
            'p': '_parameters',
            'k': '_constants',
        }
        self._total_derivative_prototypes = {
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
        self._multiplication_prototypes = {
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
                'sign': -1,
            },
            'xBdot': {
                'x': 'JB',
                'y': 'xB',
                'sign': -1,
            },
        }

        self._lock_total_derivative = False

    def import_from_sbml_importer(self, si):
        """Imports a model specification from a amici.SBMLImporter instance.

        Arguments:
            si: imported SBML model @type amici.SbmlImporter


        Returns:


        Raises:

        """

        symbols = copy.copy(si.symbols)

        # setting these equations prevents native equation generation
        self._eqs['dxdotdw'] = si.stoichiometricMatrix
        self._eqs['w'] = si.fluxVector
        self._syms['w'] = sp.DenseMatrix(
            [sp.Symbol(f'flux_r{idx}') for idx in range(len(si.fluxVector))]
        )
        self._eqs['dxdotdx'] = sp.zeros(si.stoichiometricMatrix.shape[0])
        symbols['species']['dt'] = \
            si.stoichiometricMatrix * self.sym('w')

        for symbol in [s for s in symbols if s != 'my']:
            # transform dict of lists into a list of dicts
            protos = [dict(zip(symbols[symbol], t))
                      for t in zip(*symbols[symbol].values())]
            for proto in protos:
                self.add_component(symbol_to_type[symbol](**proto))

        self.generateBasicVariables()

    def add_component(self, component):
        """Adds a new ModelQuantity to the model.

        Arguments:
            component: model quantity to be added @type ModelQuantity

        Returns:

        Raises:

        """
        for comp_type in [Observable, Expression, Parameter, Constant, State,
                          LogLikelihood, SigmaY]:
            if isinstance(component, comp_type):
                getattr(self, f'_{type(component).__name__.lower()}s').append(
                    component
                )
                return
        Exception(f'Invalid component type {type(component)}')

    def nx(self):
        """Number of states.

        Arguments:

        Returns:
        number of state variable symbols

        Raises:

        """
        return len(self.sym('x'))

    def ny(self):
        """Number of Observables.

        Arguments:

        Returns:
        number of observable symbols

        Raises:

        """
        return len(self.sym('y'))

    def nk(self):
        """Number of Constants.

        Arguments:

        Returns:
        number of constant symbols

        Raises:

        """
        return len(self.sym('k'))

    def np(self):
        """Number of Parameters.

        Arguments:

        Returns:
        number of parameter symbols

        Raises:

        """
        return len(self.sym('p'))

    def sym(self, name):
        """Returns (and constructs if necessary) the identifiers for a symbolic
        entity.

        Arguments:
            name: name of the symbolic variable @type str

        Returns:
        containing the symbolic identifiers @type symengine.DenseMatrix

        Raises:

        """
        if name not in self._syms:
            self._generateSymbol(name)
        return self._syms[name]

    def sparsesym(self, name):
        """Returns (and constructs if necessary) the sparsified identifiers for a
        sparsified symbolic variable.

        Arguments:
            name: name of the symbolic variable @type str

        Returns:
        linearized symengine.DenseMatrix containing the symbolic identifiers

        Raises:

        """
        if name not in sparse_functions:
            raise Exception(f'{name} is not marked as sparse')
        if name not in self._sparsesyms:
            self._generateSparseSymbol(name)
        return self._sparsesyms[name]

    def eq(self, name):
        """Returns (and constructs if necessary) the formulas for a symbolic
        entity.

        Arguments:
            name: name of the symbolic variable @type str

        Returns:
        symengine.DenseMatrix containing the symbolic identifiers

        Raises:

        """
        if name not in self._eqs:
            self._computeEquation(name)
        return self._eqs[name]

    def sparseeq(self, name):
        """Returns (and constructs if necessary) the sparsified formulas for a
        sparsified symbolic variable.

        Arguments:
            name: name of the symbolic variable @type str

        Returns:
        linearized symengine.DenseMatrix containing the symbolic formulas

        Raises:

        """
        if name not in sparse_functions:
            raise Exception(f'{name} is not marked as sparse')
        if name not in self._sparseeqs:
            self._generateSparseSymbol(name)
        return self._sparseeqs[name]

    def colptr(self, name):
        """Returns (and constructs if necessary) the column pointers for
        a sparsified symbolic variable.

        Arguments:
            name: name of the symbolic variable @type str

        Returns:
        symengine.DenseMatrix containing the column pointers

        Raises:

        """
        if name not in sparse_functions:
            raise Exception(f'{name} is not marked as sparse')
        if name not in self._sparseeqs:
            self._generateSparseSymbol(name)
        return self._colptrs[name]

    def rowval(self, name):
        """Returns (and constructs if necessary) the row values for a sparsified
        symbolic variable.

        Arguments:
            name: name of the symbolic variable @type str

        Returns:
        symengine.DenseMatrix containing the row values

        Raises:

        """
        if name not in sparse_functions:
            raise Exception(f'{name} is not marked as sparse')
        if name not in self._sparseeqs:
            self._generateSparseSymbol(name)
        return self._rowvals[name]

    def val(self, name):
        """Returns (and constructs if necessary) the numeric values of a symbolic
        entity

        Arguments:
            name: name of the symbolic variable @type str

        Returns:
        list containing the numeric values

        Raises:

        """
        if name not in self._vals:
            self._generateValue(name)
        return self._vals[name]

    def name(self, name):
        """Returns (and constructs if necessary) the names of a symbolic variable

        Arguments:
            name: name of the symbolic variable @type str

        Returns:
        list of names

        Raises:

        """
        if name not in self._names:
            self._generateName(name)
        return self._names[name]

    def _generateSymbol(self, name):
        """Generates the symbolic identifiers for a symbolic variable

        Arguments:
            name: name of the symbolic variable @type str

        Returns:

        Raises:

        """
        if name in self._variable_prototype:
            component = self._variable_prototype[name]
            self._syms[name] = sp.DenseMatrix(
                [comp._identifier for comp in getattr(self, component)]
            )
            if name == 'y':
                self._syms['my'] = sp.DenseMatrix(
                    [sp.Symbol(f'm{comp._identifier}')
                    for comp in getattr(self, component)]
                )
            return
        elif name in sparse_functions:
            self._generateSparseSymbol(name)
            return
        elif name in self._symboldim_funs:
            length = self._symboldim_funs[name]()
        elif name in sensi_functions:
            length = self.eq(name).shape[0]
        else:
            length = len(self.eq(name))

        self._syms[name] = sp.DenseMatrix([
            sp.Symbol(f'{name}{i}') for i in range(length)
        ])

    def generateBasicVariables(self):
        """Generates the symbolic identifiers for all variables in
        ODEModel.variable_prototype

        Arguments:

        Returns:

        Raises:

        """
        for var in self._variable_prototype:
            if var not in self._syms:
                self._generateSymbol(var)

    def _generateSparseSymbol(self, name):
        """Generates the sparse symbolic identifiers, symbolic identifiers,
        sparse equations, column pointers and row values for a symbolic variable

        Arguments:
            name: name of the symbolic variable @type str

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

    def _computeEquation(self, name):
        """computes the symbolic formula for a symbolic variable

        Arguments:
            name: name of the symbolic variable @type str

        Returns:

        Raises:

        """
        match_deriv = re.match(r'd([\w]+)d([a-z]+)', name)

        if name in self._equation_prototype:
            self._equationFromComponent(name, self._equation_prototype[name])

        elif name in self._total_derivative_prototypes:
            args = self._total_derivative_prototypes[name]
            args['name'] = name
            self._totalDerivative(**args)

        elif name in self._multiplication_prototypes:
            args = self._multiplication_prototypes[name]
            args['name'] = name
            self._multiplication(**args)

        elif name == 'xdot':
            self._eqs['xdot'] = sp.DenseMatrix(
                [comp._dt for comp in self._states]
            )

        elif name == 'sx0':
            self._derivative(name[1:], 'p', name=name)

        elif name == 'sx0_fixedParameters':
            # deltax = -x+x0_fixedParameters if x0_fixedParameters>0 else 0
            # deltasx = -sx+dx0_fixedParametersdx*sx+dx0_fixedParametersdp
            # if x0_fixedParameters>0 else 0
            # sx0_fixedParameters = sx+deltasx =
            # dx0_fixedParametersdx*sx+dx0_fixedParametersdp
            self._eqs[name] = \
                self.eq('x0_fixedParameters').jacobian(self.sym('p'))

            for ip in range(self._eqs[name].shape[1]):
                self._eqs[name][:,ip] += \
                    self.eq('x0_fixedParameters').jacobian(self.sym('x')) \
                    * self.sym('sx0') \

            for index, formula in enumerate(
                    self.eq('x0_fixedParameters')
            ):
                if formula == 0 or formula == 0.0:
                    self._eqs[name][index, :] = \
                        sp.zeros(1, self._eqs[name].shape[1])


        elif name == 'JB':
            self._eqs[name] = self.eq('J').transpose()

        elif name == 'JDiag':
            self._eqs[name] = getSymbolicDiagonal(self.eq('J'))

        elif name == 'x0_fixedParameters':
            k = self.sym('k')
            self._eqs[name] = sp.DenseMatrix([
                eq
                # check if the equation contains constants
                if any([sym in eq.free_symbols for sym in k])
                # if not set to zero
                else 0.0
                for eq in self.eq('x0')
            ])

        elif name in ['JSparse', 'JSparseB']:
            self._eqs[name] = self.eq(name.replace('Sparse', ''))

        elif match_deriv:
            self._derivative(match_deriv.group(1), match_deriv.group(2))

        else:
            raise Exception(f'Unknown equation {name}')

        if name in ['Jy', 'dydx']:
            # do not transpose if we compute the partial derivative as part of
            # a total derivative
            if not self._lock_total_derivative:
                self._eqs[name] = self._eqs[name].transpose()

    def symNames(self):
        """Returns a list of names of generated symbolic variables

        Arguments:

        Returns:
        list of names

        Raises:

        """
        return list(self._syms.keys())

    def _derivative(self, eq, var, name=None):
        """Creates a new symbolic variable according to a derivative

        Arguments:
            eq: name of the symbolic variable that defines the formula
            @type str

            var: name of the symbolic variable that defines the identifiers
            with respect to which the derivatives are to be computed @type str

            name: name of resulting symbolic variable, default is d{eq}d{var}
            @type str


        Returns:

        Raises:

        """
        if not name:
            name = f'd{eq}d{var}'

        # automatically detect chainrule
        if var_in_function_signature(eq, 'w') \
                and not self._lock_total_derivative \
                and var is not 'w'\
                and self.sym('w').size:
            self._lock_total_derivative = True
            self._totalDerivative(name, eq, 'w', var)
            self._lock_total_derivative = False
            return

        # partial derivative
        if self.eq(eq).size and self.sym(var).size:
            if eq == 'Jy':
                eq = self.eq(eq).transpose()
            else:
                eq = self.eq(eq)
            self._eqs[name] = eq.jacobian(self.sym(var))
        else:
            self._eqs[name] = sp.DenseMatrix([])

    def _totalDerivative(self, name, eq, chainvar, var,
                         dydx_name=None, dxdz_name=None):
        """Creates a new symbolic variable according to a total derivative
        using the chain rule

        Arguments:
            name: name of resulting symbolic variable @type str

            eq: name of the symbolic variable that defines the formula @type str

            chainvar: name of the symbolic variable that defines the
            identifiers whith respect to which the chain rule is applied
            @type str

            var: name of the symbolic variable that defines the identifiers
            whith respect to which the derivatives are to be computed @type str

            dydxname: defines the name of the symbolic variable that
            defines the derivative of the `eq` with respect to `chainvar`,
            default is d{eq}d{chainvar} @type str

            dxdzname: defines the name of the symbolic variable that
            defines the derivative of the `chainvar` with respect to `var`,
            default is d{chainvar}d{var} @type str

        Returns:

        Raises:

        """
        # compute total derivative according to chainrule
        # Dydz = dydx*dxdz + dydz
        variables = dict()
        variables['dydx'] = dict()
        variables['dxdz'] = dict()
        variables['dydz'] = dict()

        variables['dydx']['name'] = f'd{eq}d{chainvar}' \
            if dydx_name is None else dydx_name
        variables['dxdz']['name'] = f'd{chainvar}d{var}' \
            if dxdz_name is None else dxdz_name
        variables['dydz']['name'] = f'd{eq}d{var}'

        # if the symbol appears in the function signature,
        # we can use the respective symbol instead of the equation
        for var in variables:
            varname = variables[var]["name"]
            if var_in_function_signature(name, varname):
                variables[var]['sym'] = self.sym(varname)
            else:
                variables[var]['sym'] = self.eq(varname)

        self._eqs[name] = sp.DenseMatrix([])

        if variables['dydz']['sym'].size:
            self._eqs[name] += variables['dydz']['sym']

        if variables['dydx']['sym'].size and variables['dxdz']['sym'].size:
            self._eqs[name] += \
                variables['dydx']['sym'] * variables['dxdz']['sym']



    def _multiplication(self, name, x, y,
                        transpose_x=False, sign=1):
        """Creates a new symbolic variable according to a multiplication

        Arguments:
            name: name of resulting symbolic variable, default is d{eq}d{var}
             @type str

            x: name of the symbolic variable that defines the first factor
            @type str

            y: name of the symbolic variable that defines the second factor
            @type str

            transpose_x: indicates whether the first factor should be
            transposed before multiplication @type bool

            sign: defines the sign of the product, should be +1 or -1 @type int


        Returns:

        Raises:

        """
        if sign not in [-1, 1]:
            raise TypeError(f'sign must be +1 or -1, was {sign}')

        variables = dict()
        for varname in [x, y]:
            if var_in_function_signature(name, varname):
                variables[varname] = self.sym(varname)
            else:
                variables[varname] = self.eq(varname)

        if transpose_x:
            xx = variables[x].transpose()
        else:
            xx = variables[x]

        self._eqs[name] = sign * xx * variables[y]

    def _equationFromComponent(self, name, component):
        """Generates the formulas of a symbolic variable from the attributes

        Arguments:
            name: name of resulting symbolic variable @type str

            component: name of the attribute @type str

        Returns:

        Raises:

        """
        self._eqs[name] = sp.DenseMatrix(
            [comp._value for comp in getattr(self, component)]
        )

    def _generateValue(self, name):
        """Generates the numeric values of a symbolic variable from value
        prototypes

        Arguments:
            name: name of resulting symbolic variable @type str

        Returns:

        Raises:

        """
        if name in self._value_prototype:
            component = self._value_prototype[name]
        else:
            raise Exception(f'No values for {name}')

        self._vals[name] = [comp._value for comp in getattr(self, component)]

    def _generateName(self, name):
        """Generates the names of a symbolic variable from variable prototypes or
        equation prototypes

        Arguments:
            name: name of resulting symbolic variable @type str

        Returns:

        Raises:

        """
        if name in self._variable_prototype:
            component = self._variable_prototype[name]
        elif name in self._equation_prototype:
            component = self._equation_prototype[name]
        else:
            raise Exception(f'No names for {name}')

        self._names[name] = [comp._name for comp in getattr(self, component)]


class ODEExporter:
    """The ODEExporter class generates AMICI C++ files for ODE model as
    defined in symbolic expressions.

    Attributes:

        model: ODE definition @type ODEModel

        outdir: see sbml_import.setPaths() @type str

        verbose: more verbose output if True @type bool

        assume_pow_positivity: if set to true, a special pow function is
        used to avoid problems with state variables that may become negative
        due to numerical errors @type bool

        compiler: distutils/setuptools compiler selection to build the
        python extension @type str

        codeprinter: allows export of symbolic variables as C++ code

        functions: carries C++ function signatures and other specifications
        @type dict

        modelName: name of the model that will be used for compilation
        @type str

        modelPath: path to the generated model specific files @type str

        modelSwigPath: path to the generated swig files @type str

        allow_reinit_fixpar_initcond: indicates whether reinitialization of
        initial states depending on fixedParmeters is allowed for this model
        @type bool
    """

    def __init__(
            self,
            ode_model,
            outdir=None,
            verbose=False,
            assume_pow_positivity=False,
            compiler=None
    ):
        """Generate AMICI C++ files for the ODE provided to the constructor.

        Arguments:
            ode_model: ODE definition @type ODEModel

            outdir: see sbml_import.setPaths() @type str

            verbose: more verbose output if True @type bool

            assume_pow_positivity: if set to true, a special pow function is
            used to avoid problems with state variables that may become
            negative due to numerical errors @type bool

            compiler: distutils/setuptools compiler selection to build the
            python extension @type str

        Returns:

        Raises:

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

        # To only generate a subset of functions, apply subselection here
        self.functions = copy.deepcopy(functions)

        self.allow_reinit_fixpar_initcond = False

    def generateModelCode(self):
        """Generates the native C++ code for the loaded model

        Arguments:

        Returns:

        Raises:

        """
        self._prepareModelFolder()
        self._generateCCode()

    def compileModel(self):
        """Compiles the generated code it into a simulatable module

        Arguments:

        Returns:

        Raises:

        """
        self._compileCCode(compiler=self.compiler, verbose=self.verbose)

    def _prepareModelFolder(self):
        """Remove all files from the model folder.

        Arguments:

        Returns:

        Raises:

        """
        for file in os.listdir(self.modelPath):
            file_path = os.path.join(self.modelPath, file)
            if os.path.isfile(file_path):
                os.remove(file_path)

    def _generateCCode(self):
        """Create C++ code files for the model based on ODEExporter.model

        Arguments:

        Returns:

        Raises:

        """
        for function in self.functions.keys():
            self._writeFunctionFile(function)

        for name in self.model.symNames():
            self._writeIndexFiles(name)

        self._writeWrapfunctionsCPP()
        self._writeWrapfunctionsHeader()
        self._writeModelHeader()
        self._writeCMakeFile()
        self._writeSwigFiles()
        self._writeModuleSetup()

        shutil.copy(os.path.join(amiciSrcPath, 'main.template.cpp'),
                    os.path.join(self.modelPath, 'main.cpp'))


    def _compileCCode(self, verbose=False, compiler=None):
        """Compile the generated model code

        Arguments:
            verbose: Make model compilation verbose @type bool

            compiler: distutils/setuptools compiler selection to build the
            python extension @type str

        Returns:

        Raises:

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

    def _writeIndexFiles(self, name):
        """Write index file for a symbolic array.

        Arguments:
            name: key in self.model._syms for which the respective file should
            be written @type str

        Returns:

        Raises:


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

    def _writeFunctionFile(self, function):
        """Generate equations and write the C++ code for the function `function`.

        Arguments:
            function: name of the function to be written (see self.functions)
            @type str

        Returns:

        Raises:

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
            '#include <cmath>',
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
        body = self._getFunctionBody(function, symbol)
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

    def _getFunctionBody(self, function, symbol):
        """Generate C++ code for body of function `function`.

        Arguments:
        function: name of the function to be written (see self.functions)
        @type str

        symbol: symbolic defintion of the function body @type str
        symengine.DenseMatrix

        Returns:

        Raises:

        """
        lines = []

        if symbol.size == 0:
            return lines

        if function == 'sx0_fixedParameters':
            # here we only want to overwrite values where x0_fixedParameters
            # was applied
            lines.append(' ' * 4 + 'switch(ip) {')
            for ipar in range(self.model.np()):
                lines.append(' ' * 8 + f'case {ipar}:')
                for index, formula in enumerate(
                        self.model.eq('x0_fixedParameters')
                ):
                    if formula != 0 and formula != 0.0:
                        lines.append(' ' * 12 + f'{function}[{index}] = '
                                                f'{symbol[index, ipar]};')
                lines.append(' ' * 12 + 'break;')
            lines.append('}')
        elif function in sensi_functions:
            lines.append(' '*4 + 'switch(ip) {')
            for ipar in range(self.model.np()):
                lines.append(' ' * 8 + f'case {ipar}:')
                lines += self._getSymLines(symbol[:, ipar], function, 12)
                lines.append(' ' * 12 + 'break;')
            lines.append('}')
        elif function in multiobs_functions:
            lines.append(' '*4 + 'switch(iy) {')
            for iobs in range(self.model.ny()):
                lines.append(' ' * 8 + f'case {iobs}:')
                lines += self._getSymLines(symbol[:, iobs], function, 12)
                lines.append(' ' * 12 + 'break;')
            lines.append('}')
        else:
            if function in ['JSparse', 'JSparseB']:
                rowVals = self.model.rowval(function)
                colPtrs = self.model.colptr(function)
                lines += self._getSparseSymLines(
                    symbol, rowVals, colPtrs, function, 4
                )
            else:
                lines += self._getSymLines(symbol, function, 4)

        return [line for line in lines if line]

    def _writeWrapfunctionsCPP(self):
        """Write model-specific 'wrapper' file (wrapfunctions.cpp).

        Arguments:

        Returns:

        Raises:

        """
        templateData = {'MODELNAME': self.modelName}
        applyTemplate(
            os.path.join(amiciSrcPath, 'wrapfunctions.template.cpp'),
            os.path.join(self.modelPath, 'wrapfunctions.cpp'),
            templateData
        )

    def _writeWrapfunctionsHeader(self):
        """Write model-specific header file (wrapfunctions.h).

        Arguments:

        Returns:

        Raises:

        """
        templateData = {'MODELNAME': str(self.modelName)}
        applyTemplate(
            os.path.join(amiciSrcPath, 'wrapfunctions.ODE_template.h'),
            os.path.join(self.modelPath, 'wrapfunctions.h'),
            templateData
        )

    def _writeModelHeader(self):
        """Write model-specific header file (MODELNAME.h).

        Arguments:

        Returns:

        Raises:

        """

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
                self._getSymbolNameInitializerList('p'),
            'STATE_NAMES_INITIALIZER_LIST':
                self._getSymbolNameInitializerList('x'),
            'FIXED_PARAMETER_NAMES_INITIALIZER_LIST':
                self._getSymbolNameInitializerList('k'),
            'OBSERVABLE_NAMES_INITIALIZER_LIST':
                self._getSymbolNameInitializerList('y'),
            'PARAMETER_IDS_INITIALIZER_LIST':
                self._getSymbolIDInitializerList('p'),
            'STATE_IDS_INITIALIZER_LIST':
                self._getSymbolIDInitializerList('x'),
            'FIXED_PARAMETER_IDS_INITIALIZER_LIST':
                self._getSymbolIDInitializerList('k'),
            'OBSERVABLE_IDS_INITIALIZER_LIST':
                self._getSymbolIDInitializerList('y'),
            'REINIT_FIXPAR_INITCOND':
                'true' if self.allow_reinit_fixpar_initcond else
                'false',
        }
        applyTemplate(
            os.path.join(amiciSrcPath,'model_header.ODE_template.h'),
            os.path.join(self.modelPath,f'{self.modelName}.h'),
            templateData
        )

    def _getSymbolNameInitializerList(self, name):
        """Get SBML name initializer list for vector of names for the given model
        entity

        Arguments:
            name: any key present in self.model._syms @type str

        Returns:
        Template initializer list of names

        Raises:

        """
        return '\n'.join(
            [f'"{symbol}",' for symbol in self.model.name(name)]
        )

    def _getSymbolIDInitializerList(self, name):
        """Get C++ initializer list for vector of names for the given model entity

        Arguments:
            name: any key present in self.model._syms @type str

        Returns:
        Template initializer list of ids

        Raises:

        """
        return '\n'.join(
            [f'"{symbol}",' for symbol in self.model.sym(name)]
        )

    def _writeCMakeFile(self):
        """Write CMake CMakeLists.txt file for this model.

        Arguments:

        Returns:

        Raises:

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

    def _writeSwigFiles(self):
        """Write SWIG interface files for this model.

        Arguments:

        Returns:

        Raises:

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

    def _writeModuleSetup(self):
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

        applyTemplate(
            os.path.join(amiciModulePath, '__init__.template.py'),
            os.path.join(self.modelPath, self.modelName, '__init__.py'),
            templateData
        )

    def _getSymLines(self, symbols, variable, indentLevel):
        """Generate C++ code for assigning symbolic terms in symbols to C++ array
        `variable`.

        Arguments:
            symbols: vectors of symbolic terms @type list

            variable: name of the C++ array to assign to @type str

            indentLevel: indentation level (number of leading blanks) @type int

        Returns:
            C++ code as list of lines

        Raises:

        """

        lines = [' ' * indentLevel + f'{variable}[{index}] = '
                                     f'{self._printWithException(math)};'
                 if not (math == 0 or math == 0.0)
                 else ''
                 for index, math in enumerate(symbols)]

        try:
            lines.remove('')
        except:
            pass

        return lines

    def _getSparseSymLines(
            self, symbolList, RowVals, ColPtrs, variable, indentLevel
    ):
        """Generate C++ code for assigning sparse symbolic matrix to a C++ array
        `variable`.

        Arguments:
            symbolList: symbolic terms @type int

            RowVals: row indices of each nonzero entry (see CVODES SlsMat
            documentation for details) @type list

            ColPtrs: indices of the first column entries (see CVODES SlsMat
            documentation for details) @type list

            variable: name of the C++ array to assign to @type str

            indentLevel: indentation level (number of leading blanks) @type int

        Returns:
        C++ code as list of lines

        Raises:

        """
        lines = [
            ' ' * indentLevel + f'{variable}->indexvals[{index}] = '
                                f'{self._printWithException(math)};'
            for index, math in enumerate(RowVals)
        ]

        lines.extend(
            [' ' * indentLevel + f'{variable}->indexptrs[{index}] = '
                                 f'{self._printWithException(math)};'
             for index, math in enumerate(ColPtrs)]
        )
        lines.extend(
            [' ' * indentLevel + f'{variable}->data[{index}] = '
                                 f'{self._printWithException(math)};'
             for index, math in enumerate(symbolList)]
        )

        return lines

    def _printWithException(self, math):
        """Generate C++ code for a symbolic expression

        Arguments:
            math: symbolic expression @type symengine.Basic

        Returns:
        C++ code for the specified expression

        Raises:
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
        """Set output paths for the model and create if necessary

        Arguments:
            output_dir: relative or absolute path where the generated model
            code is to be placed. will be created if does not exists. @type str

        Returns:

        Raises:

        """
        self.modelPath = os.path.abspath(output_dir)
        self.modelSwigPath = os.path.join(self.modelPath, 'swig')

        for directory in [self.modelPath, self.modelSwigPath]:
            if not os.path.exists(directory):
                os.makedirs(directory)

    def setName(self, modelName):
        """Sets the model name

        Arguments:
            modelName: name of the model (must only contain valid filename
            characters) @type str

        Returns:

        Raises:

        """
        self.modelName = modelName

def getSymbolicDiagonal(matrix):
    """Get symbolic matrix with diagonal of matrix `matrix`.

    Arguments:
        matrix: Matrix from which to return the diagonal
        @type symengine.DenseMatrix

    Returns:
    A Symbolic matrix with the diagonal of `matrix`.

    Raises:
    Exception: The provided matrix was not square
    """
    if not matrix.cols == matrix.rows:
        raise Exception('Provided matrix is not square!')

    diagonal = [matrix[index,index] for index in range(matrix.cols)]

    return sp.DenseMatrix(diagonal)

class TemplateAmici(Template):
    """Template format used in AMICI (see string.template for more details).

    Attributes:
        delimiter: delimiter that identifies template variables @type str

    """
    delimiter = 'TPL_'


def applyTemplate(sourceFile,targetFile,templateData):
    """Load source file, apply template substitution as provided in templateData and save as targetFile.

    Arguments:
        sourceFile: relative or absolute path to template file @type str

        targetFile: relative or absolute path to output file @type str

        templateData: template keywords to substitute (key is template
        variable without TemplateAmici.delimiter) @type dict

    Returns:

    Raises:

    """
    with open(sourceFile) as filein:
        src = TemplateAmici(filein.read())
    result = src.safe_substitute(templateData)
    with open(targetFile, 'w') as fileout:
        fileout.write(result)
