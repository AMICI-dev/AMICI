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

from symengine.printing import CCodePrinter
from string import Template

from . import amiciSwigPath, amiciSrcPath, amiciModulePath

"""
prototype for generated C++ functions, keys are the names of functions

signature: str
    defines the argument part of the function signature, input variables
    should have a const flag

assume_pow_positivity: bool
    identifies the functions on which assume_pow_positivity will have an
    effect when specified during model generation. generally these are
    functions that are used for solving the ODE, where negative values may
    negatively affect convergence of the integration algorithm

sparse: bool
    specifies whether the result of this function will be stored in sparse
    format. sparse format means that the function will only return an array of
    nonzero values and not a full matrix.
"""
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
            ' const realtype *h)',
    },
    'dydp': {
        'signature':
            '(double *dydp, const realtype t, const realtype *x,'
            ' const realtype *p, const realtype *k, const realtype *h,'
            ' const int ip)',
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
    if 'const int ip' in functions[function]['signature']
    and function is not 'sxdot'
]

multiobs_functions = [
    function for function in functions
    if 'const int iy' in functions[function]['signature']
]

def var_in_function_signature(name, varname):
    """
    Checks if the values for a symbolic variable is passed in the signature
    of a function

    Arguments:
    ----------
    name: str
        name of the function

    varname: str
        name of the symbolic variable

    Returns:
    ----------

    Raises:
    ----------

    """
    varname = varname.replace('sparse', '')
    return name in functions \
           and re.search(
                    f'const (realtype|double) \*{varname}[0]*[,)]+',
                    functions[name]['signature']
                )

class ModelQuantity:
    """
    Base class for model components

    Attributes:
    ----------

    _identifier: sympy.Symbol
        unique identifier of the quantity

    _name: str
        individual name of the quantity (does not need to be unique)

    _value: sympy.Symbol or float
        either formula (of the initial value) or numeric value

    """
    def __init__(self, identifier, name, value):
        """
        Create a new ModelQuantity instance.

        Arguments:
        ----------
        identifier: sympy.Symbol
            unique identifier of the quantity

        name: str
            individual name of the quantity (does not need to be unique)

        value: symengine.Basic or float
            either formula, numeric value or initial value

        Returns:
        ----------
        ModelQuantity instance

        Raises:
        ----------
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

        if isinstance(value, sp.RealNumber):
            value = float(value)
        if not isinstance(value, sp.Basic) and not isinstance(value, float):
            raise TypeError(f'value must be sympy.Symbol or float, was '
                            f'{type(value)}')
        self._value = value

    def __repr__(self):
        return str(self._identifier)


class State(ModelQuantity):
    """
    A State variable defines an entity that evolves with time according to the
    provided time derivative, abbreviated by `x`

    Attributes:
    ----------
    _identifier: sympy.Symbol
        unique identifier of the quantity

    _name: str
        individual name of the quantity (does not need to be unique)

    _value: symengine.Basic
        initial value

    _dt: symengine.Basic
        time derivative

    """
    def __init__(self, identifier, name, value, dt):
        """
        Create a new State instance. Extends ModelQuantity.__init__ by dt

        Arguments:
        ----------
        identifier: sympy.Symbol
            unique identifier of the quantity

        name: str
            individual name of the quantity (does not need to be unique)

        value: symengine.Basic
            initial value

        dt: sympy.Symbol
            time derivative

        Returns:
        ----------
        ModelQuantity instance

        Raises:
        ----------
        TypeError:
            is thrown if input types do not match documented types
        """
        super(State, self).__init__(identifier, name, value)
        if not isinstance(dt, sp.Basic):
            raise TypeError(f'dt must be sympy.Symbol, was '
                            f'{type(dt)}')
        self._dt = dt


class Observable(ModelQuantity):
    """
    An Observable links model simulations to experimental measurements,
    abbreviated by `y`

    Attributes:
    ----------
    _identifier: sympy.Symbol
        unique identifier of the quantity

    _name: str
        individual name of the quantity (does not need to be unique)

    _value: symengine.Basic
        formula

    """
    pass


class SigmaY(ModelQuantity):
    """
    A Standard Deviation SigmaY rescales the distance between simulations and
    measurements when computing residuals, abbreviated by `sigmay`

    Attributes:
    ----------
    _identifier: sympy.Symbol
        unique identifier of the quantity

    _name: str
        individual name of the quantity (does not need to be unique)

    _value: symengine.Basic
        formula

    """
    pass


class Expression(ModelQuantity):
    """
    An Expressions is a recurring elements in symbolic formulas. Specifying
    this may yield more compact expression which may lead to substantially
    shorter model compilation times, but may also reduce model simulation time,
    abbreviated by `w`

    Attributes:
    ----------
    _identifier: sympy.Symbol
        unique identifier of the quantity

    _name: str
        individual name of the quantity (does not need to be unique)

    _value: symengine.Basic
        formula

    """
    pass


class Parameter(ModelQuantity):
    """
    A Parameter is a free variable in the model with respect to which
    sensitivites may be computed, abbreviated by `p`

    Attributes:
    ----------
    _identifier: sympy.Symbol
        unique identifier of the quantity

    _name: str
        individual name of the quantity (does not need to be unique)

    _value: float
        numeric value

    """
    pass


class Constant(ModelQuantity):
    """
    A Constant is a fixed variable in the model with respect to which
    sensitivites cannot be computed, abbreviated by `k`

    Attributes:
    ----------
    _identifier: sympy.Symbol
        unique identifier of the quantity

    _name: str
        individual name of the quantity (does not need to be unique)

    _value: float
        numeric value

    """
    pass


class LogLikelihood(ModelQuantity):
    """
    A LogLikelihood defines the distance between measurements and
    experiments for a particular observable. The final LogLikelihood value
    in the simulation will be the sum of all specified LogLikelihood
    instances evaluated at all timepoints, abbreviated by `Jy`

    Attributes:
    ----------
    _identifier: sympy.Symbol
        unique identifier of the quantity

    _name: str
        individual name of the quantity (does not need to be unique)

    _value: symengine.Basic
        formula

    """
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
    """
    An ODEModel defines an Ordinay Differential Equation as set of
    ModelQuantities. This class provides general purpose interfaces to
    compute arbitrary symbolic derivatives that are necessary for model
    simulation or sensitivity computation

    Attributes:
    ----------
    _states: list
        list of State instances

    _observables: list
        list of Observable instances

    _sigmays: list
        list of SigmaY instances

    _parameters: list
        list of Parameter instances

    _loglikelihoods: list
        list of LogLikelihood instances

    _expressions: list
        list of Expression instances

    _symboldim_funs: dict
        define functions that compute model dimensions, these are functions
        as the underlying symbolic expressions have not been populated at
        compile time

    _eqs: dict
        carries symbolic formulas of the symbolic variables of the model

    _sparseeq: dict
        carries linear list of all symbolic formulas for sparsified
        variables

    _vals: dict
        carries numeric values of symbolic identifiers of the symbolic
        variables of the model

    _names: dict
        carries names of symbolic identifiers of the symbolic variables of the
        model

    _syms: dict
        carries symbolic identifiers of the symbolic variables of the model

    _sparsesyms: dict
        carries linear list of all symbolic identifiers for sparsified
        variables

    _colptrs: dict
        carries column pointers for sparsified variables
        see SlsMat definition in CVODES for more details about ColPtrs

    _rowvals: dict
        carries row values for sparsified variables
        see SlsMat definition in CVODES for more details about RowVals

    _equation_prototype: dict
        defines the attribute from which an equation should be generated via
        list comprehension (see generateEquation)

    _variable_prototype: dict
        defines the attribute from which a variable should be generated via
        list comprehension (see generateSymbol)

    _value_prototype: dict
        defines the attribute from which a value should be generated via
        list comprehension (see generateValue)

    _total_derivative_prototypes: dict
        defines how a total derivative equation is computed for an equation,
        key defines the name and values should be arguments for
        ODEModel.totalDerivative

    _multiplication_prototypes: dict
        defines how a multiplication equation is computed for an equation,
        key defines the name and values should be arguments for
        ODEModel.multiplication

    _lock_total_derivative: bool
        set this to true when computing a total derivative from a partial
        derivative call to enforce a partial derivative in the next recursion.
        prevents infinite recursion

    """
    def __init__(self):
        """
        Create a new ODEModel instance.

        Arguments:
        ----------

        Returns:
        ----------

        Raises:
        ----------

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

    def import_from_sbml_importer(self, si, flux_as_expressions=True):
        """
        Imports a model specification from a amici.SBMLImporter instance.

        Arguments:
        ----------
        si: amici.SBMLImporter
            imported SBML model

        flux_as_expressions: bool
            defines whether fluxes should be used as Expressions

        Returns:
        ----------


        Raises:
        ----------

        """

        self.symbols = copy.copy(si.symbols)
        if flux_as_expressions:
            # setting these equations prevents native equation generation
            self._eqs['dxdotdw'] = si.stoichiometricMatrix
            self._eqs['w'] = si.fluxVector
            self.symbols['species']['dt'] = \
                si.stoichiometricMatrix * self.sym('w')
        else:
            self.symbols['species']['dt'] = \
                si.stoichiometricMatrix * si.fluxVector

        for symbol in [s for s in self.symbols if s != 'my']:
            # transform dict of lists into a list of dicts
            protos = [dict(zip(self.symbols[symbol], t))
                      for t in zip(*self.symbols[symbol].values())]
            for proto in protos:
                self.add_component(symbol_to_type[symbol](**proto))

        self._generateBasicVariables()

    def add_component(self, component):
        """
        Adds a new ModelQuantity to the model.

        Arguments:
        ----------
        component: ModelQuantity
            model quantity to be added

        Returns:
        ----------

        Raises:
        ----------

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
        """
        Number of states.

        Arguments:
        ----------

        Returns:
        ----------
        number of state variable symbols

        Raises:
        ----------

        """
        return len(self.sym('x'))

    def ny(self):
        """
        Number of Observables.

        Arguments:
        ----------

        Returns:
        ----------
        number of observable symbols

        Raises:
        ----------

        """
        return len(self.sym('y'))

    def nk(self):
        """
        Number of Constants.

        Arguments:
        ----------

        Returns:
        ----------
        number of constant symbols

        Raises:
        ----------

        """
        return len(self.sym('k'))

    def np(self):
        """
        Number of Parameters.

        Arguments:
        ----------

        Returns:
        ----------
        number of parameter symbols

        Raises:
        ----------

        """
        return len(self.sym('p'))

    def sym(self, name):
        """
        Returns (and constructs if necessary) the identifiers for a symbolic
        entity.

        Arguments:
        ----------
        name: str
            name of the symbolic variable

        Returns:
        ----------
        symengine.DenseMatrix containing the symbolic identifiers

        Raises:
        ----------

        """
        if name not in self._syms:
            self._generateSymbol(name)
        return self._syms[name]

    def sparsesym(self, name):
        """
        Returns (and constructs if necessary) the sparsified identifiers for a
        sparsified symbolic variable.

        Arguments:
        ----------
        name: str
            name of the symbolic variable

        Returns:
        ----------
        linearized symengine.DenseMatrix containing the symbolic identifiers

        Raises:
        ----------

        """
        if name not in sparse_functions:
            raise Exception(f'{name} is not marked as sparse')
        if name not in self._sparsesyms:
            self._generateSparseSymbol(name)
        return self._sparsesyms[name]

    def eq(self, name):
        """
        Returns (and constructs if necessary) the formulas for a symbolic
        entity.

        Arguments:
        ----------
        name: str
            name of the symbolic variable

        Returns:
        ----------
        symengine.DenseMatrix containing the symbolic identifiers

        Raises:
        ----------

        """
        if name not in self._eqs:
            self._computeEquation(name)
        return self._eqs[name]

    def sparseeq(self, name):
        """
        Returns (and constructs if necessary) the sparsified formulas for a
        sparsified symbolic variable.

        Arguments:
        ----------
        name: str
            name of the symbolic variable

        Returns:
        ----------
        linearized symengine.DenseMatrix containing the symbolic formulas

        Raises:
        ----------

        """
        if name not in sparse_functions:
            raise Exception(f'{name} is not marked as sparse')
        if name not in self._sparseeqs:
            self._generateSparseSymbol(name)
        return self._sparseeqs[name]

    def colptr(self, name):
        """
        Returns (and constructs if necessary) the column pointers for
        a sparsified symbolic variable.

        Arguments:
        ----------
        name: str
            name of the symbolic variable

        Returns:
        ----------
        symengine.DenseMatrix containing the column pointers

        Raises:
        ----------

        """
        if name not in sparse_functions:
            raise Exception(f'{name} is not marked as sparse')
        if name not in self._sparseeqs:
            self._generateSparseSymbol(name)
        return self._colptrs[name]

    def rowval(self, name):
        """
        Returns (and constructs if necessary) the row values for a sparsified
        symbolic variable.

        Arguments:
        ----------
        name: str
            name of the symbolic variable

        Returns:
        ----------
        symengine.DenseMatrix containing the row values

        Raises:
        ----------

        """
        if name not in sparse_functions:
            raise Exception(f'{name} is not marked as sparse')
        if name not in self._sparseeqs:
            self._generateSparseSymbol(name)
        return self._rowvals[name]

    def val(self, name):
        """
        Returns (and constructs if necessary) the numeric values of a symbolic
        entity

        Arguments:
        ----------
        name: str
            name of the symbolic variable

        Returns:
        ----------
        list containing the numeric values

        Raises:
        ----------

        """
        if name not in self._vals:
            self._generateValue(name)
        return self._vals[name]

    def name(self, name):
        """
        Returns (and constructs if necessary) the names of a symbolic variable

        Arguments:
        ----------
        name: str
            name of the symbolic variable

        Returns:
        ----------
        list of names

        Raises:
        ----------

        """
        if name not in self._names:
            self._generateName(name)
        return self._names[name]

    def _generateSymbol(self, name):
        """
        Generates the symbolic identifiers for a symbolic variable

        Arguments:
        ----------
        name: str
            name of the symbolic variable

        Returns:
        ----------

        Raises:
        ----------

        """
        if name in self._variable_prototype:
            component = self._variable_prototype[name]
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

    def _generateBasicVariables(self):
        """
        Generates the symbolic identifiers for all variables in
        ODEModel.variable_prototype

        Arguments:
        ----------

        Returns:
        ----------

        Raises:
        ----------

        """
        for var in self._variable_prototype:
            self._generateSymbol(var)

    def _generateSparseSymbol(self, name):
        """
        Generates the sparse symbolic identifiers, symbolic identifiers,
        sparse equations, column pointers and row values for a symbolic variable

        Arguments:
        ----------
        name: str
            name of the symbolic variable

        Returns:
        ----------

        Raises:
        ----------

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
        """
        computes the symbolic formula for a symbolic variable

        Arguments:
        ----------
        name: str
            name of the symbolic variable

        Returns:
        ----------

        Raises:
        ----------

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
                [comp.dt for comp in self._states]
            )

        elif name in ['sx0', 'sx0_fixedParameters']:
            self._derivative(name[1:], 'p', name=name)

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
            self._eqs[name] = self._eqs[name].transpose()

    def symNames(self):
        """
        Returns a list of names of generated symbolic variables

        Arguments:
        ----------

        Returns:
        ----------
        list of names

        Raises:
        ----------

        """
        return list(self._syms.keys())

    def _derivative(self, eq, var, name=None):
        """
        Creates a new symbolic variable according to a derivative

        Arguments:
        ----------
        eq: str
            name of the symbolic variable that defines the formula

        var: str
            name of the symbolic variable that defines the identifiers whith
            respect to which the derivatives are to be computed

        name: str
            name of resulting symbolic variable, default is d{eq}d{var}


        Returns:
        ----------

        Raises:
        ----------

        """
        if not name:
            name = f'd{eq}d{var}'

        # automatically detect chainrule
        if var_in_function_signature(eq, 'w') and \
                not self._lock_total_derivative:
            self._lock_total_derivative = True
            self._totalDerivative(name, eq, 'w', var)
            self._lock_total_derivative = False
            return

        # partial derivative
        self._eqs[name] = self.eq(eq).jacobian(self.sym(var))

    def _totalDerivative(self, name, eq, chainvar, var,
                         dydx_name=None, dxdz_name=None):
        """
        Creates a new symbolic variable according to a total derivative
        using the chain rule

        Arguments:
        ----------
        name: str
            name of resulting symbolic variable

        eq: str
            name of the symbolic variable that defines the formula

        chainvar: str
            name of the symbolic variable that defines the identifiers whith
            respect to which the chain rule is applied

        var: str
            name of the symbolic variable that defines the identifiers whith
            respect to which the derivatives are to be computed

        dydx_name: str
            defines the name of the symbolic variable that defines the
            derivative of the `eq` with respect to `chainvar`, default is
            d{eq}d{chainvar}

        dxdz_name: str
            defines the name of the symbolic variable that defines the
            derivative of the `chainvar` with respect to `var`, default is
            d{chainvar}d{var}


        Returns:
        ----------

        Raises:
        ----------

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

        self._eqs[name] = \
            variables['dydx']['sym'] * variables['dxdz']['sym'] + variables['dydz']['sym']

    def _multiplication(self, name, x, y,
                        transpose_x=False, sign=1):
        """
        Creates a new symbolic variable according to a multiplication

        Arguments:
        ----------
        name: str
            name of resulting symbolic variable, default is d{eq}d{var}

        x: str
            name of the symbolic variable that defines the first factor

        y: str
            name of the symbolic variable that defines the second factor

        transpose_x: bool
            indicates whether the first factor should be transposed before
            multiplication

        sign: int
            defines the sign of the product, should be +1 or -1


        Returns:
        ----------

        Raises:
        ----------

        """
        if sign not in [-1, 1]:
            raise TypeError(f'sign must be +1 or -1, was {sign}')

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

    def _equationFromComponent(self, name, component):
        """
        Generates the formulas of a symbolic variable from the attributes

        Arguments:
        ----------
        name: str
            name of resulting symbolic variable

        component: str
            name of the attribute


        Returns:
        ----------

        Raises:
        ----------

        """
        self._eqs[name] = sp.DenseMatrix(
            [comp.value for comp in getattr(self, component)]
        )

    def _generateValue(self, name):
        """
        Generates the numeric values of a symbolic variable from value
        prototypes

        Arguments:
        ----------
        name: str
            name of resulting symbolic variable

        Returns:
        ----------

        Raises:
        ----------

        """
        if name in self._value_prototype:
            component = self._value_prototype[name]
        else:
            raise Exception(f'No values for {name}')

        self._vals[name] = [comp.value for comp in getattr(self, component)]

    def _generateName(self, name):
        """
        Generates the names of a symbolic variable from variable prototypes or
        equation prototypes

        Arguments:
        ----------
        name: str
            name of resulting symbolic variable

        Returns:
        ----------

        Raises:
        ----------

        """
        if name in self._variable_prototype:
            component = self._variable_prototype[name]
        elif name in self._equation_prototype:
            component = self._equation_prototype[name]
        else:
            raise Exception(f'No names for {name}')

        self._names[name] = [comp.name for comp in getattr(self, component)]


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

    modelName: str
        name of the model that will be used for compilation

    modelPath: str
        path to the generated model specific files

    modelSwigPath: str
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

        output_dir: str
            see sbml_import.setPaths()

        verbose: bool
            more verbose output if True

        assume_pow_positivity: bool
            if set to true, a special pow function is used to avoid
            problems with state variables that may become negative due
            to numerical errors

        compiler: str
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

        # To only generate a subset of functions, apply subselection here
        self.functions = copy.deepcopy(functions)

        self.allow_reinit_fixpar_initcond = False

    def generateModelCode(self):
        """
        Generates the native C++ code for the loaded model

        Arguments:
        ----------

        Returns:
        ----------

        Raises:
        ----------

        """
        self._prepareModelFolder()
        self._generateCCode()

    def compileModel(self):
        """
        Compiles the generated code it into a simulatable module

        Arguments:
        ----------

        Returns:
        ----------

        Raises:
        ----------

        """
        self._compileCCode(compiler=self.compiler, verbose=self.verbose)

    def _prepareModelFolder(self):
        """
        Remove all files from the model folder.

        Arguments:
        ----------

        Returns:
        ----------

        Raises:
        ----------

        """
        for file in os.listdir(self.modelPath):
            file_path = os.path.join(self.modelPath, file)
            if os.path.isfile(file_path):
                os.remove(file_path)

    def _generateCCode(self):
        """
        Create C++ code files for the model based on ODEExporter.ODEModel.

        Arguments:
        ----------

        Returns:
        ----------

        Raises:
        ----------

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
        """
        Compile the generated model code

        Arguments:
        ----------
        verbose: bool
            Make model compilation verbose

        compiler: str
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

    def _writeIndexFiles(self, name):
        """
        Write index file for a symbolic array.

        Arguments:
        ----------
        name: str
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

    def _writeFunctionFile(self, function):
        """
        Generate equations and write the C++ code for the function `function`.

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
        """
        Generate C++ code for body of function `function`.

        Arguments:
        ----------
        function: str
            name of the function to be written (see self.functions)

        symbol: symengine.DenseMatrix
            symbolic defintion of the function body

        Returns:
        ----------

        Raises:
        ----------

        """
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
                        lines.append(' ' * 12 + f'{function}[{index}] = 0.0;')
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

    def _writeWrapfunctionsHeader(self):
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

    def _writeModelHeader(self):
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
        """
        Get SBML name initializer list for vector of names for the given model
        entity

        Arguments:
        ----------
        name: str
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

    def _getSymbolIDInitializerList(self, name):
        """
        Get C++ initializer list for vector of names for the given model entity

        Arguments:
        ----------
        name: str
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

    def _writeCMakeFile(self):
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

    def _writeSwigFiles(self):
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

    def _writeModuleSetup(self):
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

    def _getSymLines(self, symbols, variable, indentLevel):
        """
        Generate C++ code for assigning symbolic terms in symbols to C++ array
        `variable`.

        Arguments:
        ----------
        symbols: list
            vectors of symbolic terms

        variable: str
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

        variable: str
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
        output_dir: str
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
        modelName: str
            name of the model (must only contain valid filename characters)

        Returns:
        ----------

        Raises:
        ----------

        """
        self.modelName = modelName

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
    delimiter: str
        delimiter that identifies template variables

    """
    delimiter = 'TPL_'


def applyTemplate(sourceFile,targetFile,templateData):
    """Load source file, apply template substitution as provided in templateData and save as targetFile.

    Arguments:
    ----------
    sourceFile: str
        relative or absolute path to template file

    targetFile: str
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
