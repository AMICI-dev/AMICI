""" @package amici.ode_export The C++ ODE export module for python
"""
#!/usr/bin/env python3

import sympy as sp
import re
import shutil
import subprocess
import sys
import os
import copy
import numbers
import itertools
try:
    import pysb
except ImportError:
    ## pysb import dummy
    pysb = None


from typing import Callable, Optional
from string import Template
import sympy.printing.ccode as ccode
from sympy.matrices.immutable import ImmutableDenseMatrix
from sympy.matrices.dense import MutableDenseMatrix

from . import (
    amiciSwigPath, amiciSrcPath, amiciModulePath, __version__, __commit__
)

## Template for model simulation main.cpp file
CXX_MAIN_TEMPLATE_FILE = os.path.join(amiciSrcPath, 'main.template.cpp')
## Template for model/swig/CMakeLists.txt
SWIG_CMAKE_TEMPLATE_FILE = os.path.join(amiciSwigPath,
                                        'CMakeLists_model.cmake')
## Template for model/CMakeLists.txt
MODEL_CMAKE_TEMPLATE_FILE = os.path.join(amiciSrcPath,
                                         'CMakeLists.template.cmake')

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
            '(realtype *J, const realtype t, const realtype *x, '
            'const realtype *p, const realtype *k, const realtype *h, '
            'const realtype *w, const realtype *dwdx)',
        'flags': ['assume_pow_positivity']
    },
    'JB': {
        'signature':
            '(realtype *JB, const realtype t, const realtype *x, '
            'const realtype *p, const realtype *k, const realtype *h, '
            'const realtype *xB, const realtype *w, const realtype *dwdx)',
        'flags': ['assume_pow_positivity']
    },
    'JDiag': {
        'signature':
            '(realtype *JDiag, const realtype t, const realtype *x, '
            'const realtype *p, const realtype *k, const realtype *h, '
            'const realtype *w, const realtype *dwdx)',
        'flags': ['assume_pow_positivity']
    },
    'JSparse': {
        'signature':
            '(realtype *JSparse, const realtype t, const realtype *x, '
            'const realtype *p, const realtype *k, const realtype *h, '
            'const realtype *w, const realtype *dwdx)',
        'flags': ['assume_pow_positivity', 'sparse']
    },
    'JSparseB': {
        'signature':
            '(realtype *JSparseB, const realtype t, const realtype *x, '
            'const realtype *p, const realtype *k, const realtype *h, '
            'const realtype *xB, const realtype *w, const realtype *dwdx)',
        'flags': ['assume_pow_positivity', 'sparse']
    },
    'Jy': {
        'signature':
            '(realtype *Jy, const int iy, const realtype *p, '
            'const realtype *k, const realtype *y, const realtype *sigmay, '
            'const realtype *my)',
    },
    'dJydsigmay': {
        'signature':
            '(realtype *dJydsigmay, const int iy, const realtype *p, '
            'const realtype *k, const realtype *y, const realtype *sigmay, '
            'const realtype *my)',
    },
    'dJydy': {
        'signature':
            '(realtype *dJydy, const int iy, const realtype *p, '
            'const realtype *k, const realtype *y, '
            'const realtype *sigmay, const realtype *my)',
        'flags': ['sparse']
    },
    'dwdp': {
        'signature':
            '(realtype *dwdp, const realtype t, const realtype *x, '
            'const realtype *p, const realtype *k, const realtype *h, '
            'const realtype *w, const realtype *tcl, const realtype *dtcldp)',
        'flags': ['assume_pow_positivity', 'sparse']
    },
    'dwdx': {
        'signature':
            '(realtype *dwdx, const realtype t, const realtype *x, '
            'const realtype *p, const realtype *k, const realtype *h, '
            'const realtype *w, const realtype *tcl)',
        'flags': ['assume_pow_positivity', 'sparse']
    },
    'dxdotdw': {
        'signature':
            '(realtype *dxdotdw, const realtype t, const realtype *x, '
            'const realtype *p, const realtype *k, const realtype *h, '
            'const realtype *w)',
        'flags': ['assume_pow_positivity', 'sparse']
    },
    'dxdotdp_explicit': {
        'signature':
            '(realtype *dxdotdp_explicit, const realtype t, '
            'const realtype *x, const realtype *p, '
            'const realtype *k, const realtype *h, '
            'const realtype *w)',
        'flags': ['assume_pow_positivity', 'sparse']
    },
    'dxdotdp_implicit': {
        'signature':
            '(realtype *dxdotdp_implicit, const realtype t, '
            'const realtype *x, const realtype *p, const realtype *k, '
            'const realtype *h, const int ip, const realtype *w)',
        'flags': ['assume_pow_positivity', 'sparse', 'dont_generate_body']
    },
    'dydx': {
        'signature':
            '(realtype *dydx, const realtype t, const realtype *x, '
            'const realtype *p, const realtype *k, const realtype *h, '
            'const realtype *w, const realtype *dwdx)',
    },
    'dydp': {
        'signature':
            '(realtype *dydp, const realtype t, const realtype *x, '
            'const realtype *p, const realtype *k, const realtype *h, '
            'const int ip, const realtype *w, const realtype *dwdp)',
    },
    'dsigmaydp': {
        'signature':
            '(realtype *dsigmaydp, const realtype t, const realtype *p, '
            'const realtype *k, const int ip)',
    },
    'sigmay': {
        'signature':
            '(realtype *sigmay, const realtype t, const realtype *p, '
            'const realtype *k)',
    },
    'w': {
        'signature':
            '(realtype *w, const realtype t, const realtype *x, '
            'const realtype *p, const realtype *k, '
            'const realtype *h, const realtype *tcl)',
        'flags': ['assume_pow_positivity']
    },
    'x0': {
        'signature':
            '(realtype *x0, const realtype t, const realtype *p, '
            'const realtype *k)',
    },
    'x0_fixedParameters': {
        'signature':
            '(realtype *x0_fixedParameters, const realtype t, '
            'const realtype *p, const realtype *k)',
    },
    'sx0': {
        'signature':
            '(realtype *sx0, const realtype t,const realtype *x, '
            'const realtype *p, const realtype *k, const int ip)',
    },
    'sx0_fixedParameters': {
        'signature':
            '(realtype *sx0_fixedParameters, const realtype t, '
            'const realtype *x0, const realtype *p, const realtype *k, '
            'const int ip)',
    },
    'xdot': {
        'signature':
            '(realtype *xdot, const realtype t, const realtype *x, '
            'const realtype *p, const realtype *k, const realtype *h, '
            'const realtype *w)',
        'flags': ['assume_pow_positivity']
    },
    'y': {
        'signature':
            '(realtype *y, const realtype t, const realtype *x, '
            'const realtype *p, const realtype *k, '
            'const realtype *h, const realtype *w)',
    },
    'x_rdata': {
        'signature':
            '(realtype *x_rdata, const realtype *x, const realtype *tcl)',
    },
    'total_cl': {
        'signature':
            '(realtype *total_cl, const realtype *x_rdata)',
    },
    'x_solver': {
        'signature':
            '(realtype *x_solver, const realtype *x_rdata)',
    },
}

## list of sparse functions
sparse_functions = [
    function for function in functions
    if 'sparse' in functions[function].get('flags', [])
]
## list of nobody functions
nobody_functions = [
    function for function in functions
    if 'dont_generate_body' in functions[function].get('flags', [])
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

        if isinstance(value, sp.RealNumber) \
                or isinstance(value, numbers.Number):
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

    def get_id(self):
        """ModelQuantity identifier

        Arguments:

        Returns:
        identifier of the ModelQuantity

        Raises:

        """
        return self._identifier

    def get_name(self):
        """ModelQuantity name

        Arguments:

        Returns:
        name of the ModelQuantity

        Raises:

        """
        return self._name

    def get_val(self):
        """ModelQuantity value

        Arguments:

        Returns:
        value of the ModelQuantity

        Raises:

        """
        return self._value


class State(ModelQuantity):
    """A State variable defines an entity that evolves with time according to
    the provided time derivative, abbreviated by `x`

    Attributes:
        conservation_law: algebraic formula that allows computation of this
        species according to a conservation law

    """

    conservation_law = None

    def __init__(self, identifier, name, value, dt):
        """Create a new State instance. Extends ModelQuantity.__init__ by dt

        Arguments:
            identifier: unique identifier of the state @type sympy.Symbol

            name: individual name of the state (does not need to be unique)
            @type str

            value: initial value @type symengine.Basic

            dt: time derivative @type symengine.Basic

        Raises:
        TypeError:
            is thrown if input types do not match documented types
        """
        super(State, self).__init__(identifier, name, value)
        if not isinstance(dt, sp.Basic):
            raise TypeError(f'dt must have type sympy.Basic, was '
                            f'{type(dt)}')

        self._dt = dt
        self.conservation_law = None

    def set_conservation_law(self, law):
        """Sets the conservation law of a state. If the a conservation law
        is set, the respective state will be replaced by an algebraic
        formula according to the respective conservation law.

        Arguments:
            law: linear sum of states that if added to this state remain
            constant over time

        Returns:

        Raises:
        TypeError:
            is thrown if law type does not match documented type
        """
        if not isinstance(law, sp.Basic):
            raise TypeError(f'conservation law must have type sympy.Basic, '
                            f'was {type(law)}')

        self.conservation_law = law

    def set_dt(self, dt):
        """Sets the time derivative

        Arguments:
            dt: time derivative @type symengine.Basic

        Returns:

        Raises:
        TypeError:
            is thrown if dt type does not match documented type
        """
        if not isinstance(dt, sp.Basic):
            raise TypeError(f'time derivative must have type sympy.Basic, '
                            f'was {type(dt)}')
        self._dt = dt

    def get_dt(self):
        """Sets the time derivative

        Arguments:

        Returns:
        time derivative @type symengine.Basic

        Raises:
        """
        return self._dt


class ConservationLaw(ModelQuantity):
    """ A conservation law defines the absolute the total amount of a
    (weighted) sum of states

    """
    def __init__(self, identifier, name, value):
        """Create a new ConservationLaw instance.

        Arguments:
            identifier: unique identifier of the ConservationLaw
            @type sympy.Symbol

            name: individual name of the ConservationLaw (does not need to be
            unique) @type str

            value: formula (sum of states) @type symengine.Basic

        Raises:
        TypeError:
            is thrown if input types do not match documented types
        """
        super(ConservationLaw, self).__init__(identifier, name, value)


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

        Raises:
        TypeError:
            is thrown if input types do not match documented types
        """
        super(Observable, self).__init__(identifier, name, value)


class SigmaY(ModelQuantity):
    """A Standard Deviation SigmaY rescales the distance between simulations
    and measurements when computing residuals, abbreviated by `sigmay`

    """
    def __init__(self, identifier, name, value):
        """Create a new Standard Deviation instance.

        Arguments:
            identifier: unique identifier of the Standard Deviation
            @type sympy.Symbol

            name: individual name of the Standard Deviation (does not need to
            be unique) @type str

            value: formula @type symengine.Basic

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

        _conservationlaws: list of ConservationLaw instances @type list

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

        _strippedsyms: carries symbolic identifiers that were stripped of
        additional class information @type dict

        _sparsesyms: carries linear list of all symbolic identifiers for
        sparsified variables  @type dict

        _colptrs: carries column pointers for sparsified variables. See
        SUNMatrixContent_Sparse definition in <sunmatrix/sunmatrix_sparse.h>
        @type dict

        _rowvals: carries row values for sparsified variables. See
        SUNMatrixContent_Sparse definition in <sunmatrix/sunmatrix_sparse.h>
        @type dict

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

        _simplify: If not None, this function will be used to simplify symbolic
         derivative expressions. Receives sympy expressions as only argument.
         To apply multiple simplifications, wrap them in a lambda expression.
         NOTE: This does currently not work with PySB symbols.
    """

    def __init__(self, simplify: Optional[Callable] = sp.powsimp):
        """Create a new ODEModel instance.

        Arguments:
            simplify: see ODEModel._simplify

        Raises:

        """
        self._states = []
        self._observables = []
        self._sigmays = []
        self._parameters = []
        self._constants = []
        self._loglikelihoods = []
        self._expressions = []
        self._conservationlaws = []
        self._symboldim_funs = {
            'sx': self.nx_solver,
            'v': self.nx_solver,
            'vB': self.nx_solver,
            'xB': self.nx_solver,
            'sigmay': self.ny,
        }
        self._eqs = dict()
        self._sparseeqs = dict()
        self._vals = dict()
        self._names = dict()
        self._syms = dict()
        self._strippedsyms = dict()
        self._sparsesyms = dict()
        self._colptrs = dict()
        self._rowvals = dict()

        self._equation_prototype = {
            'total_cl': '_conservationlaws',
            'x0': '_states',
            'y': '_observables',
            'Jy': '_loglikelihoods',
            'w': '_expressions',
            'sigmay': '_sigmays',
        }
        self._variable_prototype = {
            'tcl': '_conservationlaws',
            'x_rdata': '_states',
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
                'chainvars': ['w'],
                'var': 'x',
            },
            'sxdot': {
                'eq': 'xdot',
                'chainvars': ['x'],
                'var': 'p',
                'dydx_name': 'JSparse',
                'dxdz_name': 'sx',
            },
            'sx_rdata': {
                'eq': 'x_rdata',
                'chainvars': ['x'],
                'var': 'p',
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
            'xBdot': {
                'x': 'JB',
                'y': 'xB',
            },
            'dxdotdp_implicit': {
                'x': 'dxdotdw',
                'y': 'dwdp',
            },
        }

        self._lock_total_derivative = False
        self._simplify = simplify

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
        self._syms['w'] = sp.Matrix(
            [sp.Symbol(f'flux_r{idx}', real=True)
             for idx in range(len(si.fluxVector))]
        )
        self._eqs['dxdotdx'] = sp.zeros(si.stoichiometricMatrix.shape[0])
        if len(si.stoichiometricMatrix):
            symbols['species']['dt'] = \
                si.stoichiometricMatrix * self.sym('w')
        else:
            symbols['species']['dt'] = sp.zeros(
                *symbols['species']['identifier'].shape
            )

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
            Exception: invalid component type

        """
        for comp_type in [Observable, Expression, Parameter, Constant, State,
                          LogLikelihood, SigmaY, ConservationLaw]:
            if isinstance(component, comp_type):
                getattr(self, f'_{type(component).__name__.lower()}s').append(
                    component
                )
                return
        Exception(f'Invalid component type {type(component)}')

    def add_conservation_law(self, state, total_abundance, state_expr,
                             abundance_expr):
        """Adds a new conservation law to the model.

        Arguments:
            state: symbolic identifier of the state that should be replaced by
            the conservation law
            total_abundance: symbolic identifier of the total abundance
            state_expr: symbolic algebraic formula that replaces the the state
            abundance_expr: symbolic algebraic formula that computes the
            total abundance

        Returns:

        Raises:
            Exception if the specified state does not correspond to a valid
            state in the model

        """
        try:
            ix = [
                s.get_id()
                for s in self._states
            ].index(state)
        except ValueError:
            raise Exception(f'Specified state {state} was not found in the '
                            f'model states.')

        state_id = self._states[ix].get_id()

        self.add_component(
            Expression(state_id, str(state_id), state_expr)
        )

        self.add_component(
            ConservationLaw(
                total_abundance,
                f'total_{state_id}',
                abundance_expr
            )
        )

        self._states[ix].set_conservation_law(state_expr)

    def nx_rdata(self):
        """Number of states.

        Arguments:

        Returns:
        number of state variable symbols

        Raises:

        """
        return len(self.sym('x_rdata'))

    def nx_solver(self):
        """Number of states after applying conservation laws.

        Arguments:

        Returns:
        number of state variable symbols

        Raises:

        """
        return len(self.sym('x'))

    def ncl(self):
        """Number of conservation laws.

        Arguments:

        Returns:
        number of conservation laws

        Raises:

        """
        return self.nx_rdata()-self.nx_solver()

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

    def sym(self, name, stripped=False):
        """Returns (and constructs if necessary) the identifiers for a symbolic
        entity.

        Arguments:
            name: name of the symbolic variable @type str
            stripped: (optional) should additional class information be
            stripped from the symbolic variables? @type bool

        Returns:
        containing the symbolic identifiers @type symengine.DenseMatrix

        Raises:
        ValueError if stripped symbols not available
        """
        if name not in self._syms:
            self._generateSymbol(name)

        if stripped:
            if name not in self._variable_prototype:
                raise ValueError('Stripped symbols only available for '
                                 'variables from variable prototypes')
            return self._strippedsyms[name]
        else:
            return self._syms[name]

    def sparsesym(self, name):
        """Returns (and constructs if necessary) the sparsified identifiers for
        a sparsified symbolic variable.

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
            self._compute_equation(name)
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

    def colptrs(self, name):
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

    def rowvals(self, name):
        """Returns (and constructs if necessary) the row values for a
        sparsified symbolic variable.

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
        """Returns (and constructs if necessary) the numeric values of a
        symbolic entity

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
        """Returns (and constructs if necessary) the names of a symbolic
        variable

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
            self._syms[name] = sp.Matrix([
                comp.get_id()
                for comp in getattr(self, component)
            ])
            self._strippedsyms[name] = sp.Matrix([
                sp.Symbol(comp.get_name())
                for comp in getattr(self, component)
            ])
            if name == 'y':
                self._syms['my'] = sp.Matrix([
                    sp.Symbol(f'm{strip_pysb(comp.get_id())}', real=True)
                    for comp in getattr(self, component)
                ])
            return
        elif name == 'x':
            self._syms[name] = sp.Matrix([
                state.get_id()
                for state in self._states
                if state.conservation_law is None
            ])
            return
        elif name == 'dtcldp':
            self._syms[name] = sp.Matrix([
                [
                    sp.Symbol(f's{strip_pysb(tcl.get_id())}__'
                              f'{strip_pysb(par.get_id())}',
                              real=True)
                    for par in self._parameters
                ]
                for tcl in self._conservationlaws
            ])
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

        self._syms[name] = sp.Matrix([
            sp.Symbol(f'{name}{i}', real=True) for i in range(length)
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

    def get_appearance_counts(self, idxs):
        """Counts how often a state appears in the time derivative of
        another state and expressions for a subset of states

        Arguments:
            idxs: list of state indices for which counts are to be computed

        Returns:
            list of counts for the states ordered according to the provided
            indices

        Raises:

        """
        free_symbols_dt = list(itertools.chain.from_iterable(
            [
                str(symbol)
                for symbol in state.get_dt().free_symbols
            ]
            for state in self._states
        ))

        free_symbols_expr = list(itertools.chain.from_iterable(
            [
                str(symbol)
                for symbol in expr.get_val().free_symbols
            ]
            for expr in self._expressions
        ))

        return [
            free_symbols_dt.count(str(self._states[idx].get_id()))
            +
            free_symbols_expr.count(str(self._states[idx].get_id()))
            for idx in idxs
        ]

    def _generateSparseSymbol(self, name):
        """Generates the sparse symbolic identifiers, symbolic identifiers,
        sparse equations, column pointers and row values for a symbolic
        variable

        Arguments:
            name: name of the symbolic variable @type str

        Returns:

        Raises:

        """
        matrix = self.eq(name)

        if name == 'dJydy':
            # One entry per y-slice
            self._colptrs[name] = []
            self._rowvals[name] = []
            self._sparseeqs[name] = []
            self._sparsesyms[name] = []
            self._syms[name] = []
            base_index = 0
            for iy in range(self.ny()):
                symbolColPtrs, symbolRowVals, sparseList, symbolList, \
                    sparseMatrix = csc_matrix(matrix[iy, :], name,
                                              base_index=base_index)
                base_index += len(symbolList)
                self._colptrs[name].append(symbolColPtrs)
                self._rowvals[name].append(symbolRowVals)
                self._sparseeqs[name].append(sparseList)
                self._sparsesyms[name].append(symbolList)
                self._syms[name].append(sparseMatrix)
        else:
            symbolColPtrs, symbolRowVals, sparseList, symbolList, \
                sparseMatrix = csc_matrix(
                    matrix, name, pattern_only=name in nobody_functions
                )

            self._colptrs[name] = symbolColPtrs
            self._rowvals[name] = symbolRowVals
            self._sparseeqs[name] = sparseList
            self._sparsesyms[name] = symbolList
            self._syms[name] = sparseMatrix

    def _compute_equation(self, name):
        """computes the symbolic formula for a symbolic variable

        Arguments:
            name: name of the symbolic variable @type str

        Returns:

        Raises:

        """
        match_deriv = re.match(r'd([\w_]+)d([a-z_]+)', name)

        if name in self._equation_prototype:
            self._equationFromComponent(name, self._equation_prototype[name])

        elif name in self._total_derivative_prototypes:
            args = self._total_derivative_prototypes[name]
            args['name'] = name
            self._lock_total_derivative = True
            self._total_derivative(**args)
            self._lock_total_derivative = False

        elif name in self._multiplication_prototypes:
            args = self._multiplication_prototypes[name]
            args['name'] = name
            self._multiplication(**args)

        elif name == 'xdot':
            self._eqs[name] = sp.Matrix([
                s._dt for s in self._states
                if s.conservation_law is None
            ])

        elif name == 'x_rdata':
            self._eqs[name] = sp.Matrix([
                state.get_id()
                if state.conservation_law is None
                else state.conservation_law
                for state in self._states
            ])

        elif name == 'x_solver':
            self._eqs[name] = sp.Matrix([
                state.get_id()
                for state in self._states
                if state.conservation_law is None
            ])

        elif name == 'sx_solver':
            self._eqs[name] = sp.Matrix([
                self.sym('sx_rdata')[ix]
                for ix, state in enumerate(self._states)
                if state.conservation_law is None
            ])

        elif name == 'sx0':
            self._derivative(name[1:], 'p', name=name)

        elif name == 'sx0_fixedParameters':
            # deltax = -x+x0_fixedParameters if x0_fixedParameters>0 else 0
            # deltasx = -sx+dx0_fixedParametersdx*sx+dx0_fixedParametersdp
            # if x0_fixedParameters>0 else 0
            # sx0_fixedParameters = sx+deltasx =
            # dx0_fixedParametersdx*sx+dx0_fixedParametersdp
            if len(self.sym('p')):
                self._eqs[name] = \
                    self.eq('x0_fixedParameters').jacobian(self.sym('p'))
            else:
                self._eqs[name] = sp.zeros(
                    len(self.eq('x0_fixedParameters')),
                    len(self.sym('p'))
                )

            if len(self.sym('x')):
                dx0_fixedParametersdx = \
                    self.eq('x0_fixedParameters').jacobian(self.sym('x'))
            else:
                dx0_fixedParametersdx = sp.zeros(
                    len(self.eq('x0_fixedParameters')),
                    len(self.sym('x'))
                )

            if dx0_fixedParametersdx.is_zero is not True:
                for ip in range(self._eqs[name].shape[1]):
                    self._eqs[name][:,ip] += \
                        dx0_fixedParametersdx \
                        * self.sym('sx0') \

            for index, formula in enumerate(self.eq('x0_fixedParameters')):
                if formula == 0 or formula == 0.0:
                    # sp.simplify returns ImmutableDenseMatrix, if we need to
                    # change them, they need to be made mutable
                    if isinstance(self._eqs[name], ImmutableDenseMatrix):
                        self._eqs[name] = MutableDenseMatrix(self._eqs[name])
                    self._eqs[name][index, :] = \
                        sp.zeros(1, self._eqs[name].shape[1])

        elif name == 'JB':
            self._eqs[name] = -self.eq('J').transpose()

        elif name == 'JDiag':
            self._eqs[name] = getSymbolicDiagonal(self.eq('J'))

        elif name == 'x0_fixedParameters':
            k = self.sym('k')
            self._eqs[name] = sp.Matrix([
                eq
                # check if the equation contains constants
                if any([sym in eq.free_symbols for sym in k])
                # if not set to zero
                else 0.0
                for eq in self.eq('x0')
            ])

        elif name in ['JSparse', 'JSparseB']:
            self._eqs[name] = self.eq(name.replace('Sparse', ''))

        elif name == 'dtotal_cldx_rdata':
            # not correctly parsed in regex
            self._derivative('total_cl', 'x_rdata')

        elif name == 'dtcldx':
            # this is always zero
            self._eqs[name] = sp.zeros(self.ncl(), self.nx_solver())

        elif name == 'dtcldp':
            # force symbols
            self._eqs[name] = self.sym(name)

        elif name == 'dxdotdp_explicit':
            # force symbols
            self._derivative('xdot', 'p', name=name)

        elif match_deriv:
            self._derivative(match_deriv.group(1), match_deriv.group(2))

        else:
            raise Exception(f'Unknown equation {name}')

        if name in ['Jy', 'dydx']:
            # do not transpose if we compute the partial derivative as part of
            # a total derivative
            if not self._lock_total_derivative:
                self._eqs[name] = self._eqs[name].transpose()

        if self._simplify:
            self._eqs[name] = self._simplify(self._eqs[name])

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
        chainvars = []
        ignore_chainrule = {
            ('xdot', 'p'): 'w'  # has generic implementation in c++ code
        }
        for cv in ['w', 'tcl']:
            if var_in_function_signature(eq, cv) \
                    and not self._lock_total_derivative \
                    and var is not cv \
                    and min(self.sym(cv).shape) \
                    and (
                            (eq, var) not in ignore_chainrule
                            or ignore_chainrule[(eq, var)] != cv
                    ):
                chainvars.append(cv)

        if len(chainvars):
            self._lock_total_derivative = True
            self._total_derivative(name, eq, chainvars, var)
            self._lock_total_derivative = False
            return

        # this is the basic requirement check
        needs_stripped_symbols = eq == 'xdot' and var != 'x'

        # partial derivative
        if eq == 'Jy':
            eq = self.eq(eq).transpose()
        else:
            eq = self.eq(eq)

        if pysb is not None and needs_stripped_symbols:
            needs_stripped_symbols = not any(
                isinstance(sym, pysb.Component)
                for sym in eq.free_symbols
            )

        # now check whether we are working with energy_modeling branch
        # where pysb class info is already stripped
        # TODO: fixme as soon as energy_modeling made it to the main pysb
        #  branch
        sym_var = self.sym(var, needs_stripped_symbols)

        if min(eq.shape) and min(sym_var.shape) \
                and eq.is_zero is not True and sym_var.is_zero is not True:
            self._eqs[name] = eq.jacobian(sym_var)
        else:
            self._eqs[name] = sp.zeros(eq.shape[0], self.sym(var).shape[0])

    def _total_derivative(self, name, eq, chainvars, var,
                          dydx_name=None, dxdz_name=None):
        """Creates a new symbolic variable according to a total derivative
        using the chain rule

        Arguments:
            name: name of resulting symbolic variable @type str

            eq: name of the symbolic variable that defines the formula
            @type str

            chainvars: names of the symbolic variable that define the
            identifiers with respect to which the chain rules are applied
            @type list

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

        # initialize with partial derivative dydz without chain rule
        self._eqs[name] = copy.deepcopy(
            self.sym_or_eq(name, f'd{eq}d{var}')
        )

        for chainvar in chainvars:
            if dydx_name is None:
                dydx_name = f'd{eq}d{chainvar}'
            if dxdz_name is None:
                dxdz_name = f'd{chainvar}d{var}'

            dydx = self.sym_or_eq(name, dydx_name)
            dxdz = self.sym_or_eq(name, dxdz_name)
            # Save time for for large models if one multiplicand is zero,
            # which is not checked for by sympy
            if dydx.is_zero is not True and dxdz.is_zero is not True:
                if dxdz.shape[1] == 1 and \
                        self._eqs[name].shape[1] != dxdz.shape[1]:
                    for iz in range(self._eqs[name].shape[1]):
                        self._eqs[name][:, iz] += dydx * dxdz
                else:
                    self._eqs[name] += dydx * dxdz

    def sym_or_eq(self, name, varname):
        """Returns symbols or equations depending on whether a given
        variable appears in the function signature or not.

        Arguments:
            name: name of function for which the signature should be checked
            varname: name of the variable which should be contained in the
            function signature

        Returns:
        the variable symbols if the variable is part of the signature and
        the variable equations otherwise.

        Raises:

        """
        if var_in_function_signature(name, varname):
            return self.sym(varname)
        else:
            return self.eq(varname)

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

        yy = variables[y]

        if not xx.shape[0] or not yy.shape[1] or xx.is_zero is True or \
                yy.is_zero is True:
            self._eqs[name] = sp.zeros(xx.shape[0], yy.shape[1])
        else:
            self._eqs[name] = sign * xx * yy

    def _equationFromComponent(self, name, component):
        """Generates the formulas of a symbolic variable from the attributes

        Arguments:
            name: name of resulting symbolic variable @type str

            component: name of the attribute @type str

        Returns:

        Raises:

        """
        self._eqs[name] = sp.Matrix(
            [comp._value for comp in getattr(self, component)]
        )
        # flatten conservation laws in expressions
        if name == 'w':
            self._eqs[name] = self._eqs[name].subs(
                self.get_conservation_laws()
            )

    def get_conservation_laws(self):
        """ Returns a list of states with conservation law set

        Arguments:

        Returns:
        list of state identifiers

        Raises:

        """
        return [
            (state.get_id(), state.conservation_law)
            for state in self._states
            if state.conservation_law is not None
        ]

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

    def state_has_fixed_parameter_initial_condition(self, ix):
        """Checks whether the state at specified index has a fixed parameter
        initial condition

        Arguments:
            ix: state index

        Returns:
            boolean indicating if any of the initial condition free
            variables is contained in the model constants

        Raises:

        """
        ic = self._states[ix].get_val()
        if not isinstance(ic, sp.Basic):
            return False
        return any([
            fp in [c.get_id() for c in self._constants]
            for fp in ic.free_symbols
        ])

    def state_has_conservation_law(self, ix):
        """Checks whether the state at specified index has a conservation
        law set

        Arguments:
            ix: state index

        Returns:
            boolean indicating if conservation_law is not None

        Raises:

        """
        return self._states[ix].conservation_law is not None

    def state_is_constant(self, ix):
        """Checks whether the temporal derivative of the state is zero

        Arguments:
            ix: state index

        Returns:
            boolean indicating if constant over time

        Raises:

        """
        return self._states[ix].get_dt() == 0.0



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

        functions: carries C++ function signatures and other specifications
        @type dict

        modelName: name of the model that will be used for compilation
        @type str

        modelPath: path to the generated model specific files @type str

        modelSwigPath: path to the generated swig files @type str

        allow_reinit_fixpar_initcond: indicates whether reinitialization of
        initial states depending on fixedParameters is allowed for this model
        @type bool
    """

    def __init__(
            self,
            ode_model,
            outdir=None,
            verbose=False,
            assume_pow_positivity=False,
            compiler=None,
            allow_reinit_fixpar_initcond=True
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

            allow_reinit_fixpar_initcond: see ODEExporter

        Raises:

        """

        self.outdir = outdir
        self.verbose = verbose
        self.assume_pow_positivity = assume_pow_positivity
        self.compiler = compiler

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

        self.allow_reinit_fixpar_initcond = allow_reinit_fixpar_initcond

    def generateModelCode(self):
        """Generates the native C++ code for the loaded model and a Matlab
        script that can be run to compile a mex file from the C++ code

        Arguments:

        Returns:

        Raises:

        """
        self._prepareModelFolder()
        self._generateCCode()
        self._generateMCode()

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
            if 'dont_generate_body' not in \
                    self.functions[function].get('flags', []):
                self._writeFunctionFile(function)
            if function in sparse_functions:
                self._write_function_index(function, 'colptrs')
                self._write_function_index(function, 'rowvals')

        for name in self.model.symNames():
            self._writeIndexFiles(name)

        self._writeWrapfunctionsCPP()
        self._writeWrapfunctionsHeader()
        self._writeModelHeader()
        self._writeCMakeFile()
        self._writeSwigFiles()
        self._writeModuleSetup()

        shutil.copy(CXX_MAIN_TEMPLATE_FILE,
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

    def _generateMCode(self):
        """Create a Matlab script for compiling code files to a mex file

        Arguments:

        Returns:

        Raises:

        """
        # creating the code lines for the Matlab compile script
        lines = []

        # Events are not yet implemented. Once this is done, the variable nz
        # will have to be replaced by "self.model.nz()"
        nz = 0

        # Second order code is not yet implemented. Once this is done,
        # those variables will have to be replaced by
        # "self.model.<var>true()", or the corresponding "model.self.o2flag"
        nxtrue_rdata = self.model.nx_rdata()
        nytrue = self.model.ny()
        o2flag = 0

        # a preliminary comment
        lines.append('% This compile script was automatically created from Python SBML import.')
        lines.append('% If mex compiler is set up within MATLAB, it can be run from MATLAB ')
        lines.append('% in order to compile a mex-file from the Python generated C++ files.')
        lines.append('')

        # write the actual compiling code
        lines.append('''modelName = '{model_name}';'''.format(
            model_name=self.modelName))
        lines.append('''amimodel.compileAndLinkModel(modelName, '', [], [], [], []);''')
        lines.append('''amimodel.generateMatlabWrapper({nx}, {ny}, {np}, {nk}, {nz}, {o2flag}, [], ...
            ['simulate_' modelName '.m'], modelName, 'lin', 1, 1);'''.format(
            nx=nxtrue_rdata,
            ny=nytrue,
            np=self.model.np(),
            nk=self.model.nk(),
            nz=nz,
            o2flag=o2flag
            ))

        # write compile script (for mex)
        with open(os.path.join(self.modelPath, 'compileMexFile.m'), 'w') as fileout:
            fileout.write('\n'.join(lines))

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
                symbols = self.model.sym(name).T
        else:
            raise Exception('Unknown symbolic array')

        for index, symbol in enumerate(symbols):
            symbol_name = strip_pysb(symbol)
            lines.append(
                f'#define {symbol_name} {name}[{index}]'
            )

        with open(os.path.join(self.modelPath, f'{name}.h'), 'w') as fileout:
            fileout.write('\n'.join(lines))

    def _writeFunctionFile(self, function):
        """Generate equations and write the C++ code for the function
        `function`.

        Arguments:
            function: name of the function to be written (see self.functions)
            @type str

        Returns:

        Raises:

        """

        # first generate the equations to make sure we have everything we
        # need in subsequent steps
        if function in sparse_functions:
            symbol = self.model.sparseeq(function)
        elif not self.allow_reinit_fixpar_initcond \
                and function == 'sx0_fixedParameters':
            # Not required. Will create empty function body.
            symbol = sp.Matrix()
        else:
            symbol = self.model.eq(function)

        # function header
        lines = [
            '#include "amici/symbolic_functions.h"',
            '#include "amici/defines.h"',
            '#include "sundials/sundials_types.h"',
            '#include <cmath>',
            ''
        ]

        # function signature
        signature = self.functions[function]['signature']

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
        if self.assume_pow_positivity and 'assume_pow_positivity' \
                in self.functions[function].get('flags', []):
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

    def _write_function_index(self, function, indextype):
        """Generate equations and write the C++ code for the function
        `function`.

        Arguments:
            function: name of the function to be written (see self.functions)
            @type str
            indextype: type of index {'colptrs', 'rowvals'}

        Returns:

        Raises:

        """

        if indextype == 'colptrs':
            values = self.model.colptrs(function)
        elif indextype == 'rowvals':
            values = self.model.rowvals(function)
        else:
            raise ValueError('Invalid value for type, must be colptr or '
                             'rowval')

        # function signature
        if function in multiobs_functions:
            signature = f'(sunindextype *{indextype}, int index)'
        else:
            signature = f'(sunindextype *{indextype})'

        lines = [
            '#include "sundials/sundials_types.h"',
            '',
            f'void {function}_{indextype}_{self.modelName}{signature}{{',
        ]
        if function in multiobs_functions:
            # list of index vectors
            cases = {switch_case: [' ' * 4 + f'{indextype}[{index}] = {value};'
                     for index, value in enumerate(idx_vector)]
                         for switch_case, idx_vector in enumerate(values)}
            lines.extend(getSwitchStatement('index', cases, 1))
        else:
            # single index vector
            lines.extend(
                [' ' * 4 + f'{indextype}[{index}] = {value};'
                 for index, value in enumerate(values)]
            )
        lines.append('}')
        with open(os.path.join(
                self.modelPath,
                f'{self.modelName}_{function}_{indextype}.cpp'
        ), 'w') as fileout:
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

        if len(symbol) == 0 or (isinstance(symbol, sp.Matrix)
                                and min(symbol.shape) == 0):
            return lines

        if not self.allow_reinit_fixpar_initcond \
                and function in ['sx0_fixedParameters', 'x0_fixedParameters']:
            return lines

        if function == 'sx0_fixedParameters':
            # here we only want to overwrite values where x0_fixedParameters
            # was applied
            cases = dict()
            for ipar in range(self.model.np()):
                expressions = []
                for index, formula in enumerate(
                        self.model.eq('x0_fixedParameters')
                ):
                    if formula != 0 and formula != 0.0:
                        expressions.append(f'{function}[{index}] = '
                                                f'{symbol[index, ipar]};')
                cases[ipar] = expressions
            lines.extend(getSwitchStatement('ip', cases, 1))

        elif function in sensi_functions:
            cases = {ipar : self._getSymLines(symbol[:, ipar], function, 0)
                     for ipar in range(self.model.np())}
            lines.extend(getSwitchStatement('ip', cases, 1))

        elif function in multiobs_functions:
            if function == 'dJydy':
                cases = {iobs: self._getSymLines(symbol[iobs], function, 0)
                         for iobs in range(self.model.ny())}
            else:
                cases = {iobs : self._getSymLines(symbol[:, iobs], function, 0)
                         for iobs in range(self.model.ny())}
            lines.extend(getSwitchStatement('iy', cases, 1))

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

        tplData = {
            'MODELNAME': str(self.modelName),
            'NX_RDATA': str(self.model.nx_rdata()),
            'NXTRUE_RDATA': str(self.model.nx_rdata()),
            'NX_SOLVER': str(self.model.nx_solver()),
            'NXTRUE_SOLVER': str(self.model.nx_solver()),
            'NY': str(self.model.ny()),
            'NYTRUE': str(self.model.ny()),
            'NZ': '0',
            'NZTRUE': '0',
            'NEVENT': '0',
            'NOBJECTIVE': '1',
            'NW': str(len(self.model.sym('w'))),
            'NDWDP': str(len(self.model.sparsesym('dwdp'))),
            'NDWDX': str(len(self.model.sparsesym('dwdx'))),
            'NDXDOTDW': str(len(self.model.sparsesym('dxdotdw'))),
            'NDXDOTDP_EXPLICIT': str(len(self.model.sparsesym(
                'dxdotdp_explicit'))),
            'NDXDOTDP_IMPLICIT': str(len(self.model.sparsesym(
                'dxdotdp_implicit'))),
            'NDJYDY': 'std::vector<int>{%s}'
                      % ','.join(str(len(x))
                                 for x in self.model.sparsesym('dJydy')),
            'NNZ': str(len(self.model.sparsesym('JSparse'))),
            'UBW': str(self.model.nx_solver()),
            'LBW': str(self.model.nx_solver()),
            'NP': str(self.model.np()),
            'NK': str(self.model.nk()),
            'O2MODE': 'amici::SecondOrderMode::none',
            'PARAMETERS': str(self.model.val('p'))[1:-1],
            'FIXED_PARAMETERS': str(self.model.val('k'))[1:-1],
            'PARAMETER_NAMES_INITIALIZER_LIST':
                self._getSymbolNameInitializerList('p'),
            'STATE_NAMES_INITIALIZER_LIST':
                self._getSymbolNameInitializerList('x_rdata'),
            'FIXED_PARAMETER_NAMES_INITIALIZER_LIST':
                self._getSymbolNameInitializerList('k'),
            'OBSERVABLE_NAMES_INITIALIZER_LIST':
                self._getSymbolNameInitializerList('y'),
            'PARAMETER_IDS_INITIALIZER_LIST':
                self._getSymbolIDInitializerList('p'),
            'STATE_IDS_INITIALIZER_LIST':
                self._getSymbolIDInitializerList('x_rdata'),
            'FIXED_PARAMETER_IDS_INITIALIZER_LIST':
                self._getSymbolIDInitializerList('k'),
            'OBSERVABLE_IDS_INITIALIZER_LIST':
                self._getSymbolIDInitializerList('y'),
            'REINIT_FIXPAR_INITCOND':
                'true' if self.allow_reinit_fixpar_initcond else
                'false',
            'AMICI_VERSION_STRING':  __version__,
            'AMICI_COMMIT_STRING': __commit__,
        }

        for fun in [
            'w', 'dwdp', 'dwdx', 'x_rdata', 'x_solver', 'total_cl', 'dxdotdw',
            'dxdotdp_explicit', 'dxdotdp_implicit', 'JSparse', 'JSparseB',
            'dJydy'
        ]:
            tplData[f'{fun.upper()}_DEF'] = \
                get_function_extern_declaration(fun, self.modelName)
            tplData[f'{fun.upper()}_IMPL'] = \
                get_model_override_implementation(fun, self.modelName)
            if fun in sparse_functions:
                tplData[f'{fun.upper()}_COLPTRS_DEF'] = \
                    get_sunindex_extern_declaration(fun, self.modelName, 'colptrs')
                tplData[f'{fun.upper()}_COLPTRS_IMPL'] = \
                    get_sunindex_override_implementation(fun, self.modelName, 'colptrs')
                tplData[f'{fun.upper()}_ROWVALS_DEF'] = \
                    get_sunindex_extern_declaration(fun, self.modelName, 'rowvals')
                tplData[f'{fun.upper()}_ROWVALS_IMPL'] = \
                    get_sunindex_override_implementation(fun, self.modelName, 'rowvals')

        if self.model.nx_solver() == self.model.nx_rdata():
            tplData['X_RDATA_DEF'] = ''
            tplData['X_RDATA_IMPL'] = ''

        applyTemplate(
            os.path.join(amiciSrcPath, 'model_header.ODE_template.h'),
            os.path.join(self.modelPath, f'{self.modelName}.h'),
            tplData
        )

    def _getSymbolNameInitializerList(self, name):
        """Get SBML name initializer list for vector of names for the given
        model entity

        Arguments:
            name: any key present in self.model._syms @type str

        Returns:
        Template initializer list of names

        Raises:

        """
        return '\n'.join(
            [
                f'"{strip_pysb(symbol)}",'
                for symbol in self.model.name(name)
            ]
        )

    def _getSymbolIDInitializerList(self, name):
        """Get C++ initializer list for vector of names for the given model
        entity

        Arguments:
            name: any key present in self.model._syms @type str

        Returns:
        Template initializer list of ids

        Raises:

        """
        return '\n'.join(
            [
                f'"{strip_pysb(symbol)}",'
                for symbol in self.model.sym(name)
            ]
        )

    def _writeCMakeFile(self):
        """Write CMake CMakeLists.txt file for this model.

        Arguments:

        Returns:

        Raises:

        """

        sources = [self.modelName + '_' + function + '.cpp '
                   for function in self.functions.keys()
                   if self.functions[function].get('body', None) is not None]

        # add extra source files for sparse matrices
        for function in sparse_functions:
            sources.append(self.modelName + '_' + function
                           + '_colptrs.cpp')
            sources.append(self.modelName + '_' + function
                           + '_rowvals.cpp ')

        templateData = {'MODELNAME': self.modelName,
                        'SOURCES': '\n'.join(sources),
                        'AMICI_VERSION': __version__}
        applyTemplate(
            MODEL_CMAKE_TEMPLATE_FILE,
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
        shutil.copy(SWIG_CMAKE_TEMPLATE_FILE,
                    os.path.join(self.modelSwigPath, 'CMakeLists.txt'))

    def _writeModuleSetup(self):
        """Create a distutils setup.py file for compile the model module.

        Arguments:

        Returns:

        Raises:

        """

        templateData = {'MODELNAME': self.modelName,
                        'AMICI_VERSION': __version__,
                        'PACKAGE_VERSION': '0.1.0'}
        applyTemplate(os.path.join(amiciModulePath, 'setup.template.py'),
                      os.path.join(self.modelPath, 'setup.py'), templateData)
        applyTemplate(os.path.join(amiciModulePath, 'MANIFEST.template.in'),
                      os.path.join(self.modelPath, 'MANIFEST.in'), {})
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

        return [' ' * indentLevel + f'{variable}[{index}] = '
                                    f'{self._printWithException(math)};'
                for index, math in enumerate(symbols)
                if not (math == 0 or math == 0.0)]

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
            ret = ccode(math)
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

    return sp.Matrix(diagonal)


class TemplateAmici(Template):
    """Template format used in AMICI (see string.template for more details).

    Attributes:
        delimiter: delimiter that identifies template variables @type str

    """
    delimiter = 'TPL_'


def applyTemplate(sourceFile,targetFile,templateData):
    """Load source file, apply template substitution as provided in
    templateData and save as targetFile.

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


def strip_pysb(symbol):
    """Strips pysb info from a pysb.Component object

    Arguments:
        symbol: symbolic expression @type sympy.Basic

    Returns:
    stripped sympy.Basic

    Raises:

    """
    # strip pysb type and transform into a flat sympy.Basic.
    # this ensures that the pysb type specific __repr__ is used when converting
    # to string
    if pysb and isinstance(symbol, pysb.Component):
        return sp.Symbol(symbol.name, real=True)
    else:
        # in this case we will use sympy specific transform anyways
        return symbol


def get_function_extern_declaration(fun, name):
    """Constructs the extern function declaration for a given function

    Arguments:
        fun: function name @type str
        name: model name @type str

    Returns:
    c++ function definition string

    Raises:

    """
    return \
        f'extern void {fun}_{name}{functions[fun]["signature"]};'


def get_sunindex_extern_declaration(fun, name, indextype):
    """Constructs the function declaration for an index function of a given
    function

    Arguments:
        fun: function name @type str
        name: model name @type str
        indextype: index function {'colptrs', 'rowvals'} @type str

    Returns:
    c++ function declaration string

    Raises:

    """
    index_arg = ', int index' if fun in multiobs_functions else ''
    return \
        f'extern void {fun}_{indextype}_{name}(sunindextype *{indextype}{index_arg});'


def get_model_override_implementation(fun, name):
    """Constructs amici::Model::* override implementation for a given function

    Arguments:
        fun: function name @type str
        name: model name @type str

    Returns:
    c++ function implementation string

    Raises:

    """
    return \
        'virtual void f{fun}{signature} override {{\n' \
        '{ind8}{fun}_{name}{eval_signature};\n' \
        '{ind4}}}\n'.format(
            ind4=' '*4,
            ind8=' '*8,
            fun=fun,
            name=name,
            signature=functions[fun]["signature"],
            eval_signature=remove_typedefs(functions[fun]["signature"])
        )


def get_sunindex_override_implementation(fun, name, indextype):
    """Constructs the amici::Model:: function implementation for an index
    function of a given function

    Arguments:
        fun: function name @type str
        name: model name @type str
        indextype: index function {'colptrs', 'rowvals'} @type str

    Returns:
    c++ function implementation string

    Raises:

    """
    index_arg = ', int index' if fun in multiobs_functions else ''
    index_arg_eval = ', index' if fun in multiobs_functions else ''

    return \
        'virtual void f{fun}_{indextype}{signature} override {{\n' \
        '{ind8}{fun}_{indextype}_{name}{eval_signature};\n' \
        '{ind4}}}\n'.format(
            ind4=' '*4,
            ind8=' '*8,
            fun=fun,
            indextype=indextype,
            name=name,
            signature=f'(sunindextype *{indextype}{index_arg})',
            eval_signature=f'({indextype}{index_arg_eval})',
        )


def remove_typedefs(signature):
    """Strips typedef info from a function signature

    Arguments:
        signature: function signature @type str

    Returns:
    string that can be used to construct function calls with the same
    variable names and ordering as in the function signature

    Raises:

    """
    # remove * pefix for pointers (pointer must always be removed before
    # values otherwise we will inadvertently dereference values,
    # same applies for const specifications)
    #
    # always add whitespace after type definition for cosmetic reasons
    typedefs = [
        'const realtype *',
        'const double *',
        'const realtype ',
        'double *',
        'realtype *',
        'const int ',
        'int ',
        'SUNMatrixContent_Sparse ',
    ]

    for typedef in typedefs:
        signature = signature.replace(typedef, ' ')

    return signature


def getSwitchStatement(condition, cases,
                  indentation_level=0,
                  indentation_step=' ' * 4):
    """Generate code for switch statement

    Arguments:
        condition: Condition for switch @type str

        cases: Cases as dict with expressions as keys and statement as
        list of strings @type dict

        indentation_level: indentation level

        indentation_step: indentation whitespace per level

    Returns:
    Code for switch expression as list of strings

    Raises:

    """
    lines = list()

    if not cases:
        return lines

    for expression, statements in cases.items():
        if statements:
            lines.append((indentation_level + 1) * indentation_step
                         + f'case {expression}:')
            for statement in statements:
                lines.append((indentation_level + 2) * indentation_step
                             + statement)
            lines.append((indentation_level + 2) * indentation_step + 'break;')

    if lines:
        lines.insert(0, indentation_level * indentation_step
                     + f'switch({condition}) {{')
        lines.append(indentation_level * indentation_step + '}')

    return lines


def csc_matrix(matrix, name, base_index=0, pattern_only=False):
    """Generates the sparse symbolic identifiers, symbolic identifiers,
    sparse matrix, column pointers and row values for a symbolic
    variable

    Arguments:
        matrix: dense matrix to be sparsified @type sp.Matrix
        name: name of the symbolic variable @type str
        base_index: index for first symbol name, defaults to 0
        pattern_only: flag for computing sparsity pattern without whole matrix

    Returns:
        symbolColPtrs, symbolRowVals, sparseList, symbolList, sparseMatrix
    Raises:

    """
    idx = 0
    symbol_name_idx = base_index

    nrows, ncols = matrix.shape

    if not pattern_only:
        sparseMatrix = sp.zeros(nrows, ncols)
    symbolList = []
    sparseList = []
    symbolColPtrs = []
    symbolRowVals = []

    for col in range(0, ncols):
        symbolColPtrs.append(idx)
        for row in range(0, nrows):
            if not (matrix[row, col] == 0):
                symbolRowVals.append(row)
                idx += 1
                symbolName = f'{name}{symbol_name_idx}'
                symbolList.append(symbolName)
                symbol_name_idx += 1
                if not pattern_only:
                    sparseMatrix[row, col] = sp.Symbol(symbolName, real=True)
                    sparseList.append(matrix[row, col])

    if idx == 0:
        symbolColPtrs = []  # avoid bad memory access for empty matrices
    else:
        symbolColPtrs.append(idx)

    if not pattern_only:
        sparseList = sp.Matrix(sparseList)
    else:
        sparseMatrix = None

    return symbolColPtrs, symbolRowVals, sparseList, symbolList, sparseMatrix
