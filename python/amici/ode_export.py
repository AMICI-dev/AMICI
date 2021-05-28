"""
C++ Export
----------
This module provides all necessary functionality specify an ODE model and
generate executable C++ simulation code. The user generally won't have to
directly call any function from this module as this will be done by
:py:func:`amici.pysb_import.pysb2amici`,
:py:func:`amici.sbml_import.SbmlImporter.sbml2amici` and
:py:func:`amici.petab_import.import_model`
"""
import sympy as sp
import numpy as np
import re
import shutil
import subprocess
import sys
import os
import copy
import numbers
import logging
import itertools
import contextlib

try:
    import pysb
except ImportError:
    pysb = None

from typing import (
    Callable, Optional, Union, List, Dict, Tuple, SupportsFloat, Sequence,
    Set, Any
)
from string import Template
from sympy.printing import cxxcode
from sympy.printing.cxx import _CXXCodePrinterBase
from sympy.matrices.immutable import ImmutableDenseMatrix
from sympy.matrices.dense import MutableDenseMatrix
from sympy.logic.boolalg import BooleanAtom
from itertools import chain


from . import (
    amiciSwigPath, amiciSrcPath, amiciModulePath, __version__, __commit__,
    sbml_import
)
from .logging import get_logger, log_execution_time, set_log_level
from .constants import SymbolId
from .import_utils import smart_subs_dict, toposort_symbols

# Template for model simulation main.cpp file
CXX_MAIN_TEMPLATE_FILE = os.path.join(amiciSrcPath, 'main.template.cpp')
# Template for model/swig/CMakeLists.txt
SWIG_CMAKE_TEMPLATE_FILE = os.path.join(amiciSwigPath,
                                        'CMakeLists_model.cmake')
# Template for model/CMakeLists.txt
MODEL_CMAKE_TEMPLATE_FILE = os.path.join(amiciSrcPath,
                                         'CMakeLists.template.cmake')

# prototype for generated C++ functions, keys are the names of functions
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
    'Jy': {
        'signature':
            '(realtype *Jy, const int iy, const realtype *p, '
            'const realtype *k, const realtype *y, const realtype *sigmay, '
            'const realtype *my)',
    },
    'dJydsigma': {
        'signature':
            '(realtype *dJydsigma, const int iy, const realtype *p, '
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
    'root': {
        'signature':
            '(realtype *root, const realtype t, const realtype *x, '
            'const realtype *p, const realtype *k, const realtype *h)'
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
    'dwdw': {
        'signature':
            '(realtype *dwdw, const realtype t, const realtype *x, '
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
    'dxdotdx_explicit': {
        'signature':
            '(realtype *dxdotdx_explicit, const realtype t, '
            'const realtype *x, const realtype *p, const realtype *k, '
            'const realtype *h, const realtype *w)',
        'flags': ['assume_pow_positivity', 'sparse']
    },
    'dxdotdp_explicit': {
        'signature':
            '(realtype *dxdotdp_explicit, const realtype t, '
            'const realtype *x, const realtype *p, const realtype *k, '
            'const realtype *h, const realtype *w)',
        'flags': ['assume_pow_positivity', 'sparse']
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
            'const int ip, const realtype *w, const realtype *dtcldp)',
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
    'sroot': {
        'signature':
            '(realtype *stau, const realtype t, const realtype *x, '
            'const realtype *p, const realtype *k, const realtype *h, '
            'const realtype *sx, const int ip, const int ie)',
        'flags': ['dont_generate_body']
    },
    'drootdt': {
        'signature': '()',
        'flags': ['dont_generate_body']
    },
    'drootdt_total': {
        'signature': '()',
        'flags': ['dont_generate_body']
    },
    'drootdp': {
        'signature': '()',
        'flags': ['dont_generate_body']
    },
    'drootdx': {
        'signature': '()',
        'flags': ['dont_generate_body']
    },
    'stau': {
        'signature':
            '(realtype *stau, const realtype t, const realtype *x, '
            'const realtype *p, const realtype *k, const realtype *h, '
            'const realtype *sx, const int ip, const int ie)'
    },
    'deltax': {
        'signature':
            '(double *deltax, const realtype t, const realtype *x, '
            'const realtype *p, const realtype *k, const realtype *h, '
            'const int ie, const realtype *xdot, const realtype *xdot_old)'
    },
    'ddeltaxdx': {
        'signature': '()',
        'flags': ['dont_generate_body']
    },
    'ddeltaxdt': {
        'signature': '()',
        'flags': ['dont_generate_body']
    },
    'ddeltaxdp': {
        'signature': '()',
        'flags': ['dont_generate_body']
    },
    'deltasx': {
        'signature':
            '(realtype *deltasx, const realtype t, const realtype *x, '
            'const realtype *p, const realtype *k, const realtype *h, '
            'const realtype *w, const int ip, const int ie, '
            'const realtype *xdot, const realtype *xdot_old, '
            'const realtype *sx, const realtype *stau)'
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
            'const realtype *p, const realtype *k, '
            'gsl::span<const int> reinitialization_state_idxs)',
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
            'const int ip, gsl::span<const int> reinitialization_state_idxs)',
    },
    'xdot': {
        'signature':
            '(realtype *xdot, const realtype t, const realtype *x, '
            'const realtype *p, const realtype *k, const realtype *h, '
            'const realtype *w)',
        'flags': ['assume_pow_positivity']
    },
    'xdot_old': {
        'signature': '()',
        'flags': ['dont_generate_body'],
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

# list of sparse functions
sparse_functions = [
    function for function in functions
    if 'sparse' in functions[function].get('flags', [])
]
# list of nobody functions
nobody_functions = [
    function for function in functions
    if 'dont_generate_body' in functions[function].get('flags', [])
]
# list of sensitivity functions
sensi_functions = [
    function for function in functions
    if 'const int ip' in functions[function]['signature']
]
# list of sensitivity functions
sparse_sensi_functions = [
    function for function in functions
    if 'const int ip' not in functions[function]['signature']
    and function.endswith('dp') or function.endswith('dp_explicit')
]
# list of event functions
event_functions = [
    function for function in functions
    if 'const int ie' in functions[function]['signature'] and
        'const int ip' not in functions[function]['signature']
]
event_sensi_functions = [
    function for function in functions
    if 'const int ie' in functions[function]['signature'] and
       'const int ip' in functions[function]['signature']
]
# list of multiobs functions
multiobs_functions = [
    function for function in functions
    if 'const int iy' in functions[function]['signature']
]
# list of equations that have ids which may not be unique
non_unique_id_symbols = [
    'x_rdata', 'y'
]

# custom c++ function replacements
CUSTOM_FUNCTIONS = [
    {'sympy': 'polygamma',
     'c++': 'boost::math::polygamma',
     'include': '#include <boost/math/special_functions/polygamma.hpp>',
     'build_hint': 'Using polygamma requires libboost-math header files.'
     },
    {'sympy': 'Heaviside',
     'c++': 'amici::heaviside'},
    {'sympy': 'DiracDelta',
     'c++': 'amici::dirac'}
]

# python log manager
logger = get_logger(__name__, logging.ERROR)


def var_in_function_signature(name: str, varname: str) -> bool:
    """
    Checks if the values for a symbolic variable is passed in the signature
    of a function

    :param name:
        name of the function
    :param varname:
        name of the symbolic variable

    :return:
        boolean indicating whether the variable occurs in the function
        signature
    """
    return name in functions \
        and re.search(
            rf'const (realtype|double) \*{varname}[0]*[,)]+',
            functions[name]['signature']
        )


class ModelQuantity:
    """
    Base class for model components
    """
    def __init__(self,
                 identifier: sp.Symbol,
                 name: str,
                 value: Union[SupportsFloat, numbers.Number, sp.Expr]):
        """
        Create a new ModelQuantity instance.

        :param identifier:
            unique identifier of the quantity

        :param name:
            individual name of the quantity (does not need to be unique)

        :param value:
            either formula, numeric value or initial value
        """

        if not isinstance(identifier, sp.Symbol):
            raise TypeError(f'identifier must be sympy.Symbol, was '
                            f'{type(identifier)}')
        self._identifier: sp.Symbol = identifier

        if not isinstance(name, str):
            raise TypeError(f'name must be str, was {type(name)}')

        self._name: str = name

        self._value: sp.Expr = cast_to_sym(value, 'value')

    def __repr__(self) -> str:
        """
        Representation of the ModelQuantity object

        :return:
            string representation of the ModelQuantity
        """
        return str(self._identifier)

    def get_id(self) -> sp.Symbol:
        """
        ModelQuantity identifier

        :return:
            identifier of the ModelQuantity
        """
        return self._identifier

    def get_name(self) -> str:
        """
        ModelQuantity name

        :return:
            name of the ModelQuantity
        """
        return self._name

    def get_val(self) -> sp.Expr:
        """
        ModelQuantity value

        :return:
            value of the ModelQuantity
        """
        return self._value

    def set_val(self, val: sp.Expr):
        """
        Set ModelQuantity value

        :return:
            value of the ModelQuantity
        """
        self._value = cast_to_sym(val, 'value')


class State(ModelQuantity):
    """
    A State variable defines an entity that evolves with time according to
    the provided time derivative, abbreviated by `x`

    :ivar _conservation_law:
        algebraic formula that allows computation of this
        state according to a conservation law

    :ivar _dt:
        algebraic formula that defines the temporal derivative of this state

    """

    _dt: Union[sp.Expr, None] = None
    _conservation_law: Union[sp.Expr, None] = None

    def __init__(self,
                 identifier: sp.Symbol,
                 name: str,
                 init: sp.Expr,
                 dt: sp.Expr):
        """
        Create a new State instance. Extends :meth:`ModelQuantity.__init__`
        by dt

        :param identifier:
            unique identifier of the state

        :param name:
            individual name of the state (does not need to be unique)

        :param init:
            initial value

        :param dt:
            time derivative
        """
        super(State, self).__init__(identifier, name, init)
        self._dt = cast_to_sym(dt, 'dt')
        self._conservation_law = None

    def set_conservation_law(self,
                             law: sp.Expr) -> None:
        """
        Sets the conservation law of a state. If the a conservation law
        is set, the respective state will be replaced by an algebraic
        formula according to the respective conservation law.

        :param law:
            linear sum of states that if added to this state remain
            constant over time
        """
        if not isinstance(law, sp.Expr):
            raise TypeError(f'conservation law must have type sympy.Expr, '
                            f'was {type(law)}')

        self._conservation_law = law

    def set_dt(self,
               dt: sp.Expr) -> None:
        """
        Sets the time derivative

        :param dt:
            time derivative
        """
        self._dt = cast_to_sym(dt, 'dt')

    def get_dt(self) -> sp.Expr:
        """
        Gets the time derivative

        :return:
            time derivative
        """
        return self._dt

    def get_free_symbols(self) -> Set[sp.Symbol]:
        """
        Gets the set of free symbols in time derivative and initial conditions

        :return:
            free symbols
        """
        return self._dt.free_symbols.union(self._value.free_symbols)


class ConservationLaw(ModelQuantity):
    """
    A conservation law defines the absolute the total amount of a
    (weighted) sum of states

    """
    def __init__(self,
                 identifier: sp.Symbol,
                 name: str,
                 value: sp.Expr):
        """
        Create a new ConservationLaw instance.

        :param identifier:
            unique identifier of the ConservationLaw

        :param name:
            individual name of the ConservationLaw (does not need  to be
            unique)

        :param value: formula (sum of states)
        """
        super(ConservationLaw, self).__init__(identifier, name, value)


class Observable(ModelQuantity):
    """
    An Observable links model simulations to experimental measurements,
    abbreviated by `y`

    :ivar _measurement_symbol:
        sympy symbol used in the objective function to represent
        measurements to this observable
    """

    _measurement_symbol: Union[sp.Symbol, None] = None

    def __init__(self,
                 identifier: sp.Symbol,
                 name: str,
                 value: sp.Expr,
                 measurement_symbol: Optional[sp.Symbol] = None):
        """
        Create a new Observable instance.

        :param identifier:
            unique identifier of the Observable

        :param name:
            individual name of the Observable (does not need to be unique)

        :param value:
            formula
        """
        super(Observable, self).__init__(identifier, name, value)
        self._measurement_symbol = measurement_symbol

    def get_measurement_symbol(self) -> sp.Symbol:
        if self._measurement_symbol is None:
            self._measurement_symbol = generate_measurement_symbol(
                self.get_id()
            )

        return self._measurement_symbol


class SigmaY(ModelQuantity):
    """
    A Standard Deviation SigmaY rescales the distance between simulations
    and measurements when computing residuals or objective functions,
    abbreviated by `sigmay`
    """
    def __init__(self,
                 identifier: sp.Symbol,
                 name: str,
                 value: sp.Expr):
        """
        Create a new Standard Deviation instance.

        :param identifier:
            unique identifier of the Standard Deviation

        :param name:
            individual name of the Standard Deviation (does not need to
            be unique)

        :param value:
            formula
        """
        super(SigmaY, self).__init__(identifier, name, value)


class Expression(ModelQuantity):
    """
    An Expressions is a recurring elements in symbolic formulas. Specifying
    this may yield more compact expression which may lead to substantially
    shorter model compilation times, but may also reduce model simulation time,
    abbreviated by `w`
    """
    def __init__(self,
                 identifier: sp.Symbol,
                 name: str,
                 value: sp.Expr):
        """
        Create a new Expression instance.

        :param identifier:
            unique identifier of the Expression

        :param name:
            individual name of the Expression (does not need to be unique)

        :param value:
            formula
        """
        super(Expression, self).__init__(identifier, name, value)


class Parameter(ModelQuantity):
    """
    A Parameter is a free variable in the model with respect to which
    sensitivities may be computed, abbreviated by `p`
    """

    def __init__(self,
                 identifier: sp.Symbol,
                 name: str,
                 value: numbers.Number):
        """
        Create a new Expression instance.

        :param identifier:
            unique identifier of the Parameter

        :param name:
            individual name of the Parameter (does not need to be
            unique)

        :param value:
            numeric value
        """
        super(Parameter, self).__init__(identifier, name, value)


class Constant(ModelQuantity):
    """
    A Constant is a fixed variable in the model with respect to which
    sensitivities cannot be computed, abbreviated by `k`
    """

    def __init__(self,
                 identifier: sp.Symbol,
                 name: str,
                 value: numbers.Number):
        """
        Create a new Expression instance.

        :param identifier:
            unique identifier of the Constant

        :param name:
            individual name of the Constant (does not need to be unique)

        :param value:
            numeric value
        """
        super(Constant, self).__init__(identifier, name, value)


class LogLikelihood(ModelQuantity):
    """
    A LogLikelihood defines the distance between measurements and
    experiments for a particular observable. The final LogLikelihood value
    in the simulation will be the sum of all specified LogLikelihood
    instances evaluated at all timepoints, abbreviated by `Jy`
    """

    def __init__(self,
                 identifier: sp.Symbol,
                 name: str,
                 value: sp.Expr):
        """
        Create a new Expression instance.

        :param identifier:
            unique identifier of the LogLikelihood

        :param name:
            individual name of the LogLikelihood (does not need to be
             unique)

        :param value:
            formula
        """
        super(LogLikelihood, self).__init__(identifier, name, value)


class Event(ModelQuantity):
    """
    An Event defines either a SBML event or a root of the argument of a
    Heaviside function. The Heaviside functions will be tracked via the
    vector `h` during simulation and are needed to inform the ODE solver about
    a discontinuity in either the right hand side or the states themselves,
    causing a reinitialization of the solver.
    """

    def __init__(self,
                 identifier: sp.Symbol,
                 name: str,
                 value: sp.Expr,
                 state_update: Union[sp.Expr, None],
                 event_observable: Union[sp.Expr, None]):
        """
        Create a new Event instance.

        :param identifier:
            unique identifier of the Event

        :param name:
            individual name of the Event (does not need to be unique)

        :param value:
            formula for the root / trigger function

        :param state_update:
            formula for the bolus function (None for Heaviside functions,
            zero vector for events without bolus)

        :param event_observable:
            formula a potential observable linked to the event
            (None for Heaviside functions, empty events without observable)
        """
        super(Event, self).__init__(identifier, name, value)
        # add the Event specific components
        self._state_update = state_update
        self._observable = event_observable

    def __eq__(self, other):
        """
        Check equality of events at the level of trigger/root functions, as we
        need to collect unique root functions for roots.cpp
        """
        return self.get_val() == other.get_val()


# defines the type of some attributes in ODEModel
symbol_to_type = {
    SymbolId.SPECIES: State,
    SymbolId.PARAMETER: Parameter,
    SymbolId.FIXED_PARAMETER: Constant,
    SymbolId.OBSERVABLE: Observable,
    SymbolId.SIGMAY: SigmaY,
    SymbolId.LLHY: LogLikelihood,
    SymbolId.EXPRESSION: Expression,
    SymbolId.EVENT: Event
}


@log_execution_time('running smart_jacobian', logger)
def smart_jacobian(eq: sp.MutableDenseMatrix,
                   sym_var: sp.MutableDenseMatrix) -> sp.MutableDenseMatrix:
    """
    Wrapper around symbolic jacobian with some additional checks that reduce
    computation time for large matrices

    :param eq:
        equation
    :param sym_var:
        differentiation variable
    :return:
        jacobian of eq wrt sym_var
    """
    if min(eq.shape) and min(sym_var.shape) \
            and not smart_is_zero_matrix(eq) \
            and not smart_is_zero_matrix(sym_var) \
            and not sym_var.free_symbols.isdisjoint(eq.free_symbols):
        return eq.jacobian(sym_var)
    return sp.zeros(eq.shape[0], sym_var.shape[0])


@log_execution_time('running smart_multiply', logger)
def smart_multiply(x: Union[sp.MutableDenseMatrix, sp.MutableSparseMatrix],
                   y: sp.MutableDenseMatrix
                   ) -> Union[sp.MutableDenseMatrix, sp.MutableSparseMatrix]:
    """
    Wrapper around symbolic multiplication with some additional checks that
    reduce computation time for large matrices

    :param x:
        educt 1
    :param y:
        educt 2
    :return:
        product
    """
    if not x.shape[0] or not y.shape[1] or smart_is_zero_matrix(x) or \
            smart_is_zero_matrix(y):
        return sp.zeros(x.shape[0], y.shape[1])
    return x.multiply(y)


def smart_is_zero_matrix(x: Union[sp.MutableDenseMatrix,
                                  sp.MutableSparseMatrix]) -> bool:
    """A faster implementation of sympy's is_zero_matrix

    Avoids repeated indexer type checks and double iteration to distinguish
    False/None. Found to be about 100x faster for large matrices.

    :param x: Matrix to check
    """

    if isinstance(x, sp.MutableDenseMatrix):
        nonzero = any(xx.is_zero is not True for xx in x._mat)
    else:
        nonzero = x.nnz() > 0

    return not nonzero


class ODEModel:
    """
    Defines an Ordinary Differential Equation as set of ModelQuantities.
    This class provides general purpose interfaces to ompute arbitrary
    symbolic derivatives that are necessary for model simulation or
    sensitivity computation

    :ivar _states:
        list of state variables

    :ivar _observables:
        list of observables

    :ivar _sigmays:
        list of sigmays

    :ivar _parameters:
        list of parameters

    :ivar _loglikelihoods:
        list of loglikelihoods

    :ivar _expressions:
        list of expressions instances

    :ivar _conservationlaws:
        list of conservation laws

    :ivar _symboldim_funs:
        define functions that compute model dimensions, these
        are functions as the underlying symbolic expressions have not been
        populated at compile time

    :ivar _eqs:
        carries symbolic formulas of the symbolic variables of the model

    :ivar _sparseeqs:
        carries linear list of all symbolic formulas for sparsified
        variables

    :ivar _vals:
        carries numeric values of symbolic identifiers of the symbolic
        variables of the model

    :ivar _names:
        carries names of symbolic identifiers of the symbolic variables
        of the model

    :ivar _syms:
        carries symbolic identifiers of the symbolic variables of the
        model

    :ivar _strippedsyms:
        carries symbolic identifiers that were stripped of additional class
        information

    :ivar _sparsesyms:
        carries linear list of all symbolic identifiers for sparsified
        variables

    :ivar _colptrs:
        carries column pointers for sparsified variables. See
        SUNMatrixContent_Sparse definition in <sunmatrix/sunmatrix_sparse.h>

    :ivar _rowvals:
        carries row values for sparsified variables. See
        SUNMatrixContent_Sparse definition in <sunmatrix/sunmatrix_sparse.h>

    :ivar _equation_prototype:
        defines the attribute from which an equation should be generated via
        list comprehension (see :meth:`ODEModel._generate_equation`)

    :ivar _variable_prototype:
        defines the attribute from which a variable should be generated via
        list comprehension (see :meth:`ODEModel._generate_symbol`)

    :ivar _value_prototype:
        defines the attribute from which a value should be generated via
        list comprehension (see :meth:`ODEModel._generate_value`)

    :ivar _total_derivative_prototypes:
        defines how a total derivative equation is computed for an equation,
        key defines the name and values should be arguments for
        ODEModel.totalDerivative()

    :ivar _lock_total_derivative:
        add chainvariables to this set when computing total derivative from
        a partial derivative call to enforce a partial derivative in the
        next recursion. prevents infinite recursion

    :ivar _simplify:
        If not None, this function will be used to simplify symbolic
        derivative expressions. Receives sympy expressions as only argument.
        To apply multiple simplifications, wrap them in a lambda expression.

    :ivar _x0_fixedParameters_idx:
        Index list of subset of states for which x0_fixedParameters was
        computed

    :ivar _w_recursion_depth:
        recursion depth in w, quantified as nilpotency of dwdw

    :ivar _has_quadratic_nllh:
        whether all observables have a gaussian noise model, i.e. whether
        res and FIM make sense.
    """

    def __init__(self, verbose: Optional[Union[bool, int]] = False,
                 simplify: Optional[Callable] = sp.powsimp):
        """
        Create a new ODEModel instance.

        :param verbose:
            verbosity level for logging, True/False default to
            ``logging.DEBUG``/``logging.ERROR``

        :param simplify:
            see :meth:`ODEModel._simplify`
        """
        self._states: List[State] = []
        self._observables: List[Observable] = []
        self._sigmays: List[SigmaY] = []
        self._parameters: List[Parameter] = []
        self._constants: List[Constant] = []
        self._loglikelihoods: List[LogLikelihood] = []
        self._expressions: List[Expression] = []
        self._conservationlaws: List[ConservationLaw] = []
        self._events: List[Event] = []
        self._symboldim_funs: Dict[str, Callable[[], int]] = {
            'sx': self.num_states_solver,
            'v': self.num_states_solver,
            'vB': self.num_states_solver,
            'xB': self.num_states_solver,
            'sigmay': self.num_obs,
        }
        self._eqs: Dict[str, Union[sp.Matrix, List[sp.Matrix]]] = dict()
        self._sparseeqs: Dict[str, Union[sp.Matrix, List[sp.Matrix]]] = dict()
        self._vals: Dict[str, List[float]] = dict()
        self._names: Dict[str, List[str]] = dict()
        self._syms: Dict[str, Union[sp.Matrix, List[sp.Matrix]]] = dict()
        self._strippedsyms: Dict[str, sp.Matrix] = dict()
        self._sparsesyms: Dict[str, Union[List[str], List[List[str]]]] = dict()
        self._colptrs: Dict[str, Union[List[int], List[List[int]]]] = dict()
        self._rowvals: Dict[str, Union[List[int], List[List[int]]]] = dict()

        self._equation_prototype: Dict[str, str] = {
            'total_cl': '_conservationlaws',
            'x0': '_states',
            'y': '_observables',
            'Jy': '_loglikelihoods',
            'w': '_expressions',
            'root': '_events',
            'sigmay': '_sigmays'
        }
        self._variable_prototype: Dict[str, str] = {
            'tcl': '_conservationlaws',
            'x_rdata': '_states',
            'y': '_observables',
            'p': '_parameters',
            'k': '_constants',
            'w': '_expressions',
            'sigmay': '_sigmays',
            'h': '_events'
        }
        self._value_prototype: Dict[str, str] = {
            'p': '_parameters',
            'k': '_constants',
        }
        self._total_derivative_prototypes: \
            Dict[str, Dict[str, Union[str, List[str]]]] = {
                'sx_rdata': {
                    'eq': 'x_rdata',
                    'chainvars': ['x'],
                    'var': 'p',
                    'dxdz_name': 'sx',
                },
                'sroot': {
                    'eq': 'root',
                    'chainvars': ['x'],
                    'var': 'p',
                    'dxdz_name': 'sx',
                }
            }

        self._lock_total_derivative: List[str] = list()
        self._simplify: Callable = simplify
        self._x0_fixedParameters_idx: Union[None, Sequence[int]]
        self._w_recursion_depth: int = 0
        self._has_quadratic_nllh: bool = True
        set_log_level(logger, verbose)

    @log_execution_time('importing SbmlImporter', logger)
    def import_from_sbml_importer(self,
                                  si: 'sbml_import.SbmlImporter',
                                  compute_cls: Optional[bool] = True) -> None:
        """
        Imports a model specification from a
        :class:`amici.sbml_import.SbmlImporter`
        instance.

        :param si:
            imported SBML model
        """

        # get symbolic expression from SBML importers
        symbols = copy.copy(si.symbols)
        nexpr = len(symbols[SymbolId.EXPRESSION])

        # assemble fluxes and add them as expressions to the model
        fluxes = []
        for ir, flux in enumerate(si.flux_vector):
            flux_id = generate_flux_symbol(ir)
            fluxes.append(flux_id)
        nr = len(fluxes)

        # correct time derivatives for compartment changes

        dxdotdw_updates = []

        def transform_dxdt_to_concentration(species_id, dxdt):
            """
            Produces the appropriate expression for the first derivative of a
            species with respect to time, for species that reside in
            compartments with a constant volume, or a volume that is defined by
            an assignment or rate rule.

            :param species_id:
                The identifier of the species (generated in "sbml_import.py").

            :param dxdt:
                The element-wise product of the row in the stoichiometric
                matrix that corresponds to the species (row x_index) and the
                flux (kinetic laws) vector. Ignored in the case of rate rules.
            """
            # The derivation of the below return expressions can be found in
            # the documentation. They are found by rearranging
            # $\frac{d}{dt} (vx) = Sw$ for $\frac{dx}{dt}$, where $v$ is the
            # vector of species compartment volumes, $x$ is the vector of
            # species concentrations, $S$ is the stoichiometric matrix, and $w$
            # is the flux vector. The conditional below handles the cases of
            # species in (i) compartments with a rate rule, (ii) compartments
            # with an assignment rule, and (iii) compartments with a constant
            # volume, respectively.
            species = si.symbols[SymbolId.SPECIES][species_id]

            comp = species['compartment']
            x_index = species['index']
            if comp in si.symbols[SymbolId.SPECIES]:
                dv_dt = si.symbols[SymbolId.SPECIES][comp]['dt']
                xdot = (dxdt - dv_dt * species_id) / comp
                dxdotdw_updates.extend(
                    (x_index, w_index, xdot.diff(r_flux))
                    for w_index, r_flux in enumerate(fluxes)
                )
                return xdot
            elif comp in si.compartment_assignment_rules:
                v = si.compartment_assignment_rules[comp]

                # we need to flatten out assignments in the compartment in
                # order to ensure that we catch all species dependencies
                v = smart_subs_dict(v, si.symbols[SymbolId.EXPRESSION],
                                    'value')
                dv_dt = v.diff(si.amici_time_symbol)
                # we may end up with a time derivative of the compartment
                # volume due to parameter rate rules
                comp_rate_vars = [p for p in v.free_symbols
                                  if p in si.symbols[SymbolId.SPECIES]]
                for var in comp_rate_vars:
                    dv_dt += \
                        v.diff(var) * si.symbols[SymbolId.SPECIES][var]['dt']
                dv_dx = v.diff(species_id)
                xdot = (dxdt - dv_dt * species_id) / (dv_dx * species_id + v)
                dxdotdw_updates.extend(
                    (x_index, w_index, xdot.diff(r_flux))
                    for w_index, r_flux in enumerate(fluxes)
                )
                return xdot
            else:
                v = si.compartments[comp]

                if v == 1.0:
                    return dxdt

                dxdotdw_updates.extend(
                    (x_index, w_index,
                     si.stoichiometric_matrix[x_index, w_index] / v)
                    for w_index in range(si.stoichiometric_matrix.shape[1])
                    if si.stoichiometric_matrix[x_index, w_index] != 0
                )

                return dxdt / v

        # create dynamics without respecting conservation laws first
        dxdt = smart_multiply(si.stoichiometric_matrix,
                              MutableDenseMatrix(fluxes))
        for ix, ((species_id, species), formula) in enumerate(zip(
                symbols[SymbolId.SPECIES].items(),
                dxdt
        )):
            assert ix == species['index']  # check that no reordering occurred
            # rate rules and amount species don't need to be updated
            if 'dt' in species:
                continue
            if species['amount']:
                species['dt'] = formula
            else:
                species['dt'] = transform_dxdt_to_concentration(species_id,
                                                                formula)

        # create all basic components of the ODE model and add them.
        for symbol_name in symbols:
            # transform dict of lists into a list of dicts
            args = ['name', 'identifier']

            if symbol_name == SymbolId.SPECIES:
                args += ['dt', 'init']
            else:
                args += ['value']
            if symbol_name == SymbolId.EVENT:
                args += ['state_update', 'event_observable']

            protos = [
                {
                    'identifier': var_id,
                    **{k: v for k, v in var.items() if k in args}
                }
                for var_id, var in symbols[symbol_name].items()
            ]

            for proto in protos:
                self.add_component(symbol_to_type[symbol_name](**proto))

        # add fluxes as expressions, this needs to happen after base
        # expressions from symbols have been parsed
        for flux_id, flux in zip(fluxes, si.flux_vector):
            self.add_component(Expression(
                identifier=flux_id,
                name=str(flux_id),
                value=flux
            ))

        # process conservation laws
        if compute_cls:
            dxdotdw_updates = si.process_conservation_laws(self,
                                                           dxdotdw_updates)

        nx_solver = si.stoichiometric_matrix.shape[0]
        nw = len(self._expressions)
        ncl = nw - nr - nexpr

        # set derivatives of xdot, if applicable. We do this as we can save
        # a substantial amount of computations by exploiting the structure
        # of the right hand side.
        # the tricky part is that the expressions w do not only contain the
        # flux entries, but also assignment rules and conservation laws.
        # assignment rules are added before the fluxes and
        # _process_conservation_laws is called after the fluxes,
        # but conservation law expressions are inserted at the beginning
        # of the self.eq['w']. Accordingly we concatenate a zero matrix (for
        # rule assignments and conservation laws) with the stoichiometric
        # matrix and then apply the necessary updates from
        # transform_dxdt_to_concentration

        if not any(s in [e.get_id() for e in self._expressions]
                   for s in si.stoichiometric_matrix.free_symbols):
            self._eqs['dxdotdw'] = sp.zeros(nx_solver, ncl + nexpr).row_join(
                si.stoichiometric_matrix
            )
            for ix, iw, val in dxdotdw_updates:
                # offset update according to concatenated zero matrix
                self._eqs['dxdotdw'][ix, ncl + nexpr + iw] = val

        # fill in 'self._sym' based on prototypes and components in ode_model
        self.generate_basic_variables(from_sbml=True)
        self._has_quadratic_nllh = all(
            llh['dist'] in ['normal', 'lin-normal']
            for llh in si.symbols[SymbolId.LLHY].values()
        )

    def add_component(self, component: ModelQuantity,
                      insert_first: Optional[bool] = False) -> None:
        """
        Adds a new ModelQuantity to the model.

        :param component:
            model quantity to be added

        :param insert_first:
            whether to add quantity first or last, relevant when components
            may refer to other components of the same type.
        """
        for comp_type in [Observable, Expression, Parameter, Constant, State,
                          LogLikelihood, SigmaY, ConservationLaw, Event]:
            if isinstance(component, comp_type):
                component_list = getattr(
                    self, f'_{type(component).__name__.lower()}s'
                )
                if insert_first:
                    component_list.insert(0, component)
                else:
                    component_list.append(component)
                return

        raise ValueError(f'Invalid component type {type(component)}')

    def add_conservation_law(self,
                             state: sp.Symbol,
                             total_abundance: sp.Symbol,
                             state_expr: sp.Expr,
                             abundance_expr: sp.Expr) -> None:
        """
        Adds a new conservation law to the model. A conservation law is defined
        by the conserved quantity T = sum_i(a_i * x_i), where a_i are
        coefficients and x_i are different state variables.

        :param state:
            symbolic identifier of the state that should be replaced by
            the conservation law (x_j)

        :param total_abundance:
            symbolic identifier of the total abundance (T/a_j)

        :param state_expr:
            symbolic algebraic formula that replaces the the state. This is
            used to compute the numeric value of of `state` during simulations.
            x_j = T/a_j - sum_i≠j(a_i * x_i)/a_j

        :param abundance_expr:
            symbolic algebraic formula that computes the value of the
            conserved quantity. This is used to update the numeric value for
            `total_abundance` after (re-)initialization.
            T/a_j = sum_i≠j(a_i * x_i)/a_j + x_j
        """
        try:
            ix = [
                s.get_id()
                for s in self._states
            ].index(state)
        except ValueError:
            raise ValueError(f'Specified state {state} was not found in the '
                             f'model states.')

        state_id = self._states[ix].get_id()

        self.add_component(
            Expression(state_id, str(state_id), state_expr),
            insert_first=True
        )

        self.add_component(
            ConservationLaw(
                total_abundance,
                f'total_{state_id}',
                abundance_expr
            )
        )

        self._states[ix].set_conservation_law(state_expr)

    def num_states_rdata(self) -> int:
        """
        Number of states.

        :return:
            number of state variable symbols
        """
        return len(self.sym('x_rdata'))

    def num_states_solver(self) -> int:
        """
        Number of states after applying conservation laws.

        :return:
            number of state variable symbols
        """
        return len(self.sym('x'))

    def num_cons_law(self) -> int:
        """
        Number of conservation laws.

        :return:
            number of conservation laws
        """
        return self.num_states_rdata() - self.num_states_solver()

    def num_state_reinits(self) -> int:
        """
        Number of solver states which would be reinitialized after
        preequilibration

        :return:
            number of state variable symbols with reinitialization
        """
        reinit_states = self.eq('x0_fixedParameters')
        solver_states = self.eq('x_solver')
        return sum([1 for ix in reinit_states if ix in solver_states])

    def num_obs(self) -> int:
        """
        Number of Observables.

        :return:
            number of observable symbols
        """
        return len(self.sym('y'))

    def num_const(self) -> int:
        """
        Number of Constants.

        :return:
            number of constant symbols
        """
        return len(self.sym('k'))

    def num_par(self) -> int:
        """
        Number of Parameters.

        :return:
            number of parameter symbols
        """
        return len(self.sym('p'))

    def num_expr(self) -> int:
        """
        Number of Expressions.

        :return:
            number of expression symbols
        """
        return len(self.sym('w'))

    def num_events(self) -> int:
        """
        Number of Events.

        :return:
            number of event symbols (length of the root vector in AMICI)
        """
        return len(self.sym('h'))

    def sym(self,
            name: str,
            stripped: Optional[bool] = False) -> sp.Matrix:
        """
        Returns (and constructs if necessary) the identifiers for a symbolic
        entity.

        :param name:
            name of the symbolic variable
        :param stripped:
            should additional class information be stripped from the
            symbolic variables (default=False)

        :return:
            matrix of symbolic identifiers
        """
        if name not in self._syms:
            self._generate_symbol(name)

        if stripped and name in self._variable_prototype:
            return self._strippedsyms[name]
        else:
            return self._syms[name]

    def sparsesym(self, name: str, force_generate: bool = True) -> List[str]:
        """
        Returns (and constructs if necessary) the sparsified identifiers for
        a sparsified symbolic variable.

        :param name:
            name of the symbolic variable

        :param force_generate:
            whether the symbols should be generated if not available

        :return:
            linearized Matrix containing the symbolic identifiers

        """
        if name not in sparse_functions:
            raise ValueError(f'{name} is not marked as sparse')
        if name not in self._sparsesyms and force_generate:
            self._generate_sparse_symbol(name)
        return self._sparsesyms.get(name, [])

    def eq(self, name: str) -> sp.Matrix:
        """
        Returns (and constructs if necessary) the formulas for a symbolic
        entity.

        :param name:
            name of the symbolic variable

        :return:
            matrix of symbolic formulas
        """

        if name not in self._eqs:
            dec = log_execution_time(f'computing {name}', logger)
            dec(self._compute_equation)(name)
        return self._eqs[name]

    def sparseeq(self, name) -> sp.Matrix:
        """
        Returns (and constructs if necessary) the sparsified formulas for a
        sparsified symbolic variable.

        :param name:
            name of the symbolic variable

        :return:
            linearized matrix containing the symbolic formulas

        """
        if name not in sparse_functions:
            raise ValueError(f'{name} is not marked as sparse')
        if name not in self._sparseeqs:
            self._generate_sparse_symbol(name)
        return self._sparseeqs[name]

    def colptrs(self, name: str) -> Union[List[sp.Number],
                                          List[List[sp.Number]]]:
        """
        Returns (and constructs if necessary) the column pointers for
        a sparsified symbolic variable.

        :param name:
            name of the symbolic variable

        :return:
            list containing the column pointers

        """
        if name not in sparse_functions:
            raise ValueError(f'{name} is not marked as sparse')
        if name not in self._sparseeqs:
            self._generate_sparse_symbol(name)
        return self._colptrs[name]

    def rowvals(self, name: str) -> Union[List[sp.Number],
                                          List[List[sp.Number]]]:
        """
        Returns (and constructs if necessary) the row values for a
        sparsified symbolic variable.

        :param name:
            name of the symbolic variable

        :return:
            list containing the row values

        """
        if name not in sparse_functions:
            raise ValueError(f'{name} is not marked as sparse')
        if name not in self._sparseeqs:
            self._generate_sparse_symbol(name)
        return self._rowvals[name]

    def val(self, name: str) -> List[float]:
        """
        Returns (and constructs if necessary) the numeric values of a
        symbolic entity

        :param name:
            name of the symbolic variable

        :return:
            list containing the numeric values

        """
        if name not in self._vals:
            self._generate_value(name)
        return self._vals[name]

    def name(self, name: str) -> List[str]:
        """
        Returns (and constructs if necessary) the names of a symbolic
        variable

        :param name:
            name of the symbolic variable

        :return:
            list of names

        """
        if name not in self._names:
            self._generate_name(name)
        return self._names[name]

    def free_symbols(self) -> Set[sp.Basic]:
        """
        Returns list of free symbols that appear in ODE rhs and initial
        conditions.
        """
        return set(chain.from_iterable(
            state.get_free_symbols()
            for state in self._states
        ))

    def _generate_symbol(self, name: str, *, from_sbml: bool = False) -> None:
        """
        Generates the symbolic identifiers for a symbolic variable

        :param name:
            name of the symbolic variable

        """
        if name in self._variable_prototype:
            component = self._variable_prototype[name]
            self._syms[name] = sp.Matrix([
                comp.get_id()
                for comp in getattr(self, component)
            ])
            # this gives us access to the "stripped" symbols that were
            # generated by pysb (if compiling a pysb model). To ensure
            # correctness of derivatives, the same assumptions as in pysb
            # have to be used (currently no assumptions)
            # NB if we are compiling a SBML model, then it will be the same
            # as the "non-stripped" in order to preserve assumptions
            self._strippedsyms[name] = self._syms[name] if from_sbml \
                else sp.Matrix([
                    sp.Symbol(comp.get_name())
                    for comp in getattr(self, component)
                ])
            if name == 'y':
                self._syms['my'] = sp.Matrix([
                    comp.get_measurement_symbol()
                    for comp in getattr(self, component)
                ])
            return
        elif name == 'x':
            self._syms[name] = sp.Matrix([
                state.get_id()
                for state in self._states
                if state._conservation_law is None
            ])
            return
        elif name == 'sx0':
            self._syms[name] = sp.Matrix([
                f's{state.get_id()}_0'
                for state in self._states
                if state._conservation_law is None
            ])
            return
        elif name == 'dtcldp':
            # check, whether the CL consists of only one state. Then,
            # sensitivities drop out, otherwise generate symbols
            self._syms[name] = sp.Matrix([
                [sp.Symbol(f's{strip_pysb(tcl.get_id())}__'
                           f'{strip_pysb(par.get_id())}', real=True)
                    for par in self._parameters]
                if self.conservation_law_has_multispecies(tcl)
                else [0] * self.num_par()
                for tcl in self._conservationlaws
            ])
            return
        elif name == 'xdot_old':
            length = len(self.eq('xdot'))
        elif name in sparse_functions:
            self._generate_sparse_symbol(name)
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

    def generate_basic_variables(self, *, from_sbml: bool = False) -> None:
        """
        Generates the symbolic identifiers for all variables in
        ODEModel.variable_prototype

        :param from_sbml:
            whether the model is generated from SBML
        """
        # We need to process events and Heaviside functions in the ODE Model,
        # before adding it to ODEExporter
        self.parse_events()

        for var in self._variable_prototype:
            if var not in self._syms:
                self._generate_symbol(var, from_sbml=from_sbml)

        self._generate_symbol('x', from_sbml=from_sbml)

    def parse_events(self) -> None:
        """
        This functions checks the right hand side for roots of Heaviside
        functions or events, collects the roots, removes redundant roots,
        and replaces the formulae of the found roots by identifiers of AMICI's
        Heaviside function implementation in the right hand side
        """
        # Track all roots functions in the right hand side
        roots = copy.deepcopy(self._events)
        for state in self._states:
            state.set_dt(self._process_heavisides(state.get_dt(), roots))

        for expr in self._expressions:
            expr.set_val(self._process_heavisides(expr.get_val(), roots))

        # remove all possible Heavisides from roots, which may arise from
        # the substitution of `'w'` in `_collect_heaviside_roots`
        for root in roots:
            root.set_val(self._process_heavisides(root.get_val(), roots))

        # Now add the found roots to the model components
        for root in roots:
            # skip roots of SBML events, as these have already been added
            if root in self._events:
                continue
            # add roots of heaviside functions
            self.add_component(root)

    def get_appearance_counts(self, idxs: List[int]) -> List[int]:
        """
        Counts how often a state appears in the time derivative of
        another state and expressions for a subset of states

        :param idxs:
            list of state indices for which counts are to be computed

        :return:
            list of counts for the states ordered according to the provided
            indices

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

    def _generate_sparse_symbol(self, name: str) -> None:
        """
        Generates the sparse symbolic identifiers, symbolic identifiers,
        sparse equations, column pointers and row values for a symbolic
        variable

        :param name:
            name of the symbolic variable

        """
        matrix = self.eq(name)
        match_deriv = re.match(r'd([\w]+)d([a-z]+)', name)
        if match_deriv:
            rownames = self.sym(match_deriv.group(1))
            colnames = self.sym(match_deriv.group(2))

        if name == 'dJydy':
            # One entry per y-slice
            self._colptrs[name] = []
            self._rowvals[name] = []
            self._sparseeqs[name] = []
            self._sparsesyms[name] = []
            self._syms[name] = []
            for iy in range(self.num_obs()):
                symbol_col_ptrs, symbol_row_vals, sparse_list, symbol_list, \
                    sparse_matrix = csc_matrix(matrix[iy, :],
                                               rownames=rownames,
                                               colnames=colnames,
                                               identifier=iy)
                self._colptrs[name].append(symbol_col_ptrs)
                self._rowvals[name].append(symbol_row_vals)
                self._sparseeqs[name].append(sparse_list)
                self._sparsesyms[name].append(symbol_list)
                self._syms[name].append(sparse_matrix)
        else:
            symbol_col_ptrs, symbol_row_vals, sparse_list, symbol_list, \
                sparse_matrix = csc_matrix(
                    matrix, rownames=rownames, colnames=colnames,
                    pattern_only=name in nobody_functions
                )

            self._colptrs[name] = symbol_col_ptrs
            self._rowvals[name] = symbol_row_vals
            self._sparseeqs[name] = sparse_list
            self._sparsesyms[name] = symbol_list
            self._syms[name] = sparse_matrix

    def _compute_equation(self, name: str) -> None:
        """
        computes the symbolic formula for a symbolic variable

        :param name:
            name of the symbolic variable

        """
        # replacement ensures that we don't have to adapt name in abstract
        # model and keep backwards compatibility with matlab
        match_deriv = re.match(r'd([\w_]+)d([a-z_]+)',
                               name.replace('dJydsigma', 'dJydsigmay'))
        time_symbol = sp.Matrix([symbol_with_assumptions('t')])

        if name in self._equation_prototype:
            self._equation_from_component(name, self._equation_prototype[name])

        elif name in self._total_derivative_prototypes:
            args = self._total_derivative_prototypes[name]
            args['name'] = name
            self._lock_total_derivative += args['chainvars']
            self._total_derivative(**args)
            for cv in args['chainvars']:
                self._lock_total_derivative.remove(cv)

        elif name == 'xdot':
            self._eqs[name] = sp.Matrix([
                s.get_dt() for s in self._states
                if s._conservation_law is None
            ])

        elif name == 'x_rdata':
            self._eqs[name] = sp.Matrix([
                state.get_id()
                if state._conservation_law is None
                else state._conservation_law
                for state in self._states
            ])

        elif name == 'x_solver':
            self._eqs[name] = sp.Matrix([
                state.get_id()
                for state in self._states
                if state._conservation_law is None
            ])

        elif name == 'sx_solver':
            self._eqs[name] = sp.Matrix([
                self.sym('sx_rdata')[ix]
                for ix, state in enumerate(self._states)
                if state._conservation_law is None
            ])

        elif name == 'sx0':
            self._derivative(name[1:], 'p', name=name)

        elif name == 'sx0_fixedParameters':
            # deltax = -x+x0_fixedParameters if x0_fixedParameters>0 else 0
            # deltasx = -sx+dx0_fixed_parametersdx*sx+dx0_fixedParametersdp
            # if x0_fixedParameters>0 else 0
            # sx0_fixedParameters = sx+deltasx =
            # dx0_fixed_parametersdx*sx+dx0_fixedParametersdp
            self._eqs[name] = smart_jacobian(
                self.eq('x0_fixedParameters'), self.sym('p')
            )

            dx0_fixed_parametersdx = smart_jacobian(
                self.eq('x0_fixedParameters'), self.sym('x')
            )

            if not smart_is_zero_matrix(dx0_fixed_parametersdx):
                if isinstance(self._eqs[name], ImmutableDenseMatrix):
                    self._eqs[name] = MutableDenseMatrix(self._eqs[name])
                for ip in range(self._eqs[name].shape[1]):
                    self._eqs[name][:, ip] += smart_multiply(
                        dx0_fixed_parametersdx, self.sym('sx0')
                    )

        elif name == 'x0_fixedParameters':
            k = self.sym('k')
            self._x0_fixedParameters_idx = [
                ix
                for ix, eq in enumerate(self.eq('x0'))
                if any([sym in eq.free_symbols for sym in k])
            ]
            eq = self.eq('x0')
            self._eqs[name] = sp.Matrix([eq[ix] for ix in
                                         self._x0_fixedParameters_idx])

        elif name == 'dtotal_cldx_rdata':
            # not correctly parsed in regex
            self._derivative('total_cl', 'x_rdata')

        elif name == 'dtcldx':
            # this is always zero
            self._eqs[name] = \
                sp.zeros(self.num_cons_law(), self.num_states_solver())

        elif name == 'dtcldp':
            # force symbols
            self._eqs[name] = self.sym(name)

        elif name == 'dxdotdx_explicit':
            # force symbols
            self._derivative('xdot', 'x', name=name)

        elif name == 'dxdotdp_explicit':
            # force symbols
            self._derivative('xdot', 'p', name=name)

        elif name == 'drootdt':
            self._eqs[name] = smart_jacobian(self.eq('root'), time_symbol)

        elif name == 'drootdt_total':
            # backsubstitution of optimized right hand side terms into RHS
            # calling subs() is costly. Due to looping over events though, the
            # following lines are only evaluated if a model has events
            w_sorted = \
                toposort_symbols(dict(zip(self._syms['w'], self._eqs['w'])))
            tmp_xdot = smart_subs_dict(self._eqs['xdot'], w_sorted)
            self._eqs[name] = (
                smart_multiply(self.eq('drootdx'), tmp_xdot)
                + self.eq('drootdt')
            )

        elif name == 'deltax':
            # fill boluses for Heaviside functions, as empty state updates
            # would cause problems when writing the function file later
            event_eqs = []
            for event in self._events:
                if event._state_update is None:
                    event_eqs.append(sp.zeros(self.num_states_solver(), 1))
                else:
                    event_eqs.append(event._state_update)

            self._eqs[name] = event_eqs

        elif name == 'ddeltaxdx':
            self._eqs[name] = [
                smart_jacobian(self.eq('deltax')[ie], self.sym('x'))
                for ie in range(self.num_events())
            ]

        elif name == 'ddeltaxdt':
            self._eqs[name] = [
                smart_jacobian(self.eq('deltax')[ie], time_symbol)
                for ie in range(self.num_events())
            ]

        elif name == 'ddeltaxdp':
            self._eqs[name] = [
                smart_jacobian(self.eq('deltax')[ie], self.sym('p'))
                for ie in range(self.num_events())
            ]

        elif name == 'stau':
            self._eqs[name] = [
                -self.eq('sroot')[ie, :] / self.eq('drootdt_total')[ie]
                for ie in range(self.num_events())
            ]

        elif name == 'deltasx':
            event_eqs = []
            for ie, event in enumerate(self._events):
                if event._state_update is not None:
                    # ====== chain rule for the state variables ===============
                    # get xdot with expressions back-substituted
                    tmp_eq = smart_multiply(
                        (self.sym('xdot_old') - self.sym('xdot')),
                        self.eq('stau')[ie])
                    # construct an enhanced state sensitivity, which accounts
                    # for the time point sensitivity as well
                    tmp_dxdp = self.sym('sx') * sp.ones(1, self.num_par())
                    tmp_dxdp += smart_multiply(self.sym('xdot'),
                                               self.eq('stau')[ie])
                    tmp_eq += smart_multiply(self.eq('ddeltaxdx')[ie],
                                             tmp_dxdp)
                    # ====== chain rule for the time point ====================
                    tmp_eq += smart_multiply(self.eq('ddeltaxdt')[ie],
                                             self.eq('stau')[ie])
                    # ====== partial derivative for the parameters ============
                    tmp_eq += self.eq('ddeltaxdp')[ie]
                else:
                    tmp_eq = smart_multiply(
                        (self.eq('xdot_old') - self.eq('xdot')),
                        self.eq('stau')[ie])

                event_eqs.append(tmp_eq)

            self._eqs[name] = event_eqs

        elif name == 'xdot_old':
            # force symbols
            self._eqs[name] = self.sym(name)

        elif match_deriv:
            self._derivative(match_deriv.group(1), match_deriv.group(2), name)

        else:
            raise ValueError(f'Unknown equation {name}')

        if name == 'root':
            # Events are processed after the ODE model has been set up.
            # Equations are there, but symbols for roots must be added
            self.sym('h')

        if name in ['Jy', 'dydx']:
            # do not transpose if we compute the partial derivative as part of
            # a total derivative
            if not len(self._lock_total_derivative):
                self._eqs[name] = self._eqs[name].transpose()

        if self._simplify:
            dec = log_execution_time(f'simplifying {name}', logger)
            if isinstance(self._eqs[name], list):
                self._eqs[name] = [dec(sub_eq.applyfunc)(self._simplify)
                                   for sub_eq in self._eqs[name]]
            else:
                self._eqs[name] = \
                    dec(self._eqs[name].applyfunc)(self._simplify)

    def sym_names(self) -> List[str]:
        """
        Returns a list of names of generated symbolic variables


        :return:
            list of names

        """
        return list(self._syms.keys())

    def _derivative(self, eq: str, var: str, name: str = None) -> None:
        """
        Creates a new symbolic variable according to a derivative

        :param eq:
            name of the symbolic variable that defines the formula

        :param var:
            name of the symbolic variable that defines the identifiers
            with respect to which the derivatives are to be computed

        :param name:
            name of resulting symbolic variable, default is d{eq}d{var}
        """
        if not name:
            name = f'd{eq}d{var}'

        # automatically detect chainrule
        chainvars = []
        ignore_chainrule = {
            ('xdot', 'p'): 'w',  # has generic implementation in c++ code
            ('xdot', 'x'): 'w',  # has generic implementation in c++ code
            ('w', 'w'): 'tcl',   # dtcldw = 0
            ('w', 'x'): 'tcl',   # dtcldx = 0
        }
        for cv in ['w', 'tcl']:
            if var_in_function_signature(eq, cv) \
                    and cv not in self._lock_total_derivative \
                    and var is not cv \
                    and min(self.sym(cv).shape) \
                    and (
                            (eq, var) not in ignore_chainrule
                            or ignore_chainrule[(eq, var)] != cv
                    ):
                chainvars.append(cv)

        if len(chainvars):
            self._lock_total_derivative += chainvars
            self._total_derivative(name, eq, chainvars, var)
            for cv in chainvars:
                self._lock_total_derivative.remove(cv)
            return

        # this is the basic requirement check
        needs_stripped_symbols = eq == 'xdot' and var != 'x'

        # partial derivative
        if eq == 'Jy':
            sym_eq = self.eq(eq).transpose()
        else:
            sym_eq = self.eq(eq)

        if pysb is not None and needs_stripped_symbols:
            needs_stripped_symbols = not any(
                isinstance(sym, pysb.Component)
                for sym in sym_eq.free_symbols
            )

        # now check whether we are working with energy_modeling branch
        # where pysb class info is already stripped
        # TODO: fixme as soon as energy_modeling made it to the main pysb
        #  branch
        sym_var = self.sym(var, needs_stripped_symbols)

        derivative = smart_jacobian(sym_eq, sym_var)

        self._eqs[name] = derivative

        # compute recursion depth based on nilpotency of jacobian. computing
        # nilpotency can be done more efficiently on numerical sparsity pattern
        if name == 'dwdw':
            nonzeros = np.asarray(
                derivative.applyfunc(lambda x: int(not x.is_zero))
            ).astype(np.int64)
            if max(nonzeros.shape):
                while nonzeros.max():
                    nonzeros = nonzeros.dot(nonzeros)
                    self._w_recursion_depth += 1
                    if self._w_recursion_depth > len(sym_eq):
                        raise RuntimeError(
                            'dwdw is not nilpotent. Something, somewhere went '
                            'terribly wrong. Please file a bug report at '
                            'https://github.com/AMICI-dev/AMICI/issues and '
                            'attach this model.'
                        )

        if name == 'dydw' and not smart_is_zero_matrix(derivative):
            dwdw = self.eq('dwdw')
            # h(k) = d{eq}dw*dwdw^k* (k=1)
            h = smart_multiply(derivative, dwdw)
            while not smart_is_zero_matrix(h):
                self._eqs[name] += h
                # h(k+1) = d{eq}dw*dwdw^(k+1) = h(k)*dwdw
                h = smart_multiply(h, dwdw)

    def _total_derivative(self, name: str, eq: str, chainvars: List[str],
                          var: str, dydx_name: str = None,
                          dxdz_name: str = None) -> None:
        """
        Creates a new symbolic variable according to a total derivative
        using the chain rule

        :param name:
            name of resulting symbolic variable

        :param eq:
            name of the symbolic variable that defines the formula

        :param chainvars:
            names of the symbolic variable that define the
            identifiers with respect to which the chain rules are applied

        :param var:
            name of the symbolic variable that defines the identifiers
            whith respect to which the derivatives are to be computed

        :param dydx_name:
            defines the name of the symbolic variable that
            defines the derivative of the `eq` with respect to `chainvar`,
            default is d{eq}d{chainvar}

        :param dxdz_name:
            defines the name of the symbolic variable that
            defines the derivative of the `chainvar` with respect to `var`,
            default is d{chainvar}d{var}

        """
        # compute total derivative according to chainrule
        # Dydz = dydx*dxdz + dydz

        # initialize with partial derivative dydz without chain rule
        self._eqs[name] = self.sym_or_eq(name, f'd{eq}d{var}')
        if not isinstance(self._eqs[name], sp.Symbol):
            # if not a Symbol, create a copy using sympy API
            # NB deepcopy does not work safely, see sympy issue #7672
            self._eqs[name] = self._eqs[name].copy()

        for chainvar in chainvars:
            if dydx_name is None:
                dydx_name = f'd{eq}d{chainvar}'
            if dxdz_name is None:
                dxdz_name = f'd{chainvar}d{var}'

            dydx = self.sym_or_eq(name, dydx_name)
            dxdz = self.sym_or_eq(name, dxdz_name)
            # Save time for for large models if one multiplicand is zero,
            # which is not checked for by sympy
            if not smart_is_zero_matrix(dydx) and not \
                    smart_is_zero_matrix(dxdz):
                if dxdz.shape[1] == 1 and \
                        self._eqs[name].shape[1] != dxdz.shape[1]:
                    for iz in range(self._eqs[name].shape[1]):
                        self._eqs[name][:, iz] += smart_multiply(dydx, dxdz)
                else:
                    self._eqs[name] += smart_multiply(dydx, dxdz)

    def sym_or_eq(self, name: str, varname: str) -> sp.Matrix:
        """
        Returns symbols or equations depending on whether a given
        variable appears in the function signature or not.

        :param name:
            name of function for which the signature should be checked

        :param varname:
            name of the variable which should be contained in the
            function signature

        :return:
            the variable symbols if the variable is part of the signature and
            the variable equations otherwise.

        """
        # dwdx and dwdp will be dynamically computed and their ordering
        # within a column may differ from the initialization of symbols here,
        # so those are not safe to use. Not removing them from signature as
        # this would break backwards compatibility.
        if var_in_function_signature(name, varname) \
                and varname not in ['dwdx', 'dwdp']:
            return self.sym(varname)
        else:
            return self.eq(varname)

    def _multiplication(self, name: str, x: str, y: str,
                        transpose_x: Optional[bool] = False,
                        sign: Optional[int] = 1):
        """
        Creates a new symbolic variable according to a multiplication

        :param name:
            name of resulting symbolic variable, default is d{eq}d{var}

        :param x:
            name of the symbolic variable that defines the first factor

        :param y:
            name of the symbolic variable that defines the second factor

        :param transpose_x:
            indicates whether the first factor should be
            transposed before multiplication

        :param sign:
            defines the sign of the product, should be +1 or -1


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

        self._eqs[name] = sign * smart_multiply(xx, yy)

    def _equation_from_component(self, name: str, component: str) -> None:
        """
        Generates the formulas of a symbolic variable from the attributes

        :param name:
            name of resulting symbolic variable

        :param component:
            name of the attribute

        """
        self._eqs[name] = sp.Matrix(
            [comp.get_val() for comp in getattr(self, component)]
        )

    def get_conservation_laws(self) -> List[Tuple[sp.Symbol, sp.Basic]]:
        """ Returns a list of states with conservation law set


        :return:
            list of state identifiers

        """
        return [
            (state.get_id(), state._conservation_law)
            for state in self._states
            if state._conservation_law is not None
        ]

    def _generate_value(self, name: str) -> None:
        """
        Generates the numeric values of a symbolic variable from value
        prototypes

        :param name:
            name of resulting symbolic variable

        """
        if name in self._value_prototype:
            component = self._value_prototype[name]
        else:
            raise ValueError(f'No values for {name}')

        self._vals[name] = [comp.get_val()
                            for comp in getattr(self, component)]

    def _generate_name(self, name: str) -> None:
        """
        Generates the names of a symbolic variable from variable prototypes or
        equation prototypes

        :param name:
            name of resulting symbolic variable

        """
        if name in self._variable_prototype:
            component = self._variable_prototype[name]
        elif name in self._equation_prototype:
            component = self._equation_prototype[name]
        else:
            raise ValueError(f'No names for {name}')

        self._names[name] = [comp.get_name()
                             for comp in getattr(self, component)]

    def state_has_fixed_parameter_initial_condition(self, ix: int) -> bool:
        """
        Checks whether the state at specified index has a fixed parameter
        initial condition

        :param ix:
            state index

        :return:
            boolean indicating if any of the initial condition free
            variables is contained in the model constants

        """
        ic = self._states[ix].get_val()
        if not isinstance(ic, sp.Basic):
            return False
        return any([
            fp in [c.get_id() for c in self._constants]
            for fp in ic.free_symbols
        ])

    def state_has_conservation_law(self, ix: int) -> bool:
        """
        Checks whether the state at specified index has a conservation
        law set

        :param ix:
            state index

        :return:
            boolean indicating if conservation_law is not None

        """
        return self._states[ix]._conservation_law is not None

    def state_is_constant(self, ix: int) -> bool:
        """
        Checks whether the temporal derivative of the state is zero

        :param ix:
            state index

        :return:
            boolean indicating if constant over time

        """
        return self._states[ix].get_dt() == 0.0

    def conservation_law_has_multispecies(self,
                                          tcl: ConservationLaw) -> bool:
        """
        Checks whether a conservation law has multiple species or it just
        defines one constant species

        :param tcl:
            conservation law

        :return:
            boolean indicating if conservation_law is not None

        """
        state_set = set(self.sym('x_rdata'))
        n_species = len(state_set.intersection(tcl.get_val().free_symbols))
        return n_species > 1

    def _expr_is_time_dependent(self, expr: sp.Expr) -> bool:
        """Determine whether an expression is time-dependent.

        :param expr:
            The expression.

        :returns:
            Whether the expression is time-dependent.
        """
        # `expr.free_symbols` will be different to `self._states.keys()`, so
        # it's easier to compare as `str`.
        expr_syms = {str(sym) for sym in expr.free_symbols}

        # Check if the time variable is in the expression.
        if 't' in expr_syms:
            return True

        # Check if any time-dependent states are in the expression.
        state_syms = [str(sym) for sym in self._states]
        for state in expr_syms.intersection(state_syms):
            if not self.state_is_constant(state_syms.index(state)):
                return True

        return False

    def _get_unique_root(
            self,
            root_found: sp.Expr,
            roots: List[Event],
    ) -> sp.Symbol:
        """
        Collects roots of Heaviside functions and events and stores them in
        the roots list. It checks for redundancy to not store symbolically
        equivalent root functions more than once.

        :param root_found:
            equation of the root function
        :param roots:
            list of already known root functions with identifier

        :returns:
            unique identifier for root, or `None` if the root is not
            time-dependent
        """

        if not self._expr_is_time_dependent(root_found):
            return None

        for root in roots:
            if sp.simplify(root_found - root.get_val()) == 0:
                return root.get_id()

        # create an event for a new root function
        root_symstr = f'Heaviside_{len(roots)}'
        roots.append(Event(
            identifier=sp.Symbol(root_symstr),
            name=root_symstr,
            value=root_found,
            state_update=None,
            event_observable=None
        ))
        return roots[-1].get_id()

    def _collect_heaviside_roots(
            self,
            args: Sequence[sp.Expr],
    ) -> List[sp.Expr]:
        """
        Recursively checks an expression for the occurrence of Heaviside
        functions and return all roots found

        :param args:
            args attribute of the expanded expression

        :returns:
            root functions that were extracted from Heaviside function
            arguments
        """
        root_funs = []
        for arg in args:
            if arg.func == sp.Heaviside:
                root_funs.append(arg.args[0])
            elif arg.has(sp.Heaviside):
                root_funs.extend(self._collect_heaviside_roots(arg.args))

        # substitute 'w' expressions into root expressions now, to avoid
        # rewriting '{model_name}_root.cpp' and '{model_name}_stau.cpp' headers
        # to include 'w.h'
        w_sorted = toposort_symbols(dict(zip(
            [expr.get_id()  for expr in self._expressions],
            [expr.get_val() for expr in self._expressions],
        )))
        root_funs = [
            r.subs(w_sorted)
            for r in root_funs
        ]

        return root_funs

    def _process_heavisides(
            self,
            dxdt: sp.Expr,
            roots: List[Event],
    ) -> sp.Expr:
        """
        Parses the RHS of a state variable, checks for Heaviside functions,
        collects unique roots functions that can be tracked by SUNDIALS and
        replaces Heaviside Functions by amici helper variables that will be
        updated based on SUNDIALS root tracking.

        :param dxdt:
            right hand side of state variable
        :param roots:
            list of known root functions with identifier

        :returns:
            dxdt with Heaviside functions replaced by amici helper variables
        """

        # expanding the rhs will in general help to collect the same
        # heaviside function
        dt_expanded = dxdt.expand()
        # track all the old Heaviside expressions in tmp_roots_old
        # replace them later by the new expressions
        heavisides = []
        # run through the expression tree and get the roots
        tmp_roots_old = self._collect_heaviside_roots(dt_expanded.args)
        for tmp_old in tmp_roots_old:
            # we want unique identifiers for the roots
            tmp_new = self._get_unique_root(tmp_old, roots)
            # `tmp_new` is None if the root is not time-dependent.
            if tmp_new is None:
                continue
            # For Heavisides, we need to add the negative function as well
            self._get_unique_root(sp.sympify(- tmp_old), roots)
            heavisides.append((sp.Heaviside(tmp_old), tmp_new))

        if heavisides:
            # only apply subs if necessary
            for heaviside_sympy, heaviside_amici in heavisides:
                dxdt = dxdt.subs(heaviside_sympy, heaviside_amici)

        return dxdt


def _print_with_exception(math: sp.Expr) -> str:
    """
    Generate C++ code for a symbolic expression

    :param math:
        symbolic expression

    :return:
        C++ code for the specified expression
    """
    # get list of custom replacements
    user_functions = {fun['sympy']: fun['c++'] for fun in CUSTOM_FUNCTIONS}

    try:
        # Required until https://github.com/sympy/sympy/pull/20558 is released
        with _monkeypatched(_CXXCodePrinterBase, '_print_Max',
                            _custom_print_max),\
                _monkeypatched(_CXXCodePrinterBase, '_print_Min',
                               _custom_print_min):
            ret = cxxcode(math, standard='c++11',
                          user_functions=user_functions)
        ret = re.sub(r'(^|\W)M_PI(\W|$)', r'\1amici::pi\2', ret)
        return ret
    except TypeError as e:
        raise ValueError(
            f'Encountered unsupported function in expression "{math}": '
            f'{e}!'
        )


def _get_sym_lines_array(equations: sp.Matrix,
                         variable: str,
                         indent_level: int) -> List[str]:
    """
    Generate C++ code for assigning symbolic terms in symbols to C++ array
    `variable`.

    :param equations:
        vectors of symbolic expressions

    :param variable:
        name of the C++ array to assign to

    :param indent_level:
        indentation level (number of leading blanks)

    :return:
        C++ code as list of lines

    """

    return [' ' * indent_level + f'{variable}[{index}] = '
                                 f'{_print_with_exception(math)};'
            for index, math in enumerate(equations)
            if not (math == 0 or math == 0.0)]


def _get_sym_lines_symbols(symbols: sp.Matrix,
                           equations: sp.Matrix,
                           variable: str,
                           indent_level: int) -> List[str]:
    """
    Generate C++ code for where array elements are directly replaced with
    their corresponding macro symbol

    :param symbols:
        vectors of symbols that equations are assigned to

    :param equations:
        vectors of expressions

    :param variable:
        name of the C++ array to assign to, only used in comments

    :param indent_level:
        indentation level (number of leading blanks)

    :return:
        C++ code as list of lines

    """

    return [f'{" " * indent_level}{sym} = {_print_with_exception(math)};'
            f'  // {variable}[{index}]'.replace('\n',
                                                '\n' + ' ' * indent_level)
            for index, (sym, math) in enumerate(zip(symbols, equations))
            if not (math == 0 or math == 0.0)]


class ODEExporter:
    """
    The ODEExporter class generates AMICI C++ files for ODE model as
    defined in symbolic expressions.

    :ivar model:
        ODE definition

    :ivar outdir:
        see :meth:`amici.ode_export.ODEExporter.set_paths`

    :ivar verbose:
        more verbose output if True

    :ivar assume_pow_positivity:
        if set to true, a special pow function is
        used to avoid problems with state variables that may become negative
        due to numerical errors

        compiler: distutils/setuptools compiler selection to build the
        python extension

    :ivar functions:
        carries C++ function signatures and other specifications

    :ivar model_name:
        name of the model that will be used for compilation

    :ivar model_path:
        path to the generated model specific files

    :ivar model_swig_path:
        path to the generated swig files

    :ivar allow_reinit_fixpar_initcond:
        indicates whether reinitialization of
        initial states depending on fixedParameters is allowed for this model

    :ivar _build_hints:
        If the given model uses special functions, this set contains hints for
        model building.

    :ivar generate_sensitivity_code:
        Specifies whether code for sensitivity computation is to be generated
    """

    def __init__(
            self,
            ode_model: ODEModel,
            outdir: Optional[str] = None,
            verbose: Optional[Union[bool, int]] = False,
            assume_pow_positivity: Optional[bool] = False,
            compiler: Optional[str] = None,
            allow_reinit_fixpar_initcond: Optional[bool] = True,
            generate_sensitivity_code: Optional[bool] = True
    ):
        """
        Generate AMICI C++ files for the ODE provided to the constructor.

        :param ode_model:
            ODE definition

        :param outdir:
            see :meth:`amici.ode_export.ODEExporter.set_paths`

        :param verbose:
            verbosity level for logging, True/False default to
            logging.Error/logging.DEBUG

        :param assume_pow_positivity:
            if set to true, a special pow function is
            used to avoid problems with state variables that may become
            negative due to numerical errors

        :param compiler: distutils/setuptools compiler selection to build the
            python extension

        :param allow_reinit_fixpar_initcond:
            see :class:`amici.ode_export.ODEExporter`

        :param generate_sensitivity_code specifies whether code required for
            sensitivity computation will be generated
        """
        set_log_level(logger, verbose)

        self.outdir: str = outdir
        self.verbose: bool = logger.getEffectiveLevel() <= logging.DEBUG
        self.assume_pow_positivity: bool = assume_pow_positivity
        self.compiler: str = compiler

        self.model_name: str = 'model'
        output_dir = os.path.join(os.getcwd(),
                                  f'amici-{self.model_name}')
        self.model_path: str = os.path.abspath(output_dir)
        self.model_swig_path: str = os.path.join(self.model_path, 'swig')

        # Signatures and properties of generated model functions (see
        # include/amici/model.h for details)
        self.model: ODEModel = ode_model

        # To only generate a subset of functions, apply subselection here
        self.functions: Dict[str, Dict[str, Union[str, List[str]]]] = \
            copy.deepcopy(functions)

        self.allow_reinit_fixpar_initcond: bool = allow_reinit_fixpar_initcond
        self._build_hints = set()
        self.generate_sensitivity_code: bool = generate_sensitivity_code

    @log_execution_time('generating cpp code', logger)
    def generate_model_code(self) -> None:
        """
        Generates the native C++ code for the loaded model and a Matlab
        script that can be run to compile a mex file from the C++ code


        """
        with _monkeypatched(sp.Pow, '_eval_derivative',
                            _custom_pow_eval_derivative):

            self._prepare_model_folder()
            self._generate_c_code()
            self._generate_m_code()

    @log_execution_time('compiling cpp code', logger)
    def compile_model(self) -> None:
        """
        Compiles the generated code it into a simulatable module


        """
        self._compile_c_code(compiler=self.compiler,
                             verbose=self.verbose)

    def _prepare_model_folder(self) -> None:
        """
        Remove all files from the model folder.
        """
        for file in os.listdir(self.model_path):
            file_path = os.path.join(self.model_path, file)
            if os.path.isfile(file_path):
                os.remove(file_path)

    def _generate_c_code(self) -> None:
        """
        Create C++ code files for the model based on ODEExporter.model
        """
        for function in self.functions.keys():
            if function in sensi_functions + sparse_sensi_functions and \
                    not self.generate_sensitivity_code:
                continue

            if 'dont_generate_body' not in \
                    self.functions[function].get('flags', []):
                dec = log_execution_time(f'writing {function}.cpp', logger)
                dec(self._write_function_file)(function)
            if function in sparse_functions \
                    and 'body' in self.functions[function]:
                self._write_function_index(function, 'colptrs')
                self._write_function_index(function, 'rowvals')

        for name in self.model.sym_names():
            # only generate for those that have nontrivial implementation,
            # check for both basic variables (not in functions) and function
            # computed values
            if (name in self.functions and
                'body' not in self.functions[name] and
                name not in nobody_functions) or \
               (name not in self.functions and
                    len(self.model.sym(name)) == 0):
                continue
            self._write_index_files(name)

        self._write_wrapfunctions_cpp()
        self._write_wrapfunctions_header()
        self._write_model_header_cpp()
        self._write_c_make_file()
        self._write_swig_files()
        self._write_module_setup()

        shutil.copy(CXX_MAIN_TEMPLATE_FILE,
                    os.path.join(self.model_path, 'main.cpp'))

    def _compile_c_code(self,
                        verbose: Optional[Union[bool, int]] = False,
                        compiler: Optional[str] = None) -> None:
        """
        Compile the generated model code

        :param verbose:
            Make model compilation verbose

        :param compiler:
            distutils/setuptools compiler selection to build the python
            extension

        """

        # setup.py assumes it is run from within the model directory
        module_dir = self.model_path
        script_args = [sys.executable, os.path.join(module_dir, 'setup.py')]

        if verbose:
            script_args.append('--verbose')
        else:
            script_args.append('--quiet')

        script_args.extend(['build_ext', f'--build-lib={module_dir}'])

        if compiler is not None:
            script_args.extend([f'--compiler={compiler}'])

        # distutils.core.run_setup looks nicer, but does not let us check the
        # result easily
        try:
            result = subprocess.run(script_args,
                                    cwd=module_dir,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT,
                                    check=True)
        except subprocess.CalledProcessError as e:
            print(e.output.decode('utf-8'))
            print("Failed building the model extension.")
            if self._build_hints:
                print("Note:")
                print('\n'.join(self._build_hints))
            raise

        if verbose:
            print(result.stdout.decode('utf-8'))

    def _generate_m_code(self) -> None:
        """
        Create a Matlab script for compiling code files to a mex file
        """
        # creating the code lines for the Matlab compile script
        lines = []

        # Events are not yet implemented. Once this is done, the variable nz
        # will have to be replaced by "self.model.nz()"
        nz = 0

        # Second order code is not yet implemented. Once this is done,
        # those variables will have to be replaced by
        # "self.model.<var>true()", or the corresponding "model.self.o2flag"
        nxtrue_rdata = self.model.num_states_rdata()
        nytrue = self.model.num_obs()
        o2flag = 0

        # a preliminary comment
        lines.append('% This compile script was automatically created from'
                     ' Python SBML import.')
        lines.append('% If mex compiler is set up within MATLAB, it can be run'
                     ' from MATLAB ')
        lines.append('% in order to compile a mex-file from the Python'
                     ' generated C++ files.')
        lines.append('')

        # write the actual compiling code
        lines.append(f"modelName = '{self.model_name}';")
        lines.append("amimodel.compileAndLinkModel"
                     "(modelName, '', [], [], [], []);")
        lines.append(f"amimodel.generateMatlabWrapper({nxtrue_rdata}, "
                     f"{nytrue}, {self.model.num_par()}, "
                     f"{self.model.num_const()}, {nz}, {o2flag}, ...\n    [], "
                     "['simulate_' modelName '.m'], modelName, ...\n"
                     "    'lin', 1, 1);")

        # write compile script (for mex)
        compile_script = os.path.join(self.model_path, 'compileMexFile.m')
        with open(compile_script, 'w') as fileout:
            fileout.write('\n'.join(lines))

    def _write_index_files(self, name: str) -> None:
        """
        Write index file for a symbolic array.

        :param name:
            key in self.model._syms for which the respective file should
            be written

        """
        lines = []
        if name in self.model.sym_names():
            if name in sparse_functions:
                symbols = self.model.sparsesym(name)
            else:
                symbols = self.model.sym(name).T
            # flatten multiobs
            if isinstance(next(iter(symbols), None), list):
                symbols = [symbol for obs in symbols for symbol in obs]
        else:
            raise ValueError(f'Unknown symbolic array: {name}')

        for index, symbol in enumerate(symbols):
            symbol_name = strip_pysb(symbol)
            if str(symbol) == '0':
                continue
            if str(symbol_name) == '':
                raise ValueError(f'{name} contains a symbol called ""')
            lines.append(
                f'#define {symbol_name} {name}[{index}]'
            )

        filename = os.path.join(self.model_path, f'{self.model_name}_{name}.h')
        with open(filename, 'w') as fileout:
            fileout.write('\n'.join(lines))

    def _write_function_file(self, function: str) -> None:
        """
        Generate equations and write the C++ code for the function
        `function`.

        :param function:
            name of the function to be written (see self.functions)
        """

        # first generate the equations to make sure we have everything we
        # need in subsequent steps
        if function in sparse_functions:
            equations = self.model.sparseeq(function)
        elif not self.allow_reinit_fixpar_initcond \
                and function == 'sx0_fixedParameters':
            # Not required. Will create empty function body.
            equations = sp.Matrix()
        else:
            equations = self.model.eq(function)

        # function header
        lines = [
            '#include "amici/symbolic_functions.h"',
            '#include "amici/defines.h"',
            '#include "sundials/sundials_types.h"',
            '',
            '#include <gsl/gsl-lite.hpp>',
            '#include <array>',
        ]

        # function signature
        signature = self.functions[function]['signature']

        lines.append('')

        # extract symbols that need definitions from signature
        # don't add includes for files that won't be generated.
        # Unfortunately we cannot check for `self.functions[sym]['body']`
        # here since it may not have been generated yet.
        for match in re.findall(
                fr'const (realtype|double) \*([\w]+)[0]*[,\)]+', signature
        ):
            sym = match[1]
            if sym not in self.model.sym_names():
                continue

            if sym in sparse_functions:
                iszero = smart_is_zero_matrix(self.model.sparseeq(sym))
            elif sym in self.functions:
                iszero = smart_is_zero_matrix(self.model.eq(sym))
            else:
                iszero = len(self.model.sym(sym)) == 0

            if iszero:
                continue

            lines.append(f'#include "{self.model_name}_{sym}.h"')

        # include return symbols
        if function in self.model.sym_names() and \
                function not in non_unique_id_symbols:
            lines.append(f'#include "{self.model_name}_{function}.h"')

        lines.extend([
            '',
            'namespace amici {',
            f'namespace model_{self.model_name} {{',
            '',
        ])

        lines.append(f'void {function}_{self.model_name}{signature}{{')

        # function body
        body = self._get_function_body(function, equations)
        if self.assume_pow_positivity and 'assume_pow_positivity' \
                in self.functions[function].get('flags', []):
            body = [re.sub(r'(^|\W)std::pow\(', r'\1amici::pos_pow(', line)
                    for line in body]
            # execute this twice to catch cases where the ending ( would be the
            # starting (^|\W) for the following match
            body = [re.sub(r'(^|\W)std::pow\(', r'\1amici::pos_pow(', line)
                    for line in body]

        if body:
            self.functions[function]['body'] = body
        else:
            return
        lines += body
        lines.extend([
            '}',
            '',
            f'}} // namespace model_{self.model_name}',
            '} // namespace amici\n',
        ])

        # check custom functions
        for fun in CUSTOM_FUNCTIONS:
            if 'include' in fun and any(fun['c++'] in line for line in lines):
                if 'build_hint' in fun:
                    self._build_hints.add(fun['build_hint'])
                lines.insert(0, fun['include'])

        # if not body is None:
        with open(os.path.join(
                self.model_path, f'{self.model_name}_{function}.cpp'), 'w'
        ) as fileout:
            fileout.write('\n'.join(lines))

    def _write_function_index(self, function: str, indextype: str) -> None:
        """
        Generate equations and write the C++ code for the function
        `function`.

        :param function:
            name of the function to be written (see self.functions)

        :param indextype:
            type of index {'colptrs', 'rowvals'}

        """

        if indextype == 'colptrs':
            values = self.model.colptrs(function)
            setter = 'indexptrs'
        elif indextype == 'rowvals':
            values = self.model.rowvals(function)
            setter = 'indexvals'
        else:
            raise ValueError('Invalid value for indextype, must be colptrs or '
                             f'rowvals: {indextype}')

        # function signature
        if function in multiobs_functions:
            signature = f'(SUNMatrixWrapper &{function}, int index)'
        else:
            signature = f'(SUNMatrixWrapper &{function})'

        lines = [
            '#include "amici/sundials_matrix_wrapper.h"',
            '#include "sundials/sundials_types.h"',
            '',
            '#include <array>',
            '#include <algorithm>',
            '',
            'namespace amici {',
            f'namespace model_{self.model_name} {{',
            '',
        ]

        # Generate static array with indices
        if len(values):
            static_array_name = f"{function}_{indextype}_{self.model_name}_"
            if function in multiobs_functions:
                # list of index vectors
                lines.append(
                    "static constexpr std::array<std::array<sunindextype, "
                    f"{len(values[0])}>, {len(values)}> "
                    f"{static_array_name} = {{{{"
                )
                lines.extend(['    {'
                              + ', '.join(map(str, index_vector)) + '}, '
                              for index_vector in values])
                lines.append("}};")
            else:
                # single index vector
                lines.append("static constexpr std::array<sunindextype, "
                             f"{len(values)}> {static_array_name} = {{")
                lines.append('    ' + ', '.join(map(str, values)))
                lines.append("};")

        lines.extend([
            '',
            f'void {function}_{indextype}_{self.model_name}{signature}{{',
        ])

        if len(values):
            if function in multiobs_functions:
                lines.append(
                    f"    {function}.set_{setter}"
                    f"(gsl::make_span({static_array_name}[index]));"
                )
            else:
                lines.append(
                    f"    {function}.set_{setter}"
                    f"(gsl::make_span({static_array_name}));"
                )

        lines.extend([
            '}'
            '',
            f'}} // namespace model_{self.model_name}',
            '} // namespace amici\n',
        ])

        filename = f'{self.model_name}_{function}_{indextype}.cpp'
        filename = os.path.join(self.model_path, filename)

        with open(filename, 'w') as fileout:
            fileout.write('\n'.join(lines))

    def _get_function_body(self,
                           function: str,
                           equations: sp.Matrix) -> List[str]:
        """
        Generate C++ code for body of function `function`.

        :param function:
            name of the function to be written (see self.functions)

        :param equations:
            symbolic definition of the function body

        :return:
            generated C++ code

        """

        lines = []

        if (
                len(equations) == 0
                or (
                    isinstance(equations, (sp.Matrix, sp.ImmutableDenseMatrix))
                    and min(equations.shape) == 0
                )
        ):
            # dJydy is a list
            return lines

        if not self.allow_reinit_fixpar_initcond \
                and function in ['sx0_fixedParameters', 'x0_fixedParameters']:
            return lines

        if function == 'sx0_fixedParameters':
            # here we only want to overwrite values where x0_fixedParameters
            # was applied

            lines.extend([
                # Keep list of indices of fixed parameters occurring in x0
                "    static const std::array<int, "
                + str(len(self.model._x0_fixedParameters_idx))
                + "> _x0_fixedParameters_idxs = {",
                "        "
                + ', '.join(str(x)
                            for x in self.model._x0_fixedParameters_idx),
                "    };",
                "",
                # Set all parameters that are to be reset to 0, so that the
                #  switch statement below only needs to handle non-zero entries
                #  (which usually reduces file size and speeds up
                #  compilation significantly).
                "    for(auto idx: reinitialization_state_idxs) {",
                "        if(std::find(_x0_fixedParameters_idxs.cbegin(), "
                "_x0_fixedParameters_idxs.cend(), idx) != "
                "_x0_fixedParameters_idxs.cend())\n"
                "            sx0_fixedParameters[idx] = 0.0;",
                "    }"])

            cases = dict()
            for ipar in range(self.model.num_par()):
                expressions = []
                for index, formula in zip(
                        self.model._x0_fixedParameters_idx,
                        equations[:, ipar]
                ):
                    if not formula.is_zero:
                        expressions.extend([
                            f'if(std::find('
                            'reinitialization_state_idxs.cbegin(), '
                            f'reinitialization_state_idxs.cend(), {index}) != '
                            'reinitialization_state_idxs.cend())',
                            f'    {function}[{index}] = '
                            f'{_print_with_exception(formula)};'
                        ])
                cases[ipar] = expressions
            lines.extend(get_switch_statement('ip', cases, 1))

        elif function == 'x0_fixedParameters':
            for index, formula in zip(
                    self.model._x0_fixedParameters_idx,
                    equations
            ):
                lines.append(
                    f'    if(std::find(reinitialization_state_idxs.cbegin(), '
                    f'reinitialization_state_idxs.cend(), {index}) != '
                    'reinitialization_state_idxs.cend())\n        '
                    f'{function}[{index}] = '
                    f'{_print_with_exception(formula)};')

        elif function in event_functions:
            cases = {ie: _get_sym_lines_array(equations[ie], function, 0)
                     for ie in range(self.model.num_events())
                     if not smart_is_zero_matrix(equations[ie])}
            lines.extend(get_switch_statement('ie', cases, 1))

        elif function in event_sensi_functions:
            outer_cases = {}
            for ie, inner_equations in enumerate(equations):
                inner_lines = []
                inner_cases = {
                    ipar: _get_sym_lines_array(inner_equations[:, ipar],
                                               function, 0)
                    for ipar in range(self.model.num_par())
                    if not smart_is_zero_matrix(inner_equations[:, ipar])}
                inner_lines.extend(get_switch_statement(
                    'ip', inner_cases, 0))
                outer_cases[ie] = copy.copy(inner_lines)
            lines.extend(get_switch_statement('ie', outer_cases, 1))

        elif function in sensi_functions:
            cases = {ipar: _get_sym_lines_array(equations[:, ipar], function,
                                                0)
                     for ipar in range(self.model.num_par())
                     if not smart_is_zero_matrix(equations[:, ipar])}
            lines.extend(get_switch_statement('ip', cases, 1))

        elif function in multiobs_functions:
            if function == 'dJydy':
                cases = {iobs: _get_sym_lines_array(equations[iobs], function,
                                                    0)
                         for iobs in range(self.model.num_obs())
                         if not smart_is_zero_matrix(equations[iobs])}
            else:
                cases = {
                    iobs: _get_sym_lines_array(equations[:, iobs], function, 0)
                    for iobs in range(self.model.num_obs())
                    if not smart_is_zero_matrix(equations[:, iobs])
                }
            lines.extend(get_switch_statement('iy', cases, 1))

        elif function in self.model.sym_names() \
                and function not in non_unique_id_symbols:
            if function in sparse_functions:
                symbols = self.model.sparsesym(function)
            else:
                symbols = self.model.sym(function, stripped=True)
            lines += _get_sym_lines_symbols(symbols, equations, function, 4)

        else:
            lines += _get_sym_lines_array(equations, function, 4)

        return [line for line in lines if line]

    def _write_wrapfunctions_cpp(self) -> None:
        """
        Write model-specific 'wrapper' file (wrapfunctions.cpp).
        """
        template_data = {'MODELNAME': self.model_name}
        apply_template(
            os.path.join(amiciSrcPath, 'wrapfunctions.template.cpp'),
            os.path.join(self.model_path, 'wrapfunctions.cpp'),
            template_data
        )

    def _write_wrapfunctions_header(self) -> None:
        """
        Write model-specific header file (wrapfunctions.h).
        """
        template_data = {'MODELNAME': str(self.model_name)}
        apply_template(
            os.path.join(amiciSrcPath, 'wrapfunctions.ODE_template.h'),
            os.path.join(self.model_path, 'wrapfunctions.h'),
            template_data
        )

    def _write_model_header_cpp(self) -> None:
        """
        Write model-specific header and cpp file (MODELNAME.{h,cpp}).
        """

        tpl_data = {
            'MODELNAME': str(self.model_name),
            'NX_RDATA': str(self.model.num_states_rdata()),
            'NXTRUE_RDATA': str(self.model.num_states_rdata()),
            'NX_SOLVER': str(self.model.num_states_solver()),
            'NXTRUE_SOLVER': str(self.model.num_states_solver()),
            'NX_SOLVER_REINIT': str(self.model.num_state_reinits()),
            'NY': str(self.model.num_obs()),
            'NYTRUE': str(self.model.num_obs()),
            'NZ': '0',
            'NZTRUE': '0',
            'NEVENT': str(self.model.num_events()),
            'NOBJECTIVE': '1',
            'NW': str(len(self.model.sym('w'))),
            'NDWDP': str(len(self.model.sparsesym(
                'dwdp', force_generate=self.generate_sensitivity_code
            ))),
            'NDWDX': str(len(self.model.sparsesym('dwdx'))),
            'NDWDW': str(len(self.model.sparsesym('dwdw'))),
            'NDXDOTDW': str(len(self.model.sparsesym('dxdotdw'))),
            'NDXDOTDP_EXPLICIT': str(len(self.model.sparsesym(
                'dxdotdp_explicit',
                force_generate=self.generate_sensitivity_code
            ))),
            'NDXDOTDX_EXPLICIT': str(len(self.model.sparsesym(
                'dxdotdx_explicit'))),
            'NDJYDY': 'std::vector<int>{%s}'
                      % ','.join(str(len(x))
                                 for x in self.model.sparsesym('dJydy')),
            'UBW': str(self.model.num_states_solver()),
            'LBW': str(self.model.num_states_solver()),
            'NP': str(self.model.num_par()),
            'NK': str(self.model.num_const()),
            'O2MODE': 'amici::SecondOrderMode::none',
            # using cxxcode ensures proper handling of nan/inf
            'PARAMETERS': _print_with_exception(self.model.val('p'))[1:-1],
            'FIXED_PARAMETERS': _print_with_exception(self.model.val('k'))[
                                1:-1],
            'PARAMETER_NAMES_INITIALIZER_LIST':
                self._get_symbol_name_initializer_list('p'),
            'STATE_NAMES_INITIALIZER_LIST':
                self._get_symbol_name_initializer_list('x_rdata'),
            'FIXED_PARAMETER_NAMES_INITIALIZER_LIST':
                self._get_symbol_name_initializer_list('k'),
            'OBSERVABLE_NAMES_INITIALIZER_LIST':
                self._get_symbol_name_initializer_list('y'),
            'EXPRESSION_NAMES_INITIALIZER_LIST':
                self._get_symbol_name_initializer_list('w'),
            'PARAMETER_IDS_INITIALIZER_LIST':
                self._get_symbol_id_initializer_list('p'),
            'STATE_IDS_INITIALIZER_LIST':
                self._get_symbol_id_initializer_list('x_rdata'),
            'FIXED_PARAMETER_IDS_INITIALIZER_LIST':
                self._get_symbol_id_initializer_list('k'),
            'OBSERVABLE_IDS_INITIALIZER_LIST':
                self._get_symbol_id_initializer_list('y'),
            'EXPRESSION_IDS_INITIALIZER_LIST':
                self._get_symbol_id_initializer_list('w'),
            'REINIT_FIXPAR_INITCOND':
                'true' if self.allow_reinit_fixpar_initcond else
                'false',
            'AMICI_VERSION_STRING':  __version__,
            'AMICI_COMMIT_STRING': __commit__,
            'W_RECURSION_DEPTH': self.model._w_recursion_depth,
            'QUADRATIC_LLH': 'true'
                if self.model._has_quadratic_nllh else 'false',
        }

        for fun, fundef in self.functions.items():
            if fun in nobody_functions:
                continue

            if 'body' not in fundef:
                tpl_data[f'{fun.upper()}_DEF'] = ''

                if fun in sensi_functions + sparse_sensi_functions and \
                        not self.generate_sensitivity_code:
                    impl = ''
                else:
                    impl = get_model_override_implementation(
                        fun, self.model_name, nobody=True
                    )

                tpl_data[f'{fun.upper()}_IMPL'] = impl

                if fun in sparse_functions:
                    for indexfield in ['colptrs', 'rowvals']:
                        if fun in sparse_sensi_functions and \
                                not self.generate_sensitivity_code:
                            impl = ''
                        else:
                            impl = get_sunindex_override_implementation(
                                fun, self.model_name, indexfield, nobody=True
                            )
                        tpl_data[f'{fun.upper()}_{indexfield.upper()}_DEF'] \
                            = ''
                        tpl_data[f'{fun.upper()}_{indexfield.upper()}_IMPL'] \
                            = impl

                continue

            tpl_data[f'{fun.upper()}_DEF'] = \
                get_function_extern_declaration(fun, self.model_name)
            tpl_data[f'{fun.upper()}_IMPL'] = \
                get_model_override_implementation(fun, self.model_name)
            if fun in sparse_functions:
                tpl_data[f'{fun.upper()}_COLPTRS_DEF'] = \
                    get_sunindex_extern_declaration(fun, self.model_name,
                                                    'colptrs')
                tpl_data[f'{fun.upper()}_COLPTRS_IMPL'] = \
                    get_sunindex_override_implementation(fun, self.model_name,
                                                         'colptrs')
                tpl_data[f'{fun.upper()}_ROWVALS_DEF'] = \
                    get_sunindex_extern_declaration(fun, self.model_name,
                                                    'rowvals')
                tpl_data[f'{fun.upper()}_ROWVALS_IMPL'] = \
                    get_sunindex_override_implementation(fun, self.model_name,
                                                         'rowvals')

        if self.model.num_states_solver() == self.model.num_states_rdata():
            tpl_data['X_RDATA_DEF'] = ''
            tpl_data['X_RDATA_IMPL'] = ''

        apply_template(
            os.path.join(amiciSrcPath, 'model_header.ODE_template.h'),
            os.path.join(self.model_path, f'{self.model_name}.h'),
            tpl_data
        )

        apply_template(
            os.path.join(amiciSrcPath, 'model.ODE_template.cpp'),
            os.path.join(self.model_path, f'{self.model_name}.cpp'),
            tpl_data
        )

    def _get_symbol_name_initializer_list(self, name: str) -> str:
        """
        Get SBML name initializer list for vector of names for the given
        model entity

        :param name:
            any key present in self.model._syms

        :return:
            Template initializer list of names
        """
        return '\n'.join(
            [
                f'"{symbol}", // {name}[{idx}]'
                for idx, symbol in enumerate(self.model.name(name))
            ]
        )

    def _get_symbol_id_initializer_list(self, name: str) -> str:
        """
        Get C++ initializer list for vector of names for the given model
        entity

        :param name:
            any key present in self.model._syms

        :return:
            Template initializer list of ids
        """
        return '\n'.join(
            [
                f'"{strip_pysb(symbol)}", // {name}[{idx}]'
                for idx, symbol in enumerate(self.model.sym(name))
            ]
        )

    def _write_c_make_file(self):
        """
        Write CMake CMakeLists.txt file for this model.
        """

        sources = [
            f + ' ' for f in os.listdir(self.model_path)
            if f.endswith('.cpp') and f != 'main.cpp'
        ]

        template_data = {'MODELNAME': self.model_name,
                         'SOURCES': '\n'.join(sources),
                         'AMICI_VERSION': __version__}
        apply_template(
            MODEL_CMAKE_TEMPLATE_FILE,
            os.path.join(self.model_path, 'CMakeLists.txt'),
            template_data
        )

    def _write_swig_files(self) -> None:
        """
        Write SWIG interface files for this model.
        """
        if not os.path.exists(self.model_swig_path):
            os.makedirs(self.model_swig_path)
        template_data = {'MODELNAME': self.model_name}
        apply_template(
            os.path.join(amiciSwigPath, 'modelname.template.i'),
            os.path.join(self.model_swig_path, self.model_name + '.i'),
            template_data
        )
        shutil.copy(SWIG_CMAKE_TEMPLATE_FILE,
                    os.path.join(self.model_swig_path, 'CMakeLists.txt'))

    def _write_module_setup(self) -> None:
        """
        Create a distutils setup.py file for compile the model module.
        """

        template_data = {'MODELNAME': self.model_name,
                         'AMICI_VERSION': __version__,
                         'PACKAGE_VERSION': '0.1.0'}
        apply_template(os.path.join(amiciModulePath, 'setup.template.py'),
                       os.path.join(self.model_path, 'setup.py'),
                       template_data)
        apply_template(os.path.join(amiciModulePath, 'MANIFEST.template.in'),
                       os.path.join(self.model_path, 'MANIFEST.in'), {})
        # write __init__.py for the model module
        if not os.path.exists(os.path.join(self.model_path, self.model_name)):
            os.makedirs(os.path.join(self.model_path, self.model_name))

        apply_template(
            os.path.join(amiciModulePath, '__init__.template.py'),
            os.path.join(self.model_path, self.model_name, '__init__.py'),
            template_data
        )

    def set_paths(self, output_dir: str) -> None:
        """
        Set output paths for the model and create if necessary

        :param output_dir:
            relative or absolute path where the generated model
            code is to be placed. will be created if does not exists.

        """
        self.model_path = os.path.abspath(output_dir)
        self.model_swig_path = os.path.join(self.model_path, 'swig')

        for directory in [self.model_path, self.model_swig_path]:
            if not os.path.exists(directory):
                os.makedirs(directory)

    def set_name(self, model_name: str) -> None:
        """
        Sets the model name

        :param model_name:
            name of the model (may only contain upper and lower case letters,
            digits and underscores, and must not start with a digit)

        """
        if not is_valid_identifier(model_name):
            raise ValueError(
                f"'{model_name}' is not a valid model name. "
                "Model name may only contain upper and lower case letters, "
                "digits and underscores, and must not start with a digit.")

        self.model_name = model_name


class TemplateAmici(Template):
    """
    Template format used in AMICI (see string.template for more details).

    :ivar delimiter:
        delimiter that identifies template variables

    """
    delimiter = 'TPL_'


def apply_template(source_file: str,
                   target_file: str,
                   template_data: Dict[str, str]) -> None:
    """
    Load source file, apply template substitution as provided in
    templateData and save as targetFile.

    :param source_file:
        relative or absolute path to template file

    :param target_file:
        relative or absolute path to output file

    :param template_data:
        template keywords to substitute (key is template
        variable without :attr:`TemplateAmici.delimiter`)

    """
    with open(source_file) as filein:
        src = TemplateAmici(filein.read())
    result = src.safe_substitute(template_data)
    with open(target_file, 'w') as fileout:
        fileout.write(result)


def strip_pysb(symbol: sp.Basic) -> sp.Basic:
    """
    Strips pysb info from a :class:`pysb.Component` object

    :param symbol:
        symbolic expression

    :return:
        stripped expression

    """
    # strip pysb type and transform into a flat sympy.Symbol.
    # this ensures that the pysb type specific __repr__ is used when converting
    # to string
    if pysb and isinstance(symbol, pysb.Component):
        return sp.Symbol(symbol.name, real=True)
    else:
        # in this case we will use sympy specific transform anyways
        return symbol


def get_function_extern_declaration(fun: str, name: str) -> str:
    """
    Constructs the extern function declaration for a given function

    :param fun:
        function name
    :param name:
        model name

    :return:
        c++ function definition string

    """
    return \
        f'extern void {fun}_{name}{functions[fun]["signature"]};'


def get_sunindex_extern_declaration(fun: str, name: str,
                                    indextype: str) -> str:
    """
    Constructs the function declaration for an index function of a given
    function

    :param fun:
        function name

    :param name:
        model name

    :param indextype:
        index function {'colptrs', 'rowvals'}

    :return:
        c++ function declaration string

    """
    index_arg = ', int index' if fun in multiobs_functions else ''
    return \
        f'extern void {fun}_{indextype}_{name}' \
        f'(SUNMatrixWrapper &{indextype}{index_arg});'


def get_model_override_implementation(fun: str, name: str,
                                      nobody: bool = False) -> str:
    """
    Constructs amici::Model::* override implementation for a given function

    :param fun:
        function name

    :param name:
        model name

    :param nobody:
        whether the function has a nontrivial implementation

    :return:
        c++ function implementation string

    """
    impl = 'virtual void f{fun}{signature} override {{'

    if nobody:
        impl += '}}\n'
    else:
        impl += '\n{ind8}{fun}_{name}{eval_signature};\n{ind4}}}\n'

    return impl.format(
            ind4=' '*4,
            ind8=' '*8,
            fun=fun,
            name=name,
            signature=functions[fun]["signature"],
            eval_signature=remove_typedefs(functions[fun]["signature"])
        )


def get_sunindex_override_implementation(fun: str, name: str,
                                         indextype: str,
                                         nobody: bool = False) -> str:
    """
    Constructs the amici::Model:: function implementation for an index
    function of a given function

    :param fun:
        function name

    :param name:
        model name

    :param indextype:
        index function {'colptrs', 'rowvals'}

    :param nobody:
        whether the corresponding function has a nontrivial implementation

    :return:
        c++ function implementation string

    """
    index_arg = ', int index' if fun in multiobs_functions else ''
    index_arg_eval = ', index' if fun in multiobs_functions else ''

    impl = 'virtual void f{fun}_{indextype}{signature} override {{'

    if nobody:
        impl += '}}\n'
    else:
        impl += '{ind8}{fun}_{indextype}_{name}{eval_signature};\n{ind4}}}\n'

    return impl.format(
            ind4=' '*4,
            ind8=' '*8,
            fun=fun,
            indextype=indextype,
            name=name,
            signature=f'(SUNMatrixWrapper &{indextype}{index_arg})',
            eval_signature=f'({indextype}{index_arg_eval})',
        )


def remove_typedefs(signature: str) -> str:
    """
    Strips typedef info from a function signature

    :param signature:
        function signature

    :return:
        string that can be used to construct function calls with the same
        variable names and ordering as in the function signature
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
        'gsl::span<const int>'
    ]

    for typedef in typedefs:
        signature = signature.replace(typedef, '')

    return signature


def get_switch_statement(condition: str, cases: Dict[int, List[str]],
                         indentation_level: Optional[int] = 0,
                         indentation_step: Optional[str] = ' ' * 4):
    """
    Generate code for switch statement

    :param condition:
        Condition for switch

    :param cases:
        Cases as dict with expressions as keys and statement as
        list of strings

    :param indentation_level:
        indentation level

    :param indentation_step:
        indentation whitespace per level

    :return:
        Code for switch expression as list of strings

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


def csc_matrix(matrix: sp.Matrix,
               rownames: List[sp.Symbol],
               colnames: List[sp.Symbol],
               identifier: Optional[int] = 0,
               pattern_only: Optional[bool] = False) -> Tuple[
    List[int], List[int], sp.Matrix, List[str], sp.Matrix
]:
    """
    Generates the sparse symbolic identifiers, symbolic identifiers,
    sparse matrix, column pointers and row values for a symbolic
    variable

    :param matrix:
        dense matrix to be sparsified

    :param rownames:
        ids of the variable of which the derivative is computed (assuming
        matrix is the jacobian)

    :param colnames:
        ids of the variable with respect to which the derivative is computed
        (assuming matrix is the jacobian)

    :param identifier:
        additional identifier that gets appended to symbol names to
        ensure their uniqueness in outer loops

    :param pattern_only:
        flag for computing sparsity pattern without whole matrix

    :return:
        symbol_col_ptrs, symbol_row_vals, sparse_list, symbol_list,
        sparse_matrix

    """
    idx = 0

    nrows, ncols = matrix.shape

    if not pattern_only:
        sparse_matrix = sp.zeros(nrows, ncols)
    symbol_list = []
    sparse_list = []
    symbol_col_ptrs = []
    symbol_row_vals = []

    for col in range(0, ncols):
        symbol_col_ptrs.append(idx)
        for row in range(0, nrows):
            if matrix[row, col] == 0:
                continue

            symbol_row_vals.append(row)
            idx += 1
            symbol_name = f'd{_print_with_exception(rownames[row])}' \
                          f'_d{_print_with_exception(colnames[col])}'
            if identifier:
                symbol_name += f'_{identifier}'
            symbol_list.append(symbol_name)
            if pattern_only:
                continue

            sparse_matrix[row, col] = sp.Symbol(symbol_name, real=True)
            sparse_list.append(matrix[row, col])

    if idx == 0:
        symbol_col_ptrs = []  # avoid bad memory access for empty matrices
    else:
        symbol_col_ptrs.append(idx)

    if pattern_only:
        sparse_matrix = None
    else:
        sparse_list = sp.Matrix(sparse_list)

    return symbol_col_ptrs, symbol_row_vals, sparse_list, symbol_list, \
        sparse_matrix


def is_valid_identifier(x: str) -> bool:
    """
    Check whether `x` is a valid identifier for conditions, parameters,
    observables... . Identifiers may only contain upper and lower case letters,
    digits and underscores, and must not start with a digit.

    :param x:
        string to check

    :return:
        ``True`` if valid, ``False`` otherwise
    """

    return re.match(r'^[a-zA-Z_]\w*$', x) is not None


def generate_measurement_symbol(observable_id: Union[str, sp.Symbol]):
    """
    Generates the appropriate measurement symbol for the provided observable

    :param observable_id:
        symbol (or string representation) of the observable

    :return:
        symbol for the corresponding measurement
    """
    if not isinstance(observable_id, str):
        observable_id = strip_pysb(observable_id)
    return symbol_with_assumptions(f'm{observable_id}')


def generate_flux_symbol(reaction_index: int) -> sp.Symbol:
    """
    Generate identifier symbol for a reaction flux.
    This function will always return the same unique python object for a
    given entity.

    :param reaction_index:
        index of the reaction to which the flux corresponds
    :return:
        identifier symbol
    """
    return symbol_with_assumptions(f'flux_r{reaction_index}')


def symbol_with_assumptions(name: str):
    """
    Central function to create symbols with consistent, canonical assumptions

    :param name:
        name of the symbol

    :return:
        symbol with canonical assumptions
    """
    return sp.Symbol(name, real=True)


def cast_to_sym(value: Union[SupportsFloat, sp.Expr, BooleanAtom],
                input_name: str) -> sp.Expr:
    """
    Typecasts the value to sympy.Float if possible, and ensures the
    value is a symbolic expression.

    :param value:
        value to be cast

    :param input_name:
        name of input variable

    :return:
        typecast value
    """
    if isinstance(value, (sp.RealNumber, numbers.Number)):
        value = sp.Float(float(value))
    elif isinstance(value, BooleanAtom):
        value = sp.Float(float(bool(value)))

    if not isinstance(value, sp.Expr):
        raise TypeError(f"Couldn't cast {input_name} to sympy.Expr, was "
                        f"{type(value)}")

    return value


@contextlib.contextmanager
def _monkeypatched(obj: object, name: str, patch: Any):
    """
    Temporarily monkeypatches an object.

    :param obj:
        object to be patched

    :param name:
        name of the attribute to be patched

    :param patch:
        patched value

    """
    pre_patched_value = getattr(obj, name)
    setattr(obj, name, patch)
    try:
        yield object
    finally:
        setattr(obj, name, pre_patched_value)


def _custom_pow_eval_derivative(self, s):
    """
    Custom Pow derivative that removes a removeable singularity for
    self.base == 0 and self.base.diff(s) == 0. This function is intended to
    be monkeypatched into sp.Pow._eval_derivative.

    :param self:
        sp.Pow class

    :param s:
        variable with respect to which the derivative will be computed

    """
    dbase = self.base.diff(s)
    dexp = self.exp.diff(s)
    part1 = sp.Pow(self.base, self.exp - 1) * self.exp * dbase
    part2 = self * dexp * sp.log(self.base)
    if self.base.is_nonzero or dbase.is_nonzero or part2.is_zero:
        # first piece never applies or is zero anyways
        return part1 + part2

    return part1 + sp.Piecewise(
        (self.base, sp.And(sp.Eq(self.base, 0), sp.Eq(dbase, 0))),
        (part2, True)
    )


def _custom_print_max(self, expr):
    """
    Custom Max printing function, see https://github.com/sympy/sympy/pull/20558
    """
    from sympy import Max
    if len(expr.args) == 1:
        return self._print(expr.args[0])
    return "%smax(%s, %s)" % (self._ns, self._print(expr.args[0]),
                              self._print(Max(*expr.args[1:])))


def _custom_print_min(self, expr):
    """
    Custom Min printing function, see https://github.com/sympy/sympy/pull/20558
    """
    from sympy import Min
    if len(expr.args) == 1:
        return self._print(expr.args[0])
    return "%smin(%s, %s)" % (self._ns, self._print(expr.args[0]),
                              self._print(Min(*expr.args[1:])))
