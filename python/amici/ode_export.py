"""
C++ Export
----------
This module provides all necessary functionality specify an ODE model and
generate executable C++ simulation code. The user generally won't have to
directly call any function from this module as this will be done by
:func:`amici.pysb_import.pysb2amici`,
:meth:`amici.sbml_import.SbmlImporter.sbml2amici` and
:func:`amici.petab_import.import_model`
"""
import sympy as sp
import re
import shutil
import subprocess
import sys
import os
import copy
import numbers
import logging
import itertools

try:
    import pysb
except ImportError:
    pysb = None

from typing import (
    Callable, Optional, Union, List, Dict, Tuple, SupportsFloat, Sequence,
    Set
)
from string import Template
import sympy.printing.cxxcode as cxxcode
from sympy.matrices.immutable import ImmutableDenseMatrix
from sympy.matrices.dense import MutableDenseMatrix
from itertools import chain

from . import (
    amiciSwigPath, amiciSrcPath, amiciModulePath, __version__, __commit__,
    sbml_import
)
from .logging import get_logger, log_execution_time, set_log_level

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
    and function != 'sxdot'
]
# list of multiobs functions
multiobs_functions = [
    function for function in functions
    if 'const int iy' in functions[function]['signature']
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
                 value: Union[SupportsFloat, numbers.Number, sp.Basic]):
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

    def get_val(self) -> Union[float, sp.Basic]:
        """
        ModelQuantity value

        :return:
            value of the ModelQuantity
        """
        return self._value


class State(ModelQuantity):
    """
    A State variable defines an entity that evolves with time according to
    the provided time derivative, abbreviated by `x`

    :ivar conservation_law:
        algebraic formula that allows computation of this
        species according to a conservation law

    """

    conservation_law: Union[sp.Basic, None] = None

    def __init__(self,
                 identifier: sp.Symbol,
                 name: str,
                 value: sp.Basic,
                 dt: sp.Basic):
        """
        Create a new State instance. Extends :meth:`ModelQuantity.__init__`
        by dt

        :param identifier:
            unique identifier of the state

        :param name:
            individual name of the state (does not need to be unique)

        :param value:
            initial value

        :param dt:
            time derivative
        """
        super(State, self).__init__(identifier, name, value)
        if not isinstance(dt, sp.Basic):
            raise TypeError(f'dt must have type sympy.Basic, was '
                            f'{type(dt)}')

        self._dt = dt
        self.conservation_law = None

    def set_conservation_law(self,
                             law: sp.Basic) -> None:
        """
        Sets the conservation law of a state. If the a conservation law
        is set, the respective state will be replaced by an algebraic
        formula according to the respective conservation law.

        :param law:
            linear sum of states that if added to this state remain
            constant over time
        """
        if not isinstance(law, sp.Basic):
            raise TypeError(f'conservation law must have type sympy.Basic, '
                            f'was {type(law)}')

        self.conservation_law = law

    def set_dt(self,
               dt: sp.Basic) -> None:
        """
        Sets the time derivative

        :param dt:
            time derivative
        """
        if not isinstance(dt, sp.Basic):
            raise TypeError(f'time derivative must have type sympy.Basic, '
                            f'was {type(dt)}')
        self._dt = dt

    def get_dt(self) -> sp.Basic:
        """
        Gets the time derivative

        :return:
            time derivative
        """
        return self._dt

    def get_free_symbols(self) -> Set[sp.Basic]:
        """
        Gets the set of free symbols in time derivative and inital conditions

        :return:
            free symbols
        """
        symbols = self._dt.free_symbols
        if isinstance(self._value, sp.Basic):
            symbols = symbols.union(self._value.free_symbols)

        return symbols


class ConservationLaw(ModelQuantity):
    """
    A conservation law defines the absolute the total amount of a
    (weighted) sum of states

    """
    def __init__(self,
                 identifier: sp.Symbol,
                 name: str,
                 value: sp.Basic):
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
    """
    def __init__(self,
                 identifier: sp.Symbol,
                 name: str,
                 value: sp.Basic):
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


class SigmaY(ModelQuantity):
    """
    A Standard Deviation SigmaY rescales the distance between simulations
    and measurements when computing residuals or objective functions,
    abbreviated by `sigmay`
    """
    def __init__(self,
                 identifier: sp.Symbol,
                 name: str,
                 value: sp.Basic):
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
                 value: sp.Basic):
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
    sensitivites may be computed, abbreviated by `p`
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
    sensitivites cannot be computed, abbreviated by `k`
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
                 value: sp.Basic):
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


# defines the type of some attributes in ODEModel
symbol_to_type = {
    'species': State,
    'parameter': Parameter,
    'fixed_parameter': Constant,
    'observable': Observable,
    'sigmay': SigmaY,
    'llhy': LogLikelihood,
    'expression': Expression,
}


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


def smart_multiply(x: sp.MutableDenseMatrix,
                   y: sp.MutableDenseMatrix) -> sp.MutableDenseMatrix:
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
    return x * y


def smart_is_zero_matrix(x: sp.MutableDenseMatrix) -> bool:
    """A faster implementation of sympy's is_zero_matrix

    Avoids repeated indexer type checks and double iteration to distinguish
    False/None. Found to be about 100x faster for large matrices.

    :param x: Matrix to check
    """

    return not any(xx.is_zero is not True for xx in x._mat)


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

    :ivar _multiplication_prototypes:
        defines how a multiplication equation is computed for an equation,
        key defines the name and values should be
        arguments for ODEModel.multiplication()

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
    """

    def __init__(self, simplify: Optional[Callable] = sp.powsimp):
        """
        Create a new ODEModel instance.

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
        self._symboldim_funs: Dict[str, Callable[[], int]] = {
            'sx': self.nx_solver,
            'v': self.nx_solver,
            'vB': self.nx_solver,
            'xB': self.nx_solver,
            'sigmay': self.ny,
        }
        self._eqs: Dict[str, sp.Matrix] = dict()
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
            'sigmay': '_sigmays',
        }
        self._variable_prototype: Dict[str, str] = {
            'tcl': '_conservationlaws',
            'x_rdata': '_states',
            'y': '_observables',
            'p': '_parameters',
            'k': '_constants',
            'w': '_expressions',
            'sigmay': '_sigmays'
        }
        self._value_prototype: Dict[str, str] = {
            'p': '_parameters',
            'k': '_constants',
        }
        self._total_derivative_prototypes: \
            Dict[str, Dict[str, Union[str, List[str]]]] = {
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
        self._multiplication_prototypes: Dict[str, Dict[str, str]] = {
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

        self._lock_total_derivative: List[str] = list()
        self._simplify: Callable = simplify
        self._x0_fixedParameters_idx: Union[None, Sequence[int]]

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

        # assemble fluxes and add them as expressions to the model
        fluxes = []
        for ir, flux in enumerate(si.flux_vector):
            flux_id = sp.Symbol(f'flux_r{ir}', real=True)
            self.add_component(Expression(
                identifier=flux_id,
                name=str(flux),
                value=flux
            ))
            fluxes.append(flux_id)
        nr = len(fluxes)

        # correct time derivatives for compartment changes

        dxdotdw_updates = []
        def dx_dt(x_index, x_Sw):
            '''
            Produces the appropriate expression for the first derivative of a
            species with respect to time, for species that reside in
            compartments with a constant volume, or a volume that is defined by
            an assignment or rate rule.

            :param x_index:
                The index (not identifier) of the species in the variables
                (generated in "sbml_import.py") that describe the model.

            :param x_Sw:
                The element-wise product of the row in the stoichiometric
                matrix that corresponds to the species (row x_index) and the
                flux (kinetic laws) vector. Ignored in the case of rate rules.
            '''
            x_id = symbols['species']['identifier'][x_index]

            # Rate rules specify dx_dt.
            # Note that the rate rule of species may describe amount, not
            # concentration.
            if x_id in si.compartment_rate_rules:
                return si.compartment_rate_rules[x_id]
            elif x_id in si.species_rate_rules:
                return si.species_rate_rules[x_id]

            # The derivation of the below return expressions can be found in
            # the documentation. They are found by rearranging
            # $\frac{d}{dt} (vx) = Sw$ for $\frac{dx}{dt}$, where $v$ is the
            # vector of species compartment volumes, $x$ is the vector of
            # species concentrations, $S$ is the stoichiometric matrix, and $w$
            # is the flux vector. The conditional below handles the cases of
            # species in (i) compartments with a rate rule, (ii) compartments
            # with an assignment rule, and (iii) compartments with a constant
            # volume, respectively.
            v_name = si.species_compartment[x_index]
            if v_name in si.compartment_rate_rules:
                dv_dt = si.compartment_rate_rules[v_name]
                xdot = (x_Sw - dv_dt*x_id)/v_name
                for w_index, flux in enumerate(fluxes):
                    dxdotdw_updates.append((x_index, w_index, xdot.diff(flux)))
                return xdot
            elif v_name in si.compartment_assignment_rules:
                v = si.compartment_assignment_rules[v_name]
                dv_dt = v.diff(si.amici_time_symbol)
                dv_dx = v.diff(x_id)
                xdot = (x_Sw - dv_dt*x_id)/(dv_dx*x_id + v)
                for w_index, flux in enumerate(fluxes):
                    dxdotdw_updates.append((x_index, w_index, xdot.diff(flux)))
                return xdot
            else:
                v = si.compartment_volume[list(si.compartment_symbols).index(
                    si.species_compartment[x_index])]
                for w_index, flux in enumerate(fluxes):
                    if si.stoichiometric_matrix[x_index, w_index] != 0:
                        dxdotdw_updates.append((x_index, w_index,
                            si.stoichiometric_matrix[x_index, w_index] / v))
                return x_Sw/v

        # create dynamics without respecting conservation laws first
        Sw = smart_multiply(MutableDenseMatrix(si.stoichiometric_matrix),
                            MutableDenseMatrix(fluxes))
        symbols['species']['dt'] = sp.Matrix([Sw.row(x_index).applyfunc(
            lambda x_Sw: dx_dt(x_index, x_Sw))
            for x_index in range(Sw.rows)])

        # create all basic components of the ODE model and add them.
        for symbol in [s for s in symbols if s != 'my']:
            # transform dict of lists into a list of dicts
            protos = [dict(zip(symbols[symbol], t))
                      for t in zip(*symbols[symbol].values())]
            for proto in protos:
                self.add_component(symbol_to_type[symbol](**proto))

        # process conservation laws
        if compute_cls:
            dxdotdw_updates = si.process_conservation_laws(self,
                                                           dxdotdw_updates)

        # set derivatives of xdot, this circumvents regular computation. we
        # do this as we can save a substantial amount of computations by
        # knowing the right solutions here
        nx_solver = si.stoichiometric_matrix.shape[0]
        nw = len(self._expressions)
        # append zero rows for conservation law `w`s, note that
        # _process_conservation_laws is called after the fluxes are added as
        # expressions, if this ordering needs to be changed, this will have
        # to be adapted.
        self._eqs['dxdotdw'] = si.stoichiometric_matrix.row_join(
            sp.zeros(nx_solver, nw-nr)
        )
        for ix, iw, val in dxdotdw_updates:
            self._eqs['dxdotdw'][ix, iw] = val

        # fill in 'self._sym' based on prototypes and components in ode_model
        self.generate_basic_variables()

    def add_component(self, component: ModelQuantity) -> None:
        """
        Adds a new ModelQuantity to the model.

        :param component:
            model quantity to be added
        """
        for comp_type in [Observable, Expression, Parameter, Constant, State,
                          LogLikelihood, SigmaY, ConservationLaw]:
            if isinstance(component, comp_type):
                getattr(self, f'_{type(component).__name__.lower()}s').append(
                    component
                )
                return

        raise ValueError(f'Invalid component type {type(component)}')

    def add_conservation_law(self,
                             state: sp.Symbol,
                             total_abundance: sp.Symbol,
                             state_expr: sp.Basic,
                             abundance_expr: sp.Basic) -> None:
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

    def nx_rdata(self) -> int:
        """
        Number of states.

        :return:
            number of state variable symbols
        """
        return len(self.sym('x_rdata'))

    def nx_solver(self) -> int:
        """
        Number of states after applying conservation laws.

        :return:
            number of state variable symbols
        """
        return len(self.sym('x'))

    def ncl(self) -> int:
        """
        Number of conservation laws.

        :return:
            number of conservation laws
        """
        return self.nx_rdata()-self.nx_solver()

    def nx_solver_reinit(self) -> int:
        """
        Number of solver states which would be reinitialized after
        preequilibraiton

        :return:
            number of state variable symbols with reinitialization
        """
        reinit_states = self.eq('x0_fixedParameters')
        solver_states = self.eq('x_solver')
        return sum([1 for ix in reinit_states if ix in solver_states])

    def ny(self) -> int:
        """
        Number of Observables.

        :return:
            number of observable symbols
        """
        return len(self.sym('y'))

    def nk(self) -> int:
        """
        Number of Constants.

        :return:
            number of constant symbols
        """
        return len(self.sym('k'))

    def np(self) -> int:
        """
        Number of Parameters.

        :return:
            number of parameter symbols
        """
        return len(self.sym('p'))

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

        if stripped:
            if name not in self._variable_prototype:
                raise ValueError('Stripped symbols only available for '
                                 'variables from variable prototypes')
            return self._strippedsyms[name]
        else:
            return self._syms[name]

    def sparsesym(self, name: str) -> List[str]:
        """
        Returns (and constructs if necessary) the sparsified identifiers for
        a sparsified symbolic variable.

        :param name:
            name of the symbolic variable

        :return:
            linearized Matrix containing the symbolic identifiers

        """
        if name not in sparse_functions:
            raise ValueError(f'{name} is not marked as sparse')
        if name not in self._sparsesyms:
            self._generate_sparse_symbol(name)
        return self._sparsesyms[name]

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

    def _generate_symbol(self, name: str) -> None:
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
        elif name == 'sx0':
            self._syms[name] = sp.Matrix([
                f's{state.get_id()}_0'
                for state in self._states
                if state.conservation_law is None
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
                else [0] * self.np()
                for tcl in self._conservationlaws
            ])
            return
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

    def generate_basic_variables(self) -> None:
        """
        Generates the symbolic identifiers for all variables in
        ODEModel.variable_prototype
        """
        for var in self._variable_prototype:
            if var not in self._syms:
                self._generate_symbol(var)

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

        if name == 'dJydy':
            # One entry per y-slice
            self._colptrs[name] = []
            self._rowvals[name] = []
            self._sparseeqs[name] = []
            self._sparsesyms[name] = []
            self._syms[name] = []
            base_index = 0
            for iy in range(self.ny()):
                symbol_col_ptrs, symbol_row_vals, sparse_list, symbol_list, \
                    sparse_matrix = csc_matrix(matrix[iy, :], name,
                                               base_index=base_index)
                base_index += len(symbol_list)
                self._colptrs[name].append(symbol_col_ptrs)
                self._rowvals[name].append(symbol_row_vals)
                self._sparseeqs[name].append(sparse_list)
                self._sparsesyms[name].append(symbol_list)
                self._syms[name].append(sparse_matrix)
        else:
            symbol_col_ptrs, symbol_row_vals, sparse_list, symbol_list, \
                sparse_matrix = csc_matrix(
                    matrix, name, pattern_only=name in nobody_functions
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
        match_deriv = re.match(r'd([\w_]+)d([a-z_]+)', name)

        if name in self._equation_prototype:
            self._equation_from_component(name, self._equation_prototype[name])

        elif name in self._total_derivative_prototypes:
            args = self._total_derivative_prototypes[name]
            args['name'] = name
            self._lock_total_derivative += args['chainvars']
            self._total_derivative(**args)
            for cv in args['chainvars']:
                self._lock_total_derivative.remove(cv)

        elif name in self._multiplication_prototypes:
            args = self._multiplication_prototypes[name]
            args['name'] = name
            self._multiplication(**args)

        elif name == 'xdot':
            self._eqs[name] = sp.Matrix([
                s.get_dt() for s in self._states
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

        elif name == 'JB':
            self._eqs[name] = -self.eq('J').transpose()

        elif name == 'JDiag':
            self._eqs[name] = get_symbolic_diagonal(self.eq('J'))

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
            raise ValueError(f'Unknown equation {name}')

        if name in ['Jy', 'dydx']:
            # do not transpose if we compute the partial derivative as part of
            # a total derivative
            if not len(self._lock_total_derivative):
                self._eqs[name] = self._eqs[name].transpose()

        if self._simplify:
            self._eqs[name] = self._eqs[name].applyfunc(self._simplify)

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
            ('xdot', 'p'): 'w'  # has generic implementation in c++ code
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

        self._eqs[name] = smart_jacobian(eq, sym_var)

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
        if var_in_function_signature(name, varname):
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
        # flatten conservation laws in expressions
        if name == 'w':
            self._eqs[name] = self._eqs[name].subs(
                self.get_conservation_laws()
            )

    def get_conservation_laws(self) -> List[Tuple[sp.Symbol, sp.Basic]]:
        """ Returns a list of states with conservation law set


        :return:
            list of state identifiers

        """
        return [
            (state.get_id(), state.conservation_law)
            for state in self._states
            if state.conservation_law is not None
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
        return self._states[ix].conservation_law is not None

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


def _print_with_exception(math: sp.Basic) -> str:
    """
    Generate C++ code for a symbolic expression

    :param math:
        symbolic expression

    :return:
        C++ code for the specified expression
    """
    try:
        ret = cxxcode(math, standard='c++11')
        ret = re.sub(r'(^|\W)M_PI(\W|$)', r'\1amici::pi\2', ret)
        return ret
    except TypeError as e:
        raise ValueError(
            f'Encountered unsupported function in expression "{math}": '
            f'{e}!'
        )


def _get_sym_lines(symbols: sp.Matrix,
                   variable: str,
                   indent_level: int) -> List[str]:
    """
    Generate C++ code for assigning symbolic terms in symbols to C++ array
    `variable`.

    :param symbols:
        vectors of symbolic terms

    :param variable:
        name of the C++ array to assign to

    :param indent_level:
        indentation level (number of leading blanks)

    :return:
        C++ code as list of lines

    """

    return [' ' * indent_level + f'{variable}[{index}] = '
                                 f'{_print_with_exception(math)};'
            for index, math in enumerate(symbols)
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

    """

    def __init__(
            self,
            ode_model: ODEModel,
            outdir: Optional[str] = None,
            verbose: Optional[Union[bool, int]] = False,
            assume_pow_positivity: Optional[bool] = False,
            compiler: Optional[str] = None,
            allow_reinit_fixpar_initcond: Optional[bool] = True
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

    @log_execution_time('generating cpp code', logger)
    def generate_model_code(self) -> None:
        """
        Generates the native C++ code for the loaded model and a Matlab
        script that can be run to compile a mex file from the C++ code


        """
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
            if 'dont_generate_body' not in \
                    self.functions[function].get('flags', []):
                dec = log_execution_time(f'writing {function}.cpp', logger)
                dec(self._write_function_file)(function)
            if function in sparse_functions:
                self._write_function_index(function, 'colptrs')
                self._write_function_index(function, 'rowvals')

        for name in self.model.sym_names():
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
        nxtrue_rdata = self.model.nx_rdata()
        nytrue = self.model.ny()
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
                     f"{nytrue}, {self.model.np()}, {self.model.nk()}, "
                     f"{nz}, {o2flag}, ...\n    [], "
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
        else:
            raise ValueError(f'Unknown symbolic array: {name}')

        for index, symbol in enumerate(symbols):
            symbol_name = strip_pysb(symbol)
            if str(symbol) == '0':
                continue
            lines.append(
                f'#define {symbol_name} {name}[{index}]'
            )

        with open(os.path.join(self.model_path, f'{name}.h'), 'w') as fileout:
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
        ]

        # function signature
        signature = self.functions[function]['signature']

        lines.append('')

        for sym in self.model.sym_names():
            # added |double for data
            # added '[0]*' for initial conditions
            if re.search(
                    fr'const (realtype|double) \*{sym}[0]*[,)]+', signature
            ):
                lines.append(f'#include "{sym}.h"')

        lines.extend([
            '',
            'namespace amici {',
            f'namespace model_{self.model_name} {{',
            '',
        ])

        lines.append(f'void {function}_{self.model_name}{signature}{{')

        # function body
        body = self._get_function_body(function, symbol)
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
        lines.extend([
            '}',
            '',
            '} // namespace amici',
            f'}} // namespace model_{self.model_name}',
        ])
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
        elif indextype == 'rowvals':
            values = self.model.rowvals(function)
        else:
            raise ValueError('Invalid value for indextype, must be colptrs or '
                             f'rowvals: {indextype}')

        # function signature
        if function in multiobs_functions:
            signature = f'(sunindextype *{indextype}, int index)'
        else:
            signature = f'(sunindextype *{indextype})'

        lines = [
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
                lines.append("static constexpr std::array<std::array<int, "
                             f"{len(values[0])}>, {len(values)}> "
                             f"{static_array_name} = {{{{")
                lines.extend(['    {'
                              + ', '.join(map(str, index_vector)) + '}, '
                              for index_vector in values])
                lines.append("}};")
            else:
                # single index vector
                lines.append("static constexpr std::array<int, "
                             f"{len(values)}> {static_array_name} = {{")
                lines.append('    ' + ', '.join(map(str, values)))
                lines.append("};")

        lines.extend([
            '',
            f'void {function}_{indextype}_{self.model_name}{signature}{{',
        ])

        if len(values):
            if function in multiobs_functions:
                lines.append(f"    std::copy({static_array_name}[index]"
                             f".begin(), {static_array_name}[index].end(), "
                             f"{indextype});")
            else:
                lines.append(f"    std::copy({static_array_name}.begin(), "
                             f"{static_array_name}.end(), {indextype});")

        lines.extend([
            '}'
            '',
            '} // namespace amici',
            f'}} // namespace model_{self.model_name}',
        ])

        filename = f'{self.model_name}_{function}_{indextype}.cpp'
        filename = os.path.join(self.model_path, filename)

        with open(filename, 'w') as fileout:
            fileout.write('\n'.join(lines))

    def _get_function_body(self,
                           function: str,
                           symbol: sp.Matrix) -> List[str]:
        """
        Generate C++ code for body of function `function`.

        :param function:
            name of the function to be written (see self.functions)

        :param symbol:
            symbolic defintion of the function body

        :return:
            generated C++ code

        """

        lines = []

        if len(symbol) == 0 or (isinstance(symbol, (sp.Matrix,
                                                    sp.ImmutableDenseMatrix))
                                and min(symbol.shape) == 0):
            # dJydy is a list
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
                for index, formula in zip(
                        self.model._x0_fixedParameters_idx,
                        symbol[:, ipar]
                ):
                    expressions.append(f'{function}[{index}] = '
                                       f'{_print_with_exception(formula)};')
                cases[ipar] = expressions
            lines.extend(get_switch_statement('ip', cases, 1))

        elif function == 'x0_fixedParameters':
            for index, formula in zip(
                    self.model._x0_fixedParameters_idx,
                    symbol
            ):
                lines.append(f'{function}[{index}] = '
                             f'{_print_with_exception(formula)};')

        elif function in sensi_functions:
            cases = {ipar: _get_sym_lines(symbol[:, ipar], function, 0)
                     for ipar in range(self.model.np())
                     if not smart_is_zero_matrix(symbol[:, ipar])}
            lines.extend(get_switch_statement('ip', cases, 1))

        elif function in multiobs_functions:
            if function == 'dJydy':
                cases = {iobs: _get_sym_lines(symbol[iobs], function, 0)
                         for iobs in range(self.model.ny())
                         if not smart_is_zero_matrix(symbol[iobs])}
            else:
                cases = {iobs: _get_sym_lines(symbol[:, iobs], function, 0)
                         for iobs in range(self.model.ny())
                         if not smart_is_zero_matrix(symbol[:, iobs])}
            lines.extend(get_switch_statement('iy', cases, 1))

        else:
            lines += _get_sym_lines(symbol, function, 4)

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
            'NX_RDATA': str(self.model.nx_rdata()),
            'NXTRUE_RDATA': str(self.model.nx_rdata()),
            'NX_SOLVER': str(self.model.nx_solver()),
            'NXTRUE_SOLVER': str(self.model.nx_solver()),
            'NX_SOLVER_REINIT': str(self.model.nx_solver_reinit()),
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
                self._get_symbol_name_initializer_list('p'),
            'STATE_NAMES_INITIALIZER_LIST':
                self._get_symbol_name_initializer_list('x_rdata'),
            'FIXED_PARAMETER_NAMES_INITIALIZER_LIST':
                self._get_symbol_name_initializer_list('k'),
            'OBSERVABLE_NAMES_INITIALIZER_LIST':
                self._get_symbol_name_initializer_list('y'),
            'PARAMETER_IDS_INITIALIZER_LIST':
                self._get_symbol_id_initializer_list('p'),
            'STATE_IDS_INITIALIZER_LIST':
                self._get_symbol_id_initializer_list('x_rdata'),
            'FIXED_PARAMETER_IDS_INITIALIZER_LIST':
                self._get_symbol_id_initializer_list('k'),
            'OBSERVABLE_IDS_INITIALIZER_LIST':
                self._get_symbol_id_initializer_list('y'),
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

        if self.model.nx_solver() == self.model.nx_rdata():
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
                f'"{symbol}",'
                for symbol in self.model.name(name)
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
                f'"{strip_pysb(symbol)}",'
                for symbol in self.model.sym(name)
            ]
        )

    def _write_c_make_file(self):
        """
        Write CMake CMakeLists.txt file for this model.
        """

        sources = [self.model_name + '_' + function + '.cpp '
                   for function in self.functions.keys()
                   if self.functions[function].get('body', None) is not None]

        # add extra source files for sparse matrices
        for function in sparse_functions:
            sources.append(self.model_name + '_' + function
                           + '_colptrs.cpp')
            sources.append(self.model_name + '_' + function
                           + '_rowvals.cpp ')

        sources.append(f'{self.model_name}.cpp')

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


def get_symbolic_diagonal(matrix: sp.Matrix) -> sp.Matrix:
    """
    Get symbolic matrix with diagonal of matrix `matrix`.

    :param matrix:
        Matrix from which to return the diagonal

    :return:
        A Symbolic matrix with the diagonal of `matrix`.
    """
    if not matrix.cols == matrix.rows:
        raise ValueError('Provided matrix is not square!')

    diagonal = [matrix[index, index] for index in range(matrix.cols)]

    return sp.Matrix(diagonal)


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
        f'(sunindextype *{indextype}{index_arg});'


def get_model_override_implementation(fun: str, name: str) -> str:
    """
    Constructs amici::Model::* override implementation for a given function

    :param fun:
        function name

    :param name:
        model name

    :return:
        c++ function implementation string

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


def get_sunindex_override_implementation(fun: str, name: str,
                                         indextype: str) -> str:
    """
    Constructs the amici::Model:: function implementation for an index
    function of a given function

    :param fun:
        function name

    :param name:
        model name

    :param indextype:
        index function {'colptrs', 'rowvals'}

    :return:
        c++ function implementation string

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
    ]

    for typedef in typedefs:
        signature = signature.replace(typedef, ' ')

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
               name: str,
               base_index: Optional[int] = 0,
               pattern_only: Optional[bool] = False) -> Tuple[
    List[int], List[int], sp.Matrix, List[str], sp.Matrix
]:
    """
    Generates the sparse symbolic identifiers, symbolic identifiers,
    sparse matrix, column pointers and row values for a symbolic
    variable

    :param matrix:
        dense matrix to be sparsified

    :param name:
        name of the symbolic variable

    :param base_index:
        index for first symbol name, defaults to 0

    :param pattern_only:
        flag for computing sparsity pattern without whole matrix

    :return:
        symbol_col_ptrs, symbol_row_vals, sparse_list, symbol_list,
        sparse_matrix

    """
    idx = 0
    symbol_name_idx = base_index

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
            symbol_name = f'{name}{symbol_name_idx}'
            symbol_list.append(symbol_name)
            symbol_name_idx += 1
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
    """Check whether `x` is a valid identifier

    Check whether `x` is a valid identifier for conditions, parameters,
    observables... . Identifiers may only contain upper and lower case letters,
    digits and underscores, and must not start with a digit.

    Arguments:
        x: string to check

    Returns:
        ``True`` if valid, ``False`` otherwise
    """

    return re.match(r'^[a-zA-Z_]\w*$', x) is not None
