"""
C++ Export
----------
This module provides all necessary functionality specify an DE model and
generate executable C++ simulation code. The user generally won't have to
directly call any function from this module as this will be done by
:py:func:`amici.pysb_import.pysb2amici`,
:py:func:`amici.sbml_import.SbmlImporter.sbml2amici` and
:py:func:`amici.petab_import.import_model`.
"""
import contextlib
import copy
import itertools
import logging
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from itertools import chain, starmap
from pathlib import Path
from string import Template
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Literal,
    Optional,
    Union,
)
from collections.abc import Sequence

import numpy as np
import sympy as sp
from sympy.matrices.dense import MutableDenseMatrix
from sympy.matrices.immutable import ImmutableDenseMatrix

from . import (
    __commit__,
    __version__,
    amiciModulePath,
    amiciSrcPath,
    amiciSwigPath,
    splines,
)
from .constants import SymbolId
from .cxxcodeprinter import AmiciCxxCodePrinter, get_switch_statement
from .de_model import *
from .import_utils import (
    ObservableTransformation,
    SBMLException,
    amici_time_symbol,
    generate_flux_symbol,
    smart_subs_dict,
    strip_pysb,
    toposort_symbols,
    unique_preserve_order,
)
from .logging import get_logger, log_execution_time, set_log_level

if TYPE_CHECKING:
    from . import sbml_import


# Template for model simulation main.cpp file
CXX_MAIN_TEMPLATE_FILE = os.path.join(amiciSrcPath, "main.template.cpp")
# Template for model/swig/CMakeLists.txt
SWIG_CMAKE_TEMPLATE_FILE = os.path.join(
    amiciSwigPath, "CMakeLists_model.cmake"
)
# Template for model/CMakeLists.txt
MODEL_CMAKE_TEMPLATE_FILE = os.path.join(
    amiciSrcPath, "CMakeLists.template.cmake"
)

IDENTIFIER_PATTERN = re.compile(r"^[a-zA-Z_]\w*$")
DERIVATIVE_PATTERN = re.compile(r"^d(x_rdata|xdot|\w+?)d(\w+?)(?:_explicit)?$")


@dataclass
class _FunctionInfo:
    """Information on a model-specific generated C++ function

    :ivar ode_arguments: argument list of the ODE function. input variables should be
        ``const``.
    :ivar dae_arguments: argument list of the DAE function, if different from ODE
        function. input variables should be ``const``.
    :ivar return_type: the return type of the function
    :ivar assume_pow_positivity:
        identifies the functions on which ``assume_pow_positivity`` will have
        an effect when specified during model generation. generally these are
        functions that are used for solving the ODE, where negative values may
        negatively affect convergence of the integration algorithm
    :ivar sparse:
        specifies whether the result of this function will be stored in sparse
        format. sparse format means that the function will only return an
        array of nonzero values and not a full matrix.
    :ivar generate_body:
        indicates whether a model-specific implementation is to be generated
    :ivar body:
        the actual function body. will be filled later
    """

    ode_arguments: str = ""
    dae_arguments: str = ""
    return_type: str = "void"
    assume_pow_positivity: bool = False
    sparse: bool = False
    generate_body: bool = True
    body: str = ""

    def arguments(self, ode: bool = True) -> str:
        """Get the arguments for the ODE or DAE function"""
        if ode or not self.dae_arguments:
            return self.ode_arguments
        return self.dae_arguments


# Information on a model-specific generated C++ function
# prototype for generated C++ functions, keys are the names of functions
functions = {
    "Jy": _FunctionInfo(
        "realtype *Jy, const int iy, const realtype *p, "
        "const realtype *k, const realtype *y, const realtype *sigmay, "
        "const realtype *my"
    ),
    "dJydsigma": _FunctionInfo(
        "realtype *dJydsigma, const int iy, const realtype *p, "
        "const realtype *k, const realtype *y, const realtype *sigmay, "
        "const realtype *my"
    ),
    "dJydy": _FunctionInfo(
        "realtype *dJydy, const int iy, const realtype *p, "
        "const realtype *k, const realtype *y, "
        "const realtype *sigmay, const realtype *my",
        sparse=True,
    ),
    "Jz": _FunctionInfo(
        "realtype *Jz, const int iz, const realtype *p, const realtype *k, "
        "const realtype *z, const realtype *sigmaz, const realtype *mz"
    ),
    "dJzdsigma": _FunctionInfo(
        "realtype *dJzdsigma, const int iz, const realtype *p, "
        "const realtype *k, const realtype *z, const realtype *sigmaz, "
        "const realtype *mz"
    ),
    "dJzdz": _FunctionInfo(
        "realtype *dJzdz, const int iz, const realtype *p, "
        "const realtype *k, const realtype *z, const realtype *sigmaz, "
        "const double *mz",
    ),
    "Jrz": _FunctionInfo(
        "realtype *Jrz, const int iz, const realtype *p, "
        "const realtype *k, const realtype *rz, const realtype *sigmaz"
    ),
    "dJrzdsigma": _FunctionInfo(
        "realtype *dJrzdsigma, const int iz, const realtype *p, "
        "const realtype *k, const realtype *rz, const realtype *sigmaz"
    ),
    "dJrzdz": _FunctionInfo(
        "realtype *dJrzdz, const int iz, const realtype *p, "
        "const realtype *k, const realtype *rz, const realtype *sigmaz",
    ),
    "root": _FunctionInfo(
        "realtype *root, const realtype t, const realtype *x, "
        "const realtype *p, const realtype *k, const realtype *h, "
        "const realtype *tcl"
    ),
    "dwdp": _FunctionInfo(
        "realtype *dwdp, const realtype t, const realtype *x, "
        "const realtype *p, const realtype *k, const realtype *h, "
        "const realtype *w, const realtype *tcl, const realtype *dtcldp, "
        "const realtype *spl, const realtype *sspl",
        assume_pow_positivity=True,
        sparse=True,
    ),
    "dwdx": _FunctionInfo(
        "realtype *dwdx, const realtype t, const realtype *x, "
        "const realtype *p, const realtype *k, const realtype *h, "
        "const realtype *w, const realtype *tcl, const realtype *spl",
        assume_pow_positivity=True,
        sparse=True,
    ),
    "create_splines": _FunctionInfo(
        "const realtype *p, const realtype *k",
        return_type="std::vector<HermiteSpline>",
    ),
    "spl": _FunctionInfo(generate_body=False),
    "sspl": _FunctionInfo(generate_body=False),
    "spline_values": _FunctionInfo(
        "const realtype *p, const realtype *k", generate_body=False
    ),
    "spline_slopes": _FunctionInfo(
        "const realtype *p, const realtype *k", generate_body=False
    ),
    "dspline_valuesdp": _FunctionInfo(
        "realtype *dspline_valuesdp, const realtype *p, const realtype *k, const int ip"
    ),
    "dspline_slopesdp": _FunctionInfo(
        "realtype *dspline_slopesdp, const realtype *p, const realtype *k, const int ip"
    ),
    "dwdw": _FunctionInfo(
        "realtype *dwdw, const realtype t, const realtype *x, "
        "const realtype *p, const realtype *k, const realtype *h, "
        "const realtype *w, const realtype *tcl",
        assume_pow_positivity=True,
        sparse=True,
    ),
    "dxdotdw": _FunctionInfo(
        "realtype *dxdotdw, const realtype t, const realtype *x, "
        "const realtype *p, const realtype *k, const realtype *h, "
        "const realtype *w",
        "realtype *dxdotdw, const realtype t, const realtype *x, "
        "const realtype *p, const realtype *k, const realtype *h, "
        "const realtype *dx, const realtype *w",
        assume_pow_positivity=True,
        sparse=True,
    ),
    "dxdotdx_explicit": _FunctionInfo(
        "realtype *dxdotdx_explicit, const realtype t, "
        "const realtype *x, const realtype *p, const realtype *k, "
        "const realtype *h, const realtype *w",
        "realtype *dxdotdx_explicit, const realtype t, "
        "const realtype *x, const realtype *p, const realtype *k, "
        "const realtype *h, const realtype *dx, const realtype *w",
        assume_pow_positivity=True,
        sparse=True,
    ),
    "dxdotdp_explicit": _FunctionInfo(
        "realtype *dxdotdp_explicit, const realtype t, "
        "const realtype *x, const realtype *p, const realtype *k, "
        "const realtype *h, const realtype *w",
        "realtype *dxdotdp_explicit, const realtype t, "
        "const realtype *x, const realtype *p, const realtype *k, "
        "const realtype *h, const realtype *dx, const realtype *w",
        assume_pow_positivity=True,
        sparse=True,
    ),
    "dydx": _FunctionInfo(
        "realtype *dydx, const realtype t, const realtype *x, "
        "const realtype *p, const realtype *k, const realtype *h, "
        "const realtype *w, const realtype *dwdx",
    ),
    "dydp": _FunctionInfo(
        "realtype *dydp, const realtype t, const realtype *x, "
        "const realtype *p, const realtype *k, const realtype *h, "
        "const int ip, const realtype *w, const realtype *tcl, "
        "const realtype *dtcldp, const realtype *spl, const realtype *sspl"
    ),
    "dzdx": _FunctionInfo(
        "realtype *dzdx, const int ie, const realtype t, "
        "const realtype *x, const realtype *p, const realtype *k, "
        "const realtype *h",
    ),
    "dzdp": _FunctionInfo(
        "realtype *dzdp, const int ie, const realtype t, "
        "const realtype *x, const realtype *p, const realtype *k, "
        "const realtype *h, const int ip",
    ),
    "drzdx": _FunctionInfo(
        "realtype *drzdx, const int ie, const realtype t, "
        "const realtype *x, const realtype *p, const realtype *k, "
        "const realtype *h",
    ),
    "drzdp": _FunctionInfo(
        "realtype *drzdp, const int ie, const realtype t, "
        "const realtype *x, const realtype *p, const realtype *k, "
        "const realtype *h, const int ip",
    ),
    "dsigmaydy": _FunctionInfo(
        "realtype *dsigmaydy, const realtype t, const realtype *p, "
        "const realtype *k, const realtype *y"
    ),
    "dsigmaydp": _FunctionInfo(
        "realtype *dsigmaydp, const realtype t, const realtype *p, "
        "const realtype *k, const realtype *y, const int ip",
    ),
    "sigmay": _FunctionInfo(
        "realtype *sigmay, const realtype t, const realtype *p, "
        "const realtype *k, const realtype *y",
    ),
    "dsigmazdp": _FunctionInfo(
        "realtype *dsigmazdp, const realtype t, const realtype *p,"
        " const realtype *k, const int ip",
    ),
    "sigmaz": _FunctionInfo(
        "realtype *sigmaz, const realtype t, const realtype *p, "
        "const realtype *k",
    ),
    "sroot": _FunctionInfo(
        "realtype *stau, const realtype t, const realtype *x, "
        "const realtype *p, const realtype *k, const realtype *h, "
        "const realtype *sx, const int ip, const int ie, "
        "const realtype *tcl",
        generate_body=False,
    ),
    "drootdt": _FunctionInfo(generate_body=False),
    "drootdt_total": _FunctionInfo(generate_body=False),
    "drootdp": _FunctionInfo(generate_body=False),
    "drootdx": _FunctionInfo(generate_body=False),
    "stau": _FunctionInfo(
        "realtype *stau, const realtype t, const realtype *x, "
        "const realtype *p, const realtype *k, const realtype *h, "
        "const realtype *tcl, const realtype *sx, const int ip, "
        "const int ie"
    ),
    "deltax": _FunctionInfo(
        "double *deltax, const realtype t, const realtype *x, "
        "const realtype *p, const realtype *k, const realtype *h, "
        "const int ie, const realtype *xdot, const realtype *xdot_old"
    ),
    "ddeltaxdx": _FunctionInfo(generate_body=False),
    "ddeltaxdt": _FunctionInfo(generate_body=False),
    "ddeltaxdp": _FunctionInfo(generate_body=False),
    "deltasx": _FunctionInfo(
        "realtype *deltasx, const realtype t, const realtype *x, "
        "const realtype *p, const realtype *k, const realtype *h, "
        "const realtype *w, const int ip, const int ie, "
        "const realtype *xdot, const realtype *xdot_old, "
        "const realtype *sx, const realtype *stau, const realtype *tcl"
    ),
    "w": _FunctionInfo(
        "realtype *w, const realtype t, const realtype *x, "
        "const realtype *p, const realtype *k, "
        "const realtype *h, const realtype *tcl, const realtype *spl",
        assume_pow_positivity=True,
    ),
    "x0": _FunctionInfo(
        "realtype *x0, const realtype t, const realtype *p, "
        "const realtype *k"
    ),
    "x0_fixedParameters": _FunctionInfo(
        "realtype *x0_fixedParameters, const realtype t, "
        "const realtype *p, const realtype *k, "
        "gsl::span<const int> reinitialization_state_idxs",
    ),
    "sx0": _FunctionInfo(
        "realtype *sx0, const realtype t, const realtype *x, "
        "const realtype *p, const realtype *k, const int ip",
    ),
    "sx0_fixedParameters": _FunctionInfo(
        "realtype *sx0_fixedParameters, const realtype t, "
        "const realtype *x0, const realtype *p, const realtype *k, "
        "const int ip, gsl::span<const int> reinitialization_state_idxs",
    ),
    "xdot": _FunctionInfo(
        "realtype *xdot, const realtype t, const realtype *x, "
        "const realtype *p, const realtype *k, const realtype *h, "
        "const realtype *w",
        "realtype *xdot, const realtype t, const realtype *x, "
        "const realtype *p, const realtype *k, const realtype *h, "
        "const realtype *dx, const realtype *w",
        assume_pow_positivity=True,
    ),
    "xdot_old": _FunctionInfo(generate_body=False),
    "y": _FunctionInfo(
        "realtype *y, const realtype t, const realtype *x, "
        "const realtype *p, const realtype *k, "
        "const realtype *h, const realtype *w",
    ),
    "x_rdata": _FunctionInfo(
        "realtype *x_rdata, const realtype *x, const realtype *tcl, "
        "const realtype *p, const realtype *k"
    ),
    "total_cl": _FunctionInfo(
        "realtype *total_cl, const realtype *x_rdata, "
        "const realtype *p, const realtype *k"
    ),
    "dtotal_cldp": _FunctionInfo(
        "realtype *dtotal_cldp, const realtype *x_rdata, "
        "const realtype *p, const realtype *k, const int ip"
    ),
    "dtotal_cldx_rdata": _FunctionInfo(
        "realtype *dtotal_cldx_rdata, const realtype *x_rdata, "
        "const realtype *p, const realtype *k, const realtype *tcl",
        sparse=True,
    ),
    "x_solver": _FunctionInfo("realtype *x_solver, const realtype *x_rdata"),
    "dx_rdatadx_solver": _FunctionInfo(
        "realtype *dx_rdatadx_solver, const realtype *x, "
        "const realtype *tcl, const realtype *p, const realtype *k",
        sparse=True,
    ),
    "dx_rdatadp": _FunctionInfo(
        "realtype *dx_rdatadp, const realtype *x, "
        "const realtype *tcl, const realtype *p, const realtype *k, "
        "const int ip"
    ),
    "dx_rdatadtcl": _FunctionInfo(
        "realtype *dx_rdatadtcl, const realtype *x, "
        "const realtype *tcl, const realtype *p, const realtype *k",
        sparse=True,
    ),
    "z": _FunctionInfo(
        "realtype *z, const int ie, const realtype t, const realtype *x, "
        "const realtype *p, const realtype *k, const realtype *h"
    ),
    "rz": _FunctionInfo(
        "realtype *rz, const int ie, const realtype t, const realtype *x, "
        "const realtype *p, const realtype *k, const realtype *h"
    ),
}

# list of sparse functions
sparse_functions = [
    func_name for func_name, func_info in functions.items() if func_info.sparse
]
# list of nobody functions
nobody_functions = [
    func_name
    for func_name, func_info in functions.items()
    if not func_info.generate_body
]
# list of sensitivity functions
sensi_functions = [
    func_name
    for func_name, func_info in functions.items()
    if "const int ip" in func_info.arguments()
]
# list of sensitivity functions
sparse_sensi_functions = [
    func_name
    for func_name, func_info in functions.items()
    if "const int ip" not in func_info.arguments()
    and func_name.endswith("dp")
    or func_name.endswith("dp_explicit")
]
# list of event functions
event_functions = [
    func_name
    for func_name, func_info in functions.items()
    if "const int ie" in func_info.arguments()
    and "const int ip" not in func_info.arguments()
]
event_sensi_functions = [
    func_name
    for func_name, func_info in functions.items()
    if "const int ie" in func_info.arguments()
    and "const int ip" in func_info.arguments()
]
# list of multiobs functions
multiobs_functions = [
    func_name
    for func_name, func_info in functions.items()
    if "const int iy" in func_info.arguments()
    or "const int iz" in func_info.arguments()
]
# list of equations that have ids which may not be unique
non_unique_id_symbols = ["x_rdata", "y"]

# custom c++ function replacements
CUSTOM_FUNCTIONS = [
    {
        "sympy": "polygamma",
        "c++": "boost::math::polygamma",
        "include": "#include <boost/math/special_functions/polygamma.hpp>",
        "build_hint": "Using polygamma requires libboost-math header files.",
    },
    {"sympy": "Heaviside", "c++": "amici::heaviside"},
    {"sympy": "DiracDelta", "c++": "amici::dirac"},
]

# python log manager
logger = get_logger(__name__, logging.ERROR)


def var_in_function_signature(name: str, varname: str, ode: bool) -> bool:
    """
    Checks if the values for a symbolic variable is passed in the signature
    of a function

    :param name:
        name of the function
    :param varname:
        name of the symbolic variable
    :param ode:
        whether to check the ODE or DAE signature

    :return:
        boolean indicating whether the variable occurs in the function
        signature
    """
    return name in functions and re.search(
        rf"const (realtype|double) \*{varname}[0]*(,|$)+",
        functions[name].arguments(ode=ode),
    )


# defines the type of some attributes in DEModel
symbol_to_type = {
    SymbolId.SPECIES: DifferentialState,
    SymbolId.ALGEBRAIC_STATE: AlgebraicState,
    SymbolId.ALGEBRAIC_EQUATION: AlgebraicEquation,
    SymbolId.PARAMETER: Parameter,
    SymbolId.FIXED_PARAMETER: Constant,
    SymbolId.OBSERVABLE: Observable,
    SymbolId.EVENT_OBSERVABLE: EventObservable,
    SymbolId.SIGMAY: SigmaY,
    SymbolId.SIGMAZ: SigmaZ,
    SymbolId.LLHY: LogLikelihoodY,
    SymbolId.LLHZ: LogLikelihoodZ,
    SymbolId.LLHRZ: LogLikelihoodRZ,
    SymbolId.EXPRESSION: Expression,
    SymbolId.EVENT: Event,
}


@log_execution_time("running smart_jacobian", logger)
def smart_jacobian(
    eq: sp.MutableDenseMatrix, sym_var: sp.MutableDenseMatrix
) -> sp.MutableSparseMatrix:
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
    nrow = eq.shape[0]
    ncol = sym_var.shape[0]
    if (
        not min(eq.shape)
        or not min(sym_var.shape)
        or smart_is_zero_matrix(eq)
        or smart_is_zero_matrix(sym_var)
    ):
        return sp.MutableSparseMatrix(nrow, ncol, dict())

    # preprocess sparsity pattern
    elements = (
        (i, j, a, b)
        for i, a in enumerate(eq)
        for j, b in enumerate(sym_var)
        if a.has(b)
    )

    if (n_procs := int(os.environ.get("AMICI_IMPORT_NPROCS", 1))) == 1:
        # serial
        return sp.MutableSparseMatrix(
            nrow, ncol, dict(starmap(_jacobian_element, elements))
        )

    # parallel
    from multiprocessing import get_context

    # "spawn" should avoid potential deadlocks occurring with fork
    #  see e.g. https://stackoverflow.com/a/66113051
    ctx = get_context("spawn")
    with ctx.Pool(n_procs) as p:
        mapped = p.starmap(_jacobian_element, elements)
    return sp.MutableSparseMatrix(nrow, ncol, dict(mapped))


@log_execution_time("running smart_multiply", logger)
def smart_multiply(
    x: Union[sp.MutableDenseMatrix, sp.MutableSparseMatrix],
    y: sp.MutableDenseMatrix,
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
    if (
        not x.shape[0]
        or not y.shape[1]
        or smart_is_zero_matrix(x)
        or smart_is_zero_matrix(y)
    ):
        return sp.zeros(x.shape[0], y.shape[1])
    return x.multiply(y)


def smart_is_zero_matrix(
    x: Union[sp.MutableDenseMatrix, sp.MutableSparseMatrix],
) -> bool:
    """A faster implementation of sympy's is_zero_matrix

    Avoids repeated indexer type checks and double iteration to distinguish
    False/None. Found to be about 100x faster for large matrices.

    :param x: Matrix to check
    """

    if isinstance(x, sp.MutableDenseMatrix):
        return all(xx.is_zero is True for xx in x.flat())

    if isinstance(x, list):
        return all(smart_is_zero_matrix(xx) for xx in x)

    return x.nnz() == 0


def _default_simplify(x):
    """Default simplification applied in DEModel"""
    # We need this as a free function instead of a lambda to have it picklable
    #  for parallel simplification
    return sp.powsimp(x, deep=True)


class DEModel:
    """
    Defines a Differential Equation as set of ModelQuantities.
    This class provides general purpose interfaces to compute arbitrary
    symbolic derivatives that are necessary for model simulation or
    sensitivity computation.

    :ivar _differential_states:
        list of differential state variables

    :ivar _algebraic_states:
        list of algebraic state variables

    :ivar _observables:
        list of observables

    :ivar _event_observables:
        list of event observables

    :ivar _sigma_ys:
        list of sigmas for observables

    :ivar _sigma_zs:
        list of sigmas for event observables

    :ivar _parameters:
        list of parameters

    :ivar _log_likelihood_ys:
        list of loglikelihoods for observables

    :ivar _log_likelihood_zs:
        list of loglikelihoods for event observables

    :ivar _log_likelihood_rzs:
        list of loglikelihoods for event observable regularizations

    :ivar _expressions:
        list of expressions instances

    :ivar _conservation_laws:
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
        carries the names of symbolic identifiers of the symbolic variables
        of the model

    :ivar _syms:
        carries symbolic identifiers of the symbolic variables of the
        model

    :ivar _sparsesyms:
        carries linear list of all symbolic identifiers for sparsified
        variables

    :ivar _colptrs:
        carries column pointers for sparsified variables. See
        SUNMatrixContent_Sparse definition in ``sunmatrix/sunmatrix_sparse.h``

    :ivar _rowvals:
        carries row values for sparsified variables. See
        SUNMatrixContent_Sparse definition in ``sunmatrix/sunmatrix_sparse.h``

    :ivar _equation_prototype:
        defines the attribute from which an equation should be generated via
        list comprehension (see :meth:`OEModel._generate_equation`)

    :ivar _variable_prototype:
        defines the attribute from which a variable should be generated via
        list comprehension (see :meth:`DEModel._generate_symbol`)

    :ivar _value_prototype:
        defines the attribute from which a value should be generated via
        list comprehension (see :meth:`DEModel._generate_value`)

    :ivar _total_derivative_prototypes:
        defines how a total derivative equation is computed for an equation,
        key defines the name and values should be arguments for
        :meth:`DEModel.totalDerivative`

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

    :ivar _code_printer:
        Code printer to generate C++ code

    :ivar _z2event:
        list of event indices for each event observable
    """

    def __init__(
        self,
        verbose: Optional[Union[bool, int]] = False,
        simplify: Optional[Callable] = _default_simplify,
        cache_simplify: bool = False,
    ):
        """
        Create a new DEModel instance.

        :param verbose:
            verbosity level for logging, True/False default to
            ``logging.DEBUG``/``logging.ERROR``

        :param simplify:
            see :meth:`DEModel._simplify`

        :param cache_simplify:
            Whether to cache calls to the simplify method. Can e.g. decrease
            import times for models with events.
        """
        self._differential_states: list[DifferentialState] = []
        self._algebraic_states: list[AlgebraicState] = []
        self._algebraic_equations: list[AlgebraicEquation] = []
        self._observables: list[Observable] = []
        self._event_observables: list[EventObservable] = []
        self._sigma_ys: list[SigmaY] = []
        self._sigma_zs: list[SigmaZ] = []
        self._parameters: list[Parameter] = []
        self._constants: list[Constant] = []
        self._log_likelihood_ys: list[LogLikelihoodY] = []
        self._log_likelihood_zs: list[LogLikelihoodZ] = []
        self._log_likelihood_rzs: list[LogLikelihoodRZ] = []
        self._expressions: list[Expression] = []
        self._conservation_laws: list[ConservationLaw] = []
        self._events: list[Event] = []
        self.splines = []
        self._symboldim_funs: dict[str, Callable[[], int]] = {
            "sx": self.num_states_solver,
            "v": self.num_states_solver,
            "vB": self.num_states_solver,
            "xB": self.num_states_solver,
            "sigmay": self.num_obs,
            "sigmaz": self.num_eventobs,
        }
        self._eqs: dict[
            str,
            Union[
                sp.Matrix,
                sp.SparseMatrix,
                list[Union[sp.Matrix, sp.SparseMatrix]],
            ],
        ] = dict()
        self._sparseeqs: dict[str, Union[sp.Matrix, list[sp.Matrix]]] = dict()
        self._vals: dict[str, list[sp.Expr]] = dict()
        self._names: dict[str, list[str]] = dict()
        self._syms: dict[str, Union[sp.Matrix, list[sp.Matrix]]] = dict()
        self._sparsesyms: dict[str, Union[list[str], list[list[str]]]] = dict()
        self._colptrs: dict[str, Union[list[int], list[list[int]]]] = dict()
        self._rowvals: dict[str, Union[list[int], list[list[int]]]] = dict()

        self._equation_prototype: dict[str, Callable] = {
            "total_cl": self.conservation_laws,
            "x0": self.states,
            "y": self.observables,
            "Jy": self.log_likelihood_ys,
            "Jz": self.log_likelihood_zs,
            "Jrz": self.log_likelihood_rzs,
            "w": self.expressions,
            "root": self.events,
            "sigmay": self.sigma_ys,
            "sigmaz": self.sigma_zs,
        }
        self._variable_prototype: dict[str, Callable] = {
            "tcl": self.conservation_laws,
            "x_rdata": self.states,
            "y": self.observables,
            "z": self.event_observables,
            "p": self.parameters,
            "k": self.constants,
            "w": self.expressions,
            "sigmay": self.sigma_ys,
            "sigmaz": self.sigma_zs,
            "h": self.events,
        }
        self._value_prototype: dict[str, Callable] = {
            "p": self.parameters,
            "k": self.constants,
        }
        self._total_derivative_prototypes: dict[
            str, dict[str, Union[str, list[str]]]
        ] = {
            "sroot": {
                "eq": "root",
                "chainvars": ["x"],
                "var": "p",
                "dxdz_name": "sx",
            },
        }

        self._lock_total_derivative: list[str] = list()
        self._simplify: Callable = simplify
        if cache_simplify and simplify is not None:

            def cached_simplify(
                expr: sp.Expr,
                _simplified: dict[str, sp.Expr] = {},
                _simplify: Callable = simplify,
            ) -> sp.Expr:
                """Speed up expression simplification with caching.

                NB: This can decrease model import times for models that have
                    many repeated expressions during C++ file generation.
                    For example, this can be useful for models with events.
                    However, for other models, this may increase model import
                    times.

                :param expr:
                    The SymPy expression.
                :param _simplified:
                    The cache.
                :param _simplify:
                    The simplification method.

                :return:
                    The simplified expression.
                """
                expr_str = repr(expr)
                if expr_str not in _simplified:
                    _simplified[expr_str] = _simplify(expr)
                return _simplified[expr_str]

            self._simplify = cached_simplify
        self._x0_fixedParameters_idx: Union[None, Sequence[int]]
        self._w_recursion_depth: int = 0
        self._has_quadratic_nllh: bool = True
        set_log_level(logger, verbose)

        self._code_printer = AmiciCxxCodePrinter()
        for fun in CUSTOM_FUNCTIONS:
            self._code_printer.known_functions[fun["sympy"]] = fun["c++"]

    def differential_states(self) -> list[DifferentialState]:
        """Get all differential states."""
        return self._differential_states

    def algebraic_states(self) -> list[AlgebraicState]:
        """Get all algebraic states."""
        return self._algebraic_states

    def observables(self) -> list[Observable]:
        """Get all observables."""
        return self._observables

    def parameters(self) -> list[Parameter]:
        """Get all parameters."""
        return self._parameters

    def constants(self) -> list[Constant]:
        """Get all constants."""
        return self._constants

    def expressions(self) -> list[Expression]:
        """Get all expressions."""
        return self._expressions

    def events(self) -> list[Event]:
        """Get all events."""
        return self._events

    def event_observables(self) -> list[EventObservable]:
        """Get all event observables."""
        return self._event_observables

    def sigma_ys(self) -> list[SigmaY]:
        """Get all observable sigmas."""
        return self._sigma_ys

    def sigma_zs(self) -> list[SigmaZ]:
        """Get all event observable sigmas."""
        return self._sigma_zs

    def conservation_laws(self) -> list[ConservationLaw]:
        """Get all conservation laws."""
        return self._conservation_laws

    def log_likelihood_ys(self) -> list[LogLikelihoodY]:
        """Get all observable log likelihoodss."""
        return self._log_likelihood_ys

    def log_likelihood_zs(self) -> list[LogLikelihoodZ]:
        """Get all event observable log likelihoods."""
        return self._log_likelihood_zs

    def log_likelihood_rzs(self) -> list[LogLikelihoodRZ]:
        """Get all event observable regularization log likelihoods."""
        return self._log_likelihood_rzs

    def is_ode(self) -> bool:
        """Check if model is ODE model."""
        return len(self._algebraic_equations) == 0

    def states(self) -> list[State]:
        """Get all states."""
        return self._differential_states + self._algebraic_states

    @log_execution_time("importing SbmlImporter", logger)
    def import_from_sbml_importer(
        self,
        si: "sbml_import.SbmlImporter",
        compute_cls: Optional[bool] = True,
    ) -> None:
        """
        Imports a model specification from a
        :class:`amici.sbml_import.SbmlImporter` instance.

        :param si:
            imported SBML model
        :param compute_cls:
            whether to compute conservation laws
        """

        # add splines as expressions to the model
        # saved for later substituting into the fluxes
        spline_subs = {}

        for ispl, spl in enumerate(si.splines):
            spline_expr = spl.ode_model_symbol(si)
            spline_subs[spl.sbml_id] = spline_expr
            self.add_component(
                Expression(
                    identifier=spl.sbml_id,
                    name=str(spl.sbml_id),
                    value=spline_expr,
                )
            )
        self.splines = si.splines

        # get symbolic expression from SBML importers
        symbols = copy.copy(si.symbols)

        # assemble fluxes and add them as expressions to the model
        assert len(si.flux_ids) == len(si.flux_vector)
        fluxes = [
            generate_flux_symbol(ir, name=flux_id)
            for ir, flux_id in enumerate(si.flux_ids)
        ]

        # correct time derivatives for compartment changes
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

            comp = species["compartment"]
            if comp in si.symbols[SymbolId.SPECIES]:
                dv_dt = si.symbols[SymbolId.SPECIES][comp]["dt"]
                xdot = (dxdt - dv_dt * species_id) / comp
                return xdot
            elif comp in si.compartment_assignment_rules:
                v = si.compartment_assignment_rules[comp]

                # we need to flatten out assignments in the compartment in
                # order to ensure that we catch all species dependencies
                v = smart_subs_dict(
                    v, si.symbols[SymbolId.EXPRESSION], "value"
                )
                dv_dt = v.diff(amici_time_symbol)
                # we may end up with a time derivative of the compartment
                # volume due to parameter rate rules
                comp_rate_vars = [
                    p
                    for p in v.free_symbols
                    if p in si.symbols[SymbolId.SPECIES]
                ]
                for var in comp_rate_vars:
                    dv_dt += (
                        v.diff(var) * si.symbols[SymbolId.SPECIES][var]["dt"]
                    )
                dv_dx = v.diff(species_id)
                xdot = (dxdt - dv_dt * species_id) / (dv_dx * species_id + v)
                return xdot
            elif comp in si.symbols[SymbolId.ALGEBRAIC_STATE]:
                raise SBMLException(
                    f"Species {species_id} is in a compartment {comp} that is"
                    f" defined by an algebraic equation. This is not"
                    f" supported."
                )
            else:
                v = si.compartments[comp]

                if v == 1.0:
                    return dxdt

                return dxdt / v

        # create dynamics without respecting conservation laws first
        dxdt = smart_multiply(
            si.stoichiometric_matrix, MutableDenseMatrix(fluxes)
        )
        for ix, ((species_id, species), formula) in enumerate(
            zip(symbols[SymbolId.SPECIES].items(), dxdt)
        ):
            # rate rules and amount species don't need to be updated
            if "dt" in species:
                continue
            if species["amount"]:
                species["dt"] = formula
            else:
                species["dt"] = transform_dxdt_to_concentration(
                    species_id, formula
                )

        # create all basic components of the DE model and add them.
        for symbol_name in symbols:
            # transform dict of lists into a list of dicts
            args = ["name", "identifier"]

            if symbol_name == SymbolId.SPECIES:
                args += ["dt", "init"]
            elif symbol_name == SymbolId.ALGEBRAIC_STATE:
                args += ["init"]
            else:
                args += ["value"]

            if symbol_name == SymbolId.EVENT:
                args += ["state_update", "initial_value"]
            elif symbol_name == SymbolId.OBSERVABLE:
                args += ["transformation"]
            elif symbol_name == SymbolId.EVENT_OBSERVABLE:
                args += ["event"]

            comp_kwargs = [
                {
                    "identifier": var_id,
                    **{k: v for k, v in var.items() if k in args},
                }
                for var_id, var in symbols[symbol_name].items()
            ]

            for comp_kwarg in comp_kwargs:
                self.add_component(symbol_to_type[symbol_name](**comp_kwarg))

        # add fluxes as expressions, this needs to happen after base
        # expressions from symbols have been parsed
        for flux_id, flux in zip(fluxes, si.flux_vector):
            # replace splines inside fluxes
            flux = flux.subs(spline_subs)
            self.add_component(
                Expression(identifier=flux_id, name=str(flux_id), value=flux)
            )

        # process conservation laws
        if compute_cls:
            si.process_conservation_laws(self)

        # fill in 'self._sym' based on prototypes and components in ode_model
        self.generate_basic_variables()
        self._has_quadratic_nllh = all(
            llh["dist"]
            in ["normal", "lin-normal", "log-normal", "log10-normal"]
            for llh in si.symbols[SymbolId.LLHY].values()
        )

        self._process_sbml_rate_of(
            symbols
        )  # substitute SBML-rateOf constructs

    def _process_sbml_rate_of(self, symbols) -> None:
        """Substitute any SBML-rateOf constructs in the model equations"""
        rate_of_func = sp.core.function.UndefinedFunction("rateOf")
        species_sym_to_xdot = dict(zip(self.sym("x"), self.sym("xdot")))
        species_sym_to_idx = {x: i for i, x in enumerate(self.sym("x"))}

        def get_rate(symbol: sp.Symbol):
            """Get rate of change of the given symbol"""
            nonlocal symbols

            if symbol.find(rate_of_func):
                raise SBMLException("Nesting rateOf() is not allowed.")

            # Replace all rateOf(some_species) by their respective xdot equation
            with contextlib.suppress(KeyError):
                return self._eqs["xdot"][species_sym_to_idx[symbol]]

            # For anything other than a state, rateOf(.) is 0 or invalid
            return 0

        # replace rateOf-instances in xdot by xdot symbols
        for i_state in range(len(self.eq("xdot"))):
            if rate_ofs := self._eqs["xdot"][i_state].find(rate_of_func):
                self._eqs["xdot"][i_state] = self._eqs["xdot"][i_state].subs(
                    {
                        # either the rateOf argument is a state, or it's 0
                        rate_of: species_sym_to_xdot.get(rate_of.args[0], 0)
                        for rate_of in rate_ofs
                    }
                )
        # substitute in topological order
        subs = toposort_symbols(dict(zip(self.sym("xdot"), self.eq("xdot"))))
        self._eqs["xdot"] = smart_subs_dict(self.eq("xdot"), subs)

        # replace rateOf-instances in x0 by xdot equation
        for i_state in range(len(self.eq("x0"))):
            if rate_ofs := self._eqs["x0"][i_state].find(rate_of_func):
                self._eqs["x0"][i_state] = self._eqs["x0"][i_state].subs(
                    {
                        rate_of: get_rate(rate_of.args[0])
                        for rate_of in rate_ofs
                    }
                )

        for component in chain(
            self.observables(),
            self.expressions(),
            self.events(),
            self._algebraic_equations,
        ):
            if rate_ofs := component.get_val().find(rate_of_func):
                if isinstance(component, Event):
                    # TODO froot(...) can currently not depend on `w`, so this substitution fails for non-zero rates
                    #  see, e.g., sbml test case 01293
                    raise SBMLException(
                        "AMICI does currently not support rateOf(.) inside event trigger functions."
                    )

                if isinstance(component, AlgebraicEquation):
                    # TODO IDACalcIC fails with
                    #   "The linesearch algorithm failed: step too small or too many backtracks."
                    #  see, e.g., sbml test case 01482
                    raise SBMLException(
                        "AMICI does currently not support rateOf(.) inside AlgebraicRules."
                    )

                component.set_val(
                    component.get_val().subs(
                        {
                            rate_of: get_rate(rate_of.args[0])
                            for rate_of in rate_ofs
                        }
                    )
                )

        for event in self.events():
            if event._state_update is None:
                continue

            for i_state in range(len(event._state_update)):
                if rate_ofs := event._state_update[i_state].find(rate_of_func):
                    raise SBMLException(
                        "AMICI does currently not support rateOf(.) inside event state updates."
                    )
                    # TODO here we need xdot sym, not eqs
                    # event._state_update[i_state] = event._state_update[i_state].subs(
                    #     {rate_of: get_rate(rate_of.args[0]) for rate_of in rate_ofs}
                    # )

    def add_component(
        self, component: ModelQuantity, insert_first: Optional[bool] = False
    ) -> None:
        """
        Adds a new ModelQuantity to the model.

        :param component:
            model quantity to be added

        :param insert_first:
            whether to add quantity first or last, relevant when components
            may refer to other components of the same type.
        """
        if type(component) not in {
            Observable,
            Expression,
            Parameter,
            Constant,
            DifferentialState,
            AlgebraicState,
            AlgebraicEquation,
            LogLikelihoodY,
            LogLikelihoodZ,
            LogLikelihoodRZ,
            SigmaY,
            SigmaZ,
            ConservationLaw,
            Event,
            EventObservable,
        }:
            raise ValueError(f"Invalid component type {type(component)}")

        component_list = getattr(
            self,
            "_"
            + "_".join(
                s.lower()
                for s in re.split(r"([A-Z][^A-Z]+)", type(component).__name__)
                if s
            )
            + "s",
        )
        if insert_first:
            component_list.insert(0, component)
        else:
            component_list.append(component)

    def add_conservation_law(
        self,
        state: sp.Symbol,
        total_abundance: sp.Symbol,
        coefficients: dict[sp.Symbol, sp.Expr],
    ) -> None:
        r"""
        Adds a new conservation law to the model. A conservation law is defined
        by the conserved quantity :math:`T = \sum_i(a_i * x_i)`, where
        :math:`a_i` are coefficients and :math:`x_i` are different state
        variables.

        :param state:
            symbolic identifier of the state that should be replaced by
            the conservation law (:math:`x_j`)

        :param total_abundance:
            symbolic identifier of the total abundance (:math:`T/a_j`)

        :param coefficients:
            Dictionary of coefficients {x_i: a_i}
        """
        try:
            ix = next(
                filter(
                    lambda is_s: is_s[1].get_id() == state,
                    enumerate(self._differential_states),
                )
            )[0]
        except StopIteration:
            raise ValueError(
                f"Specified state {state} was not found in the "
                f"model states."
            )

        state_id = self._differential_states[ix].get_id()

        # \sum_{i≠j}(a_i * x_i)/a_j
        target_expression = (
            sp.Add(
                *(
                    c_i * x_i
                    for x_i, c_i in coefficients.items()
                    if x_i != state
                )
            )
            / coefficients[state]
        )

        # x_j = T/a_j - \sum_{i≠j}(a_i * x_i)/a_j
        state_expr = total_abundance - target_expression

        # T/a_j = \sum_{i≠j}(a_i * x_i)/a_j + x_j
        abundance_expr = target_expression + state_id

        self.add_component(
            Expression(state_id, str(state_id), state_expr), insert_first=True
        )

        cl = ConservationLaw(
            total_abundance,
            f"total_{state_id}",
            abundance_expr,
            coefficients,
            state_id,
        )

        self.add_component(cl)
        self._differential_states[ix].set_conservation_law(cl)

    def get_observable_transformations(self) -> list[ObservableTransformation]:
        """
        List of observable transformations

        :return:
            list of transformations
        """
        return [obs.trafo for obs in self._observables]

    def num_states_rdata(self) -> int:
        """
        Number of states.

        :return:
            number of state variable symbols
        """
        return len(self.sym("x_rdata"))

    def num_states_solver(self) -> int:
        """
        Number of states after applying conservation laws.

        :return:
            number of state variable symbols
        """
        return len(self.sym("x"))

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
        reinit_states = self.eq("x0_fixedParameters")
        solver_states = self.eq("x_solver")
        return sum(ix in solver_states for ix in reinit_states)

    def num_obs(self) -> int:
        """
        Number of Observables.

        :return:
            number of observable symbols
        """
        return len(self.sym("y"))

    def num_eventobs(self) -> int:
        """
        Number of Event Observables.

        :return:
            number of event observable symbols
        """
        return len(self.sym("z"))

    def num_const(self) -> int:
        """
        Number of Constants.

        :return:
            number of constant symbols
        """
        return len(self.sym("k"))

    def num_par(self) -> int:
        """
        Number of Parameters.

        :return:
            number of parameter symbols
        """
        return len(self.sym("p"))

    def num_expr(self) -> int:
        """
        Number of Expressions.

        :return:
            number of expression symbols
        """
        return len(self.sym("w"))

    def num_events(self) -> int:
        """
        Total number of Events (those for which root-functions are added and those without).

        :return:
            number of events
        """
        return len(self.sym("h"))

    def num_events_solver(self) -> int:
        """
        Number of Events.

        :return:
            number of event symbols (length of the root vector in AMICI)
        """
        return sum(
            not event.triggers_at_fixed_timepoint() for event in self.events()
        )

    def sym(self, name: str) -> sp.Matrix:
        """
        Returns (and constructs if necessary) the identifiers for a symbolic
        entity.

        :param name:
            name of the symbolic variable

        :return:
            matrix of symbolic identifiers
        """
        if name not in self._syms:
            self._generate_symbol(name)

        return self._syms[name]

    def sparsesym(self, name: str, force_generate: bool = True) -> list[str]:
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
            raise ValueError(f"{name} is not marked as sparse")
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
            dec = log_execution_time(f"computing {name}", logger)
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
            raise ValueError(f"{name} is not marked as sparse")
        if name not in self._sparseeqs:
            self._generate_sparse_symbol(name)
        return self._sparseeqs[name]

    def colptrs(
        self, name: str
    ) -> Union[list[sp.Number], list[list[sp.Number]]]:
        """
        Returns (and constructs if necessary) the column pointers for
        a sparsified symbolic variable.

        :param name:
            name of the symbolic variable

        :return:
            list containing the column pointers
        """
        if name not in sparse_functions:
            raise ValueError(f"{name} is not marked as sparse")
        if name not in self._sparseeqs:
            self._generate_sparse_symbol(name)
        return self._colptrs[name]

    def rowvals(
        self, name: str
    ) -> Union[list[sp.Number], list[list[sp.Number]]]:
        """
        Returns (and constructs if necessary) the row values for a
        sparsified symbolic variable.

        :param name:
            name of the symbolic variable

        :return:
            list containing the row values
        """
        if name not in sparse_functions:
            raise ValueError(f"{name} is not marked as sparse")
        if name not in self._sparseeqs:
            self._generate_sparse_symbol(name)
        return self._rowvals[name]

    def val(self, name: str) -> list[sp.Number]:
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

    def name(self, name: str) -> list[str]:
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

    def free_symbols(self) -> set[sp.Basic]:
        """
        Returns list of free symbols that appear in RHS and initial
        conditions.
        """
        return set(
            chain.from_iterable(
                state.get_free_symbols()
                for state in self.states() + self.algebraic_equations()
            )
        )

    def _generate_symbol(self, name: str) -> None:
        """
        Generates the symbolic identifiers for a symbolic variable

        :param name:
            name of the symbolic variable
        """
        if name in self._variable_prototype:
            components = self._variable_prototype[name]()
            self._syms[name] = sp.Matrix(
                [comp.get_id() for comp in components]
            )
            if name == "y":
                self._syms["my"] = sp.Matrix(
                    [comp.get_measurement_symbol() for comp in components]
                )
            if name == "z":
                self._syms["mz"] = sp.Matrix(
                    [comp.get_measurement_symbol() for comp in components]
                )
                self._syms["rz"] = sp.Matrix(
                    [comp.get_regularization_symbol() for comp in components]
                )
            return
        elif name == "x":
            self._syms[name] = sp.Matrix(
                [
                    state.get_id()
                    for state in self.states()
                    if not state.has_conservation_law()
                ]
            )
            return
        elif name == "xdot":
            self._syms[name] = sp.Matrix(
                [
                    f"d{x.get_id()}dt" if self.is_ode() else f"de_{ix}"
                    for ix, x in enumerate(self._differential_states)
                    if not x.has_conservation_law()
                ]
                + [f"ae_{ix}" for ix in range(len(self._algebraic_equations))]
            )
            return
        elif name == "dx":
            self._syms[name] = sp.Matrix(
                [
                    f"d{state.get_id()}dt"
                    for state in self.states()
                    if not state.has_conservation_law()
                ]
            )
            return
        elif name == "sx0":
            self._syms[name] = sp.Matrix(
                [
                    f"s{state.get_id()}_0"
                    for state in self.states()
                    if not state.has_conservation_law()
                ]
            )
            return
        elif name == "sx_rdata":
            self._syms[name] = sp.Matrix(
                [f"sx_rdata_{i}" for i in range(len(self.states()))]
            )
            return
        elif name == "dtcldp":
            # check, whether the CL consists of only one state. Then,
            # sensitivities drop out, otherwise generate symbols
            self._syms[name] = sp.Matrix(
                [
                    [
                        sp.Symbol(
                            f"s{strip_pysb(tcl.get_id())}__"
                            f"{strip_pysb(par.get_id())}",
                            real=True,
                        )
                        for par in self._parameters
                    ]
                    if self.conservation_law_has_multispecies(tcl)
                    else [0] * self.num_par()
                    for tcl in self._conservation_laws
                ]
            )
            return
        elif name == "xdot_old":
            length = len(self.eq("xdot"))
        elif name in sparse_functions:
            self._generate_sparse_symbol(name)
            return
        elif name in self._symboldim_funs:
            length = self._symboldim_funs[name]()
        elif name == "stau":
            length = self.eq(name)[0].shape[1]
        elif name in sensi_functions:
            length = self.eq(name).shape[0]
        elif name == "spl":
            # placeholders for the numeric spline values.
            # Need to create symbols
            self._syms[name] = sp.Matrix(
                [[f"spl_{isp}" for isp in range(len(self.splines))]]
            )
            return
        elif name == "sspl":
            # placeholders for spline sensitivities. Need to create symbols
            self._syms[name] = sp.Matrix(
                [
                    [f"sspl_{isp}_{ip}" for ip in range(len(self._syms["p"]))]
                    for isp in range(len(self.splines))
                ]
            )
            return
        else:
            length = len(self.eq(name))
        self._syms[name] = sp.Matrix(
            [
                sp.Symbol(f'{name}{0 if name == "stau" else i}', real=True)
                for i in range(length)
            ]
        )

    def generate_basic_variables(self) -> None:
        """
        Generates the symbolic identifiers for all variables in
        ``DEModel._variable_prototype``
        """
        # We need to process events and Heaviside functions in the ``DEModel`,
        # before adding it to DEExporter
        self.parse_events()

        for var in self._variable_prototype:
            if var not in self._syms:
                self._generate_symbol(var)
        # symbols for spline values need to be created in addition
        for var in ["spl", "sspl"]:
            self._generate_symbol(var)

        self._generate_symbol("x")

    def parse_events(self) -> None:
        """
        This function checks the right-hand side for roots of Heaviside
        functions or events, collects the roots, removes redundant roots,
        and replaces the formulae of the found roots by identifiers of AMICI's
        Heaviside function implementation in the right-hand side
        """
        # Track all roots functions in the right-hand side
        roots = copy.deepcopy(self._events)
        for state in self._differential_states:
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

        # re-order events - first those that require root tracking, then the others
        self._events = list(
            chain(
                itertools.filterfalse(
                    Event.triggers_at_fixed_timepoint, self._events
                ),
                filter(Event.triggers_at_fixed_timepoint, self._events),
            )
        )

    def get_appearance_counts(self, idxs: list[int]) -> list[int]:
        """
        Counts how often a state appears in the time derivative of
        another state and expressions for a subset of states

        :param idxs:
            list of state indices for which counts are to be computed

        :return:
            list of counts for the states ordered according to the provided
            indices
        """
        free_symbols_dt = list(
            itertools.chain.from_iterable(
                [str(symbol) for symbol in state.get_dt().free_symbols]
                for state in self.states()
            )
        )

        free_symbols_expr = list(
            itertools.chain.from_iterable(
                [str(symbol) for symbol in expr.get_val().free_symbols]
                for expr in self._expressions
            )
        )

        return [
            free_symbols_dt.count(str(self._differential_states[idx].get_id()))
            + free_symbols_expr.count(
                str(self._differential_states[idx].get_id())
            )
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

        if match_deriv := DERIVATIVE_PATTERN.match(name):
            eq = match_deriv[1]
            var = match_deriv[2]

            rownames = self.sym(eq)
            colnames = self.sym(var)

        if name == "dJydy":
            # One entry per y-slice
            self._colptrs[name] = []
            self._rowvals[name] = []
            self._sparseeqs[name] = []
            self._sparsesyms[name] = []
            self._syms[name] = []

            for iy in range(self.num_obs()):
                (
                    symbol_col_ptrs,
                    symbol_row_vals,
                    sparse_list,
                    symbol_list,
                    sparse_matrix,
                ) = self._code_printer.csc_matrix(
                    matrix[iy, :],
                    rownames=rownames,
                    colnames=colnames,
                    identifier=iy,
                )
                self._colptrs[name].append(symbol_col_ptrs)
                self._rowvals[name].append(symbol_row_vals)
                self._sparseeqs[name].append(sparse_list)
                self._sparsesyms[name].append(symbol_list)
                self._syms[name].append(sparse_matrix)
        else:
            (
                symbol_col_ptrs,
                symbol_row_vals,
                sparse_list,
                symbol_list,
                sparse_matrix,
            ) = self._code_printer.csc_matrix(
                matrix,
                rownames=rownames,
                colnames=colnames,
                pattern_only=name in nobody_functions,
            )

            self._colptrs[name] = symbol_col_ptrs
            self._rowvals[name] = symbol_row_vals
            self._sparseeqs[name] = sparse_list
            self._sparsesyms[name] = symbol_list
            self._syms[name] = sparse_matrix

    def _compute_equation(self, name: str) -> None:
        """
        Computes the symbolic formula for a symbolic variable

        :param name:
            name of the symbolic variable
        """
        # replacement ensures that we don't have to adapt name in abstract
        # model and keep backwards compatibility with matlab
        match_deriv = DERIVATIVE_PATTERN.match(
            re.sub(r"dJ(y|z|rz)dsigma", r"dJ\1dsigma\1", name)
            .replace("sigmarz", "sigmaz")
            .replace("dJrzdz", "dJrzdrz")
        )
        time_symbol = sp.Matrix([amici_time_symbol])

        if name in self._equation_prototype:
            self._equation_from_components(
                name, self._equation_prototype[name]()
            )

        elif name in self._total_derivative_prototypes:
            args = self._total_derivative_prototypes[name]
            args["name"] = name
            self._lock_total_derivative += args["chainvars"]
            self._total_derivative(**args)
            for cv in args["chainvars"]:
                self._lock_total_derivative.remove(cv)

        elif name == "xdot":
            if self.is_ode():
                self._eqs[name] = sp.Matrix(
                    [
                        state.get_dt()
                        for state in self._differential_states
                        if not state.has_conservation_law()
                    ]
                )
            else:
                self._eqs[name] = sp.Matrix(
                    [
                        x.get_dt() - dx
                        for x, dx in zip(
                            (
                                s
                                for s in self._differential_states
                                if not s.has_conservation_law()
                            ),
                            self.sym("dx"),
                        )
                    ]
                    + [eq.get_val() for eq in self._algebraic_equations]
                )

        elif name == "x_rdata":
            self._eqs[name] = sp.Matrix(
                [state.get_x_rdata() for state in self.states()]
            )

        elif name == "x_solver":
            self._eqs[name] = sp.Matrix(
                [
                    state.get_id()
                    for state in self.states()
                    if not state.has_conservation_law()
                ]
            )

        elif name == "sx_solver":
            self._eqs[name] = sp.Matrix(
                [
                    self.sym("sx_rdata")[ix]
                    for ix, state in enumerate(self.states())
                    if not state.has_conservation_law()
                ]
            )

        elif name == "sx0":
            self._derivative(name[1:], "p", name=name)

        elif name == "sx0_fixedParameters":
            # deltax = -x+x0_fixedParameters if x0_fixedParameters>0 else 0
            # deltasx = -sx+dx0_fixed_parametersdx*sx+dx0_fixedParametersdp
            # if x0_fixedParameters>0 else 0
            # sx0_fixedParameters = sx+deltasx =
            # dx0_fixed_parametersdx*sx+dx0_fixedParametersdp
            self._eqs[name] = smart_jacobian(
                self.eq("x0_fixedParameters"), self.sym("p")
            )

            dx0_fixed_parametersdx = smart_jacobian(
                self.eq("x0_fixedParameters"), self.sym("x")
            )

            if not smart_is_zero_matrix(dx0_fixed_parametersdx):
                if isinstance(self._eqs[name], ImmutableDenseMatrix):
                    self._eqs[name] = MutableDenseMatrix(self._eqs[name])
                tmp = smart_multiply(dx0_fixed_parametersdx, self.sym("sx0"))
                for ip in range(self._eqs[name].shape[1]):
                    self._eqs[name][:, ip] += tmp

        elif name == "x0_fixedParameters":
            k = self.sym("k")
            self._x0_fixedParameters_idx = [
                ix
                for ix, eq in enumerate(self.eq("x0"))
                if any(sym in eq.free_symbols for sym in k)
            ]
            eq = self.eq("x0")
            self._eqs[name] = sp.Matrix(
                [eq[ix] for ix in self._x0_fixedParameters_idx]
            )

        elif name == "dtotal_cldx_rdata":
            x_rdata = self.sym("x_rdata")
            self._eqs[name] = sp.Matrix(
                [
                    [cl.get_ncoeff(xr) for xr in x_rdata]
                    for cl in self._conservation_laws
                ]
            )

        elif name == "dtcldx":
            # this is always zero
            self._eqs[name] = sp.zeros(
                self.num_cons_law(), self.num_states_solver()
            )

        elif name == "dtcldp":
            # force symbols
            self._eqs[name] = self.sym(name)

        elif name == "dx_rdatadx_solver":
            if self.num_cons_law():
                x_solver = self.sym("x")
                self._eqs[name] = sp.Matrix(
                    [
                        [state.get_dx_rdata_dx_solver(xs) for xs in x_solver]
                        for state in self.states()
                    ]
                )
            else:
                # so far, dx_rdatadx_solver is only required for sx_rdata
                # in case of no conservation laws, C++ code will directly use
                # sx, we don't need this
                self._eqs[name] = sp.zeros(
                    self.num_states_rdata(), self.num_states_solver()
                )

        elif name == "dx_rdatadp":
            if self.num_cons_law():
                self._eqs[name] = smart_jacobian(
                    self.eq("x_rdata"), self.sym("p")
                )
            else:
                # so far, dx_rdatadp is only required for sx_rdata
                # in case of no conservation laws, C++ code will directly use
                # sx, we don't need this
                self._eqs[name] = sp.zeros(
                    self.num_states_rdata(), self.num_par()
                )

        elif name == "dx_rdatadtcl":
            self._eqs[name] = smart_jacobian(
                self.eq("x_rdata"), self.sym("tcl")
            )

        elif name == "dxdotdx_explicit":
            # force symbols
            self._derivative("xdot", "x", name=name)

        elif name == "dxdotdp_explicit":
            # force symbols
            self._derivative("xdot", "p", name=name)

        elif name == "spl":
            self._eqs[name] = self.sym(name)

        elif name == "sspl":
            # force symbols
            self._eqs[name] = self.sym(name)

        elif name == "spline_values":
            # force symbols
            self._eqs[name] = sp.Matrix(
                [y for spline in self.splines for y in spline.values_at_nodes]
            )

        elif name == "spline_slopes":
            # force symbols
            self._eqs[name] = sp.Matrix(
                [
                    d
                    for spline in self.splines
                    for d in (
                        sp.zeros(len(spline.derivatives_at_nodes), 1)
                        if spline.derivatives_by_fd
                        else spline.derivatives_at_nodes
                    )
                ]
            )

        elif name == "drootdt":
            self._eqs[name] = smart_jacobian(self.eq("root"), time_symbol)

        elif name == "drootdt_total":
            # backsubstitution of optimized right-hand side terms into RHS
            # calling subs() is costly. Due to looping over events though, the
            # following lines are only evaluated if a model has events
            w_sorted = toposort_symbols(dict(zip(self.sym("w"), self.eq("w"))))
            tmp_xdot = smart_subs_dict(self.eq("xdot"), w_sorted)
            self._eqs[name] = self.eq("drootdt")
            if self.num_states_solver():
                self._eqs[name] += smart_multiply(self.eq("drootdx"), tmp_xdot)

        elif name == "deltax":
            # fill boluses for Heaviside functions, as empty state updates
            # would cause problems when writing the function file later
            event_eqs = []
            for event in self._events:
                if event._state_update is None:
                    event_eqs.append(sp.zeros(self.num_states_solver(), 1))
                else:
                    event_eqs.append(event._state_update)

            self._eqs[name] = event_eqs

        elif name == "z":
            event_observables = [
                sp.zeros(self.num_eventobs(), 1) for _ in self._events
            ]
            event_ids = [e.get_id() for e in self._events]
            # TODO: get rid of this stupid 1-based indexing as soon as we can
            # the matlab interface
            z2event = [
                event_ids.index(event_obs.get_event()) + 1
                for event_obs in self._event_observables
            ]
            for (iz, ie), event_obs in zip(
                enumerate(z2event), self._event_observables
            ):
                event_observables[ie - 1][iz] = event_obs.get_val()

            self._eqs[name] = event_observables
            self._z2event = z2event

        elif name in ["ddeltaxdx", "ddeltaxdp", "ddeltaxdt", "dzdp", "dzdx"]:
            if match_deriv[2] == "t":
                var = time_symbol
            else:
                var = self.sym(match_deriv[2])

            self._eqs[name] = [
                smart_jacobian(self.eq(match_deriv[1])[ie], var)
                for ie in range(self.num_events())
            ]
            if name == "dzdx":
                for ie in range(self.num_events()):
                    dtaudx = (
                        -self.eq("drootdx")[ie, :]
                        / self.eq("drootdt_total")[ie]
                    )
                    for iz in range(self.num_eventobs()):
                        if ie != self._z2event[iz] - 1:
                            continue
                        dzdt = sp.diff(self.eq("z")[ie][iz], time_symbol)
                        self._eqs[name][ie][iz, :] += dzdt * dtaudx

        elif name in ["rz", "drzdx", "drzdp"]:
            eq_events = []
            for ie in range(self.num_events()):
                val = sp.zeros(
                    self.num_eventobs(),
                    1 if name == "rz" else len(self.sym(match_deriv[2])),
                )
                # match event observables to root function
                for iz in range(self.num_eventobs()):
                    if ie == self._z2event[iz] - 1:
                        val[iz, :] = self.eq(name.replace("rz", "root"))[ie, :]
                eq_events.append(val)

            self._eqs[name] = eq_events

        elif name == "stau":
            self._eqs[name] = [
                -self.eq("sroot")[ie, :] / self.eq("drootdt_total")[ie]
                if not self.eq("drootdt_total")[ie].is_zero
                else sp.zeros(*self.eq("sroot")[ie, :].shape)
                for ie in range(self.num_events())
            ]

        elif name == "deltasx":
            if self.num_states_solver() * self.num_par() == 0:
                self._eqs[name] = []
                return

            event_eqs = []
            for ie, event in enumerate(self._events):
                tmp_eq = sp.zeros(self.num_states_solver(), self.num_par())

                # need to check if equations are zero since we are using
                # symbols
                if not smart_is_zero_matrix(
                    self.eq("stau")[ie]
                ) and not smart_is_zero_matrix(self.eq("xdot")):
                    tmp_eq += smart_multiply(
                        self.sym("xdot_old") - self.sym("xdot"),
                        self.sym("stau").T,
                    )

                # only add deltax part if there is state update
                if event._state_update is not None:
                    # partial derivative for the parameters
                    tmp_eq += self.eq("ddeltaxdp")[ie]

                    # initial part of chain rule state variables
                    tmp_dxdp = self.sym("sx") * sp.ones(1, self.num_par())

                    # need to check if equations are zero since we are using
                    # symbols
                    if not smart_is_zero_matrix(self.eq("stau")[ie]):
                        # chain rule for the time point
                        tmp_eq += smart_multiply(
                            self.eq("ddeltaxdt")[ie], self.sym("stau").T
                        )

                        # additional part of chain rule state variables
                        tmp_dxdp += smart_multiply(
                            self.sym("xdot_old"), self.sym("stau").T
                        )

                    # finish chain rule for the state variables
                    tmp_eq += smart_multiply(
                        self.eq("ddeltaxdx")[ie], tmp_dxdp
                    )

                event_eqs.append(tmp_eq)

            self._eqs[name] = event_eqs

        elif name == "xdot_old":
            # force symbols
            self._eqs[name] = self.sym(name)

        elif name == "dwdx":
            x = self.sym("x")
            self._eqs[name] = sp.Matrix(
                [
                    [-cl.get_ncoeff(xs) for xs in x]
                    # the insert first in ode_model._add_conservation_law() means
                    # that we need to reverse the order here
                    for cl in reversed(self._conservation_laws)
                ]
            ).col_join(
                smart_jacobian(self.eq("w")[self.num_cons_law() :, :], x)
            )

        elif match_deriv:
            self._derivative(match_deriv[1], match_deriv[2], name)

        else:
            raise ValueError(f"Unknown equation {name}")

        if name in ("sigmay", "sigmaz"):
            # check for states in sigma{y,z}, which is currently not supported
            syms_x = self.sym("x")
            syms_yz = self.sym(name.removeprefix("sigma"))
            xs_in_sigma = {}
            for sym_yz, eq_yz in zip(syms_yz, self._eqs[name]):
                yz_free_syms = eq_yz.free_symbols
                if tmp := {x for x in syms_x if x in yz_free_syms}:
                    xs_in_sigma[sym_yz] = tmp
            if xs_in_sigma:
                msg = ", ".join(
                    [f"{yz} depends on {xs}" for yz, xs in xs_in_sigma.items()]
                )
                raise NotImplementedError(
                    f"State-dependent observables are not supported, but {msg}."
                )

        if name == "root":
            # Events are processed after the model has been set up.
            # Equations are there, but symbols for roots must be added
            self.sym("h")

        if name in {"Jy", "dydx"}:
            # do not transpose if we compute the partial derivative as part of
            # a total derivative
            if not len(self._lock_total_derivative):
                self._eqs[name] = self._eqs[name].transpose()

        if name in {"dzdx", "drzdx"}:
            self._eqs[name] = [e.T for e in self._eqs[name]]

        if self._simplify:
            dec = log_execution_time(f"simplifying {name}", logger)
            if isinstance(self._eqs[name], list):
                self._eqs[name] = [
                    dec(_parallel_applyfunc)(sub_eq, self._simplify)
                    for sub_eq in self._eqs[name]
                ]
            else:
                self._eqs[name] = dec(_parallel_applyfunc)(
                    self._eqs[name], self._simplify
                )

    def sym_names(self) -> list[str]:
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
            name of resulting symbolic variable, default is ``d{eq}d{var}``
        """
        if not name:
            name = f"d{eq}d{var}"

        ignore_chainrule = {
            ("xdot", "p"): "w",  # has generic implementation in c++ code
            ("xdot", "x"): "w",  # has generic implementation in c++ code
            ("w", "w"): "tcl",  # dtcldw = 0
            ("w", "x"): "tcl",  # dtcldx = 0
        }
        # automatically detect chainrule
        chainvars = [
            cv
            for cv in ["w", "tcl"]
            if var_in_function_signature(eq, cv, self.is_ode())
            and cv not in self._lock_total_derivative
            and var != cv
            and min(self.sym(cv).shape)
            and (
                (eq, var) not in ignore_chainrule
                or ignore_chainrule[(eq, var)] != cv
            )
        ]
        if len(chainvars):
            self._lock_total_derivative += chainvars
            self._total_derivative(name, eq, chainvars, var)
            for cv in chainvars:
                self._lock_total_derivative.remove(cv)
            return

        # partial derivative
        sym_eq = self.eq(eq).transpose() if eq == "Jy" else self.eq(eq)

        sym_var = self.sym(var)

        derivative = smart_jacobian(sym_eq, sym_var)

        self._eqs[name] = derivative

        # compute recursion depth based on nilpotency of jacobian. computing
        # nilpotency can be done more efficiently on numerical sparsity pattern
        if name == "dwdw":
            nonzeros = np.asarray(
                derivative.applyfunc(lambda x: int(not x.is_zero))
            ).astype(np.int64)
            recursion = nonzeros.copy()
            if max(recursion.shape):
                while recursion.max():
                    recursion = recursion.dot(nonzeros)
                    self._w_recursion_depth += 1
                    if self._w_recursion_depth > len(sym_eq):
                        raise RuntimeError(
                            "dwdw is not nilpotent. Something, somewhere went "
                            "terribly wrong. Please file a bug report at "
                            "https://github.com/AMICI-dev/AMICI/issues and "
                            "attach this model."
                        )

        if name == "dydw" and not smart_is_zero_matrix(derivative):
            dwdw = self.eq("dwdw")
            # h(k) = d{eq}dw*dwdw^k* (k=1)
            h = smart_multiply(derivative, dwdw)
            while not smart_is_zero_matrix(h):
                self._eqs[name] += h
                # h(k+1) = d{eq}dw*dwdw^(k+1) = h(k)*dwdw
                h = smart_multiply(h, dwdw)

    def _total_derivative(
        self,
        name: str,
        eq: str,
        chainvars: list[str],
        var: str,
        dydx_name: str = None,
        dxdz_name: str = None,
    ) -> None:
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
            with respect to which the derivatives are to be computed

        :param dydx_name:
            defines the name of the symbolic variable that
            defines the derivative of the ``eq`` with respect to ``chainvar``,
            default is ``d{eq}d{chainvar}``

        :param dxdz_name:
            defines the name of the symbolic variable that
            defines the derivative of the ``chainvar`` with respect to ``var``,
            default is d{chainvar}d{var}
        """
        # compute total derivative according to chainrule
        # Dydz = dydx*dxdz + dydz

        # initialize with partial derivative dydz without chain rule
        self._eqs[name] = self.sym_or_eq(name, f"d{eq}d{var}")
        if not isinstance(self._eqs[name], sp.Symbol):
            # if not a Symbol, create a copy using sympy API
            # NB deepcopy does not work safely, see sympy issue #7672
            self._eqs[name] = self._eqs[name].copy()

        for chainvar in chainvars:
            if dydx_name is None:
                dydx_name = f"d{eq}d{chainvar}"
            if dxdz_name is None:
                dxdz_name = f"d{chainvar}d{var}"

            dydx = self.sym_or_eq(name, dydx_name)
            dxdz = self.sym_or_eq(name, dxdz_name)
            # Save time for large models if one multiplicand is zero,
            # which is not checked for by sympy
            if not smart_is_zero_matrix(dydx) and not smart_is_zero_matrix(
                dxdz
            ):
                dydx_times_dxdz = smart_multiply(dydx, dxdz)
                if (
                    dxdz.shape[1] == 1
                    and self._eqs[name].shape[1] != dxdz.shape[1]
                ):
                    for iz in range(self._eqs[name].shape[1]):
                        self._eqs[name][:, iz] += dydx_times_dxdz
                else:
                    self._eqs[name] += dydx_times_dxdz

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
        if var_in_function_signature(
            name, varname, self.is_ode()
        ) and varname not in [
            "dwdx",
            "dwdp",
        ]:
            return self.sym(varname)
        else:
            return self.eq(varname)

    def _multiplication(
        self,
        name: str,
        x: str,
        y: str,
        transpose_x: Optional[bool] = False,
        sign: Optional[int] = 1,
    ):
        """
        Creates a new symbolic variable according to a multiplication

        :param name:
            name of resulting symbolic variable, default is ``d{eq}d{var}``

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
            raise TypeError(f"sign must be +1 or -1, was {sign}")

        variables = {
            varname: self.sym(varname)
            if var_in_function_signature(name, varname, self.is_ode())
            else self.eq(varname)
            for varname in [x, y]
        }

        xx = variables[x].transpose() if transpose_x else variables[x]
        yy = variables[y]

        self._eqs[name] = sign * smart_multiply(xx, yy)

    def _equation_from_components(
        self, name: str, components: list[ModelQuantity]
    ) -> None:
        """
        Generates the formulas of a symbolic variable from the attributes

        :param name:
            name of resulting symbolic variable

        :param component:
            name of the attribute
        """
        self._eqs[name] = sp.Matrix([comp.get_val() for comp in components])

    def get_conservation_laws(self) -> list[tuple[sp.Symbol, sp.Expr]]:
        """Returns a list of states with conservation law set

        :return:
            list of state identifiers
        """
        return [
            (state.get_id(), state.get_x_rdata())
            for state in self.states()
            if state.has_conservation_law()
        ]

    def _generate_value(self, name: str) -> None:
        """
        Generates the numeric values of a symbolic variable from value
        prototypes

        :param name:
            name of resulting symbolic variable
        """
        if name in self._value_prototype:
            components = self._value_prototype[name]()
        else:
            raise ValueError(f"No values for {name}")

        self._vals[name] = [comp.get_val() for comp in components]

    def _generate_name(self, name: str) -> None:
        """
        Generates the names of a symbolic variable from variable prototypes or
        equation prototypes

        :param name:
            name of resulting symbolic variable
        """
        if name in self._variable_prototype:
            components = self._variable_prototype[name]()
        elif name in self._equation_prototype:
            components = self._equation_prototype[name]()
        else:
            raise ValueError(f"No names for {name}")

        self._names[name] = [comp.get_name() for comp in components]

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
        ic = self.states()[ix].get_val()
        if not isinstance(ic, sp.Basic):
            return False
        return any(
            fp in (c.get_id() for c in self._constants)
            for fp in ic.free_symbols
        )

    def state_has_conservation_law(self, ix: int) -> bool:
        """
        Checks whether the state at specified index has a conservation
        law set

        :param ix:
            state index

        :return:
            boolean indicating if conservation_law is not None
        """
        return self.states()[ix].has_conservation_law()

    def get_solver_indices(self) -> dict[int, int]:
        """
        Returns a mapping that maps rdata species indices to solver indices

        :return:
            dictionary mapping rdata species indices to solver indices
        """
        solver_index = {}
        ix_solver = 0
        for ix in range(len(self.states())):
            if self.state_has_conservation_law(ix):
                continue
            solver_index[ix] = ix_solver
            ix_solver += 1
        return solver_index

    def state_is_constant(self, ix: int) -> bool:
        """
        Checks whether the temporal derivative of the state is zero

        :param ix:
            state index

        :return:
            boolean indicating if constant over time
        """
        state = self.states()[ix]
        if isinstance(state, AlgebraicState):
            return False

        return state.get_dt() == 0.0

    def conservation_law_has_multispecies(self, tcl: ConservationLaw) -> bool:
        """
        Checks whether a conservation law has multiple species or it just
        defines one constant species

        :param tcl:
            conservation law

        :return:
            boolean indicating if conservation_law is not None
        """
        state_set = set(self.sym("x_rdata"))
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
        if "t" in expr_syms:
            return True

        # Check if any time-dependent states are in the expression.
        state_syms = [str(sym) for sym in self.states()]
        return any(
            not self.state_is_constant(state_syms.index(state))
            for state in expr_syms.intersection(state_syms)
        )

    def _get_unique_root(
        self,
        root_found: sp.Expr,
        roots: list[Event],
    ) -> Union[sp.Symbol, None]:
        """
        Collects roots of Heaviside functions and events and stores them in
        the roots list. It checks for redundancy to not store symbolically
        equivalent root functions more than once.

        :param root_found:
            equation of the root function
        :param roots:
            list of already known root functions with identifier

        :returns:
            unique identifier for root, or ``None`` if the root is not
            time-dependent
        """
        if not self._expr_is_time_dependent(root_found):
            return None

        for root in roots:
            if sp.simplify(root_found - root.get_val()) == 0:
                return root.get_id()

        # create an event for a new root function
        root_symstr = f"Heaviside_{len(roots)}"
        roots.append(
            Event(
                identifier=sp.Symbol(root_symstr),
                name=root_symstr,
                value=root_found,
                state_update=None,
            )
        )
        return roots[-1].get_id()

    def _collect_heaviside_roots(
        self,
        args: Sequence[sp.Expr],
    ) -> list[sp.Expr]:
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
        # rewriting 'root.cpp' and 'stau.cpp' headers
        # to include 'w.h'
        w_sorted = toposort_symbols(
            dict(
                zip(
                    [expr.get_id() for expr in self._expressions],
                    [expr.get_val() for expr in self._expressions],
                )
            )
        )
        root_funs = [r.subs(w_sorted) for r in root_funs]

        return root_funs

    def _process_heavisides(
        self,
        dxdt: sp.Expr,
        roots: list[Event],
    ) -> sp.Expr:
        """
        Parses the RHS of a state variable, checks for Heaviside functions,
        collects unique roots functions that can be tracked by SUNDIALS and
        replaces Heaviside Functions by amici helper variables that will be
        updated based on SUNDIALS root tracking.

        :param dxdt:
            right-hand side of state variable
        :param roots:
            list of known root functions with identifier

        :returns:
            dxdt with Heaviside functions replaced by amici helper variables
        """
        # track all the old Heaviside expressions in tmp_roots_old
        # replace them later by the new expressions
        heavisides = []
        # run through the expression tree and get the roots
        tmp_roots_old = self._collect_heaviside_roots(dxdt.args)
        for tmp_old in unique_preserve_order(tmp_roots_old):
            # we want unique identifiers for the roots
            tmp_new = self._get_unique_root(tmp_old, roots)
            # `tmp_new` is None if the root is not time-dependent.
            if tmp_new is None:
                continue
            # For Heavisides, we need to add the negative function as well
            self._get_unique_root(sp.sympify(-tmp_old), roots)
            heavisides.append((sp.Heaviside(tmp_old), tmp_new))

        if heavisides:
            # only apply subs if necessary
            for heaviside_sympy, heaviside_amici in heavisides:
                dxdt = dxdt.subs(heaviside_sympy, heaviside_amici)

        return dxdt


class DEExporter:
    """
    The DEExporter class generates AMICI C++ files for a model as
    defined in symbolic expressions.

    :ivar model:
        DE definition

    :ivar verbose:
        more verbose output if True

    :ivar assume_pow_positivity:
        if set to true, a special pow function is
        used to avoid problems with state variables that may become negative
        due to numerical errors

    :ivar compiler:
        Absolute path to the compiler executable to be used to build the Python
        extension, e.g. ``/usr/bin/clang``.

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

    .. note::
        When importing large models (several hundreds of species or
        parameters), import time can potentially be reduced by using multiple
        CPU cores. This is controlled by setting the ``AMICI_IMPORT_NPROCS``
        environment variable to the number of parallel processes that are to be
        used (default: 1). Note that for small models this may (slightly)
        increase import times.
    """

    def __init__(
        self,
        de_model: DEModel,
        outdir: Optional[Union[Path, str]] = None,
        verbose: Optional[Union[bool, int]] = False,
        assume_pow_positivity: Optional[bool] = False,
        compiler: Optional[str] = None,
        allow_reinit_fixpar_initcond: Optional[bool] = True,
        generate_sensitivity_code: Optional[bool] = True,
        model_name: Optional[str] = "model",
    ):
        """
        Generate AMICI C++ files for the DE provided to the constructor.

        :param de_model:
            DE model definition

        :param outdir:
            see :meth:`amici.de_export.DEExporter.set_paths`

        :param verbose:
            verbosity level for logging, ``True``/``False`` default to
            :data:`logging.Error`/:data:`logging.DEBUG`

        :param assume_pow_positivity:
            if set to true, a special pow function is
            used to avoid problems with state variables that may become
            negative due to numerical errors

        :param compiler: Absolute path to the compiler executable to be used
            to build the Python extension, e.g. ``/usr/bin/clang``.

        :param allow_reinit_fixpar_initcond:
            see :class:`amici.de_export.DEExporter`

        :param generate_sensitivity_code:
            specifies whether code required for sensitivity computation will be
            generated

        :param model_name:
            name of the model to be used during code generation
        """
        set_log_level(logger, verbose)

        self.verbose: bool = logger.getEffectiveLevel() <= logging.DEBUG
        self.assume_pow_positivity: bool = assume_pow_positivity
        self.compiler: str = compiler

        self.model_path: str = ""
        self.model_swig_path: str = ""

        self.set_name(model_name)
        self.set_paths(outdir)

        # Signatures and properties of generated model functions (see
        # include/amici/model.h for details)
        self.model: DEModel = de_model
        self.model._code_printer.known_functions.update(
            splines.spline_user_functions(
                self.model.splines, self._get_index("p")
            )
        )

        # To only generate a subset of functions, apply subselection here
        self.functions: dict[str, _FunctionInfo] = copy.deepcopy(functions)

        self.allow_reinit_fixpar_initcond: bool = allow_reinit_fixpar_initcond
        self._build_hints = set()
        self.generate_sensitivity_code: bool = generate_sensitivity_code

    @log_execution_time("generating cpp code", logger)
    def generate_model_code(self) -> None:
        """
        Generates the native C++ code for the loaded model and a Matlab
        script that can be run to compile a mex file from the C++ code
        """
        with _monkeypatched(
            sp.Pow, "_eval_derivative", _custom_pow_eval_derivative
        ):
            self._prepare_model_folder()
            self._generate_c_code()
            self._generate_m_code()

    @log_execution_time("compiling cpp code", logger)
    def compile_model(self) -> None:
        """
        Compiles the generated code it into a simulatable module
        """
        self._compile_c_code(compiler=self.compiler, verbose=self.verbose)

    def _prepare_model_folder(self) -> None:
        """
        Create model directory or remove all files if the output directory
        already exists.
        """
        os.makedirs(self.model_path, exist_ok=True)

        for file in os.listdir(self.model_path):
            file_path = os.path.join(self.model_path, file)
            if os.path.isfile(file_path):
                os.remove(file_path)

    def _generate_c_code(self) -> None:
        """
        Create C++ code files for the model based on
        :attribute:`DEExporter.model`.
        """
        for func_name, func_info in self.functions.items():
            if (
                func_name in sensi_functions + sparse_sensi_functions
                and not self.generate_sensitivity_code
            ):
                continue

            if func_info.generate_body:
                dec = log_execution_time(f"writing {func_name}.cpp", logger)
                dec(self._write_function_file)(func_name)

        for name in self.model.sym_names():
            # only generate for those that have nontrivial implementation,
            # check for both basic variables (not in functions) and function
            # computed values
            if (
                (
                    name in self.functions
                    and not self.functions[name].body
                    and name not in nobody_functions
                )
                or name not in self.functions
            ) and len(self.model.sym(name)) == 0:
                continue
            self._write_index_files(name)

        self._write_wrapfunctions_cpp()
        self._write_wrapfunctions_header()
        self._write_model_header_cpp()
        self._write_c_make_file()
        self._write_swig_files()
        self._write_module_setup()

        shutil.copy(
            CXX_MAIN_TEMPLATE_FILE, os.path.join(self.model_path, "main.cpp")
        )

    def _compile_c_code(
        self,
        verbose: Optional[Union[bool, int]] = False,
        compiler: Optional[str] = None,
    ) -> None:
        """
        Compile the generated model code

        :param verbose:
            Make model compilation verbose

        :param compiler:
            Absolute path to the compiler executable to be used to build the Python
            extension, e.g. ``/usr/bin/clang``.
        """
        # setup.py assumes it is run from within the model directory
        module_dir = self.model_path
        script_args = [sys.executable, os.path.join(module_dir, "setup.py")]

        if verbose:
            script_args.append("--verbose")
        else:
            script_args.append("--quiet")

        script_args.extend(
            [
                "build_ext",
                f"--build-lib={module_dir}",
                # This is generally not required, but helps to reduce the path
                # length of intermediate build files, that may easily become
                # problematic on Windows, due to its ridiculous 255-character path
                # length limit.
                f'--build-temp={Path(module_dir, "build")}',
            ]
        )

        env = os.environ.copy()
        if compiler is not None:
            # CMake will use the compiler specified in the CXX environment variable
            env["CXX"] = compiler

        # distutils.core.run_setup looks nicer, but does not let us check the
        # result easily
        try:
            result = subprocess.run(
                script_args,
                cwd=module_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=True,
                env=env,
            )
        except subprocess.CalledProcessError as e:
            print(e.output.decode("utf-8"))
            print("Failed building the model extension.")
            if self._build_hints:
                print("Note:")
                print("\n".join(self._build_hints))
            raise

        if verbose:
            print(result.stdout.decode("utf-8"))

    def _generate_m_code(self) -> None:
        """
        Create a Matlab script for compiling code files to a mex file
        """

        # Second order code is not yet implemented. Once this is done,
        # those variables will have to be replaced by
        # "self.model.<var>true()", or the corresponding "model.self.o2flag"
        nxtrue_rdata = self.model.num_states_rdata()
        nytrue = self.model.num_obs()
        nztrue = self.model.num_eventobs()
        o2flag = 0

        lines = [
            "% This compile script was automatically created from"
            " Python SBML import.",
            "% If mex compiler is set up within MATLAB, it can be run"
            " from MATLAB ",
            "% in order to compile a mex-file from the Python"
            " generated C++ files.",
            "",
            f"modelName = '{self.model_name}';",
            "amimodel.compileAndLinkModel(modelName, '', [], [], [], []);",
            f"amimodel.generateMatlabWrapper({nxtrue_rdata}, "
            f"{nytrue}, {self.model.num_par()}, "
            f"{self.model.num_const()}, {nztrue}, {o2flag}, ...",
            "    [], ['simulate_' modelName '.m'], modelName, ...",
            "    'lin', 1, 1);",
        ]

        # write compile script (for mex)
        compile_script = os.path.join(self.model_path, "compileMexFile.m")
        with open(compile_script, "w") as fileout:
            fileout.write("\n".join(lines))

    def _get_index(self, name: str) -> dict[sp.Symbol, int]:
        """
        Compute indices for a symbolic array.
        :param name:
            key in self.model._syms for which to obtain the index.
        :return:
            a dictionary of symbol/index pairs.
        """
        if name in self.model.sym_names():
            if name in sparse_functions:
                symbols = self.model.sparsesym(name)
            else:
                symbols = self.model.sym(name).T
        else:
            raise ValueError(f"Unknown symbolic array: {name}")

        return {
            strip_pysb(symbol).name: index
            for index, symbol in enumerate(symbols)
        }

    def _write_index_files(self, name: str) -> None:
        """
        Write index file for a symbolic array.

        :param name:
            key in ``self.model._syms`` for which the respective file should
            be written
        """
        if name not in self.model.sym_names():
            raise ValueError(f"Unknown symbolic array: {name}")

        symbols = (
            self.model.sparsesym(name)
            if name in sparse_functions
            else self.model.sym(name).T
        )
        if not len(symbols):
            return

        # flatten multiobs
        if isinstance(next(iter(symbols), None), list):
            symbols = [symbol for obs in symbols for symbol in obs]

        lines = []
        for index, symbol in enumerate(symbols):
            symbol_name = strip_pysb(symbol)
            if str(symbol) == "0":
                continue
            if str(symbol_name) == "":
                raise ValueError(f'{name} contains a symbol called ""')
            lines.append(f"#define {symbol_name} {name}[{index}]")
            if name == "stau":
                # we only need a single macro, as all entries have the same symbol
                break

        filename = os.path.join(self.model_path, f"{name}.h")
        with open(filename, "w") as fileout:
            fileout.write("\n".join(lines))

    def _write_function_file(self, function: str) -> None:
        """
        Generate equations and write the C++ code for the function
        ``function``.

        :param function:
            name of the function to be written (see ``self.functions``)
        """

        # first generate the equations to make sure we have everything we
        # need in subsequent steps
        if function in sparse_functions:
            equations = self.model.sparseeq(function)
        elif (
            not self.allow_reinit_fixpar_initcond
            and function == "sx0_fixedParameters"
        ):
            # Not required. Will create empty function body.
            equations = sp.Matrix()
        elif function == "create_splines":
            # nothing to do
            pass
        else:
            equations = self.model.eq(function)

        # function body
        if function == "create_splines":
            body = self._get_create_splines_body()
        else:
            body = self._get_function_body(function, equations)
        if not body:
            return

        # colptrs / rowvals for sparse matrices
        if function in sparse_functions:
            lines = self._generate_function_index(function, "colptrs")
            lines.extend(self._generate_function_index(function, "rowvals"))
            lines.append("\n\n")
        else:
            lines = []

        # function header
        lines.extend(
            [
                '#include "amici/symbolic_functions.h"',
                '#include "amici/defines.h"',
                '#include "sundials/sundials_types.h"',
                "",
                "#include <gsl/gsl-lite.hpp>",
                "#include <algorithm>",
                "",
            ]
        )
        if function == "create_splines":
            lines += [
                '#include "amici/splinefunctions.h"',
                "#include <vector>",
            ]

        func_info = self.functions[function]

        # extract symbols that need definitions from signature
        # don't add includes for files that won't be generated.
        # Unfortunately we cannot check for `self.functions[sym].body`
        # here since it may not have been generated yet.
        for sym in re.findall(
            r"const (?:realtype|double) \*([\w]+)[0]*(?:,|$)",
            func_info.arguments(self.model.is_ode()),
        ):
            if sym not in self.model.sym_names():
                continue

            if sym in sparse_functions:
                iszero = smart_is_zero_matrix(self.model.sparseeq(sym))
            elif sym in self.functions:
                iszero = smart_is_zero_matrix(self.model.eq(sym))
            else:
                iszero = len(self.model.sym(sym)) == 0

            if iszero and not (
                (sym == "y" and "Jy" in function)
                or (
                    sym == "w"
                    and "xdot" in function
                    and len(self.model.sym(sym))
                )
            ):
                continue

            lines.append(f'#include "{sym}.h"')

        # include return symbols
        if (
            function in self.model.sym_names()
            and function not in non_unique_id_symbols
        ):
            lines.append(f'#include "{function}.h"')

        lines.extend(
            [
                "",
                "namespace amici {",
                f"namespace model_{self.model_name} {{",
                "",
                f"{func_info.return_type} {function}_{self.model_name}"
                f"({func_info.arguments(self.model.is_ode())}){{",
            ]
        )

        if self.assume_pow_positivity and func_info.assume_pow_positivity:
            pow_rx = re.compile(r"(^|\W)std::pow\(")
            body = [
                # execute this twice to catch cases where the ending '(' would
                #  be the starting (^|\W) for the following match
                pow_rx.sub(
                    r"\1amici::pos_pow(",
                    pow_rx.sub(r"\1amici::pos_pow(", line),
                )
                for line in body
            ]

        self.functions[function].body = body

        lines += body
        lines.extend(
            [
                "}",
                "",
                f"}} // namespace model_{self.model_name}",
                "} // namespace amici\n",
            ]
        )

        # check custom functions
        for fun in CUSTOM_FUNCTIONS:
            if "include" in fun and any(fun["c++"] in line for line in lines):
                if "build_hint" in fun:
                    self._build_hints.add(fun["build_hint"])
                lines.insert(0, fun["include"])

        # if not body is None:
        filename = os.path.join(self.model_path, f"{function}.cpp")
        with open(filename, "w") as fileout:
            fileout.write("\n".join(lines))

    def _generate_function_index(
        self, function: str, indextype: Literal["colptrs", "rowvals"]
    ) -> list[str]:
        """
        Generate equations and C++ code for the function ``function``.

        :param function:
            name of the function to be written (see ``self.functions``)

        :param indextype:
            type of index {'colptrs', 'rowvals'}

        :returns:
            The code lines for the respective function index file
        """
        if indextype == "colptrs":
            values = self.model.colptrs(function)
            setter = "indexptrs"
        elif indextype == "rowvals":
            values = self.model.rowvals(function)
            setter = "indexvals"
        else:
            raise ValueError(
                "Invalid value for indextype, must be colptrs or "
                f"rowvals: {indextype}"
            )

        # function signature
        if function in multiobs_functions:
            signature = f"(SUNMatrixWrapper &{function}, int index)"
        else:
            signature = f"(SUNMatrixWrapper &{function})"

        lines = [
            '#include "amici/sundials_matrix_wrapper.h"',
            '#include "sundials/sundials_types.h"',
            "",
            "#include <array>",
            "#include <algorithm>",
            "",
            "namespace amici {",
            f"namespace model_{self.model_name} {{",
            "",
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
                lines.extend(
                    [
                        "    {" + ", ".join(map(str, index_vector)) + "}, "
                        for index_vector in values
                    ]
                )
                lines.append("}};")
            else:
                # single index vector
                lines.extend(
                    [
                        "static constexpr std::array<sunindextype, "
                        f"{len(values)}> {static_array_name} = {{",
                        "    " + ", ".join(map(str, values)),
                        "};",
                    ]
                )

        lines.extend(
            [
                "",
                f"void {function}_{indextype}_{self.model_name}{signature}{{",
            ]
        )

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

        lines.extend(
            [
                "}" "",
                f"}} // namespace model_{self.model_name}",
                "} // namespace amici\n",
            ]
        )

        return lines

    def _get_function_body(
        self, function: str, equations: sp.Matrix
    ) -> list[str]:
        """
        Generate C++ code for body of function ``function``.

        :param function:
            name of the function to be written (see ``self.functions``)

        :param equations:
            symbolic definition of the function body

        :return:
            generated C++ code
        """
        lines = []

        if len(equations) == 0 or (
            isinstance(equations, (sp.Matrix, sp.ImmutableDenseMatrix))
            and min(equations.shape) == 0
        ):
            # dJydy is a list
            return lines

        if not self.allow_reinit_fixpar_initcond and function in {
            "sx0_fixedParameters",
            "x0_fixedParameters",
        }:
            return lines

        if function == "sx0_fixedParameters":
            # here we only want to overwrite values where x0_fixedParameters
            # was applied

            lines.extend(
                [
                    # Keep list of indices of fixed parameters occurring in x0
                    "    static const std::array<int, "
                    + str(len(self.model._x0_fixedParameters_idx))
                    + "> _x0_fixedParameters_idxs = {",
                    "        "
                    + ", ".join(
                        str(x) for x in self.model._x0_fixedParameters_idx
                    ),
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
                    "    }",
                ]
            )

            cases = {}
            for ipar in range(self.model.num_par()):
                expressions = []
                for index, formula in zip(
                    self.model._x0_fixedParameters_idx, equations[:, ipar]
                ):
                    if not formula.is_zero:
                        expressions.extend(
                            [
                                f"if(std::find("
                                "reinitialization_state_idxs.cbegin(), "
                                f"reinitialization_state_idxs.cend(), {index}) != "
                                "reinitialization_state_idxs.cend())",
                                f"    {function}[{index}] = "
                                f"{self.model._code_printer.doprint(formula)};",
                            ]
                        )
                cases[ipar] = expressions
            lines.extend(get_switch_statement("ip", cases, 1))

        elif function == "x0_fixedParameters":
            for index, formula in zip(
                self.model._x0_fixedParameters_idx, equations
            ):
                lines.append(
                    f"    if(std::find(reinitialization_state_idxs.cbegin(), "
                    f"reinitialization_state_idxs.cend(), {index}) != "
                    "reinitialization_state_idxs.cend())\n        "
                    f"{function}[{index}] = "
                    f"{self.model._code_printer.doprint(formula)};"
                )

        elif function in event_functions:
            cases = {
                ie: self.model._code_printer._get_sym_lines_array(
                    equations[ie], function, 0
                )
                for ie in range(self.model.num_events())
                if not smart_is_zero_matrix(equations[ie])
            }
            lines.extend(get_switch_statement("ie", cases, 1))

        elif function in event_sensi_functions:
            outer_cases = {}
            for ie, inner_equations in enumerate(equations):
                inner_lines = []
                inner_cases = {
                    ipar: self.model._code_printer._get_sym_lines_array(
                        inner_equations[:, ipar], function, 0
                    )
                    for ipar in range(self.model.num_par())
                    if not smart_is_zero_matrix(inner_equations[:, ipar])
                }
                inner_lines.extend(get_switch_statement("ip", inner_cases, 0))
                outer_cases[ie] = copy.copy(inner_lines)
            lines.extend(get_switch_statement("ie", outer_cases, 1))

        elif (
            function in sensi_functions
            and equations.shape[1] == self.model.num_par()
        ):
            cases = {
                ipar: self.model._code_printer._get_sym_lines_array(
                    equations[:, ipar], function, 0
                )
                for ipar in range(self.model.num_par())
                if not smart_is_zero_matrix(equations[:, ipar])
            }
            lines.extend(get_switch_statement("ip", cases, 1))
        elif function in multiobs_functions:
            if function == "dJydy":
                cases = {
                    iobs: self.model._code_printer._get_sym_lines_array(
                        equations[iobs], function, 0
                    )
                    for iobs in range(self.model.num_obs())
                    if not smart_is_zero_matrix(equations[iobs])
                }
            else:
                cases = {
                    iobs: self.model._code_printer._get_sym_lines_array(
                        equations[:, iobs], function, 0
                    )
                    for iobs in range(equations.shape[1])
                    if not smart_is_zero_matrix(equations[:, iobs])
                }
            if function.startswith(("Jz", "dJz", "Jrz", "dJrz")):
                iterator = "iz"
            else:
                iterator = "iy"
            lines.extend(get_switch_statement(iterator, cases, 1))

        elif (
            function in self.model.sym_names()
            and function not in non_unique_id_symbols
        ):
            if function in sparse_functions:
                symbols = list(map(sp.Symbol, self.model.sparsesym(function)))
            else:
                symbols = self.model.sym(function)
            lines += self.model._code_printer._get_sym_lines_symbols(
                symbols, equations, function, 4
            )

        else:
            lines += self.model._code_printer._get_sym_lines_array(
                equations, function, 4
            )

        return [line for line in lines if line]

    def _get_create_splines_body(self):
        if not self.model.splines:
            return ["    return {};"]

        ind4 = " " * 4
        ind8 = " " * 8

        body = ["return {"]
        for ispl, spline in enumerate(self.model.splines):
            if isinstance(spline.nodes, splines.UniformGrid):
                nodes = (
                    f"{ind8}{{{spline.nodes.start}, {spline.nodes.stop}}}, "
                )
            else:
                nodes = f"{ind8}{{{', '.join(map(str, spline.nodes))}}}, "

            # vector with the node values
            values = (
                f"{ind8}{{{', '.join(map(str, spline.values_at_nodes))}}}, "
            )
            # vector with the slopes
            if spline.derivatives_by_fd:
                slopes = f"{ind8}{{}},"
            else:
                slopes = f"{ind8}{{{', '.join(map(str, spline.derivatives_at_nodes))}}},"

            body.extend(
                [
                    f"{ind4}HermiteSpline(",
                    nodes,
                    values,
                    slopes,
                ]
            )

            bc_to_cpp = {
                None: "SplineBoundaryCondition::given, ",
                "zeroderivative": "SplineBoundaryCondition::zeroDerivative, ",
                "natural": "SplineBoundaryCondition::natural, ",
                "zeroderivative+natural": "SplineBoundaryCondition::naturalZeroDerivative, ",
                "periodic": "SplineBoundaryCondition::periodic, ",
            }
            for bc in spline.bc:
                try:
                    body.append(ind8 + bc_to_cpp[bc])
                except KeyError:
                    raise ValueError(
                        f"Unknown boundary condition '{bc}' "
                        "found in spline object"
                    )
            extrapolate_to_cpp = {
                None: "SplineExtrapolation::noExtrapolation, ",
                "polynomial": "SplineExtrapolation::polynomial, ",
                "constant": "SplineExtrapolation::constant, ",
                "linear": "SplineExtrapolation::linear, ",
                "periodic": "SplineExtrapolation::periodic, ",
            }
            for extr in spline.extrapolate:
                try:
                    body.append(ind8 + extrapolate_to_cpp[extr])
                except KeyError:
                    raise ValueError(
                        f"Unknown extrapolation '{extr}' "
                        "found in spline object"
                    )
            line = ind8
            line += "true, " if spline.derivatives_by_fd else "false, "
            line += (
                "true, "
                if isinstance(spline.nodes, splines.UniformGrid)
                else "false, "
            )
            line += "true" if spline.logarithmic_parametrization else "false"
            body.append(line)
            body.append(f"{ind4}),")

        body.append("};")
        return ["    " + line for line in body]

    def _write_wrapfunctions_cpp(self) -> None:
        """
        Write model-specific 'wrapper' file (``wrapfunctions.cpp``).
        """
        template_data = {"MODELNAME": self.model_name}
        apply_template(
            os.path.join(amiciSrcPath, "wrapfunctions.template.cpp"),
            os.path.join(self.model_path, "wrapfunctions.cpp"),
            template_data,
        )

    def _write_wrapfunctions_header(self) -> None:
        """
        Write model-specific header file (``wrapfunctions.h``).
        """
        template_data = {"MODELNAME": str(self.model_name)}
        apply_template(
            os.path.join(amiciSrcPath, "wrapfunctions.template.h"),
            os.path.join(self.model_path, "wrapfunctions.h"),
            template_data,
        )

    def _write_model_header_cpp(self) -> None:
        """
        Write model-specific header and cpp file (MODELNAME.{h,cpp}).
        """
        model_type = "ODE" if self.model.is_ode() else "DAE"
        tpl_data = {
            "MODEL_TYPE_LOWER": model_type.lower(),
            "MODEL_TYPE_UPPER": model_type,
            "MODELNAME": self.model_name,
            "NX_RDATA": self.model.num_states_rdata(),
            "NXTRUE_RDATA": self.model.num_states_rdata(),
            "NX_SOLVER": self.model.num_states_solver(),
            "NXTRUE_SOLVER": self.model.num_states_solver(),
            "NX_SOLVER_REINIT": self.model.num_state_reinits(),
            "NY": self.model.num_obs(),
            "NYTRUE": self.model.num_obs(),
            "NZ": self.model.num_eventobs(),
            "NZTRUE": self.model.num_eventobs(),
            "NEVENT": self.model.num_events(),
            "NEVENT_SOLVER": self.model.num_events_solver(),
            "NOBJECTIVE": "1",
            "NSPL": len(self.model.splines),
            "NW": len(self.model.sym("w")),
            "NDWDP": len(
                self.model.sparsesym(
                    "dwdp", force_generate=self.generate_sensitivity_code
                )
            ),
            "NDWDX": len(self.model.sparsesym("dwdx")),
            "NDWDW": len(self.model.sparsesym("dwdw")),
            "NDXDOTDW": len(self.model.sparsesym("dxdotdw")),
            "NDXDOTDP_EXPLICIT": len(
                self.model.sparsesym(
                    "dxdotdp_explicit",
                    force_generate=self.generate_sensitivity_code,
                )
            ),
            "NDXDOTDX_EXPLICIT": len(self.model.sparsesym("dxdotdx_explicit")),
            "NDJYDY": "std::vector<int>{%s}"
            % ",".join(str(len(x)) for x in self.model.sparsesym("dJydy")),
            "NDXRDATADXSOLVER": len(self.model.sparsesym("dx_rdatadx_solver")),
            "NDXRDATADTCL": len(self.model.sparsesym("dx_rdatadtcl")),
            "NDTOTALCLDXRDATA": len(self.model.sparsesym("dtotal_cldx_rdata")),
            "UBW": self.model.num_states_solver(),
            "LBW": self.model.num_states_solver(),
            "NP": self.model.num_par(),
            "NK": self.model.num_const(),
            "O2MODE": "amici::SecondOrderMode::none",
            # using code printer ensures proper handling of nan/inf
            "PARAMETERS": self.model._code_printer.doprint(
                self.model.val("p")
            )[1:-1],
            "FIXED_PARAMETERS": self.model._code_printer.doprint(
                self.model.val("k")
            )[1:-1],
            "PARAMETER_NAMES_INITIALIZER_LIST": self._get_symbol_name_initializer_list(
                "p"
            ),
            "STATE_NAMES_INITIALIZER_LIST": self._get_symbol_name_initializer_list(
                "x_rdata"
            ),
            "FIXED_PARAMETER_NAMES_INITIALIZER_LIST": self._get_symbol_name_initializer_list(
                "k"
            ),
            "OBSERVABLE_NAMES_INITIALIZER_LIST": self._get_symbol_name_initializer_list(
                "y"
            ),
            "OBSERVABLE_TRAFO_INITIALIZER_LIST": "\n".join(
                f"ObservableScaling::{trafo.value}, // y[{idx}]"
                for idx, trafo in enumerate(
                    self.model.get_observable_transformations()
                )
            ),
            "EXPRESSION_NAMES_INITIALIZER_LIST": self._get_symbol_name_initializer_list(
                "w"
            ),
            "PARAMETER_IDS_INITIALIZER_LIST": self._get_symbol_id_initializer_list(
                "p"
            ),
            "STATE_IDS_INITIALIZER_LIST": self._get_symbol_id_initializer_list(
                "x_rdata"
            ),
            "FIXED_PARAMETER_IDS_INITIALIZER_LIST": self._get_symbol_id_initializer_list(
                "k"
            ),
            "OBSERVABLE_IDS_INITIALIZER_LIST": self._get_symbol_id_initializer_list(
                "y"
            ),
            "EXPRESSION_IDS_INITIALIZER_LIST": self._get_symbol_id_initializer_list(
                "w"
            ),
            "STATE_IDXS_SOLVER_INITIALIZER_LIST": ", ".join(
                str(idx)
                for idx, state in enumerate(self.model.states())
                if not state.has_conservation_law()
            ),
            "REINIT_FIXPAR_INITCOND": AmiciCxxCodePrinter.print_bool(
                self.allow_reinit_fixpar_initcond
            ),
            "AMICI_VERSION_STRING": __version__,
            "AMICI_COMMIT_STRING": __commit__,
            "W_RECURSION_DEPTH": self.model._w_recursion_depth,
            "QUADRATIC_LLH": AmiciCxxCodePrinter.print_bool(
                self.model._has_quadratic_nllh
            ),
            "ROOT_INITIAL_VALUES": ", ".join(
                map(
                    lambda event: AmiciCxxCodePrinter.print_bool(
                        event.get_initial_value()
                    ),
                    self.model.events(),
                )
            ),
            "Z2EVENT": ", ".join(map(str, self.model._z2event)),
            "STATE_INDEPENDENT_EVENTS": self._get_state_independent_event_intializer(),
            "ID": ", ".join(
                str(float(isinstance(s, DifferentialState)))
                for s in self.model.states()
                if not s.has_conservation_law()
            ),
        }

        for func_name, func_info in self.functions.items():
            if func_name in nobody_functions:
                continue

            if not func_info.body:
                tpl_data[f"{func_name.upper()}_DEF"] = ""

                if (
                    func_name in sensi_functions + sparse_sensi_functions
                    and not self.generate_sensitivity_code
                ):
                    impl = ""
                else:
                    impl = get_model_override_implementation(
                        func_name,
                        self.model_name,
                        self.model.is_ode(),
                        nobody=True,
                    )

                tpl_data[f"{func_name.upper()}_IMPL"] = impl

                if func_name in sparse_functions:
                    for indexfield in ["colptrs", "rowvals"]:
                        if (
                            func_name in sparse_sensi_functions
                            and not self.generate_sensitivity_code
                        ):
                            impl = ""
                        else:
                            impl = get_sunindex_override_implementation(
                                func_name,
                                self.model_name,
                                indexfield,
                                nobody=True,
                            )
                        tpl_data[
                            f"{func_name.upper()}_{indexfield.upper()}_DEF"
                        ] = ""
                        tpl_data[
                            f"{func_name.upper()}_{indexfield.upper()}_IMPL"
                        ] = impl
                continue

            tpl_data[
                f"{func_name.upper()}_DEF"
            ] = get_function_extern_declaration(
                func_name, self.model_name, self.model.is_ode()
            )
            tpl_data[
                f"{func_name.upper()}_IMPL"
            ] = get_model_override_implementation(
                func_name, self.model_name, self.model.is_ode()
            )
            if func_name in sparse_functions:
                tpl_data[
                    f"{func_name.upper()}_COLPTRS_DEF"
                ] = get_sunindex_extern_declaration(
                    func_name, self.model_name, "colptrs"
                )
                tpl_data[
                    f"{func_name.upper()}_COLPTRS_IMPL"
                ] = get_sunindex_override_implementation(
                    func_name, self.model_name, "colptrs"
                )
                tpl_data[
                    f"{func_name.upper()}_ROWVALS_DEF"
                ] = get_sunindex_extern_declaration(
                    func_name, self.model_name, "rowvals"
                )
                tpl_data[
                    f"{func_name.upper()}_ROWVALS_IMPL"
                ] = get_sunindex_override_implementation(
                    func_name, self.model_name, "rowvals"
                )

        if self.model.num_states_solver() == self.model.num_states_rdata():
            tpl_data["X_RDATA_DEF"] = ""
            tpl_data["X_RDATA_IMPL"] = ""

        tpl_data = {k: str(v) for k, v in tpl_data.items()}

        apply_template(
            os.path.join(amiciSrcPath, "model_header.template.h"),
            os.path.join(self.model_path, f"{self.model_name}.h"),
            tpl_data,
        )

        apply_template(
            os.path.join(amiciSrcPath, "model.template.cpp"),
            os.path.join(self.model_path, f"{self.model_name}.cpp"),
            tpl_data,
        )

    def _get_symbol_name_initializer_list(self, name: str) -> str:
        """
        Get SBML name initializer list for vector of names for the given
        model entity

        :param name:
            any key present in ``self.model._syms``

        :return:
            Template initializer list of names
        """
        return "\n".join(
            f'"{symbol}", // {name}[{idx}]'
            for idx, symbol in enumerate(self.model.name(name))
        )

    def _get_symbol_id_initializer_list(self, name: str) -> str:
        """
        Get C++ initializer list for vector of names for the given model
        entity

        :param name:
            any key present in ``self.model._syms``

        :return:
            Template initializer list of ids
        """
        return "\n".join(
            f'"{self.model._code_printer.doprint(symbol)}", // {name}[{idx}]'
            for idx, symbol in enumerate(self.model.sym(name))
        )

    def _get_state_independent_event_intializer(self) -> str:
        """Get initializer list for state independent events in amici::Model."""
        map_time_to_event_idx = {}
        for event_idx, event in enumerate(self.model.events()):
            if not event.triggers_at_fixed_timepoint():
                continue
            trigger_time = float(event.get_trigger_time())
            try:
                map_time_to_event_idx[trigger_time].append(event_idx)
            except KeyError:
                map_time_to_event_idx[trigger_time] = [event_idx]

        def vector_initializer(v):
            """std::vector initializer list with elements from `v`"""
            return f"{{{', '.join(map(str, v))}}}"

        return ", ".join(
            f"{{{trigger_time}, {vector_initializer(event_idxs)}}}"
            for trigger_time, event_idxs in map_time_to_event_idx.items()
        )

    def _write_c_make_file(self):
        """Write CMake ``CMakeLists.txt`` file for this model."""
        sources = "\n".join(
            f + " "
            for f in os.listdir(self.model_path)
            if f.endswith(
                (".cpp", ".h"),
            )
            and f != "main.cpp"
        )

        template_data = {
            "MODELNAME": self.model_name,
            "SOURCES": sources,
            "AMICI_VERSION": __version__,
        }
        apply_template(
            MODEL_CMAKE_TEMPLATE_FILE,
            Path(self.model_path, "CMakeLists.txt"),
            template_data,
        )

    def _write_swig_files(self) -> None:
        """Write SWIG interface files for this model."""
        Path(self.model_swig_path).mkdir(exist_ok=True)
        template_data = {"MODELNAME": self.model_name}
        apply_template(
            Path(amiciSwigPath, "modelname.template.i"),
            Path(self.model_swig_path, self.model_name + ".i"),
            template_data,
        )
        shutil.copy(
            SWIG_CMAKE_TEMPLATE_FILE,
            Path(self.model_swig_path, "CMakeLists.txt"),
        )

    def _write_module_setup(self) -> None:
        """
        Create a setuptools ``setup.py`` file for compile the model module.
        """

        template_data = {
            "MODELNAME": self.model_name,
            "AMICI_VERSION": __version__,
            "PACKAGE_VERSION": "0.1.0",
        }
        apply_template(
            Path(amiciModulePath, "setup.template.py"),
            Path(self.model_path, "setup.py"),
            template_data,
        )
        apply_template(
            Path(amiciModulePath, "MANIFEST.template.in"),
            Path(self.model_path, "MANIFEST.in"),
            {},
        )
        # write __init__.py for the model module
        Path(self.model_path, self.model_name).mkdir(exist_ok=True)

        apply_template(
            Path(amiciModulePath, "__init__.template.py"),
            Path(self.model_path, self.model_name, "__init__.py"),
            template_data,
        )

    def set_paths(self, output_dir: Optional[Union[str, Path]] = None) -> None:
        """
        Set output paths for the model and create if necessary

        :param output_dir:
            relative or absolute path where the generated model
            code is to be placed. If ``None``, this will default to
            ``amici-{self.model_name}`` in the current working directory.
            will be created if it does not exist.

        """
        if output_dir is None:
            output_dir = os.path.join(os.getcwd(), f"amici-{self.model_name}")

        self.model_path = os.path.abspath(output_dir)
        self.model_swig_path = os.path.join(self.model_path, "swig")

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
                "digits and underscores, and must not start with a digit."
            )

        self.model_name = model_name


class TemplateAmici(Template):
    """
    Template format used in AMICI (see :class:`string.Template` for more
    details).

    :cvar delimiter:
        delimiter that identifies template variables
    """

    delimiter = "TPL_"


def apply_template(
    source_file: Union[str, Path],
    target_file: Union[str, Path],
    template_data: dict[str, str],
) -> None:
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
    with open(target_file, "w") as fileout:
        fileout.write(result)


def get_function_extern_declaration(fun: str, name: str, ode: bool) -> str:
    """
    Constructs the extern function declaration for a given function

    :param fun:
        function name
    :param name:
        model name
    :param ode:
        whether to generate declaration for DAE or ODE

    :return:
        C++ function definition string
    """
    f = functions[fun]
    return f"extern {f.return_type} {fun}_{name}({f.arguments(ode)});"


def get_sunindex_extern_declaration(
    fun: str, name: str, indextype: str
) -> str:
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
        C++ function declaration string
    """
    index_arg = ", int index" if fun in multiobs_functions else ""
    return (
        f"extern void {fun}_{indextype}_{name}"
        f"(SUNMatrixWrapper &{indextype}{index_arg});"
    )


def get_model_override_implementation(
    fun: str, name: str, ode: bool, nobody: bool = False
) -> str:
    """
    Constructs ``amici::Model::*`` override implementation for a given function

    :param fun:
        function name

    :param name:
        model name

    :param nobody:
        whether the function has a nontrivial implementation

    :return:
        C++ function implementation string
    """
    func_info = functions[fun]
    body = (
        ""
        if nobody
        else "\n{ind8}{maybe_return}{fun}_{name}({eval_signature});{ind4}\n".format(
            ind4=" " * 4,
            ind8=" " * 8,
            maybe_return="" if func_info.return_type == "void" else "return ",
            fun=fun,
            name=name,
            eval_signature=remove_argument_types(func_info.arguments(ode)),
        )
    )
    return "{return_type} f{fun}({signature}) override {{{body}}}\n".format(
        return_type=func_info.return_type,
        fun=fun,
        signature=func_info.arguments(ode),
        body=body,
    )


def get_sunindex_override_implementation(
    fun: str, name: str, indextype: str, nobody: bool = False
) -> str:
    """
    Constructs the ``amici::Model`` function implementation for an index
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
        C++ function implementation string
    """
    index_arg = ", int index" if fun in multiobs_functions else ""
    index_arg_eval = ", index" if fun in multiobs_functions else ""

    impl = "void f{fun}_{indextype}({signature}) override {{"

    if nobody:
        impl += "}}\n"
    else:
        impl += "{ind8}{fun}_{indextype}_{name}({eval_signature});\n{ind4}}}\n"

    return impl.format(
        ind4=" " * 4,
        ind8=" " * 8,
        fun=fun,
        indextype=indextype,
        name=name,
        signature=f"SUNMatrixWrapper &{indextype}{index_arg}",
        eval_signature=f"{indextype}{index_arg_eval}",
    )


def remove_argument_types(signature: str) -> str:
    """
    Strips argument types from a function signature

    :param signature:
        function signature

    :return:
        string that can be used to construct function calls with the same
        variable names and ordering as in the function signature
    """
    # remove * prefix for pointers (pointer must always be removed before
    # values otherwise we will inadvertently dereference values,
    # same applies for const specifications)
    #
    # always add whitespace after type definition for cosmetic reasons
    known_types = [
        "const realtype *",
        "const double *",
        "const realtype ",
        "double *",
        "realtype *",
        "const int ",
        "int ",
        "SUNMatrixContent_Sparse ",
        "gsl::span<const int>",
    ]

    for type_str in known_types:
        signature = signature.replace(type_str, "")

    return signature


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

    return IDENTIFIER_PATTERN.match(x) is not None


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
    Custom Pow derivative that removes a removable singularity for
    ``self.base == 0`` and ``self.base.diff(s) == 0``. This function is
    intended to be monkeypatched into :py:method:`sympy.Pow._eval_derivative`.

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
        (part2, True),
    )


def _jacobian_element(i, j, eq_i, sym_var_j):
    """Compute a single element of a jacobian"""
    return (i, j), eq_i.diff(sym_var_j)


def _parallel_applyfunc(obj: sp.Matrix, func: Callable) -> sp.Matrix:
    """Parallel implementation of sympy's Matrix.applyfunc"""
    if (n_procs := int(os.environ.get("AMICI_IMPORT_NPROCS", 1))) == 1:
        # serial
        return obj.applyfunc(func)

    # parallel
    from multiprocessing import get_context
    from pickle import PicklingError

    from sympy.matrices.dense import DenseMatrix

    # "spawn" should avoid potential deadlocks occurring with fork
    #  see e.g. https://stackoverflow.com/a/66113051
    ctx = get_context("spawn")
    with ctx.Pool(n_procs) as p:
        try:
            if isinstance(obj, DenseMatrix):
                return obj._new(obj.rows, obj.cols, p.map(func, obj))
            elif isinstance(obj, sp.SparseMatrix):
                dok = obj.todok()
                mapped = p.map(func, dok.values())
                dok = {k: v for k, v in zip(dok.keys(), mapped) if v != 0}
                return obj._new(obj.rows, obj.cols, dok)
            else:
                raise ValueError(f"Unsupported matrix type {type(obj)}")
        except PicklingError as e:
            raise ValueError(
                f"Couldn't pickle {func}. This is likely because the argument "
                "was not a module-level function. Either rewrite the argument "
                "to a module-level function or disable parallelization by "
                "setting `AMICI_IMPORT_NPROCS=1`."
            ) from e
