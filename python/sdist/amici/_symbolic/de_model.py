"""Symbolic differential equation model."""

from __future__ import annotations

import copy
import itertools
import logging
import re
from collections.abc import Callable, Sequence
from contextlib import suppress
from itertools import chain
from operator import itemgetter
from typing import TYPE_CHECKING

import numpy as np
import sympy as sp
from sympy import ImmutableDenseMatrix, MutableDenseMatrix

from amici.exporters.sundials.cxx_functions import (
    nobody_functions,
    sensi_functions,
    sparse_functions,
    var_in_function_signature,
)
from amici.exporters.sundials.cxxcodeprinter import csc_matrix
from amici.importers.utils import (
    ObservableTransformation,
    _default_simplify,
    amici_time_symbol,
    toposort_symbols,
    unique_preserve_order,
)
from amici.logging import get_logger, log_execution_time, set_log_level

from .de_model_components import (
    AlgebraicEquation,
    AlgebraicState,
    ConservationLaw,
    DifferentialState,
    Event,
    EventObservable,
    Expression,
    FixedParameter,
    FreeParameter,
    LogLikelihood,
    LogLikelihoodRZ,
    LogLikelihoodY,
    LogLikelihoodZ,
    ModelQuantity,
    NoiseParameter,
    Observable,
    ObservableParameter,
    Sigma,
    SigmaY,
    SigmaZ,
    State,
)
from .sympy_utils import (
    _parallel_applyfunc,
    smart_is_zero_matrix,
    smart_jacobian,
    smart_multiply,
)

if TYPE_CHECKING:
    from amici.importers.sbml.splines import AbstractSpline

logger = get_logger(__name__, logging.ERROR)


DERIVATIVE_PATTERN = re.compile(r"^d(x_rdata|xdot|\w+?)d(\w+?)(?:_explicit)?$")

__all__ = [
    "DEModel",
]


class DEModel:
    """
    Defines a Differential Equation as set of ModelQuantities.

    This class provides general purpose interfaces to compute arbitrary
    symbolic derivatives that are necessary for model simulation or
    sensitivity computation.

    All occurrences of a symbolic variable with a given name must use the same
    assumptions (e.g. real, positive, etc.) throughout the model. Mixing
    different assumptions for the same variable name will result in incorrect
    derivatives and potentially other errors.

    All symbols in the model are expected to be of type sympy.Symbol. If any
    subtypes are used, they must behave identically to sympy.Symbol in all
    relevant aspects (e.g. hashing, equality testing, etc.). In particular,
    `str(symbol)` is expected to return the same value as `symbol.name`.

    :ivar _differential_states:
        differential state variables

    :ivar _algebraic_states:
        algebraic state variables

    :ivar _observables:
        observables

    :ivar _event_observables:
        event observables

    :ivar _sigma_ys:
        sigmas for observables

    :ivar _sigma_zs:
        sigmas for event observables

    :ivar _free_parameters:
        parameters included in sensitivity analysis

    :ivar _fixed_parameters:
        parameters excluded from sensitivity analysis

    :ivar _log_likelihood_ys:
        loglikelihoods for observables

    :ivar _log_likelihood_zs:
        loglikelihoods for event observables

    :ivar _log_likelihood_rzs:
        loglikelihoods for event observable regularizations

    :ivar _expressions:
        expressions instances

    :ivar _conservation_laws:
        conservation laws

    :ivar _symboldim_funs:
        define functions that compute model dimensions, these are functions as the underlying symbolic expressions have
        not been populated at compile time

    :ivar _eqs:
        symbolic formulas of the symbolic variables of the model

    :ivar _sparseeqs:
        linear list of all symbolic formulas for sparsified variables

    :ivar _vals:
        numeric values of symbolic identifiers of the symbolic variables of the model

    :ivar _names:
        the names of symbolic identifiers of the symbolic variables of the model

    :ivar _syms:
        symbolic identifiers of the symbolic variables of the model

    :ivar _sparsesyms:
        linear list of all symbols for sparsified variables

    :ivar _colptrs:
        column pointers for sparsified variables. See SUNMatrixContent_Sparse definition in
        ``sunmatrix/sunmatrix_sparse.h``

    :ivar _rowvals:
        row values for sparsified variables. See SUNMatrixContent_Sparse definition in ``sunmatrix/sunmatrix_sparse.h``

    :ivar _equation_prototype:
        attribute from which an equation should be generated via list comprehension
        (see :meth:`OEModel._generate_equation`)

    :ivar _variable_prototype:
        attribute from which a variable should be generated via list comprehension
        (see :meth:`DEModel._generate_symbol`)

    :ivar _value_prototype:
        attribute from which a value should be generated via list comprehension (see :meth:`DEModel._generate_value`)

    :ivar _total_derivative_prototypes:
        defines how a total derivative equation is computed for an equation, key defines the name and values should be
        arguments for :meth:`DEModel.totalDerivative`

    :ivar _lock_total_derivative:
        add chainvariables to this set when computing total derivative from a partial derivative call to enforce a
        partial derivative in the next recursion. prevents infinite recursion

    :ivar _simplify:
        If not None, this function will be used to simplify symbolic
        derivative expressions. Receives sympy expressions as only argument.
        To apply multiple simplifications, wrap them in a lambda expression.

    :ivar _x0_fixedParameters_idx:
        Index list of subset of states for which x0_fixedParameters was computed

    :ivar _w_recursion_depth:
        recursion depth in w, quantified as nilpotency of dwdw

    :ivar _has_quadratic_nllh:
        whether all observables have a gaussian noise model, i.e. whether res and FIM make sense.

    :ivar _static_indices:
        indices of static variables for different model entities.

    :ivar _z2event:
        event indices for each event observable
    """

    def __init__(
        self,
        verbose: bool | int | None = False,
        simplify: Callable | None = _default_simplify,
        cache_simplify: bool = False,
        hybridisation: bool = False,
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
        self._free_parameters: list[FreeParameter] = []
        self._fixed_parameters: list[FixedParameter] = []
        self._log_likelihood_ys: list[LogLikelihoodY] = []
        self._log_likelihood_zs: list[LogLikelihoodZ] = []
        self._log_likelihood_rzs: list[LogLikelihoodRZ] = []
        self._noise_parameters: list[NoiseParameter] = []
        self._observable_parameters: list[ObservableParameter] = []
        self._expressions: list[Expression] = []
        self._conservation_laws: list[ConservationLaw] = []
        self._events: list[Event] = []
        self._splines: list[AbstractSpline] = []
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
            (sp.Matrix | sp.SparseMatrix | list[sp.Matrix | sp.SparseMatrix]),
        ] = dict()
        self._sparseeqs: dict[str, sp.Matrix | list[sp.Matrix]] = dict()
        self._vals: dict[str, list[sp.Expr]] = dict()
        self._names: dict[str, list[str]] = dict()
        self._syms: dict[str, sp.Matrix | list[sp.Matrix]] = dict()
        self._sparsesyms: dict[
            str, list[sp.Symbol] | list[list[sp.Symbol]]
        ] = dict()
        self._colptrs: dict[str, list[int] | list[list[int]]] = dict()
        self._rowvals: dict[str, list[int] | list[list[int]]] = dict()

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
            "p": self.free_parameters,
            "k": self.fixed_parameters,
            "w": self.expressions,
            "sigmay": self.sigma_ys,
            "sigmaz": self.sigma_zs,
            "h": self.events,
            "np": self.noise_parameters,
            "op": self.observable_parameters,
        }
        self._value_prototype: dict[str, Callable] = {
            "p": self.free_parameters,
            "k": self.fixed_parameters,
        }
        self._total_derivative_prototypes: dict[
            str, dict[str, str | list[str]]
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
                _simplified: dict[str, sp.Expr] = {},  # noqa B006
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
        self._x0_fixedParameters_idx: None | Sequence[int]
        self._w_recursion_depth: int = 0
        self._has_quadratic_nllh: bool = True
        set_log_level(logger, verbose)

        self._static_indices: dict[str, list[int]] = {}

    def differential_states(self) -> list[DifferentialState]:
        """Get all differential states."""
        return self._differential_states

    def algebraic_states(self) -> list[AlgebraicState]:
        """Get all algebraic states."""
        return self._algebraic_states

    def algebraic_equations(self) -> list[AlgebraicEquation]:
        """Get all algebraic equations."""
        return self._algebraic_equations

    def observables(self) -> list[Observable]:
        """Get all observables."""
        return self._observables

    def free_parameters(self) -> list[FreeParameter]:
        """Get all parameters."""
        return self._free_parameters

    def fixed_parameters(self) -> list[FixedParameter]:
        """Get all constants."""
        return self._fixed_parameters

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

    def noise_parameters(self) -> list[NoiseParameter]:
        """Get all noise parameters."""
        return self._noise_parameters

    def observable_parameters(self) -> list[ObservableParameter]:
        """Get all observable parameters."""
        return self._observable_parameters

    def is_ode(self) -> bool:
        """Check if model is ODE model."""
        return len(self._algebraic_equations) == 0

    def states(self) -> list[State]:
        """Get all states."""
        return self._differential_states + self._algebraic_states

    def add_component(
        self, component: ModelQuantity, insert_first: bool | None = False
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
            FreeParameter,
            FixedParameter,
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
            NoiseParameter,
            ObservableParameter,
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
            Symbol of the state that should be replaced by
            the conservation law (:math:`x_j`)

        :param total_abundance:
            Symbol of the total abundance (:math:`T/a_j`)

        :param coefficients:
            Dictionary of coefficients {x_i: a_i}
        """
        try:
            ix = next(
                filter(
                    lambda is_s: is_s[1].get_sym() == state,
                    enumerate(self._differential_states),
                )
            )[0]
        except StopIteration:
            raise ValueError(
                f"Specified state {state} was not found in the model states."
            )

        state_id = self._differential_states[ix].get_sym()

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

    def add_spline(self, spline: AbstractSpline, spline_expr: sp.Expr) -> None:
        """Add a spline to the model.

        :param spline:
            Spline instance to be added
        :param spline_expr:
            Sympy function representation of `spline` from
            ``spline.ode_model_symbol()``.
        """
        self._splines.append(spline)
        self.add_component(
            Expression(
                symbol=spline.sbml_id,
                name=str(spline.sbml_id),
                value=spline_expr,
            )
        )

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
        Number of Events that rely on numerical root-finding.

        :return:
            number of event symbols (length of the root vector in AMICI)
        """
        # TODO(performance): we could include constant `x` here as well
        #   (dx/dt = 0 AND not target of event assignments)
        #   this will require passing `x` to `fexplicit_roots`
        # TODO(performance): cache result, and don't repeat list symbols
        # Note that`self.num_events_solver` is currently
        #  NOT the same as len(self.get_implicit_roots())
        static_syms = self._static_symbols(["k", "p", "w"])
        return sum(
            not event.has_explicit_trigger_times(static_syms)
            for event in self.events()
        )

    def sym(self, name: str) -> sp.Matrix:
        """
        Returns (and constructs if necessary) the identifiers for a symbolic
        entity.

        :param name:
            name of the symbolic variable

        :return:
            matrix of symbols
        """
        if name not in self._syms:
            self._generate_symbol(name)

        return self._syms[name]

    def sparsesym(
        self, name: str, force_generate: bool = True
    ) -> list[sp.Symbol]:
        """
        Returns (and constructs if necessary) the sparsified symbols for
        a sparsified symbolic variable.

        :param name:
            name of the symbolic variable

        :param force_generate:
            whether the symbols should be generated if not available

        :return:
            linearized Matrix containing the symbols
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

    def colptrs(self, name: str) -> list[sp.Number] | list[list[sp.Number]]:
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

    def rowvals(self, name: str) -> list[sp.Number] | list[list[sp.Number]]:
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
                state.get_free_symbols() for state in self.states()
            )
        )

    def static_indices(self, name: str) -> list[int]:
        """
        Returns the indices of static expressions in the given model entity.

        Static expressions are those that do not depend on time,
        neither directly nor indirectly.

        :param name: Name of the model entity.
        :return: List of indices of static expressions.
        """
        # already computed?
        if (res := self._static_indices.get(name)) is not None:
            return res

        if name == "w":
            dwdx = self.sym("dwdx")
            dwdw = self.sym("dwdw")
            w = self.eq("w")

            # Check for direct (via `t`) or indirect (via `x`, `h`, or splines)
            # time dependency.
            # To avoid lengthy symbolic computations, we only check if we have
            # any non-zeros in hierarchy. We currently neglect the case where
            # different hierarchy levels may cancel out. Treating a static
            # expression as dynamic in such rare cases shouldn't be a problem.
            dynamic_dependency = np.asarray(
                dwdx.applyfunc(lambda x: int(not x.is_zero))
            ).astype(np.int64)
            # to check for other time-dependence, we add a column to the dwdx
            #  matrix
            dynamic_syms = [
                # FIXME: see spline comment below
                # *self.sym("spl"),
                *self.sym("h"),
                amici_time_symbol,
            ]
            dynamic_dependency = np.hstack(
                (
                    dynamic_dependency,
                    np.array(
                        [
                            expr.has(*dynamic_syms)
                            # FIXME: the current spline implementation is a giant pita
                            #  currently, the splines occur in the form of sympy functions, e.g.:
                            #   AmiciSpline(y0, time, y0_3, y0_1)
                            #   AmiciSplineSensitivity(y0, time, y0_1, y0_3, y0_1)
                            #  until it uses the proper self.sym("spl") / self.sym("sspl")
                            #  symbols, which will require quite some refactoring,
                            #  we just do dumb string checks :|
                            or (
                                bool(self._splines)
                                and "AmiciSpline" in str(expr)
                            )
                            for expr in w
                        ]
                    )[:, np.newaxis],
                )
            )

            nonzero_dwdw = np.asarray(
                dwdw.applyfunc(lambda x: int(not x.is_zero))
            ).astype(np.int64)

            # `w` is made up an expression hierarchy. Any given entry is only
            # static if all its dependencies are static. Here, we unravel
            # the hierarchical structure of `w`.
            # If for an entry in `w`, the row sum of the intermediate products
            # is 0 across all levels, the expression is static.
            tmp = dynamic_dependency
            res = np.sum(tmp, axis=1)
            while np.any(tmp != 0):
                tmp = nonzero_dwdw.dot(tmp)
                res += np.sum(tmp, axis=1)
            self._static_indices[name] = (
                np.argwhere(res == 0).flatten().tolist()
            )

            return self._static_indices[name]

        if name in ("dwdw", "dwdx", "dwdp"):
            static_indices_w = set(self.static_indices("w"))
            dynamic_syms = [
                *(
                    sym
                    for i, sym in enumerate(self.sym("w"))
                    if i not in static_indices_w
                ),
                amici_time_symbol,
                *self.sym("x"),
                *self.sym("h"),
                # FIXME see spline comment above
                # *(self.sym("spl") if name in ("dwdw", "dwdx") else ()),
                # *(self.sym("sspl") if name == "dwdp" else ()),
            ]
            dynamic_syms = sp.Matrix(dynamic_syms)
            rowvals = self.rowvals(name)
            sparseeq = self.sparseeq(name)

            # collect the indices of static expressions of dwd* from the list
            #  of non-zeros entries of the sparse matrix
            self._static_indices[name] = [
                i
                for i, (expr, row_idx) in enumerate(
                    zip(sparseeq, rowvals, strict=True)
                )
                # derivative of a static expression is static
                if row_idx in static_indices_w
                # constant expressions
                or expr.is_Number
                # check for dependencies on non-static entities
                or (
                    # FIXME see spline comment above
                    #  (check str before diff, as diff will fail on spline functions)
                    (
                        # splines: non-static
                        not self._splines or "AmiciSpline" not in str(expr)
                    )
                    and (
                        # If the expression contains dynamic symbols, it might
                        # still be static. However, checking the derivative
                        # is currently too expensive, and we rather accept
                        # false negatives.
                        not expr.has(*dynamic_syms)
                        # or all(
                        #     expr.diff(dyn_sym).is_zero
                        #     for dyn_sym in dynamic_syms
                        # )
                    )
                )
            ]
            return self._static_indices[name]

        raise NotImplementedError(name)

    def _static_symbols(self, names: list[str]) -> set[sp.Symbol]:
        """
        Return the static symbols among the given model entities.

        E.g., `static_symbols(["p", "w"])` returns all symbols in `p` and `w`
        that do not depend on time, neither directly nor indirectly.
        """
        result = set()

        for name in names:
            if name in ("k", "p"):
                result |= set(self.sym(name))
            elif name == "w":
                static_indices = self.static_indices("w")
                if len(static_indices) == 1:
                    result.add(self.sym("w")[static_indices[0]])
                elif len(static_indices) > 1:
                    result |= set(itemgetter(*static_indices)(self.sym("w")))
            else:
                raise ValueError(name)

        return result

    def dynamic_indices(self, name: str) -> list[int]:
        """
        Return the indices of dynamic expressions in the given model entity.

        :param name: Name of the model entity.
        :return: List of indices of dynamic expressions.
        """
        static_idxs = set(self.static_indices(name))
        length = len(
            self.sparsesym(name)
            if name in sparse_functions
            else self.sym(name)
        )
        return [i for i in range(length) if i not in static_idxs]

    def _generate_symbol(self, name: str) -> None:
        """
        Generates the symbols for a symbolic variable

        :param name:
            name of the symbolic variable
        """
        if name in self._variable_prototype:
            components = self._variable_prototype[name]()
            # ensure placeholder parameters are consistently and correctly ordered
            # we want that components are ordered by their placeholder index
            if name == "op":
                components = sorted(
                    components,
                    key=lambda x: int(
                        x.get_id().replace("observableParameter", "")
                    ),
                )
            if name == "np":
                components = sorted(
                    components,
                    key=lambda x: int(
                        x.get_id().replace("noiseParameter", "")
                    ),
                )
            self._syms[name] = sp.Matrix(
                [comp.get_sym() for comp in components]
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
                    state.get_sym()
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
                            f"s{tcl.get_id()}__{par.get_id()}",
                            real=True,
                        )
                        for par in self._free_parameters
                    ]
                    if self.conservation_law_has_multispecies(tcl)
                    else [0] * self.num_par()
                    for tcl in self._conservation_laws
                ]
            )
            return
        elif name == "x_old":
            length = len(self.eq("xdot"))
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
                [[f"spl_{isp}" for isp in range(len(self._splines))]]
            )
            return
        elif name == "sspl":
            # placeholders for spline sensitivities. Need to create symbols
            self._syms[name] = sp.Matrix(
                [
                    [f"sspl_{isp}_{ip}" for ip in range(len(self._syms["p"]))]
                    for isp in range(len(self._splines))
                ]
            )
            return
        elif name == "ih":
            self._syms[name] = sp.Matrix(
                [
                    sym
                    for sym, event in zip(self.sym("h"), self._events)
                    if not event.has_explicit_trigger_times()
                ]
            )
            return
        elif name == "eh":
            self._syms[name] = sp.Matrix(
                [
                    sym
                    for sym, event in zip(self.sym("h"), self._events)
                    if event.has_explicit_trigger_times()
                ]
            )
            return
        else:
            length = len(self.eq(name))
        self._syms[name] = sp.Matrix(
            [
                sp.Symbol(f"{name}{0 if name == 'stau' else i}", real=True)
                for i in range(length)
            ]
        )

    def generate_basic_variables(self) -> None:
        """
        Generates the symbolic identifiers for all variables in
        ``DEModel._variable_prototype``
        """
        self.parse_events()
        self._reorder_events()

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

        # Now add the found roots to the model components
        for root in roots:
            # skip roots of SBML events, as these have already been added
            if root in self._events:
                continue
            # add roots of heaviside functions
            self.add_component(root)

    def _reorder_events(self) -> None:
        """
        Re-order events - first those that require root tracking,
        then the others.
        """
        # Currently, the C++ simulations relies on the order of events:
        #  those that require numerical root-finding must come first, then
        #  those with explicit trigger times that don't depend on dynamic
        #  variables.
        #  TODO: This re-ordering here is a bit ugly, because we already need
        #   to generate certain model equations to perform this ordering.
        #   Ideally, we'd split froot into explicit and implicit parts during
        #   code generation instead (as already done for jax models).
        if not self.events():
            return

        static_syms = self._static_symbols(["k", "p", "w"])

        # ensure that we don't have computed any root-related symbols/equations
        #  yet, because the re-ordering might invalidate them
        # check after `self._static_symbols` which itself generates certain
        #  equations
        if (
            generated := set(self._syms)
            | set(self._eqs)
            | set(self._sparsesyms)
            | set(self._sparseeqs)
        ) and (
            "root" in generated
            or any(
                name.startswith("droot") or name.endswith("droot")
                for name in generated
            )
        ):
            raise AssertionError(
                "This function must be called before computing any "
                "root-related symbols/equations. "
                "The following symbols/equations are already "
                f"generated: {generated}"
            )

        self._events = list(
            chain(
                itertools.filterfalse(
                    lambda e: e.has_explicit_trigger_times(static_syms),
                    self._events,
                ),
                filter(
                    lambda e: e.has_explicit_trigger_times(static_syms),
                    self._events,
                ),
            )
        )

        # regenerate after re-ordering
        with suppress(KeyError):
            del self._syms["h"]
        self.sym("h")

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
            free_symbols_dt.count(
                str(self._differential_states[idx].get_sym())
            )
            + free_symbols_expr.count(
                str(self._differential_states[idx].get_sym())
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
                ) = csc_matrix(
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
            ) = csc_matrix(
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
                            # dx contains extra elements for algebraic states
                            strict=False,
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
                    state.get_sym()
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
                [y for spline in self._splines for y in spline.values_at_nodes]
            )

        elif name == "spline_slopes":
            # force symbols
            self._eqs[name] = sp.Matrix(
                [
                    d
                    for spline in self._splines
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
            # root(t, x(t), w(t, x(t)))
            # drootdt_total = drootdt + drootdx * dxdt + drootdw * dwdt_total
            # dwdt_total = dwdt + dwdx * dxdt
            self._eqs[name] = self.eq("drootdt")
            xdot = self.eq("xdot")

            # += drootdx * dxdt
            if self.num_states_solver() and not smart_is_zero_matrix(
                drootdx_explicit := smart_jacobian(
                    self.eq("root"), self.sym("x")
                )
            ):
                self._eqs[name] += smart_multiply(drootdx_explicit, xdot)

            # += drootdw * dwdt_total
            if not smart_is_zero_matrix(drootdw := self.eq("drootdw")):
                dwdt = self.eq("dwdt")
                dwdx_explicit = smart_jacobian(self.eq("w"), self.sym("x"))
                dwdt_total = dwdt + smart_multiply(dwdx_explicit, xdot)
                self._eqs[name] += smart_multiply(drootdw, dwdt_total)

        elif name == "deltax":
            # fill boluses for Heaviside functions, as empty state updates
            # would cause problems when writing the function file later
            event_eqs = []
            for event in self._events:
                state_update = event.get_state_update(
                    x=self.sym("x"), x_old=self.sym("x_old")
                )
                if state_update is None:
                    event_eqs.append(sp.zeros(self.num_states_solver(), 1))
                else:
                    event_eqs.append(state_update)

            self._eqs[name] = event_eqs

        elif name == "z":
            event_observables = [
                sp.zeros(self.num_eventobs(), 1) for _ in self._events
            ]
            event_ids = [e.get_sym() for e in self._events]
            z2event = [
                event_ids.index(event_obs.get_event())
                for event_obs in self._event_observables
            ]
            for (iz, ie), event_obs in zip(
                enumerate(z2event), self._event_observables, strict=True
            ):
                event_observables[ie][iz] = event_obs.get_val()

            self._eqs[name] = event_observables
            self._z2event = z2event

        elif name in [
            "ddeltaxdx",
            "ddeltaxdx_old",
            "ddeltaxdp",
            "ddeltaxdt",
            "dzdp",
            "dzdx",
        ]:
            if match_deriv[2] == "t":
                var = time_symbol
            else:
                var = self.sym(match_deriv[2])

            self._eqs[name] = [
                smart_jacobian(self.eq(match_deriv[1])[ie], var)
                for ie in range(self.num_events())
            ]
            if name == "dzdx":
                dtaudx = self.eq("dtaudx")
                for ie in range(self.num_events()):
                    for iz in range(self.num_eventobs()):
                        if ie != self._z2event[iz]:
                            continue
                        dzdt = sp.diff(self.eq("z")[ie][iz], time_symbol)
                        self._eqs[name][ie][iz, :] += dzdt * -dtaudx[ie]

        elif name in ["rz", "drzdx", "drzdp"]:
            eq_events = []
            for ie in range(self.num_events()):
                val = sp.zeros(
                    self.num_eventobs(),
                    1 if name == "rz" else len(self.sym(match_deriv[2])),
                )
                # match event observables to root function
                for iz in range(self.num_eventobs()):
                    if ie == self._z2event[iz]:
                        val[iz, :] = self.eq(name.replace("rz", "root"))[ie, :]
                eq_events.append(val)

            self._eqs[name] = eq_events

        elif name == "stau":
            self._eqs[name] = [
                self.eq("sroot")[ie, :] / self.eq("drootdt_total")[ie]
                if not self.eq("drootdt_total")[ie].is_zero
                else sp.zeros(*self.eq("sroot")[ie, :].shape)
                for ie in range(self.num_events())
            ]

        elif name == "dtaudx":
            drootdx_explicit = smart_jacobian(self.eq("root"), self.sym("x"))
            drootdt_total = self.eq("drootdt_total")
            drootdw = self.eq("drootdw")
            dwdx = self.eq("dwdx")

            self._eqs[name] = [
                (
                    # drootdx + drootdw * dwdx
                    drootdx_explicit[ie, :] + drootdw[ie, :] * dwdx
                )
                / drootdt_total[ie]
                if not drootdt_total[ie].is_zero
                else sp.zeros(*drootdx_explicit[ie, :].shape)
                for ie in range(self.num_events())
            ]

        elif name == "dtaudp":
            drootdp_explicit = smart_jacobian(self.eq("root"), self.sym("p"))
            drootdt_total = self.eq("drootdt_total")
            drootdw = self.eq("drootdw")
            drootdp = self.eq("drootdp")
            dwdp = self.eq("dwdp")
            self._eqs[name] = [
                (
                    # drootdp + drootdw * dwdp
                    drootdp_explicit[ie, :] + drootdw[ie, :] * dwdp
                )
                / drootdt_total[ie]
                if not drootdt_total[ie].is_zero
                else sp.zeros(*drootdp[ie, :].shape)
                for ie in range(self.num_events())
            ]

        elif name == "deltasx":
            if (
                self.num_states_solver() * self.num_par() * self.num_events()
                == 0
            ):
                self._eqs[name] = []
                return

            xdot_is_zero = smart_is_zero_matrix(self.eq("xdot"))

            event_eqs = []
            for ie, event in enumerate(self._events):
                tmp_eq = sp.zeros(self.num_states_solver(), self.num_par())

                # need to check if equations are zero since we are using
                # symbols
                stau_is_zero = smart_is_zero_matrix(self.eq("stau")[ie])
                if not stau_is_zero and not xdot_is_zero:
                    tmp_eq += smart_multiply(
                        self.sym("xdot") - self.sym("xdot_old"),
                        self.sym("stau").T,
                    )

                # only add deltax part if there is a state update
                if event.updates_state:
                    # partial derivative for the parameters
                    tmp_eq += self.eq("ddeltaxdp")[ie]

                    # initial part of chain rule state variables
                    tmp_dxdp = self.sym("sx") * sp.ones(1, self.num_par())

                    # need to check if equations are zero since we are using
                    # symbols
                    if not stau_is_zero:
                        # chain rule for the time point
                        tmp_eq += smart_multiply(
                            self.eq("ddeltaxdt")[ie],
                            -self.sym("stau").T,
                        )

                        # additional part of chain rule state variables
                        tmp_dxdp += smart_multiply(
                            self.sym("xdot_old"),
                            -self.sym("stau").T,
                        )

                    # finish chain rule for the state variables
                    tmp_eq += smart_multiply(
                        self.eq("ddeltaxdx")[ie]
                        + self.eq("ddeltaxdx_old")[ie],
                        tmp_dxdp,
                    )
                event_eqs.append(tmp_eq)

            self._eqs[name] = event_eqs

        elif name == "deltaxB":
            event_eqs = []
            for ie, event in enumerate(self._events):
                # ==== 1st group of terms: Heaviside functions ===========
                tmp_eq = smart_multiply(
                    self.sym("xdot") - self.sym("xdot_old"),
                    self.eq("dtaudx")[ie],
                )
                if event.updates_state:
                    # ==== 2nd group of terms: Derivatives of Dirac deltas ===
                    # Part 2a: explicit time dependence of bolus function
                    tmp_eq -= smart_multiply(
                        self.eq("ddeltaxdt")[ie], self.eq("dtaudx")[ie]
                    )
                    # Part 2b: implicit time dependence of bolus function
                    tmp_eq -= smart_multiply(
                        smart_multiply(
                            self.eq("ddeltaxdx")[ie]
                            + self.eq("ddeltaxdx_old")[ie],
                            self.sym("xdot_old"),
                        ),
                        self.eq("dtaudx")[ie],
                    )
                    # ==== 3rd group of terms: Dirac deltas ==================
                    tmp_eq += (
                        self.eq("ddeltaxdx")[ie] + self.eq("ddeltaxdx_old")[ie]
                    )
                tmp_eq = smart_multiply(self.sym("xB").T, tmp_eq)
                event_eqs.append(tmp_eq)
            self._eqs[name] = event_eqs

        elif name == "deltaqB":
            event_eqs = []
            for ie, event in enumerate(self._events):
                # ==== 1st group of terms: Heaviside functions ===========
                tmp_eq = smart_multiply(
                    self.sym("xdot") - self.sym("xdot_old"),
                    self.eq("dtaudp")[ie],
                )
                if event.updates_state:
                    # ==== 2nd group of terms: Derivatives of Dirac deltas ===
                    # Part 2a: explicit time dependence of bolus function
                    tmp_eq -= smart_multiply(
                        self.eq("ddeltaxdt")[ie], self.eq("dtaudp")[ie]
                    )
                    # Part 2b: implicit time dependence of bolus function
                    tmp_eq -= smart_multiply(
                        smart_multiply(
                            self.eq("ddeltaxdx")[ie]
                            + self.eq("ddeltaxdx_old")[ie],
                            self.sym("xdot_old"),
                        ),
                        self.eq("dtaudp")[ie],
                    )
                    # ==== 3rd group of terms: Dirac deltas ==================
                    tmp_eq += self.eq("ddeltaxdp")[ie]
                event_eqs.append(smart_multiply(self.sym("xB").T, tmp_eq))
            self._eqs[name] = event_eqs

        elif name == "xdot_old":
            # force symbols
            self._eqs[name] = self.sym(name)

        elif name == "dwdx":
            if (
                expected := list(
                    map(
                        ConservationLaw.get_x_rdata,
                        reversed(self.conservation_laws()),
                    )
                )
            ) != (actual := self.eq("w")[: self.num_cons_law()]):
                raise AssertionError(
                    "Conservation laws are not at the beginning of 'w'. "
                    f"Got {actual}, expected {expected}."
                )
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

        elif name == "dwdt":
            self._eqs[name] = smart_jacobian(self.eq("w"), time_symbol)

        elif name == "iroot":
            self._eqs[name] = sp.Matrix(
                [
                    eq
                    for eq, event in zip(
                        self.eq("root"), self._events, strict=True
                    )
                    if not event.has_explicit_trigger_times()
                ]
            )

        elif name == "eroot":
            self._eqs[name] = sp.Matrix(
                [
                    eq
                    for eq, event in zip(
                        self.eq("root"), self._events, strict=True
                    )
                    if event.has_explicit_trigger_times()
                ]
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
            for i, (sym_yz, eq_sigma_yz) in enumerate(
                zip(
                    syms_yz,
                    self._eqs[name],
                    strict=True,
                )
            ):
                yz_free_syms = eq_sigma_yz.free_symbols
                if tmp := {x for x in syms_x if x in yz_free_syms}:
                    # Can we replace x symbols by an equivalent observable?
                    #  (currently, only the matching observable is supported)
                    x_to_eliminate = next(iter(tmp))
                    eq_yz = (
                        self.eq("y")[i]
                        if name == "sigmay"
                        else self.eq("z")[self._z2event[i]][i]
                    )

                    try:
                        # solve for the next best x symbol and substitute
                        #  (maybe try another one if we fail?)
                        replacement = sp.solve(
                            sp.Eq(sym_yz, eq_yz), x_to_eliminate
                        )
                    except NotImplementedError:
                        # can't solve
                        replacement = []

                    if len(replacement) == 1:
                        self._eqs[name][i] = self._eqs[name][i].subs(
                            x_to_eliminate, replacement[0]
                        )
                        if not any(
                            x in self._eqs[name][i].free_symbols
                            for x in syms_x
                        ):
                            # successfully eliminated x symbols
                            continue

                    # Report all x symbols that cannot be replaced
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
            ("root", "w"): "tcl",  # dtcldw = 0
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

        elif name in ("dydw", "drootdw") and not smart_is_zero_matrix(
            derivative
        ):
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
        if var_in_function_signature(name, varname, self.is_ode()):
            if varname in [
                "dwdx",
                "dwdp",
            ]:
                # dwdx and dwdp will be dynamically computed, and their
                #  ordering within a column may differ from the initialization
                #  of symbols here, so those are not safe to use.
                raise AssertionError()
            return self.sym(varname)
        else:
            return self.eq(varname)

    def _multiplication(
        self,
        name: str,
        x: str,
        y: str,
        transpose_x: bool | None = False,
        sign: int | None = 1,
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
            (state.get_sym(), state.get_x_rdata())
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
            fp in (c.get_sym() for c in self._fixed_parameters)
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

        return state.get_dt().is_zero

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

        This function is solely for checking whether `expr` has a discontinuity
        we need to track. Better report a false positive than miss a
        time-dependence.

        :param expr:
            The expression.

        :returns:
            Whether the expression is time-dependent.
        """
        if not (free_syms := expr.free_symbols):
            return False

        free_syms -= set(self.sym("p"))

        if not free_syms:
            return False

        free_syms -= set(self.sym("k"))

        # TODO(performance): handle static expressions,
        #  handle constant state variables,
        #  combine with other checks for time-dependence
        #   after https://github.com/AMICI-dev/AMICI/pull/3031
        return bool(free_syms)

    def _get_unique_root(
        self,
        root_found: sp.Expr,
        roots: list[Event],
    ) -> sp.Symbol | None:
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
        for root in roots:
            if (difference := (root_found - root.get_val())).is_zero or (
                self._simplify and self._simplify(difference).is_zero
            ):
                return root.get_sym()

        # create an event for a new root function
        root_symstr = f"Heaviside_{len(roots)}"
        roots.append(
            Event(
                symbol=sp.Symbol(root_symstr),
                name=root_symstr,
                value=root_found,
                assignments=None,
                use_values_from_trigger_time=True,
            )
        )

        if not self._expr_is_time_dependent(root_found):
            # Not time-dependent. Return None, but we still need to create
            # the event above
            return None

        return roots[-1].get_sym()

    def _collect_heaviside_roots(
        self,
        args: Sequence[sp.Basic],
    ) -> list[tuple[sp.Expr, sp.Expr]]:
        """
        Recursively check an expression for the occurrence of Heaviside
        functions and return all roots found.

        :param args:
            args attribute of the expanded expression

        :returns:
            List of (root function, Heaviside x0)-tuples that were extracted
            from Heaviside function arguments.
        """
        root_funs = []
        for arg in args:
            if arg.func == sp.Heaviside:
                root_funs.append(arg.args)
            elif arg.has(sp.Heaviside):
                root_funs.extend(self._collect_heaviside_roots(arg.args))

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
        tmp_roots_old = self._collect_heaviside_roots((dxdt,))
        for tmp_root_old, tmp_x0_old in unique_preserve_order(tmp_roots_old):
            # we want unique identifiers for the roots
            tmp_root_new = self._get_unique_root(tmp_root_old, roots)
            # `tmp_new` is None if the root is not time-dependent.
            if tmp_root_new is None:
                continue
            # For Heavisides, we need to add the negative function as well
            self._get_unique_root(sp.sympify(-tmp_root_old), roots)
            heavisides.append(
                (sp.Heaviside(tmp_root_old, tmp_x0_old), tmp_root_new)
            )

        if heavisides:
            # only apply subs if necessary
            for heaviside_sympy, heaviside_amici in heavisides:
                dxdt = dxdt.subs(heaviside_sympy, heaviside_amici)

        return dxdt

    @property
    def _components(self) -> list[ModelQuantity]:
        """
        Returns the components of the model

        :return:
            components of the model
        """
        return (
            self._algebraic_states
            + self._algebraic_equations
            + self._conservation_laws
            + self._fixed_parameters
            + self._differential_states
            + self._event_observables
            + self._events
            + self._expressions
            + self._log_likelihood_ys
            + self._log_likelihood_zs
            + self._log_likelihood_rzs
            + self._observables
            + self._free_parameters
            + self._sigma_ys
            + self._sigma_zs
            + self._splines
        )

    def _process_hybridization(self, hybridization: dict) -> None:
        """
        Parses the hybridization information and updates the model accordingly

        :param hybridization:
            dict representation of the hybridization information in the PEtab YAML file, see
            https://petab-sciml.readthedocs.io/latest/format.html#problem-yaml-file
        """
        added_expressions = False
        orig_obs = tuple([s.get_sym() for s in self._observables])
        for net_id, net in hybridization.items():
            if net["static"]:
                continue  # do not integrate into ODEs, handle in amici.jax.petab
            inputs = [
                comp
                for comp in self._components
                if str(comp.get_sym()) in net["input_vars"]
            ]
            # sort inputs by order in input_vars
            inputs = sorted(
                inputs,
                key=lambda comp: net["input_vars"].index(str(comp.get_sym())),
            )
            if len(inputs) != len(net["input_vars"]):
                found_vars = {str(comp.get_sym()) for comp in inputs}
                missing_vars = set(net["input_vars"]) - found_vars
                raise ValueError(
                    f"Could not find all input variables for neural network {net_id}. "
                    f"Missing variables: {sorted(missing_vars)}"
                )
            for inp in inputs:
                if isinstance(
                    inp,
                    Sigma
                    | LogLikelihood
                    | Event
                    | ConservationLaw
                    | Observable,
                ):
                    raise NotImplementedError(
                        f"{inp.get_name()} ({type(inp)}) is not supported as neural network input."
                    )

            outputs = {
                out_var: {"comp": comp, "ind": net["output_vars"][out_var]}
                for comp in self._components
                if (out_var := str(comp.get_sym())) in net["output_vars"]
                # TODO: SYNTAX NEEDS to CHANGE
                or (out_var := str(comp.get_sym()) + "_dot")
                in net["output_vars"]
            }
            if len(outputs.keys()) != len(net["output_vars"]):
                found_vars = set(outputs.keys())
                missing_vars = set(net["output_vars"]) - found_vars
                raise ValueError(
                    f"Could not find all output variables for neural network {net_id}. "
                    f"Missing variables: {sorted(missing_vars)}"
                )

            for out_var, parts in outputs.items():
                comp = parts["comp"]
                # remove output from model components
                if isinstance(comp, FreeParameter):
                    self._free_parameters.remove(comp)
                elif isinstance(comp, Expression):
                    self._expressions.remove(comp)
                elif isinstance(comp, DifferentialState):
                    pass
                else:
                    raise NotImplementedError(
                        f"{comp.get_name()} ({type(comp)}) is not supported as neural network output."
                    )

                # generate dummy Function
                out_val = sp.Function(net_id)(
                    *[input.get_sym() for input in inputs], parts["ind"]
                )

                # add to the model
                if isinstance(comp, DifferentialState):
                    ix = self._differential_states.index(comp)
                    # TODO: SYNTAX NEEDS to CHANGE
                    if out_var.endswith("_dot"):
                        self._differential_states[ix].set_dt(out_val)
                    else:
                        self._differential_states[ix].set_val(out_val)
                else:
                    self.add_component(
                        Expression(
                            symbol=comp.get_sym(),
                            name=net_id,
                            value=out_val,
                        )
                    )
                    added_expressions = True

            observables = {
                ob_var: {"comp": comp, "ind": net["observable_vars"][ob_var]}
                for comp in self._components
                if (ob_var := str(comp.get_sym())) in net["observable_vars"]
                # # TODO: SYNTAX NEEDS to CHANGE
                # or (ob_var := str(comp.get_id()) + "_dot")
                # in net["observable_vars"]
            }
            if len(observables.keys()) != len(net["observable_vars"]):
                found_vars = set(observables.keys())
                missing_vars = set(net["observable_vars"]) - found_vars
                raise ValueError(
                    f"Could not find all observable variables for neural network {net_id}. "
                    f"Missing variables: {sorted(missing_vars)}"
                )

            for ob_var, parts in observables.items():
                comp = parts["comp"]
                if isinstance(comp, Observable):
                    self._observables.remove(comp)
                else:
                    raise ValueError(
                        f"{comp.get_name()} ({type(comp)}) is not an observable."
                    )
                out_val = sp.Function(net_id)(
                    *[input.get_sym() for input in inputs], parts["ind"]
                )
                # add to the model
                self.add_component(
                    Observable(
                        symbol=comp.get_sym(),
                        name=net_id,
                        value=out_val,
                    )
                )

        new_order = [orig_obs.index(s.get_sym()) for s in self._observables]
        self._observables = [self._observables[i] for i in new_order]

        if added_expressions:
            self.toposort_expressions()

    def get_explicit_roots(self) -> set[sp.Expr]:
        """
        Returns explicit formulas for all discontinuities (events)
        that can be precomputed

        :return:
            set of symbolic roots
        """
        return {root for e in self._events for root in e.get_trigger_times()}

    def get_implicit_roots(self) -> set[sp.Expr]:
        """
        Returns implicit equations for all discontinuities (events)
        that have to be located via rootfinding

        :return:
            set of symbolic roots
        """
        return {
            e.get_val()
            for e in self._events
            if not e.has_explicit_trigger_times()
        }

    def has_algebraic_states(self) -> bool:
        """
        Checks whether the model has algebraic states

        :return:
            boolean indicating if algebraic states are present
        """
        return len(self._algebraic_states) > 0

    def has_event_assignments(self) -> bool:
        """
        Checks whether the model has event assignments

        :return:
            boolean indicating if event assignments are present
        """
        return any(event.updates_state for event in self._events)

    def toposort_expressions(
        self, reorder: bool = True
    ) -> dict[sp.Symbol, sp.Expr]:
        """
        Sort expressions in topological order.

        :param reorder:
            Whether to reorder the internal expression list (``True``) or
            just return the toposorted expressions (``False``).

        :return:
            dict of expression symbols to expressions in topological order
        """
        # NOTE: elsewhere, conservations law expressions are expected to
        #  occur before any other w expressions, so we must maintain their
        #  position.
        # toposort everything but conservation law expressions,
        #  then prepend conservation laws

        w_toposorted = toposort_symbols(
            {
                e.get_sym(): e.get_val()
                for e in self.expressions()[self.num_cons_law() :]
            }
        )

        w_toposorted = {
            e.get_sym(): e.get_val()
            for e in self.expressions()[: self.num_cons_law()]
        } | w_toposorted

        if not reorder:
            return w_toposorted

        # ensure no symbols or equations that depend on `w` have been generated
        #  yet, otherwise the re-ordering might break dependencies
        if (
            generated := set(self._syms)
            | set(self._eqs)
            | set(self._sparsesyms)
            | set(self._sparseeqs)
        ) - {"w", "p", "k", "x", "x_rdata"}:
            raise AssertionError(
                "This function must be called before computing any "
                "derivatives. The following symbols/equations are already "
                f"generated: {generated}"
            )

        old_syms = tuple(e.get_sym() for e in self.expressions())
        topo_expr_syms = tuple(w_toposorted)
        new_order = [old_syms.index(s) for s in topo_expr_syms]
        self._expressions = [self._expressions[i] for i in new_order]
        self._syms["w"] = sp.Matrix(topo_expr_syms)
        self._eqs["w"] = sp.Matrix(list(w_toposorted.values()))

        return w_toposorted
