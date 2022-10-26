"""Objects for AMICI's internal ODE model representation"""


import sympy as sp
import numbers

from typing import (
    Optional, Union, Dict, SupportsFloat, Set
)

from .import_utils import ObservableTransformation, \
    generate_measurement_symbol, generate_regularization_symbol,\
    RESERVED_SYMBOLS
from .import_utils import cast_to_sym

__all__ = [
    'ConservationLaw', 'Constant', 'Event', 'Expression', 'LogLikelihoodY',
    'LogLikelihoodZ', 'LogLikelihoodRZ', 'ModelQuantity', 'Observable',
    'Parameter', 'SigmaY', 'SigmaZ', 'State', 'EventObservable'
]


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

        if str(identifier) in RESERVED_SYMBOLS or \
                (hasattr(identifier, 'name') and
                 identifier.name in RESERVED_SYMBOLS):
            raise ValueError(f'Cannot add model quantity with name "{name}", '
                             f'please rename.')
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


class ConservationLaw(ModelQuantity):
    """
    A conservation law defines the absolute the total amount of a
    (weighted) sum of states

    """
    def __init__(self,
                 identifier: sp.Symbol,
                 name: str,
                 value: sp.Expr,
                 coefficients: Dict[sp.Symbol, sp.Expr],
                 state_id: sp.Symbol):
        """
        Create a new ConservationLaw instance.

        :param identifier:
            unique identifier of the ConservationLaw

        :param name:
            individual name of the ConservationLaw (does not need to be
            unique)

        :param value: formula (sum of states)

        :param coefficients:
            coefficients of the states in the sum

        :param state_id:
            identifier of the state that this conservation law replaces
        """
        self._state_expr: sp.Symbol = identifier - (value - state_id)
        self._coefficients: Dict[sp.Symbol, sp.Expr] = coefficients
        self._ncoeff: sp.Expr = coefficients[state_id]
        super(ConservationLaw, self).__init__(identifier, name, value)

    def get_ncoeff(self, state_id) -> Union[sp.Expr, int, float]:
        """
        Computes the normalized coefficient a_i/a_j where i is the index of
        the provided state_id and j is the index of the state that is
        replaced by this conservation law. This can be used to compute both
        dtotal_cl/dx_rdata (=ncoeff) and dx_rdata/dx_solver (=-ncoeff).

        :param state_id:
            identifier of the state

        :return: normalized coefficent of the state
        """
        return self._coefficients.get(state_id, 0.0) / self._ncoeff

    def get_x_rdata(self):
        """
        Returns the expression that allows computation of x_rdata for the state
        that this conservation law replaces.

        :return: x_rdata expression
        """
        return self._state_expr


class State(ModelQuantity):
    """
    A State variable defines an entity that evolves with time according to
    the provided time derivative, abbreviated by ``x``.

    :ivar _conservation_law:
        algebraic formula that allows computation of this
        state according to a conservation law

    :ivar _dt:
        algebraic formula that defines the temporal derivative of this state

    """
    def __init__(self,
                 identifier: sp.Symbol,
                 name: str,
                 init: sp.Expr,
                 dt: sp.Expr):
        """
        Create a new State instance. Extends :meth:`ModelQuantity.__init__`
        by ``dt``

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
        self._conservation_law: Union[ConservationLaw, None] = None

    def set_conservation_law(self, law: ConservationLaw) -> None:
        """
        Sets the conservation law of a state.

        If a conservation law is set, the respective state will be replaced by
        an algebraic formula according to the respective conservation law.

        :param law:
            linear sum of states that if added to this state remain
            constant over time
        """
        if not isinstance(law, ConservationLaw):
            raise TypeError(f'conservation law must have type ConservationLaw'
                            f', was {type(law)}')

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

    def get_free_symbols(self) -> Set[sp.Basic]:
        """
        Gets the set of free symbols in time derivative and initial conditions

        :return:
            free symbols
        """
        return self._dt.free_symbols.union(self._value.free_symbols)

    def has_conservation_law(self):
        """
        Checks whether this state has a conservation law assigned.

        :return: True if assigned, False otherwise
        """
        return self._conservation_law is not None

    def get_x_rdata(self):
        """
        Returns the expression that allows computation of x_rdata for this
        state, accounting for conservation laws.

        :return: x_rdata expression
        """
        if self._conservation_law is None:
            return self.get_id()
        else:
            return self._conservation_law.get_x_rdata()

    def get_dx_rdata_dx_solver(self, state_id):
        """
        Returns the expression that allows computation of
        ``dx_rdata_dx_solver`` for this state, accounting for conservation
        laws.

        :return: dx_rdata_dx_solver expression
        """
        if self._conservation_law is None:
            return sp.Integer(self._identifier == state_id)
        else:
            return -self._conservation_law.get_ncoeff(state_id)


class Observable(ModelQuantity):
    """
    An Observable links model simulations to experimental measurements,
    abbreviated by ``y``.

    :ivar _measurement_symbol:
        sympy symbol used in the objective function to represent
        measurements to this observable

    :ivar trafo:
        observable transformation, only applies when evaluating objective
        function or residuals
    """

    _measurement_symbol: Union[sp.Symbol, None] = None

    def __init__(self,
                 identifier: sp.Symbol,
                 name: str,
                 value: sp.Expr,
                 measurement_symbol: Optional[sp.Symbol] = None,
                 transformation: Optional[ObservableTransformation] = 'lin'):
        """
        Create a new Observable instance.

        :param identifier:
            unique identifier of the Observable

        :param name:
            individual name of the Observable (does not need to be unique)

        :param value:
            formula

        :param transformation:
            observable transformation, only applies when evaluating objective
            function or residuals
        """
        super(Observable, self).__init__(identifier, name, value)
        self._measurement_symbol = measurement_symbol
        self._regularization_symbol = None
        self.trafo = transformation

    def get_measurement_symbol(self) -> sp.Symbol:
        if self._measurement_symbol is None:
            self._measurement_symbol = generate_measurement_symbol(
                self.get_id()
            )

        return self._measurement_symbol

    def get_regularization_symbol(self) -> sp.Symbol:
        if self._regularization_symbol is None:
            self._regularization_symbol = generate_regularization_symbol(
                self.get_id()
            )

        return self._regularization_symbol


class EventObservable(Observable):
    """
    An Event Observable links model simulations to event related experimental
    measurements, abbreviated by ``z``.

    :ivar _event:
        symbolic event identifier
    """

    def __init__(self,
                 identifier: sp.Symbol,
                 name: str,
                 value: sp.Expr,
                 event: sp.Symbol,
                 measurement_symbol: Optional[sp.Symbol] = None,
                 transformation: Optional[ObservableTransformation] = 'lin',):
        """
        Create a new EventObservable instance.

        :param identifier:
            See :py:meth:`Observable.__init__`.

        :param name:
            See :py:meth:`Observable.__init__`.

        :param value:
            See :py:meth:`Observable.__init__`.

        :param transformation:
            See :py:meth:`Observable.__init__`.

        :param event:
            Symbolic identifier of the corresponding event.
        """
        super(EventObservable, self).__init__(identifier, name, value,
                                              measurement_symbol,
                                              transformation)
        self._event: sp.Symbol = event

    def get_event(self) -> sp.Symbol:
        """
        Get the symbolic identifier of the corresponding event.

        :return: symbolic identifier
        """
        return self._event


class Sigma(ModelQuantity):
    """
    A Standard Deviation Sigma rescales the distance between simulations
    and measurements when computing residuals or objective functions,
    abbreviated by ``sigma{y,z}``.
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
        if self.__class__.__name__ == "Sigma":
            raise RuntimeError(
                "This class is meant to be sub-classed, not used directly."
            )
        super(Sigma, self).__init__(identifier, name, value)


class SigmaY(Sigma):
    """
    Standard deviation for observables
    """


class SigmaZ(Sigma):
    """
    Standard deviation for event observables
    """


class Expression(ModelQuantity):
    """
    An Expression is a recurring elements in symbolic formulas. Specifying
    this may yield more compact expression which may lead to substantially
    shorter model compilation times, but may also reduce model simulation time.
    Abbreviated by ``w``.
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
    sensitivities may be computed, abbreviated by ``p``.
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
    sensitivities cannot be computed, abbreviated by ``k``.
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
    instances evaluated at all timepoints, abbreviated by ``Jy``.
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
        if self.__class__.__name__ == "LogLikelihood":
            raise RuntimeError(
                "This class is meant to be sub-classed, not used directly."
            )
        super(LogLikelihood, self).__init__(identifier, name, value)


class LogLikelihoodY(LogLikelihood):
    """
    Loglikelihood for observables
    """


class LogLikelihoodZ(LogLikelihood):
    """
    Loglikelihood for event observables
    """


class LogLikelihoodRZ(LogLikelihood):
    """
    Loglikelihood for event observables regularization
    """


class Event(ModelQuantity):
    """
    An Event defines either a SBML event or a root of the argument of a
    Heaviside function. The Heaviside functions will be tracked via the
    vector ``h`` during simulation and are needed to inform the ODE solver
    about a discontinuity in either the right-hand side or the states
    themselves, causing a reinitialization of the solver.
    """

    def __init__(self,
                 identifier: sp.Symbol,
                 name: str,
                 value: sp.Expr,
                 state_update: Union[sp.Expr, None],
                 initial_value: Optional[bool] = True):
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

        :param initial_value:
            initial boolean value of the trigger function at t0. If set to
            `False`, events may trigger at ``t==t0``, otherwise not.
        """
        super(Event, self).__init__(identifier, name, value)
        # add the Event specific components
        self._state_update = state_update
        self._initial_value = initial_value

    def get_initial_value(self) -> bool:
        """
        Return the initial value for the root function.

        :return:
            initial value formula
        """
        return self._initial_value

    def __eq__(self, other):
        """
        Check equality of events at the level of trigger/root functions, as we
        need to collect unique root functions for ``roots.cpp``
        """
        return self.get_val() == other.get_val() and \
            (self.get_initial_value() == other.get_initial_value())
