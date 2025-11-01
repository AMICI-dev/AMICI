"""Miscellaneous functions related to model import, independent of any specific
model format"""

import enum
import itertools as itt
import numbers
import sys
from collections.abc import Callable, Iterable, Sequence
from typing import (
    Any,
    SupportsFloat,
)

import sympy as sp
from sympy.functions.elementary.piecewise import ExprCondPair
from sympy.logic.boolalg import BooleanAtom
from sympy.core.relational import Relational
from toposort import toposort

RESERVED_SYMBOLS = [
    "x",
    "k",
    "p",
    "y",
    "w",
    "h",
    "t",
    "AMICI_EMPTY_BOLUS",
    "NULL",
]

try:
    import pysb
except ImportError:
    pysb = None


class SBMLException(Exception):
    pass


SymbolDef = dict[sp.Symbol, dict[str, sp.Expr] | sp.Expr]


# Monkey-patch toposort CircularDependencyError to handle non-sortable objects,
#  such as sympy objects
class CircularDependencyError(ValueError):
    def __init__(self, data):
        # Sort the data just to make the output consistent, for use in
        #  error messages.  That's convenient for doctests.
        s = "Circular dependencies exist among these items: {{{}}}".format(
            ", ".join(
                f"{key!r}:{value!r}"
                for key, value in sorted(
                    {str(k): v for k, v in data.items()}.items()
                )
            )
        )
        super().__init__(s)
        self.data = data


setattr(
    sys.modules["toposort"], "CircularDependencyError", CircularDependencyError
)

annotation_namespace = "https://github.com/AMICI-dev/AMICI"


class ObservableTransformation(str, enum.Enum):
    """
    Different modes of observable transformation.
    """

    LOG10 = "log10"
    LOG = "log"
    LIN = "lin"


class MeasurementChannel:
    """
    A measurement channel (observable) definition.

    Measurement channels define how model state and parameters are mapped to
    observables, including any associated noise models.

    Measurement channels can be time-resolved (i.e., defined over the course of
    a simulation) or event-resolved (i.e., defined at specific events during
    a simulation). Event-resolved observables are associated with a specific
    event via the `event_id` attribute.
    """

    def __init__(
        self,
        id_: str,
        name: str | None = None,
        formula: str | sp.Expr | None = None,
        noise_distribution: str | Callable[[str], str] | None = None,
        sigma: str | sp.Expr | float | None = None,
        event_id: str | None = None,
    ):
        """
        Initialize a measurement channel.

        .. note::

            When providing expressions for (event) observables and their sigmas
            as strings (see below), those will be passed to
            :func:`sympy.sympify`. The supported grammar is not well-defined.
            Note there can be issues with, for example, ``==`` or n-ary (n>2)
            comparison operators.
            Also note that operator precedence and function names may differ
            from SBML L3 formulas or PEtab math expressions.
            Passing a sympy expression directly will
            be the safer option for more complex expressions.

        .. note::

            In any math expressions passed to this function, ``time`` will
            be interpreted as the time symbol.

        :param id_:
            Unique identifier for the measurement channel.
        :param name:
            Human-readable name for the measurement channel.
        :param formula:
            Expression defining how the observable is computed from model state
            and parameters.
        :param noise_distribution:
            Noise distribution associated with the observable.
            This is usually a string identifier (e.g., 'normal', 'log-normal';
            see
            :func:`amici.import_utils.noise_distribution_to_cost_function`).
            If ``None``, a normal distribution is assumed.

            Alternatively, a callable can be passed to account for a
            custom noise model. The callable must take a single argument
            ``str_symbol``, and return a string defining the negative
            log-likelihood contribution for a single data point, using
            variables ``{str_symbol}``, ``m{str_symbol}``, and
            ``sigma{str_symbol}`` for the simulation, measurement, and
            scale parameter, respectively.
        :param sigma:
            Expression representing the scale parameter of the noise
            distribution. This can be a numeric value, a sympy expression,
            or an expression string that will be passed to
            :func:`sympy.sympify`.
        :param event_id:
            Identifier of the associated event for event-resolved observables.
            `None` for time-resolved observables.

        Example usage:
        >>> # Time-resolved observable
        >>> mc1 = MeasurementChannel(
        ...     id_="obs1",
        ...     name="Observable 1",
        ...     formula="k1 * x1 + k2",
        ...     noise_distribution="log-normal",
        ...     sigma="sigma1"
        ... )
        >>> mc1  # doctest: +NORMALIZE_WHITESPACE
        MeasurementChannel(id_='obs1', name='Observable 1', \
formula='k1 * x1 + k2', noise_distribution='log-normal', \
sigma='sigma1', event_id=None)
        >>> mc1.is_time_resolved
        True
        >>> mc1.is_event_resolved
        False
        >>> # Event-resolved observable
        >>> mc2 = MeasurementChannel(
        ...     id_="obs2",
        ...     name="Observable 2",
        ...     formula="x3",
        ...     noise_distribution="log-normal",
        ...     sigma="sigma1",
        ...     event_id="event1"
        ... )
        >>> mc2  # doctest: +NORMALIZE_WHITESPACE
        MeasurementChannel(id_='obs2', name='Observable 2', formula='x3', \
 noise_distribution='log-normal', sigma='sigma1', event_id='event1')
        >>> mc2.is_event_resolved
        True
        >>> mc2.is_time_resolved
        False
        """
        self.id = id_
        self.name = name
        self.formula = formula
        self.noise_distribution = noise_distribution or "normal"
        self.sigma = sigma
        self.event_id = event_id

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(id_={self.id!r}, name={self.name!r}, "
            f"formula={self.formula!r}, noise_distribution="
            f"{self.noise_distribution!r}, sigma={self.sigma!r}, "
            f"event_id={self.event_id!r})"
        )

    @property
    def is_time_resolved(self):
        """Whether this is a time-resolved observable."""
        return self.event_id is None

    @property
    def is_event_resolved(self) -> bool:
        """Whether this is an event-resolved observable."""
        return self.event_id is not None


def noise_distribution_to_observable_transformation(
    noise_distribution: str | Callable,
) -> ObservableTransformation:
    """
    Parse noise distribution string and extract observable transformation

    :param noise_distribution:
        see :func:`noise_distribution_to_cost_function`

    :return:
        observable transformation
    """
    if isinstance(noise_distribution, str):
        if noise_distribution.startswith("log-"):
            return ObservableTransformation.LOG
        if noise_distribution.startswith("log10-"):
            return ObservableTransformation.LOG10

    return ObservableTransformation.LIN


def noise_distribution_to_cost_function(
    noise_distribution: str | Callable,
) -> Callable[[str], str]:
    """
    Parse noise distribution string to a cost function definition amici can
    work with.

    The noise distributions listed in the following are supported. :math:`m`
    denotes the measurement, :math:`y` the simulation, and :math:`\\sigma` a
    distribution scale parameter
    (currently, AMICI only supports a single distribution parameter).

    - `'normal'`, `'lin-normal'`: A normal distribution:

      .. math::
         \\pi(m|y,\\sigma) = \\frac{1}{\\sqrt{2\\pi}\\sigma}\\
         exp\\left(-\\frac{(m-y)^2}{2\\sigma^2}\\right)

    - `'log-normal'`: A log-normal distribution (i.e. log(m) is
      normally distributed):

      .. math::
         \\pi(m|y,\\sigma) = \\frac{1}{\\sqrt{2\\pi}\\sigma m}\\
         exp\\left(-\\frac{(\\log m - \\log y)^2}{2\\sigma^2}\\right)

    - `'log10-normal'`: A log10-normal distribution (i.e. log10(m) is
      normally distributed):

      .. math::
         \\pi(m|y,\\sigma) = \\frac{1}{\\sqrt{2\\pi}\\sigma m \\log(10)}\\
         exp\\left(-\\frac{(\\log_{10} m - \\log_{10} y)^2}{2\\sigma^2}\\right)

    - `'laplace'`, `'lin-laplace'`: A laplace distribution:

      .. math::
         \\pi(m|y,\\sigma) = \\frac{1}{2\\sigma}
         \\exp\\left(-\\frac{|m-y|}{\\sigma}\\right)

    - `'log-laplace'`: A log-Laplace distribution (i.e. log(m) is Laplace
      distributed):

      .. math::
         \\pi(m|y,\\sigma) = \\frac{1}{2\\sigma m}
         \\exp\\left(-\\frac{|\\log m - \\log y|}{\\sigma}\\right)

    - `'log10-laplace'`: A log10-Laplace distribution (i.e. log10(m) is
      Laplace distributed):

      .. math::
         \\pi(m|y,\\sigma) = \\frac{1}{2\\sigma m \\log(10)}
         \\exp\\left(-\\frac{|\\log_{10} m - \\log_{10} y|}{\\sigma}\\right)

    - `'binomial'`, `'lin-binomial'`: A (continuation of a discrete) binomial
      distribution, parameterized via the success probability
      :math:`p=\\sigma`:

      .. math::
         \\pi(m|y,\\sigma) = \\operatorname{Heaviside}(y-m) \\cdot
                \\frac{\\Gamma(y+1)}{\\Gamma(m+1) \\Gamma(y-m+1)}
                \\sigma^m (1-\\sigma)^{(y-m)}

    - `'negative-binomial'`, `'lin-negative-binomial'`: A (continuation of a
      discrete) negative binomial distribution, with with `mean = y`,
      parameterized via success probability `p`:

      .. math::

         \\pi(m|y,\\sigma) = \\frac{\\Gamma(m+r)}{\\Gamma(m+1) \\Gamma(r)}
            (1-\\sigma)^m \\sigma^r

      where

      .. math::
         r = \\frac{1-\\sigma}{\\sigma} y

    The distributions above are for a single data point.
    For a collection :math:`D=\\{m_i\\}_i` of data points and corresponding
    simulations :math:`Y=\\{y_i\\}_i` and noise parameters
    :math:`\\Sigma=\\{\\sigma_i\\}_i`, AMICI assumes independence,
    i.e. the full distributions is

    .. math::
       \\pi(D|Y,\\Sigma) = \\prod_i\\pi(m_i|y_i,\\sigma_i)

    AMICI uses the logarithm :math:`\\log(\\pi(m|y,\\sigma)`.

    In addition to the above-mentioned distributions, it is also possible to
    pass a function taking a symbol string and returning a log-distribution
    string with variables '{str_symbol}', 'm{str_symbol}', 'sigma{str_symbol}'
    for y, m, sigma, respectively.

    :param noise_distribution: An identifier specifying a noise model.
        Possible values are

        {`'normal'`, `'lin-normal'`, `'log-normal'`, `'log10-normal'`,
        `'laplace'`, `'lin-laplace'`, `'log-laplace'`, `'log10-laplace'`,
        `'binomial'`, `'lin-binomial'`, `'negative-binomial'`,
        `'lin-negative-binomial'`, `<Callable>`}

        For the meaning of the values see above.

    :return: A function that takes a strSymbol and then creates a cost
        function string (negative log-likelihood) from it, which can be
        sympified.
    """

    if isinstance(noise_distribution, Callable):
        return noise_distribution

    if noise_distribution in ["normal", "lin-normal"]:
        y_string = "0.5*log(2*pi*{sigma}**2) + 0.5*(({y} - {m}) / {sigma})**2"
    elif noise_distribution == "log-normal":
        y_string = (
            "0.5*log(2*pi*{sigma}**2*{m}**2) "
            "+ 0.5*((log({y}) - log({m})) / {sigma})**2"
        )
    elif noise_distribution == "log10-normal":
        y_string = (
            "0.5*log(2*pi*{sigma}**2*{m}**2*log(10)**2) "
            "+ 0.5*((log({y}, 10) - log({m}, 10)) / {sigma})**2"
        )
    elif noise_distribution in ["laplace", "lin-laplace"]:
        y_string = "log(2*{sigma}) + Abs({y} - {m}) / {sigma}"
    elif noise_distribution == "log-laplace":
        y_string = "log(2*{sigma}*{m}) + Abs(log({y}) - log({m})) / {sigma}"
    elif noise_distribution == "log10-laplace":
        y_string = (
            "log(2*{sigma}*{m}*log(10)) "
            "+ Abs(log({y}, 10) - log({m}, 10)) / {sigma}"
        )
    elif noise_distribution in ["binomial", "lin-binomial"]:
        # Binomial noise model parameterized via success probability p
        y_string = (
            "- log(Heaviside({y} - {m})) - loggamma({y}+1) "
            "+ loggamma({m}+1) + loggamma({y}-{m}+1) "
            "- {m} * log({sigma}) - ({y} - {m}) * log(1-{sigma})"
        )
    elif noise_distribution in ["negative-binomial", "lin-negative-binomial"]:
        # Negative binomial noise model of the number of successes m
        # (data) before r=(1-sigma)/sigma * y failures occur,
        # with mean number of successes y (simulation),
        # parameterized via success probability p = sigma.
        r = "{y} * (1-{sigma}) / {sigma}"
        y_string = (
            f"- loggamma({{m}}+{r}) + loggamma({{m}}+1) "
            f"+ loggamma({r}) - {r} * log(1-{{sigma}}) "
            f"- {{m}} * log({{sigma}})"
        )
    else:
        raise ValueError(
            f"Cost identifier {noise_distribution} not recognized."
        )

    def nllh_y_string(str_symbol):
        y, m, sigma = _get_str_symbol_identifiers(str_symbol)
        return y_string.format(y=y, m=m, sigma=sigma)

    return nllh_y_string


def _get_str_symbol_identifiers(str_symbol: str) -> tuple:
    """Get identifiers for simulation, measurement, and sigma."""
    y, m, sigma = f"{str_symbol}", f"m{str_symbol}", f"sigma{str_symbol}"
    return y, m, sigma


def smart_subs_dict(
    sym: sp.Expr,
    subs: SymbolDef,
    field: str | None = None,
    reverse: bool = True,
) -> sp.Expr:
    """
    Substitutes expressions completely flattening them out. Requires
    sorting of expressions with toposort.

    :param sym:
        Symbolic expression in which expressions will be substituted

    :param subs:
        Substitutions

    :param field:
        Field of substitution expressions in subs.values(), if applicable

    :param reverse:
        Whether ordering in subs should be reversed. Note that substitution
        requires the reverse order of what is required for evaluation.

    :return:
        Substituted symbolic expression
    """
    s = [
        (eid, expr[field] if field is not None else expr)
        for eid, expr in subs.items()
    ]
    if reverse:
        s.reverse()
    for substitution in s:
        # note that substitution may change free symbols, so we have to do
        # this recursively
        if sym.has(substitution[0]):
            sym = sym.subs(*substitution)
    return sym


def smart_subs(element: sp.Expr, old: sp.Symbol, new: sp.Expr) -> sp.Expr:
    """
    Optimized substitution that checks whether anything needs to be done first

    :param element:
        substitution target

    :param old:
        to be substituted

    :param new:
        subsitution value

    :return:
        substituted expression
    """
    return element.subs(old, new) if element.has(old) else element


def toposort_symbols(
    symbols: SymbolDef, field: str | None = None
) -> SymbolDef:
    """
    Topologically sort symbol definitions according to their interdependency

    :param symbols:
        symbol definitions

    :param field:
        field of definition.values() that is used to compute interdependency

    :return:
        ordered symbol definitions
    """
    sorted_symbols = toposort(
        {
            identifier: {
                s
                for s in (
                    definition[field] if field is not None else definition
                ).free_symbols
                if s in symbols
            }
            for identifier, definition in symbols.items()
        }
    )
    return {
        s: symbols[s]
        for symbol_group in sorted_symbols
        for s in sorted(symbol_group, key=str)
    }


def _parse_special_functions(
    sym: sp.Expr,
    parameters: Iterable[sp.Symbol],
    toplevel: bool = True,
) -> sp.Expr:
    """
    Recursively checks the symbolic expression for functions which have be
    to parsed in a special way, such as piecewise functions

    :param sym:
        symbolic expressions

    :param toplevel:
        as this is called recursively, are we in the top level expression?
    """

    args = tuple(
        arg
        if arg.__class__.__name__ == "piecewise"
        and sym.__class__.__name__ == "piecewise"
        else _parse_special_functions(arg, parameters, False)
        for arg in sym.args
    )

    fun_mappings = {
        "times": sp.Mul,
        "xor": sp.Xor,
        "abs": sp.Abs,
        "min": sp.Min,
        "max": sp.Max,
        "ceil": sp.functions.ceiling,
        "floor": sp.functions.floor,
        "factorial": sp.functions.factorial,
        "arcsin": sp.functions.asin,
        "arccos": sp.functions.acos,
        "arctan": sp.functions.atan,
        "arccot": sp.functions.acot,
        "arcsec": sp.functions.asec,
        "arccsc": sp.functions.acsc,
        "arcsinh": sp.functions.asinh,
        "arccosh": sp.functions.acosh,
        "arctanh": sp.functions.atanh,
        "arccoth": sp.functions.acoth,
        "arcsech": sp.functions.asech,
        "arccsch": sp.functions.acsch,
        "lt": sp.StrictLessThan,
        "gt": sp.StrictGreaterThan,
        "geq": sp.GreaterThan,
        "leq": sp.LessThan,
    }

    if sym.__class__.__name__ in fun_mappings:
        return fun_mappings[sym.__class__.__name__](*args)

    elif sym.__class__.__name__ == "piecewise" or isinstance(
        sym, sp.Piecewise
    ):
        if isinstance(sym, sp.Piecewise):
            # this is sympy piecewise, can't be nested
            denested_args = args
        else:
            # this is sbml piecewise, can be nested
            denested_args = _denest_piecewise(args, parameters)
        return _parse_piecewise_to_heaviside(denested_args, parameters)

    if sym.__class__.__name__ == "plus" and not sym.args:
        return sp.Float(0.0)

    if isinstance(sym, (sp.Function | sp.Mul | sp.Add | sp.Pow)):
        sym._args = args

    elif toplevel and isinstance(sym, BooleanAtom):
        # Replace boolean constants by numbers so they can be differentiated
        #  must not replace in Piecewise function. Therefore, we only replace
        #  it the complete expression consists only of a Boolean value.
        sym = sp.Float(int(bool(sym)))

    return sym


def _denest_piecewise(
    args: Sequence[sp.Expr | sp.logic.boolalg.Boolean | bool],
    parameters: Sequence[sp.Symbol] | None = None,
) -> tuple[sp.Expr | sp.logic.boolalg.Boolean | bool]:
    """
    Denest piecewise functions that contain piecewise as condition

    :param args:
        Arguments to the piecewise function

    :return:
        Arguments where conditions no longer contain piecewise functions and
        the conditional dependency is flattened out
    """
    if parameters is None:
        parameters = []

    args_out = []
    for coeff, cond in grouper(args, 2, True):
        # handling of this case is explicitely disabled in
        # _parse_special_functions as keeping track of coeff/cond
        # arguments is tricky. Simpler to just parse them out here
        if coeff.__class__.__name__ == "piecewise":
            coeff = _parse_special_functions(coeff, parameters, False)

        # we can have conditions that are piecewise function
        # returning True or False
        if cond.__class__.__name__ == "piecewise":
            # this keeps track of conditional that the previous
            # piece was picked
            previous_was_picked = sp.false
            # recursively denest those first
            for sub_coeff, sub_cond in grouper(
                _denest_piecewise(cond.args, parameters), 2, True
            ):
                # flatten the individual pieces
                pick_this = sp.And(sp.Not(previous_was_picked), sub_cond)
                if sub_coeff == sp.true:
                    args_out.extend([coeff, pick_this])
                previous_was_picked = pick_this

        else:
            args_out.extend([coeff, cond])
    # cut off last condition as that's the default
    return tuple(args_out[:-1])


def _is_c1_piecewise(
    pw: sp.Piecewise, parameters: Sequence[sp.Symbol]
) -> bool:
    """Return ``True`` if ``pw`` is continuously differentiable with respect to
    ``parameters``.

    This check ensures that piecewise expressions which are already smooth are
    not transformed into events.
    """

    from sympy.calculus.util import continuous_domain

    pieces = pw.args

    # collect boundaries appearing in the conditions
    boundaries: list[tuple[sp.Symbol, sp.Expr]] = []
    for _, cond in pieces:
        if cond in (True, False, sp.true, sp.false):
            continue
        if not isinstance(cond, Relational):
            return False
        if isinstance(cond.lhs, sp.Symbol) and not cond.rhs.has(cond.lhs):
            boundaries.append((cond.lhs, cond.rhs))
        elif isinstance(cond.rhs, sp.Symbol) and not cond.lhs.has(cond.rhs):
            boundaries.append((cond.rhs, cond.lhs))
        else:
            return False

    # check that each piece and its derivatives are continuous on R
    for expr, _ in pieces:
        for var in parameters:
            try:
                if continuous_domain(expr, var, sp.S.Reals) != sp.S.Reals:
                    return False
            except NotImplementedError:
                return False
            try:
                dexpr = sp.diff(expr, var)
                if continuous_domain(dexpr, var, sp.S.Reals) != sp.S.Reals:
                    return False
            except NotImplementedError:
                return False

    # check continuity and derivative continuity at boundaries
    for (sym, boundary), (expr_left, _), (expr_right, _) in zip(
        boundaries, pieces[:-1], pieces[1:]
    ):
        if not sp.simplify(
            expr_left.subs(sym, boundary) - expr_right.subs(sym, boundary)
        ).is_zero:
            return False
        for var in parameters:
            if not sp.simplify(
                sp.diff(expr_left, var).subs(sym, boundary)
                - sp.diff(expr_right, var).subs(sym, boundary)
            ).is_zero:
                return False

    return True


def _parse_piecewise_to_heaviside(
    args: Iterable[sp.Expr], parameters: Sequence[sp.Symbol] | None = None
) -> sp.Expr:
    """
    Piecewise functions cannot be transformed into C++ right away, but AMICI
    has a special interface for Heaviside functions, so we transform them.

    :param args:
        symbolic expressions for arguments of the piecewise function
    """
    # how many condition-expression pairs will we have?
    formula = sp.Integer(0)
    not_condition = sp.Integer(1)

    if all(isinstance(arg, ExprCondPair) for arg in args):
        # sympy piecewise
        grouped_args = args
    else:
        # smbl piecewise
        grouped_args = grouper(args, 2, True)

    pw = sp.Piecewise(*grouped_args)

    if parameters is None:
        parameters = []

    if _is_c1_piecewise(pw, parameters):
        return pw

    for coeff, trigger in grouped_args:
        if isinstance(coeff, BooleanAtom):
            coeff = sp.Integer(int(bool(coeff)))

        if trigger == sp.true:
            return formula + coeff * not_condition

        if trigger == sp.false:
            continue

        tmp = _parse_heaviside_trigger(trigger)
        formula += coeff * sp.simplify(not_condition * tmp)
        not_condition *= sp.Integer(1) - tmp

    return formula


def _parse_heaviside_trigger(trigger: sp.Expr) -> sp.Expr:
    """
    Recursively translates a boolean trigger function into a real valued
    root function

    :param trigger:
    :return: real valued root function expression
    """
    if trigger.is_Relational:
        root = trigger.args[0] - trigger.args[1]
        _check_unsupported_functions(root, "sympy.Expression")

        # normalize such that we always implement <,
        # this ensures that we can correctly evaluate the condition if
        # simulation starts at H(0). This is achieved by translating
        # conditionals into Heaviside functions H that is implemented as unit
        # step with H(0) = 1
        if isinstance(trigger, sp.core.relational.StrictLessThan):
            # x < y => x - y < 0 => r < 0
            return sp.Integer(1) - sp.Heaviside(root, 1)
        if isinstance(trigger, sp.core.relational.LessThan):
            # x <= y => not(y < x) => not(y - x < 0) => not -r < 0
            return sp.Heaviside(-root, 1)
        if isinstance(trigger, sp.core.relational.StrictGreaterThan):
            # y > x => y - x < 0 => -r < 0
            return sp.Integer(1) - sp.Heaviside(-root, 1)
        if isinstance(trigger, sp.core.relational.GreaterThan):
            # y >= x => not(x < y) => not(x - y < 0) => not r < 0
            return sp.Heaviside(root, 1)

    # rewrite n-ary XOR to OR to be handled below:
    trigger = trigger.replace(sp.Xor, _xor_to_or)
    # rewrite ==, !==
    trigger = trigger.replace(sp.Eq, _eq_to_and)
    trigger = trigger.replace(sp.Ne, _ne_to_or)

    # or(x,y) = not(and(not(x),not(y))
    if isinstance(trigger, sp.Or):
        return sp.Integer(1) - sp.Mul(
            *[
                sp.Integer(1) - _parse_heaviside_trigger(arg)
                for arg in trigger.args
            ]
        )

    if isinstance(trigger, sp.And):
        return sp.Mul(*[_parse_heaviside_trigger(arg) for arg in trigger.args])

    raise RuntimeError(
        "AMICI can not parse piecewise/event trigger functions with argument "
        f"{trigger}."
    )


def _xor_to_or(*args):
    """
    Replace XOR by OR expression.

    ``xor(x, y, ...) = (x & ~y & ...) | (~x & y & ...) | ...``.

    to be used in ``trigger = trigger.replace(sp.Xor, _xor_to_or)``.
    """
    res = sp.false
    for i in range(len(args)):
        res = sp.Or(
            res,
            sp.And(
                *(arg if i == j else sp.Not(arg) for j, arg in enumerate(args))
            ),
        )
    return res.simplify()


def _eq_to_and(*args):
    """
    Replace equality expression with numerical arguments by inequalities.

    ``Eq(x, y) = (x >= y) & (x <= y)``.

    to be used in ``trigger = trigger.replace(sp.Eq, _eq_to_and)``.
    """
    x, y = args
    return (x >= y) & (x <= y)


def _ne_to_or(*args):
    """
    Replace not-equal expression with numerical arguments by inequalities.

    ``Ne(x, y) = (x > y) | (x < y)``.

    to be used in ``trigger = trigger.replace(sp.Ne, _ne_to_or)``.

    This expects x and y to be not-NaN. No model should rely on NaN semantics
    anyways: In sympy, NaNs are equal, but they don't support <, >, >=, <=.
    In SBML, NaNs are equal, but support all comparisons. In IEEE 754, NaNs
    are not equal, but support all comparisons.
    """
    x, y = args
    return (x > y) | (x < y)


def grouper(
    iterable: Iterable, n: int, fillvalue: Any = None
) -> Iterable[tuple[Any]]:
    """
    Collect data into fixed-length chunks or blocks

    grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"

    :param iterable:
        any iterable

    :param n:
        chunk length

    :param fillvalue:
        padding for last chunk if length < n

    :return: itertools.zip_longest of requested chunks
    """
    args = [iter(iterable)] * n
    return itt.zip_longest(*args, fillvalue=fillvalue)


def _check_unsupported_functions(
    sym: sp.Expr, expression_type: str, full_sym: sp.Expr | None = None
):
    """
    Recursively checks the symbolic expression for unsupported symbolic
    functions

    :param sym:
        symbolic expressions

    :param expression_type:
        type of expression, only used when throwing errors

    :param full sym:
        outermost symbolic expression in recursive checks, only used for errors
    """
    if full_sym is None:
        full_sym = sym

    # note that sp.functions.factorial, sp.functions.ceiling,
    # sp.functions.floor applied to numbers should be simplified out and
    # thus pass this test
    unsupported_functions = (
        sp.functions.factorial,
        sp.functions.ceiling,
        sp.functions.floor,
        sp.functions.tan,
        sp.functions.sec,
        sp.functions.csc,
        sp.functions.cot,
        sp.functions.asec,
        sp.functions.acsc,
        sp.functions.acot,
        sp.functions.acsch,
        sp.functions.acoth,
        sp.Mod,
        sp.core.function.UndefinedFunction,
    )

    if (
        isinstance(sym.func, unsupported_functions)
        or isinstance(sym, unsupported_functions)
    ) and getattr(sym.func, "name", "") != "rateOf":
        raise RuntimeError(
            f"Encountered unsupported expression "
            f'"{sym.func}" of type '
            f'"{type(sym.func)}" as part of a '
            f'{expression_type}: "{full_sym}"!'
        )
    for arg in list(sym.args):
        _check_unsupported_functions(arg, expression_type)


def cast_to_sym(
    value: SupportsFloat | sp.Expr | BooleanAtom, input_name: str
) -> sp.Expr:
    """
    Typecasts the value to :py:class:`sympy.Float` if possible, and ensures the
    value is a symbolic expression.

    :param value:
        value to be cast

    :param input_name:
        name of input variable

    :return:
        typecast value
    """
    if isinstance(value, (sp.RealNumber | numbers.Number)):
        value = sp.Float(float(value))
    elif isinstance(value, BooleanAtom):
        value = sp.Float(float(bool(value)))

    if not isinstance(value, sp.Expr):
        raise TypeError(
            f"Couldn't cast {input_name} to sympy.Expr, was {type(value)}"
        )

    return value


def generate_measurement_symbol(observable_id: str | sp.Symbol):
    """
    Generates the appropriate measurement symbol for the provided observable

    :param observable_id:
        symbol (or string representation) of the observable

    :return:
        symbol for the corresponding measurement
    """
    if not isinstance(observable_id, str):
        observable_id = strip_pysb(observable_id)
    return symbol_with_assumptions(f"m{observable_id}")


def generate_regularization_symbol(observable_id: str | sp.Symbol):
    """
    Generates the appropriate regularization symbol for the provided observable

    :param observable_id:
        symbol (or string representation) of the observable

    :return:
        symbol for the corresponding regularization
    """
    if not isinstance(observable_id, str):
        observable_id = strip_pysb(observable_id)
    return symbol_with_assumptions(f"r{observable_id}")


def generate_flux_symbol(
    reaction_index: int, name: str | None = None
) -> sp.Symbol:
    """
    Generate identifier symbol for a reaction flux.
    This function will always return the same unique python object for a
    given entity.

    :param reaction_index:
        index of the reaction to which the flux corresponds
    :param name:
        an optional identifier of the reaction to which the flux corresponds
    :return:
        identifier symbol
    """
    if name is not None:
        return symbol_with_assumptions(name)

    return symbol_with_assumptions(f"flux_r{reaction_index}")


def symbol_with_assumptions(name: str):
    """
    Central function to create symbols with consistent, canonical assumptions

    :param name:
        name of the symbol

    :return:
        symbol with canonical assumptions
    """
    return sp.Symbol(name, real=True)


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


def unique_preserve_order(seq: Sequence) -> list:
    """Return a list of unique elements in Sequence, keeping only the first
    occurrence of each element

    Parameters:
        seq: Sequence to prune

    Returns:
        List of unique elements in ``seq``
    """
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


sbml_time_symbol = symbol_with_assumptions("time")
amici_time_symbol = symbol_with_assumptions("t")


def _default_simplify(x):
    """Default simplification applied in DEModel"""
    # We need this as a free function instead of a lambda to have it picklable
    #  for parallel simplification
    return sp.powsimp(x, deep=True)


def contains_periodic_subexpression(expr: sp.Expr, symbol: sp.Symbol) -> bool:
    """
    Check if a sympy expression contains any periodic subexpression.

    :param expr: The sympy expression to check.
    :param symbol: The variable with respect to which periodicity is checked.
    :return: `True` if the expression contains a periodic subexpression,
        `False` otherwise.
    """
    for subexpr in expr.atoms(sp.Function):
        if sp.periodicity(subexpr, symbol) is not None:
            return True
    return False
