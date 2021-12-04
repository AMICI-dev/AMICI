"""
Splines
-------
This module provides helper functions for reading/writing splines with AMICI
annotations from/to SBML files and for adding such splines to the AMICI C++
code.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .sbml_import import SbmlImporter
    from typing import (
        Sequence,
        Union,
        Tuple,
        Optional,
        Dict,
        Any,
        List,
    )

    BClike = Union[
        None, str,
        Tuple[Union[None, str], Union[None, str]]
    ]

import collections.abc
import xml.etree.ElementTree as ET
import numpy as np
import sympy as sp
import libsbml
from abc import abstractmethod, ABC
from itertools import count
from sympy.core.parameters import evaluate
from .sbml_utils import (
    sbml_time_symbol,
    amici_time_symbol,
    pretty_xml,
    mathml2sympy,
    sbmlMathML,
    annotation_namespace,
    getSbmlUnits,
    hasParameter,
    addParameter,
    addAssignmentRule,
    SbmlException
)
from numbers import Integral


def sympify_noeval(x):
    with evaluate(False):
        return sp.sympify(x)


###############################################################################


class UniformGrid(collections.abc.Sequence):
    """
    A grid of uniformly-spaced real points, computed with rational arithmetic.

    Implements the :py:class:`collections.abc.Sequence` interface and can be
    converted to a :py:class:`numpy.ndarray` (conversion to float can be
    specified with ``dtype=float``).

    TODO: document attributes
    """

    def __init__(
            self, start, stop, step=None, *,
            length: Optional[Integral] = None,
            include_stop: bool = True
    ):
        """
        Create a new ``UniformGrid``.
        :param start:
            First point in the grid
        :param stop:
            Last point in the grid (some caveats apply, see ``include_stop``)
        :param step:
            Desired step size of the grid. Mutually exclusive with ``length``.
        :param length:
            Desired length of the grid. Mutually exclusive with ``step``.
        :param include_stop:
            Controls the behaviour when ``step`` is not ``None`` and
            ``stop - start`` is not an integer multiple of ``step``.
            If ``True`` (default), the actual endpoint of the grid will be
            larger than ``stop``. Otherwise it will be smaller.
        """
        start = sp.nsimplify(sp.sympify(start))
        stop = sp.nsimplify(sp.sympify(stop))
        if step is None:
            if length is None:
                raise ValueError("One of step/length must be specified!")
            elif not isinstance(length, Integral):
                raise TypeError("Length must be an integer!")
            elif length < 2:
                raise ValueError("Length must be at least 2!")
            step = (stop - start) / (length - 1)
        elif length is not None:
            raise ValueError("Only one of step/length can be specified!")
        else:
            step = sp.nsimplify(sp.sympify(step))

        if start > stop:
            raise ValueError(
                f'Start point {start} greater than stop point {stop}!'
            )

        if step <= 0:
            raise ValueError(f'Step size {step} must be strictly positive!')

        xx = []
        for i in count():
            x = start + i * step
            if not include_stop and x > stop:
                break
            xx.append(x)
            if include_stop and x >= stop:
                break

        self._xx = np.asarray(xx)

    @property
    def start(self):
        """First point."""
        return self._xx[0]

    @property
    def stop(self):
        """Last point."""
        return self._xx[-1]

    @property
    def step(self):
        """Distance between consecutive points."""
        if len(self._xx) > 1:
            return self._xx[1] - self._xx[0]
        return sp.core.numbers.Zero

    def __getitem__(self, i):
        return self._xx[i]

    def __len__(self):
        return len(self._xx)

    def __array__(self, dtype=None):
        if dtype is None:
            return self._xx
        return np.array(self._xx, dtype=dtype)

    def __repr__(self):
        return (f"UniformGrid(start={self.start}, stop={self.stop}, "
                f"step={self.step})")


###############################################################################


class AbstractSpline(ABC):
    """
    Base class for spline functions which can be computed efficiently
    thanks to tailored C++ implementations in AMICI.
    Inside an SBML file, such splines are implemented with
    an assignment rule containing both a symbolic piecewise formula
    for the spline (allowing compatibility with any SBML-aware software)
    and annotations which encode the necessary information for AMICI to
    recreate the spline object (allowing for fast computations when the SBML
    file is used together with AMICI).
    """

    def __init__(
            self,
            sbml_id: Union[str, sp.Symbol],
            x: Union[str, sp.Basic],
            xx: Sequence,
            yy: Sequence,
            *,
            bc: BClike = None,
            extrapolate: BClike = None,
            logarithmic_parametrization: bool = False
    ):
        """
        Base constructor for ``AbstractSpline`` objects.

        :param sbml_id:
            The SBML ID of the parameter associated to the spline
            as a string or a SymPy symbol.

        :param x:
            The point at which the spline is evaluated.
            It will be sympified.

        :param xx:
            The points at which the spline values are known.
            Currently they can only depend on constant parameters.
            These points should be strictly increasing.
            This argument will be sympified.

        :param yy:
            The spline values at each of the points in ``xx``.
            They may not depend on model species.
            This argument will be sympified.

        :param bc:
            Tuple of applied boundary conditions, one for each side of the
            spline domain. If a single boundary condition is given it will be
            applied to both sides.
            Possible boundary conditions
            (allowed values depend on the ``AbstractSpline`` subclass):
            * `None` or `'no_bc'`: boundary conditions are not needed
              for this spline object;
            * `'zeroderivative'`: first derivative set to zero;
            * `'natural'`: second derivative set to zero;
            * `'zeroderivative+natural'`: first and second derivatives
              set to zero;
            * `'periodic'`: periodic bc.

        :param extrapolate:
            Whether to extrapolate the spline outside the base interval
            defined by ``(xx[0], xx[-1])``.
            It is a tuple of extrapolation methods, one for each side of the
            base interval.
            If it is not a tuple, then the same extrapolation will be applied
            on both sides.
            Extrapolation methods supported:
            * `None` or `'no_extrapolation'`: no extrapolation should be
              performed. An exception will be raise in the C++ code
              if the spline is evaluated outside the base interval.
              In the fallback SBML symbolic expression
              `'polynomial'` extrapolation will be used.
            * `'polynomial'`: the cubic polynomial used in the nearest
              spline segment will be used.
            * `'constant'`: constant extrapolation will be used.
              Requires `'zeroderivative'` boundary condition.
              For splines which are continuous up to the second derivative,
              it requires the stricter `'zeroderivative+natural'`
              boundary condition.
            * `'linear'`: linear extrapolation will be used.
              For splines which are continuous up to the second derivative,
              this requires the `'natural'` boundary condition.
            * `'periodic'`: periodic extrapolation.
              Requires `'periodic'` boundary conditions.

        :param logarithmic_parametrization:
            Whether interpolation should be done in log-scale.
        """

        if isinstance(sbml_id, str):
            sbml_id = sp.Symbol(sbml_id)
        elif not isinstance(sbml_id, sp.Symbol):
            raise TypeError(
                'sbml_id must be either a string or a SymPy symbol, '
                f'got {sbml_id} instead!'
            )

        x = sympify_noeval(x)
        if not isinstance(x, sp.Basic):
            # It may still be e.g. a list!
            raise ValueError(f'Invalid x = {x}!')

        if not isinstance(xx, UniformGrid):
            xx = np.asarray([sympify_noeval(x) for x in xx])
        yy = np.asarray([sympify_noeval(y) for y in yy])

        if len(xx) != len(yy):
            raise ValueError(
                'Length of xx and yy must be the same '
                f'(instead len(xx) = {len(xx)} and len(yy) = {len(yy)})!'
            )

        if all(x.is_Number for x in xx) and not np.all(np.diff(xx) >= 0):
            raise ValueError('xx should be strictly increasing!')

        if (
                logarithmic_parametrization
                and all(y.is_Number for y in yy)
                and any(y <= 0 for y in yy)
        ):
            raise ValueError(
                'When interpolation is done in log-scale, '
                'yy should be strictly positive!'
            )

        bc, extrapolate = self._normalize_bc_and_extrapolate(bc, extrapolate)
        if bc == ('periodic', 'periodic') and yy[0] != yy[-1]:
            raise ValueError(
                'If the spline is to be periodic, '
                'the first and last elements of yy must be equal!'
            )

        self._sbml_id = sbml_id
        self._x = x
        self._xx = xx
        self._yy = yy
        self._bc = bc
        self._extrapolate = extrapolate
        self._logarithmic_parametrization = logarithmic_parametrization

        self._formula_cache = {}

    def _normalize_bc_and_extrapolate(self, bc, extrapolate):
        bc = AbstractSpline._normalize_bc(bc)
        return self._normalize_extrapolate(bc, extrapolate)

    @staticmethod
    def _normalize_bc(bc):
        """
        Preprocess the boundary condition `bc` to a standard form.
        """
        if not isinstance(bc, tuple):
            bc = (bc, bc)
        elif len(bc) != 2:
            raise TypeError(f'bc should be a 2-tuple, got {bc} instead!')

        bc = list(bc)

        valid_bc = (
            'periodic',
            'zeroderivative',
            'zeroderivative+natural',
            'natural',
            'no_bc',
            'auto',
            None
        )

        for i in (0, 1):
            if bc[i] not in valid_bc:
                raise ValueError(
                    f'Unsupported bc = {bc[i]}! '
                    f'The currently supported bc methods are: {valid_bc}'
                )
            elif bc[i] == 'no_bc':
                bc[i] = None

        if (bc[0] == 'periodic' or bc[1] == 'periodic') and bc[0] != bc[1]:
            raise ValueError(
                'If the bc on one side is periodic, '
                'then the bc on the other side must be periodic too!'
            )

        return tuple(bc)

    def _normalize_extrapolate(self, bc, extrapolate):
        """
        Preprocess `extrapolate` to a standard form
        and perform consistency checks
        """
        if not isinstance(extrapolate, tuple):
            extrapolate = (extrapolate, extrapolate)
        elif len(extrapolate) != 2:
            raise TypeError(
                f'extrapolate should be a 2-tuple, got {extrapolate} instead!'
            )
        extrapolate = list(extrapolate)

        if not isinstance(bc, tuple) or len(bc) != 2:
            raise TypeError(f'bc should be a 2-tuple, got {bc} instead!')
        bc = list(bc)

        valid_extrapolate = (
            'no_extrapolation',
            'constant',
            'linear',
            'polynomial',
            'periodic',
            None
        )

        for i in (0, 1):
            if extrapolate[i] not in valid_extrapolate:
                raise ValueError(
                    f'Unsupported extrapolate= {extrapolate[i]}!' +
                    ' The currently supported extrapolation methods are: ' +
                    str(valid_extrapolate)
                )

            if extrapolate[i] == 'no_extrapolation':
                extrapolate[i] = None

            if extrapolate[i] == 'periodic':
                if bc[0] == 'auto':
                    bc[0] = 'periodic'
                if bc[1] == 'auto':
                    bc[1] = 'periodic'
                if not (bc[0] == bc[1] == 'periodic'):
                    raise ValueError(
                        'The spline must satisfy periodic boundary conditions '
                        'on both sides of the base interval '
                        'in order for periodic extrapolation to be used!'
                    )

            elif extrapolate[i] == 'constant':
                assert self.smoothness > 0
                if self.smoothness == 1:
                    if bc[i] == 'auto':
                        bc[i] = 'zeroderivative'
                    elif bc[i] != 'zeroderivative':
                        raise ValueError(
                            'The spline must satisfy zero-derivative bc '
                            'in order for constant extrapolation to be used!'
                        )
                elif bc[i] == 'auto':
                    bc[i] = 'zeroderivative+natural'
                elif bc[i] != 'zeroderivative+natural':
                    raise ValueError(
                        'The spline must satisfy zero-derivative and natural'
                        ' bc in order for constant extrapolation to be used!'
                    )

            elif extrapolate[i] == 'linear':
                assert self.smoothness > 0
                if self.smoothness > 1:
                    if bc[i] == 'auto':
                        bc[i] = 'natural'
                    elif bc[i] != 'natural':
                        raise ValueError(
                            'The spline must satisfy natural bc '
                            'in order for linear extrapolation to be used!'
                        )
                elif bc[i] == 'auto':
                    bc[i] = None

            elif bc[i] == 'auto':
                bc[i] = None

        if (
                (extrapolate[0] == 'periodic' or extrapolate[1] == 'periodic')
                and extrapolate[0] != extrapolate[1]
                and extrapolate[0] is not None
                and extrapolate[1] is not None
        ):
            raise NotImplementedError(
                'At the moment if periodic extrapolation is applied '
                'to one side, the extrapolation at the other side '
                'must either be periodic or not be applied '
                '(in which case it will be periodic anyway).'
            )

        return tuple(bc), tuple(extrapolate)

    @property
    def sbml_id(self) -> sp.Symbol:
        """SBML ID of the spline parameter."""
        return self._sbml_id

    @property
    def x(self) -> sp.Basic:
        """The symbolic argument at which the spline is evaluated."""
        return self._x

    @property
    def xx(self) -> np.ndarray:
        """The points at which the spline values are known."""
        return self._xx

    @property
    def yy(self) -> np.ndarray:
        """The spline values at each of the points in ``xx``."""
        return self._yy

    @property
    def bc(self) -> Union[str, None]:
        """Boundary conditions applied to this spline."""
        return self._bc

    @property
    def extrapolate(self) -> Tuple[Union[None, str], Union[None, str]]:
        """Whether to extrapolate the spline outside the base interval."""
        return self._extrapolate

    @property
    def logarithmic_parametrization(self) -> bool:
        """Whether interpolation is done in log-scale."""
        return self._logarithmic_parametrization

    @property
    @abstractmethod
    def smoothness(self) -> str:
        """Smoothness of this spline."""
        raise NotImplementedError()

    @property
    @abstractmethod
    def method(self) -> str:
        """Spline method."""
        raise NotImplementedError()

    def check_if_valid(self, importer: SbmlImporter):
        """
        Check if the spline described by this object can be correctly
        be implemented by AMICI. E.g., check whether the formulas
        for spline grid points, values, ... contain species symbols.
        """
        # TODO this is very much a draft
        from .ode_export import SymbolId
        fixed_parameters = importer.symbols[SymbolId.FIXED_PARAMETER][
            'identifier']
        species = list(importer.symbols[SymbolId.SPECIES]['identifier'])

        for x in self.xx:
            if not x.free_symbols.issubset(fixed_parameters):
                raise ValueError(
                    'xx should only depend on constant parameters!'
                )

        for y in self.yy:
            if y.free_symbols.intersection(species):
                raise ValueError('yy should not depend on model species!')

        fixed_parameters_values = importer.symbols[SymbolId.FIXED_PARAMETER][
            'value']
        subs = dict(zip(fixed_parameters, fixed_parameters_values))
        xx_values = [sp.simplify(x.subs(subs)) for x in self.xx]
        for x in xx_values:
            assert x.is_Number
        if not np.all(np.diff(xx) >= 0):
            raise ValueError('xx should be strictly increasing!')

    def poly(self, i: int, *, x=None) -> sp.Basic:
        """
        Get the polynomial interpolant on the ``(xx[i], xx[i+1])`` interval.
        The polynomial is written in Horner form with respect to the scaled
        variable ``poly_variable(x, i)``.
        If no variable ``x`` is provided, it will default to the one given at
        initialization time.
        """
        if i < 0:
            i += len(self.xx) - 1

        if not 0 <= i < len(self.xx) - 1:
            raise ValueError(f'Interval index {i} is out of bounds!')

        if x is None:
            x = self.x

        # Compute polynomial in Horner form for the scaled variable
        t = sp.Dummy('t')
        poly = sp.poly(self._poly(t, i), wrt=t)
        poly = sp.horner(poly)

        # Replace scaled variable with its value,
        # without changing the expression form
        t_value = self._poly_variable(x, i)
        with evaluate(False):
            return poly.subs(t, t_value)

    def poly_variable(self, x, i) -> sp.Basic:
        """
        Given an evaluation point, return the value of the variable
        in which the polynomial on the ``i``-th interval is expressed.
        """
        if not 0 <= i < len(self.xx) - 1:
            raise ValueError(f'Interval index {i} is out of bounds!')
        return self._poly_variable(x, i)

    @abstractmethod
    def _poly_variable(self, x, i) -> sp.Basic:
        # This function (and not poly_variable)
        # should be implemented by the subclasses
        raise NotImplementedError()

    @abstractmethod
    def _poly(self, t, i) -> sp.Basic:
        """
        Return the symbolic expression for the spline restricted to the `i`-th
        interval as a polynomial in the scaled variable `t`.
        """
        raise NotImplementedError()

    def segment_formula(self, i: int, *, x=None) -> sp.Basic:
        """
        Return the formula for the actual value of the spline expression
        on the ``(xx[i], xx[i+1])`` interval.
        Unless logarithmic parametrization is used,
        this is equal to the interpolating polynomial.
        """
        if x is None:
            x = self.x
        poly = self.poly(i, x=x)
        if self.logarithmic_parametrization:
            return sp.exp(poly)
        return poly

    def y_scaled(self, i: Integral):
        """
        Return the values which should be interpolated by a polynomial.
        Unless logarithmic parametrization is used,
        they are equal to the values given at initialization time.
        """
        if self.logarithmic_parametrization:
            return sp.log(self.yy[i])
        return self.yy[i]

    @property
    def extrapolation_formulas(self) \
            -> Tuple[Union[None, sp.Basic], Union[None, sp.Basic]]:
        """
        Returns the extrapolation formulas on the left and right side
        of the interval ``(xx[0], xx[-1])``.
        A value of ``None`` means that no extrapolation is required.
        """
        return self._extrapolation_formulas(self.x)

    def _extrapolation_formulas(self, x, extrapolate=None):
        if extrapolate is None:
            extrLeft, extrRight = self.extrapolate
        else:
            extrLeft, extrRight = extrapolate

        if extrLeft == 'constant':
            extrLeft = self.yy[0]
        elif extrLeft == 'linear':
            dx = x - self.xx[0]
            dydx = self.derivative(self.xx[0], extrapolate=None)
            extrLeft = self.yy[0] + dydx * dx
        elif extrLeft == 'polynomial':
            extrLeft = None
        else:
            assert extrLeft is None

        if extrRight == 'constant':
            extrRight = self.yy[-1]
        elif extrRight == 'linear':
            dx = x - self.xx[-1]
            dydx = self.derivative(self.xx[-1], extrapolate=None)
            extrRight = self.yy[-1] + dydx * dx
        elif extrRight == 'polynomial':
            extrRight = None
        else:
            assert extrRight is None

        return extrLeft, extrRight

    @property
    def formula(self) -> sp.Piecewise:
        """
        Compute a symbolic piecewise formula for the spline.
        """
        return self._formula(sbml_syms=False, sbml_ops=False)

    @property
    def sbmlFormula(self) -> sp.Piecewise:
        """
        Compute a symbolic piecewise formula for the spline,
        using SBML symbol naming
        (the AMICI time symbol will be replaced with its SBML counterpart).
        """
        return self._formula(sbml_syms=True, sbml_ops=False)

    @property
    def mathmlFormula(self) -> sp.Piecewise:
        """
        Compute a symbolic piecewise formula for the spline for use inside
        a SBML assignment rule: SBML symbol naming will be used
        and operations not supported by SBML MathML will be avoided.
        """
        return self._formula(sbml_syms=True, sbml_ops=True)

    def _formula(
            self,
            *,
            x=None,
            sbml_syms=False,
            sbml_ops=False,
            cache=True,
            **kwargs
    ) -> sp.Piecewise:

        # Cache formulas in the case they are reused
        if cache:
            if 'extrapolate' in kwargs.keys():
                key = (x, sbml_syms, sbml_ops, kwargs['extrapolate'])
            else:
                key = (x, sbml_syms, sbml_ops)
            if key in self._formula_cache.keys():
                return self._formula_cache[key]

        if x is None:
            x = self.x
        if 'extrapolate' in kwargs.keys():
            _bc, extrapolate = self._normalize_extrapolate(
                self.bc, kwargs['extrapolate']
            )
            assert self.bc == _bc
        else:
            extrapolate = self.extrapolate

        pieces = []

        if extrapolate[0] == 'periodic' or extrapolate[1] == 'periodic':
            if sbml_ops:
                # NB mod is not supported in SBML
                x = sp.Symbol(self.sbml_id.name + '_x_in_fundamental_period')
                # NB we will do the parameter substitution in SBML
                #    because the formula for x will be a piecewise
                #    and sympy handles Piecewises inside other Piecewises
                #    really badly.
            else:
                x = self._to_base_interval(x)
            extrLeft, extrRight = None, None
        else:
            extrLeft, extrRight = self._extrapolation_formulas(x, extrapolate)

        if extrLeft is not None:
            pieces.append((extrLeft, x < self.xx[0]))

        for i in range(len(self.xx) - 2):
            pieces.append((self.segment_formula(i, x=x), x < self.xx[i + 1]))

        if extrRight is not None:
            pieces.append((self.segment_formula(-1, x=x), x < self.xx[-1]))
            pieces.append((extrRight, sp.sympify(True)))
        else:
            pieces.append((self.segment_formula(-1, x=x), sp.sympify(True)))

        with evaluate(False):
            if sbml_syms:
                pieces = [
                    (
                        p.subs(amici_time_symbol, sbml_time_symbol),
                        c.subs(amici_time_symbol, sbml_time_symbol)
                    )
                    for (p, c) in pieces
                ]
            formula = sp.Piecewise(*pieces)

        if cache:
            self._formula_cache[key] = formula
        return formula

    @property
    def period(self):
        """Period of a periodic spline. `None` if the spline is not periodic.
        """
        if self.bc == ('periodic', 'periodic'):
            return self.xx[-1] - self.xx[0]
        return None

    def _to_base_interval(self, x, *, with_interval_number: bool = False):
        """For periodic splines, maps the real point `x` to the reference
        period."""
        if self.bc != ('periodic', 'periodic'):
            raise ValueError(
                '_to_base_interval makes no sense with non-periodic bc'
            )

        xA = self.xx[0]
        xB = self.xx[-1]
        T = self.period
        z = xA + sp.Mod(x - xA, T)
        assert not z.is_Number or xA <= z < xB

        if with_interval_number:
            k = sp.floor((x - xA) / T)
            assert x == z + k * T
            return k, z
        return z

    def evaluate(self, x):
        """Evaluate the spline at the point `x`."""
        _x = sp.Dummy('x')
        return self._formula(x=_x, cache=False).subs(_x, x)

    def derivative(self, x, **kwargs):
        """Evaluate the spline derivative at the point `x`."""
        # NB kwargs are used to pass on extrapolate=None
        #    when called from .extrapolation_formulas()
        _x = sp.Dummy('x')
        return self._formula(x=_x, cache=False, **kwargs).diff(_x).subs(_x, x)

    def second_derivative(self, x):
        """Evaluate the spline second derivative at the point `x`."""
        _x = sp.Dummy('x')
        return self._formula(x=_x, cache=False).diff(_x).diff(_x).subs(_x, x)

    def integrate(self, x0, x1):
        """Integrate the spline between the points `x0` and `x1`."""
        x = sp.Dummy('x')
        x0, x1 = sp.sympify((x0, x1))

        if x0 > x1:
            raise ValueError('x0 > x1')

        if x0 == x1:
            return sp.sympify(0)

        if self.extrapolate != ('periodic', 'periodic'):
            return self._formula(x=x, cache=False).integrate((x, x0, x1))

        formula = self._formula(x=x, cache=False, extrapolate=None)

        T = self.period
        xA, xB = self.xx[0], self.xx[-1]
        k0, z0 = self._to_base_interval(x0, with_interval_number=True)
        k1, z1 = self._to_base_interval(x1, with_interval_number=True)

        assert k0 <= k1

        if k0 == k1:
            return formula.integrate((x, z0, z1))

        if k0 + 1 == k1:
            return formula.integrate((x, z0, xB)) + \
                   formula.integrate((x, xA, z1))

        return formula.integrate((x, z0, xB)) + \
               (k1 - k0 - 1) * formula.integrate((x, xA, xB)) + \
               formula.integrate((x, xA, z1))

    @property
    def amiciAnnotation(self) -> str:
        """
        An SBML annotation describing the spline.
        """
        annotation = f'<amici:spline xmlns:amici="{annotation_namespace}"'
        for (attr, value) in self._annotation_attributes().items():
            if isinstance(value, bool):
                value = str(value).lower()
            value = f'"{value}"'
            annotation += f' amici:{attr}={value}'
        annotation += '>'

        for (child, grandchildren) in self._annotation_children().items():
            if isinstance(grandchildren, str):
                grandchildren = [grandchildren]
            annotation += f'<amici:{child}>'
            for gc in grandchildren:
                annotation += gc
            annotation += f'</amici:{child}>'

        annotation += '</amici:spline>'

        # Check XML and prettify
        return pretty_xml(annotation)

    def _annotation_attributes(self) -> Dict[str, Any]:
        attributes = {'spline_method': self.method}

        if self.bc[0] == self.bc[1]:
            if self.bc[0] is not None:
                attributes['spline_bc'] = self.bc[0]
        else:
            bc1, bc2 = self.bc
            bc1 = 'no_bc' if bc1 is None else bc1
            bc2 = 'no_bc' if bc2 is None else bc2
            attributes['spline_bc'] = f'({bc1}, {bc2})'

        if self.extrapolate[0] == self.extrapolate[1]:
            extr = None if self.extrapolate is None else self.extrapolate[0]
        else:
            extr1, extr2 = self.extrapolate
            extr1 = 'no_extrapolation' if extr1 is None else extr1
            extr2 = 'no_extrapolation' if extr2 is None else extr2
            extr = f'({extr1}, {extr2})'
        if extr is not None:
            attributes['spline_extrapolate'] = extr

        if self.logarithmic_parametrization:
            attributes['spline_logarithmic_parametrization'] = True

        return attributes

    def _annotation_children(self) -> Dict[str, Union[str, List[str]]]:
        children = {}

        with evaluate(False):
            x = self.x.subs(amici_time_symbol, sbml_time_symbol)
        children['spline_evaluation_point'] = sbmlMathML(x)

        if isinstance(self.xx, UniformGrid):
            children['spline_uniform_grid'] = [
                sbmlMathML(self.xx.start),
                sbmlMathML(self.xx.stop),
                sbmlMathML(self.xx.step),
            ]
        else:
            for x in self.xx:
                assert amici_time_symbol not in x.free_symbols
            children['spline_grid'] = [sbmlMathML(x) for x in self.xx]

        children['spline_values'] = [sbmlMathML(y) for y in self.yy]

        return children

    def addToSbmlModel(
            self,
            model: libsbml.Model,
            *,
            auto_add: Union[bool, str] = False,
            x_nominal: Sequence[float] = None,
            y_nominal: Optional[Union[Sequence[float], float]] = None,
            y_units: Optional[str] = None,
            x_units: Optional[str] = None
    ):
        """
        Function add the spline to an SBML model using an assignment rule
        with AMICI-specific annotations.

        :param model:
            A :py:class:`libsbml.Model` to which the spline is to be added.

        :param auto_add:
            Automatically add missing parameters to the SBML model
            (defaults to `False`).
            Only used for expressions consisting in a single symbol.
            If equal to `'spline'`,
            only the parameter representing the spline will be added.

        :param x_nominal:
            Nominal values used when auto-adding parameters for `xx`.

        :param y_nominal:
            Nominal values used when auto-adding parameters for `yy`.

        :param x_units:
            Units used when auto-adding parameters for `xx`.

        :param y_units:
            Units used when auto-adding parameters for `yy`.
        """
        # Convert time from AMICI to SBML naming
        with evaluate(False):
            x = self.x.subs(amici_time_symbol, sbml_time_symbol)

        # Try to autodetermine units
        if x_units is None:
            x_units = getSbmlUnits(model, x)
            for _x in self.xx:
                if x_units is not None:
                    break
                x_units = getSbmlUnits(model, _x)
        if y_units is None:
            for _y in self.yy:
                y_units = getSbmlUnits(model, _y)
                if y_units is not None:
                    break

        # Autoadd parameters
        if auto_add is True or auto_add == 'spline':
            if not hasParameter(model, self.sbml_id):
                addParameter(model, self.sbml_id, constant=False,
                             units=y_units)

            if auto_add is True:
                if isinstance(x_nominal, collections.abc.Sequence):
                    if len(x_nominal) != len(self.xx):
                        raise ValueError(
                            'If x_nominal is a list, then it must have '
                            'the same length as the spline grid!'
                        )
                else:
                    x_nominal = len(self.xx) * [x_nominal]
                for (_x, _val) in zip(self.xx, x_nominal):
                    if _x.is_Symbol and not hasParameter(model, _x.name):
                        addParameter(model, _x.name, value=_val, units=x_units)

                if isinstance(y_nominal, collections.abc.Sequence):
                    if len(y_nominal) != len(self.yy):
                        raise ValueError(
                            'If y_nominal is a list, then it must have '
                            'the same length as the spline values!'
                        )
                else:
                    y_nominal = len(self.yy) * [y_nominal]
                for (_y, _val) in zip(self.yy, y_nominal):
                    if _y.is_Symbol and not hasParameter(model, _y.name):
                        addParameter(model, _y.name, value=_val, units=y_units)

        elif auto_add is not False:
            raise ValueError(f'Invalid value {auto_add} for auto_add!')

        # Create assignment rule for spline
        rule = addAssignmentRule(model, self.sbml_id, self.mathmlFormula)

        # Add annotation specifying spline method
        retcode = rule.setAnnotation(self.amiciAnnotation)
        if retcode != libsbml.LIBSBML_OPERATION_SUCCESS:
            raise SbmlException('Could not set SBML annotation!')

        # Create additional assignment rule for periodic extrapolation
        # NB mod is not in the subset of MathML supported by SBML
        if any(extr == 'periodic' for extr in self.extrapolate):
            parameterId = self.sbml_id.name + '_x_in_fundamental_period'
            T = self.xx[-1] - self.xx[0]
            x0 = self.xx[0]
            s = 2 * sp.pi * ((x - x0) / T - sp.sympify(1) / 4)
            k = sp.Piecewise((3, sp.cos(s) < 0), (1, True))
            formula = x0 + T * (sp.atan(sp.tan(s)) / (2 * sp.pi) + k / 4)
            assert amici_time_symbol not in formula.free_symbols
            par = addParameter(
                model, parameterId, constant=False, units=x_units
            )
            retcode = par.setAnnotation(
                f'<amici:discard xmlns:amici="{annotation_namespace}" />'
            )
            if retcode != libsbml.LIBSBML_OPERATION_SUCCESS:
                raise SbmlException('Could not set SBML annotation!')
            addAssignmentRule(model, parameterId, formula)

    # def _replace_sbml_time_with_amici_time(self):
    #     self._replace_in_all_expressions(
    #         sbml_time_symbol, amici_time_symbol
    #     )

    def _replace_in_all_expressions(self, old: sp.Symbol, new: sp.Symbol):
        if self.sbml_id == old:
            self._sbml_id = new
        self._x = self.x.subs(old, new)
        if not isinstance(self.xx, UniformGrid):
            self._xx = [x.subs(old, new) for x in self.xx]
        self._yy = [y.subs(old, new) for y in self.yy]

    @staticmethod
    def isSpline(rule: libsbml.AssignmentRule):
        """
        Determine if an SBML assignment rule (given as a
        :py:class:`libsbml.AssignmentRule` object) is an AMICI-annotated
        spline formula.
        """
        return AbstractSpline.getAnnotation(rule) is not None

    @staticmethod
    def getAnnotation(rule: libsbml.AssignmentRule):
        """
        Extract AMICI spline annotation from an SBML assignment rule
        (given as a :py:class:`libsbml.AssignmentRule` object).
        Return ``None`` if any such annotation could not be found.
        """
        if not isinstance(rule, libsbml.AssignmentRule):
            raise TypeError('Rule must be an AssignmentRule!')
        if rule.isSetAnnotation():
            annotation = ET.fromstring(rule.getAnnotationString())
            for child in annotation:
                if child.tag == f'{{{annotation_namespace}}}spline':
                    return child
        return None

    @staticmethod
    def fromAnnotation(sbmlId: sp.Symbol, annotation: ET.Element, *, locals):
        """Create a spline object from a SBML annotation.

        This function extracts annotation and children from the XML annotation
        and gives them to the ``_fromAnnotation`` function for parsing.
        Subclass behaviour should be implemented by extending
        ``_fromAnnotation``.
        However, the mapping between method strings and subclasses
        must be hard-coded into this function here (at the moment).
        """
        if annotation.tag != f'{{{annotation_namespace}}}spline':
            raise ValueError(
                'The given annotation is not an AMICI spline annotation!'
            )

        attributes = {}
        for key, value in annotation.items():
            if not key.startswith(f'{{{annotation_namespace}}}'):
                raise ValueError(
                    f'Unexpected attribute {key} inside spline annotation!'
                )
            key = key[len(annotation_namespace) + 2:]
            if value == 'true':
                value = True
            elif value == 'false':
                value = False
            attributes[key] = value

        children = {}
        for child in annotation:
            if not child.tag.startswith(f'{{{annotation_namespace}}}'):
                raise ValueError(
                    f'Unexpected node {child.tag} inside spline annotation!'
                )
            key = child.tag[len(annotation_namespace) + 2:]
            value = [
                mathml2sympy(
                    ET.tostring(gc).decode(),
                    evaluate=False,
                    locals=locals,
                    expression_type='Rule'
                )
                for gc in child
            ]
            children[key] = value

        if attributes['spline_method'] == 'cubic_hermite':
            cls = CubicHermiteSpline
        else:
            raise ValueError(
                f"Unknown spline method {attributes['spline_method']}!"
            )

        del attributes['spline_method']
        kwargs = cls._fromAnnotation(attributes, children)

        if len(attributes) != 0:
            raise ValueError(
                'Unprocessed attributes in spline annotation!\n' +
                str(attributes)
            )

        if len(children) != 0:
            raise ValueError(
                'Unprocessed children in spline annotation!\n' +
                str(children)
            )

        return cls(sbmlId, **kwargs)

    @classmethod
    def _fromAnnotation(cls, attributes, children):
        """
        Given the attributes and children of a AMICI spline annotation,
        returns the keyword arguments to be passed
        to the spline object ``__init__`` function.
        """
        kwargs = {}

        bc = attributes.pop('spline_bc', None)
        if isinstance(bc, str) and bc.startswith('('):
            if not bc.endswith(')'):
                raise ValueError(f'Ill-formatted bc {bc}!')
            bc_cmps = bc[1:-1].split(',')
            if len(bc_cmps) != 2:
                raise ValueError(f'Ill-formatted bc {bc}!')
            bc = (bc_cmps[0].strip(), bc_cmps[1].strip())
        kwargs['bc'] = bc

        extr = attributes.pop('spline_extrapolate', None)
        if isinstance(extr, str) and extr.startswith('('):
            if not extr.endswith(')'):
                raise ValueError(f'Ill-formatted extrapolation {extr}!')
            extr_cmps = extr[1:-1].split(',')
            if len(extr_cmps) != 2:
                raise ValueError(f'Ill-formatted extrapolation {extr}!')
            extr = (extr_cmps[0].strip(), extr_cmps[1].strip())
        kwargs['extrapolate'] = extr

        kwargs['logarithmic_parametrization'] = \
            attributes.pop('spline_logarithmic_parametrization', False)

        if 'spline_evaluation_point' not in children.keys():
            raise ValueError(
                "Required spline annotation 'spline_evaluation_point' missing!"
            )
        x = children.pop('spline_evaluation_point')
        if len(x) != 1:
            raise ValueError(
                "Ill-formatted spline annotation 'spline_evaluation_point' " +
                "(more than one children is present)!"
            )
        kwargs['x'] = x[0]

        if 'spline_uniform_grid' in children.keys():
            start, stop, step = children.pop('spline_uniform_grid')
            kwargs['xx'] = UniformGrid(start, stop, step)
        elif 'spline_grid' in children.keys():
            kwargs['xx'] = children.pop('spline_grid')
        else:
            raise ValueError(
                "Spline annotation requires either "
                "'spline_grid' or 'spline_uniform_grid' to be specified!"
            )

        if 'spline_values' not in children.keys():
            raise ValueError(
                "Required spline annotation 'spline_values' missing!"
            )
        kwargs['yy'] = children.pop('spline_values')

        return kwargs

    def parameters(self, importer: SbmlImporter):
        """Returns the SBML parameters used by this spline"""
        from .ode_export import SymbolId
        return self._parameters().intersection(
            set(importer.symbols[SymbolId.PARAMETER].keys())
        )

    def _parameters(self):
        parameters = set()
        for y in self.yy:
            parameters.update(y.free_symbols)
        return parameters

    def odeModelSymbol(self, importer: SbmlImporter):
        """
        Returns the `sympy` object to be used by
        :py:class:`amici.ode_export.ODEModel`.
        This expression can be differentiated and easily mapped to the C++
        code.
        """
        parameters = list(self.parameters(importer))

        class AmiciSpline(sp.Function):
            # AmiciSpline(splineId, x, *parameters)
            nargs = (len(parameters) + 2,)

            @classmethod
            def eval(cls, *args):
                return None  # means leave unevaluated

            def fdiff(self, argindex=1):
                if argindex == 1:
                    # Derivative with respect to the spline SBML ID
                    # Since the SBML ID is not a real function parameter
                    # (more like a subscript), the derivative will be zero
                    return sp.Integer(0)

                if argindex == 2:
                    class AmiciSplineDerivative(sp.Function):
                        # Spline derivative
                        # AmiciSplineDerivative(splineId, x, *parameters)
                        nargs = (len(parameters) + 2,)

                        @classmethod
                        def eval(cls, *args):
                            return None  # means leave unevaluated

                        def fdiff(self, argindex=1):
                            return NotImplementedError(
                                'Second order derivatives for spline '
                                'are not implemented yet.'
                            )

                        def _eval_is_real(self):
                            return True

                    return AmiciSplineDerivative(*self.args)

                pindex = argindex - 3
                assert 0 <= pindex < len(parameters)

                class AmiciSplineSensitivity(sp.Function):
                    # Derivative with respect to a parameter paramId
                    # AmiciSplineSensitivity(splineId, x, paramId, *parameters)
                    nargs = (len(parameters) + 3,)

                    @classmethod
                    def eval(cls, *args):
                        return None  # means leave unevaluated

                    def fdiff(self, argindex=1):
                        return NotImplementedError(
                            'Second order derivatives for spline '
                            'are not implemented yet.'
                        )

                    def _eval_is_real(self):
                        return True

                return AmiciSplineSensitivity(
                    self.args[0],
                    self.args[1],
                    parameters[pindex],
                    *self.args[2:]
                )

            def _eval_is_real(self):
                return True

        return AmiciSpline(self.sbml_id, self.x, *parameters)


def spline_user_functions(
        splines: List[AbstractSpline],
        p_index: Dict[sp.Symbol, int]
):
    """
    Custom user functions to be used in `ODEExporter`
    for linking spline expressions to C++ code.
    """
    splineIds = [spline.sbml_id.name for spline in splines]
    return {
        'AmiciSpline': [(
            lambda *args: True,
            lambda splineId, x, *p: f"spl_{splineIds.index(splineId)}"
        )],
        'AmiciSplineDerivative': [(
            lambda *args: True,
            lambda splineId, x, *p: f"dspl_{splineIds.index(splineId)}"
        )],
        'AmiciSplineSensitivity': [(
            lambda *args: True,
            lambda splineId, x, paramId, *p:
            f"sspl_{splineIds.index(splineId)}_{p_index[paramId]}"
        )],
    }


class CubicHermiteSpline(AbstractSpline):
    def __init__(
            self,
            sbml_id: Union[str, sp.Symbol],
            x: Union[str, sp.Basic],
            xx: Sequence,
            yy: Sequence,
            dd: Sequence = None,
            *,
            bc: BClike = 'auto',
            extrapolate: BClike = None,
            logarithmic_parametrization: bool = False
    ):
        """
        Constructor for `CubicHermiteSpline` objects.

        :param sbml_id:
            The SBML ID of the parameter associated to the spline
            as a string or a SymPy symbol.

        :param x:
            The point at which the spline is evaluated.
            It will be sympified.

        :param xx:
            The points at which the spline values are known.
            Currently they can only depend on constant parameters.
            These points should be strictly increasing.
            This argument will be sympified.

        :param yy:
            The spline values at each of the points in `xx`.
            They may not depend on model species.
            This argument will be sympified.

        :param dd:
            The spline derivatives at each of the points in `xx`.
            They may not depend on model species.
            This argument will be sympified.
            If not specified, it will be computed by finite differences.

        :param bc:
            Applied boundary conditions
            (see `AbstractSpline` documentation).
            If `'auto'` (the default), the boundary conditions will be
            automatically set depending on the extrapolation methods.

        :param extrapolate:
            Extrapolation method (see `AbstractSpline` documentation).

        :param logarithmic_parametrization:
            Whether interpolation should be done in log-scale.
        """

        if not isinstance(xx, UniformGrid):
            xx = np.asarray([sympify_noeval(x) for x in xx])
        yy = np.asarray([sympify_noeval(y) for y in yy])

        if len(xx) != len(yy):
            # NB this would be checked in AbstractSpline.__init__()
            #    however, we check it now so that an informative message
            #    can be printed (otherwise finite difference computation fails)
            raise ValueError(
                'Length of xx and yy must be the same '
                f'(instead len(xx) = {len(xx)} and len(yy) = {len(yy)})!'
            )

        bc, extrapolate = self._normalize_bc_and_extrapolate(bc, extrapolate)
        if bc[0] == 'zeroderivative+natural' \
                or bc[1] == 'zeroderivative+natural':
            raise ValueError("zeroderivative+natural bc not supported by "
                             "CubicHermiteSplines!")

        if dd is None:
            dd = finite_differences(xx, yy, bc)
            self._derivatives_by_fd = True
        else:
            dd = np.asarray([sympify_noeval(d) for d in dd])
            self._derivatives_by_fd = False

        if len(xx) != len(dd):
            raise ValueError(
                'Length of xx and dd must be the same '
                f'(instead len(xx) = {len(xx)} and len(dd) = {len(dd)})!'
            )

        if bc == ('periodic', 'periodic') \
                and (yy[0] != yy[-1] or dd[0] != dd[-1]):
            raise ValueError(
                'bc=periodic but given yy and dd do not satisfy '
                'periodic boundary conditions!'
            )

        super().__init__(
            sbml_id, x, xx, yy,
            bc=bc,
            extrapolate=extrapolate,
            logarithmic_parametrization=logarithmic_parametrization
        )

        self._dd = dd

    @property
    def dd(self) -> np.ndarray:
        """The spline derivatives at each of the points in `xx`."""
        return self._dd

    @property
    def smoothness(self) -> str:
        """
        Smoothness of this spline (equal to 1 for cubic Hermite splines
        since they are continuous up to the first derivative).
        """
        return 1

    @property
    def method(self) -> str:
        """Spline method (cubic Hermite spline)"""
        return 'cubic_hermite'

    @property
    def derivatives_by_fd(self):
        return self._derivatives_by_fd

    def check_if_valid(self, importer: SbmlImporter):
        """
        Check if the spline described by this object can be correctly
        be implemented by AMICI. E.g., check whether the formulas
        for spline grid points, values, ... contain species symbols.
        """
        # TODO this is very much a draft
        from .ode_export import SymbolId
        species = list(importer.symbols[SymbolId.SPECIES]['identifier'])
        for d in self.dd:
            if len(d.free_symbols.intersection(species)) != 0:
                raise ValueError('dd should not depend on model species')

        super().check_if_valid(importer)

    def d_scaled(self, i: int):
        """
        Return the derivative of the polynomial interpolant at the `i`-th
        point. Unless logarithmic parametrization is used, it is equal to the
        derivative of the spline expression.
        """
        if self.logarithmic_parametrization:
            return self.dd[i] / self.yy[i]
        return self.dd[i]

    def _poly_variable(self, x, i) -> sp.Basic:
        assert 0 <= i < len(self.xx) - 1
        dx = self.xx[i + 1] - self.xx[i]
        with evaluate(False):
            return (x - self.xx[i]) / dx

    def _poly(self, t, i) -> sp.Basic:
        """
        Return the symbolic expression for the spline restricted to the `i`-th
        interval as polynomial in the scaled variable `t`.
        """
        assert 0 <= i < len(self.xx) - 1

        dx = self.xx[i + 1] - self.xx[i]

        h00 = 2 * t ** 3 - 3 * t ** 2 + 1
        h10 = t ** 3 - 2 * t ** 2 + t
        h01 = -2 * t ** 3 + 3 * t ** 2
        h11 = t ** 3 - t ** 2

        y0 = self.y_scaled(i)
        y1 = self.y_scaled(i + 1)
        dy0 = self.d_scaled(i)
        dy1 = self.d_scaled(i + 1)

        with evaluate(False):
            return h00 * y0 + h10 * dx * dy0 + h01 * y1 + h11 * dx * dy1

    def _annotation_children(self) -> Dict[str, Union[str, List[str]]]:
        children = super()._annotation_children()
        if not self._derivatives_by_fd:
            children['spline_derivatives'] = [sbmlMathML(d) for d in self.dd]
        return children

    def _parameters(self):
        parameters = super()._parameters()
        for d in self.dd:
            parameters.update(d.free_symbols)
        return parameters

    def _replace_in_all_expressions(self, old: sp.Symbol, new: sp.Symbol):
        super()._replace_in_all_expressions(old, new)
        self._dd = [d.subs(old, new) for d in self.dd]

    @classmethod
    def _fromAnnotation(cls, attributes, children):
        kwargs = super()._fromAnnotation(attributes, children)

        if 'spline_derivatives' in children.keys():
            kwargs['dd'] = children.pop('spline_derivatives')

        return kwargs

    def __str__(self):
        s = 'HermiteCubicSpline ' + \
            f'on ({self.xx[0]}, {self.xx[-1]}) with {len(self.xx)} points'
        cmps = []
        if self.bc != (None, None):
            if self.bc == ('periodic', 'periodic'):
                cmps.append('periodic')
            else:
                cmps.append(f'bc = {self.bc}')
        if self.derivatives_by_fd:
            cmps.append('finite differences')
        if self.extrapolate != (None, None):
            cmps.append(f'extrapolate = {self.extrapolate}')
        if not cmps:
            return s
        return s + ' [' + ', '.join(cmps) + ']'


def finite_differences(xx, yy, bc):
    dd = []

    if bc[0] == 'periodic':
        fd = centeredFD(yy[-2], yy[0], yy[1], xx[-1] - xx[-2], xx[1] - xx[0])
    elif bc[0] == 'zeroderivative':
        fd = sp.Integer(0)
    elif bc[0] == 'natural':
        if len(xx) < 3:
            raise ValueError(
                'At least 3 nodes are needed '
                'for computing finite differences with natural bc!'
            )
        fd = naturalFD(yy[0], xx[1] - xx[0], yy[1], xx[2] - xx[1], yy[2])
    else:
        fd = onesidedFD(yy[0], yy[1], xx[1] - xx[0])
    dd.append(fd)

    for i in range(1, len(xx) - 1):
        dd.append(
            centeredFD(yy[i - 1], yy[i], yy[i + 1], xx[i] - xx[i - 1],
                       xx[i + 1] - xx[i])
        )

    if bc[1] == 'periodic':
        fd = dd[0]
    elif bc[1] == 'zeroderivative':
        fd = sp.Integer(0)
    elif bc[1] == 'natural':
        if len(xx) < 3:
            raise ValueError(
                'At least 3 nodes are needed '
                'for computing finite differences with natural bc!'
            )
        fd = naturalFD(yy[-1], xx[-2] - xx[-1], yy[-2], xx[-3] - xx[-2],
                       yy[-3])
    else:
        fd = onesidedFD(yy[-2], yy[-1], xx[-1] - xx[-2])
    dd.append(fd)

    return np.asarray(dd)


def onesidedFD(y0, y1, h):
    return sp.Mul(1 / h, y1 - y0, evaluate=False)


def centeredFD(ym1, y0, yp1, hm, hp):
    if hm == hp:
        return sp.Mul(1 / (2 * hm), yp1 - ym1, evaluate=False)
    else:
        return ((yp1 - y0) / hp + (y0 - ym1) / hm) / 2


def naturalFD(y0, dx1, y1, dx2, y2):
    if dx1 == dx2:
        den = 4 * dx1
        with evaluate(False):
            return (6 * y1 - 5 * y0 - y2) / den
    else:
        with evaluate(False):
            return ((y1 - y2) / dx2 - 5 * (y0 - y1) / dx1) / 4
    # Another formula, which depends on
    # y0, dx1 = x1 - x0, y1 and dy1 (derivative at x1)
    # (-dx1*dy1 - 3*y0 + 3*y1)/(2*dx1)
