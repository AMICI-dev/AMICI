"""
Splines
------------
This module provides helper functions for reading/writing splines
with AMICI annotations from/to SBML files
and for adding such splines to the AMICI C++ code.
"""

import collections

import xml.dom.minidom
import numpy as np
import sympy as sp
import libsbml

from abc import abstractmethod, ABC
from itertools import count
from typing import Sequence, Union, Tuple, Optional, Dict, Any, List
from sympy.core.parameters import evaluate

from .sbml_import import SbmlImporter
from .sbml_utils import (
    sbml_time_symbol,
    amici_time_symbol,
    sbmlMathML
)


def sympify_noeval(x):
    with evaluate(False):
        return sp.sympify(x)


################################################################################


class UniformGrid(collections.abc.Sequence):
    """
    A grid of uniformly-spaced real points, computed with rational arithmetic.
    Implements the `collections.abc.Sequence` interface
    and can be converted to a `numpy.ndarray`
    (conversion to float can be specified with `dtype=float`).
    """

    def __init__(self, start, stop, step):
        """
        Create a new `UniformGrid`.
        The `stop` attribute of the resulting object will be larger
        than the given `stop` argument in the case `stop - start` is not
        an integer multiple of `step`.
        """
        start = sp.nsimplify(sp.sympify(start))
        stop = sp.nsimplify(sp.sympify(stop))
        step = sp.nsimplify(sp.sympify(step))

        if start > stop:
            raise ValueError(
                f'start point {start} greater than stop point {stop}'
            )

        if step <= 0:
            raise ValueError(f'step size {step} must be strictly positive')

        xx = []
        for i in count():
            x = start + i * step
            xx.append(x)
            if x >= stop:
                break

        self._xx = np.asarray(xx)

    @property
    def start(self):
        "First point."
        return self._xx[0]

    @property
    def stop(self):
        "Last point."
        return self._xx[-1]

    @property
    def step(self):
        "Distance between consecutive points."
        if len(self._xx) > 1:
            return self._xx[1] - self._xx[0]
        else:
            return sp.sympify(0)

    def __getitem__(self, i):
        return self._xx[i]

    def __len__(self):
        return len(self._xx)

    def __array__(self, dtype=None):
        if dtype is None:
            return self._xx
        else:
            return np.array(self._xx, dtype=dtype)

    def __repr__(self):
        return f'UniformGrid(start={self.start}, stop={self.stop}, step={self.step})'


################################################################################


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
            sbmlId: Union[str, sp.Symbol],
            x: Union[str, sp.Basic],
            xx: Sequence,
            yy: Sequence,
            *,
            periodic: bool = False,
            extrapolate: Union[
                None, str,
                Tuple[Union[None, str], Union[None, str]]
            ] = None
        ):
        """
        Base constructor for `AbstractSpline` objects.

            :param sbmlId:
            The SBML ID of the parameter associated to the spline as a string or
            a SymPy symbol.

            :param x:
            The symbolic argument at which the spline is evaluated.
            It will be sympified.

            :param xx:
            The points at which the spline values are known.
            Currently they can only contain depend on constant parameters.
            It should be strictly increasing.
            It will be sympified.

            :param yy:
            The spline values at each of the points in `xx`.
            They may not depend on model species.
            It will be sympified.

            :param periodic:
            Whether the spline satisfies periodic boundary conditions.

            :param extrapolate:
            Whether to extrapolate the spline outside the base interval
            defined by `(xx[0], xx[-1])`.
            It is a tuple of extrapolation methods, one for each side of the
            base interval.
            `None` or `'no_extrapolation'`
            means no extrapolation is to be performed
            (evaluating the spline outside the base interval is not guaranteed
            to work).
            If it is not a tuple, then the same extrapolation will be applied
            on both sides.
        """

        if isinstance(sbmlId, str):
            sbmlId = sp.Symbol(sbmlId)
        elif not isinstance(sbmlId, sp.Symbol):
            raise TypeError(
                'sbmlId must be either a string or a SymPy symbol, '
                f'got {sbmlId} instead'
            )

        x = sympify_noeval(x)

        if not isinstance(xx, UniformGrid):
            xx = np.asarray([sympify_noeval(x) for x in xx])
        yy = np.asarray([sympify_noeval(y) for y in yy])

        if len(xx) != len(yy):
            raise ValueError(
                'length of xx and yy must be the same '
                f'(instead len(xx) = {len(xx)} and len(yy) = {len(yy)})'
            )

        if periodic and yy[0] != yy[-1]:
            raise ValueError(
                'if the spline is to be periodic, '
                'the first and last elements of yy must be equal'
            )

        if all(x.is_Number for x in xx) and not np.all(np.diff(xx) >= 0):
            raise ValueError('xx should be stricly increasing')

        extrapolate = \
            AbstractSpline._normalize_extrapolate(extrapolate, periodic)

        self._sbmlId = sbmlId
        self._x = x
        self._xx = xx
        self._yy = yy
        self._periodic = periodic
        self._extrapolate = extrapolate

    @staticmethod
    def _normalize_extrapolate(extrapolate, periodic):
        """
        Preprocess `extrapolate` to a standard form
        and perform consistency checks
        """
        if not isinstance(extrapolate, tuple):
            extrapolate = (extrapolate, extrapolate)
        if len(extrapolate) != 2:
            raise TypeError(
                f'extrapolate should be a 2-tuple, gor {extrapolate} instead'
            )
        valid_extrapolate = (
            'constant',
            'no_extrapolation',
            'periodic',
            None
        )
        for i in (0, 1):
            if extrapolate[i] not in valid_extrapolate:
                raise ValueError(
                    'currently the supported extrapolation methods are: ' +
                    str(valid_extrapolate) +
                    f' (got {extrapolate[i]} instead)'
                )
            if extrapolate[i] == 'no_extrapolation':
                extrapolate[i] = None
            if extrapolate[i] == 'periodic' and not periodic:
                raise ValueError(
                    'The spline must be satisfy periodic boundary conditions '
                    'in order for periodic extrapolation to be used.'
                )

        if (extrapolate[0] == 'periodic' or extrapolate[1] == 'periodic') and \
                extrapolate[0] != extrapolate[1] and \
                not (extrapolate[0] is None or extrapolate[1] is None):
            raise NotImplementedError(
                'At the moment if periodic extrapolation is applied '
                'to one side, the extrapolation at the other side '
                'must either be periodic or not be applied '
                '(in which case it will be periodic anyway).'
            )

        return extrapolate

    @property
    def sbmlImporter(self) -> sp.Symbol:
        return self._sbmlImporter

    @property
    def sbmlId(self) -> sp.Symbol:
        "SBML ID of the spline parameter."
        return self._sbmlId

    @property
    def x(self) -> sp.Basic:
        "The symbolic argument at which the spline is evaluated."
        return self._x

    @property
    def xx(self) -> np.ndarray:
        "The points at which the spline values are known."
        return self._xx

    @property
    def yy(self) -> np.ndarray:
        "The spline values at each of the points in `xx`."
        return self._yy

    @property
    def periodic(self) -> bool:
        "Whether periodic boundary conditions are applied."
        return self._periodic

    @property
    def extrapolate(self) -> Tuple[Union[None, str], Union[None, str]]:
        "Whether to extrapolate the spline outside the base interval."
        return self._extrapolate

    @property
    @abstractmethod
    def method(self) -> str:
        "Spline method."
        return NotImplemented

    @property
    def bc(self) -> Union[str, None]:
        """
        Return the boundary conditions applied to this spline.
        A value of `None` means that the spline type does not require
        to explicitly define boundary conditions.
        """
        return None

    def check_if_valid(self, importer: SbmlImporter):
        """
        Check if the spline described by this object can be correctly
        be implemented by AMICI. E.g., check whether the formulas
        for spline grid points, values, ... contain species symbols.
        """
        # TODO this is very much a draft

        fixed_parameters = importer.symbols['fixed_parameter']['identifier']
        species = list(importer.symbols['species']['identifier'])

        for x in self.xx:
            if not x.free_symbols.issubset(fixed_parameters):
                raise ValueError('xx should only depend on constant parameters')

        for y in self.yy:
            if len(y.free_symbols.intersection(species)) != 0:
                raise ValueError('yy should not depend on model species')

        fixed_parameters_values = importer.symbols['fixed_parameter']['value']
        subs = dict(zip(fixed_parameters, fixed_parameters_values))
        xx_values = [sp.simplify(x.subs(subs)) for x in self.xx]
        for x in xx_values:
            assert x.is_Number
        if not np.all(np.diff(xx) >= 0):
            raise ValueError('xx should be stricly increasing')

    def poly(self, i: int) -> sp.Basic:
        """
        Get the polynomial interpolant on the `(xx[i], xx[i+1])` interval
        in Horner form.
        """
        if i < 0:
            i += len(self.xx) - 1

        assert 0 <= i < len(self.xx) - 1

        poly = sp.poly(self._poly(self.x, i), wrt=self.x)
        return sp.horner(poly)

    @abstractmethod
    def _poly(self, x, i) -> sp.Basic:
        """
        Return the symbolic expression for polynomial interpolant
        on the `(xx[i], xx[i+1])` interval using `x` as the spline parameter.
        """
        return NotImplemented

    @property
    def extrapolation_formulas(self) \
            -> Tuple[Union[None, sp.Basic], Union[None, sp.Basic]]:
        """
        Returns the extrapolation formulas on the left and right side
        of the interval `(xx[0], xx[-1])`.
        A value of `None` means that no extrapolation is required.
        """
        extrLeft, extrRight = self.extrapolate

        if extrLeft == 'constant':
            extrLeft = self.yy[0]
        elif extrLeft is not None:
            raise ValueError(f'unsupported extrapolation method {extrLeft}')

        if extrRight == 'constant':
            extrRight = self.yy[-1]
        elif extrRight is not None:
            raise ValueError(f'unsupported extrapolation method {extrRight}')

        return (extrLeft, extrRight)

    @property
    def formula(self) -> sp.Piecewise:
        """
        Compute a symbolic piecewise formula for the spline.
        """
        return self._formula(sbml=False)

    @property
    def sbmlFormula(self) -> sp.Piecewise:
        """
        Compute a symbolic piecewise formula for the spline for use inside
        a SBML assignment rule (e.g., the AMICI time symbol will be replaced
        with its SBML counterpart).
        """
        return self._formula(sbml=True)

    def _formula(self, *, sbml: bool = False) -> sp.Piecewise:
        x = self.x
        pieces = []

        if self.extrapolate[0] == 'periodic' or self.extrapolate[1] == 'periodic':
            if sbml:
                # NB mod is not supported in SBML
                x = sympy.Symbol(self.sbmlId + '_x_in_fundamental_period')
                # NB we will do the parameter substitution in SBML
                #    because the formula for x will be a piecewise
                #    and sympy handles Piecewises inside other Piecewises
                #    really badly.
            else:
                x = self.xx[0] + sp.Mod(x - self.xx[0], self.xx[-1] - self.xx[0])
            extrLeft, extrRight = None, None
        else:
            extrLeft, extrRight = self.extrapolation_formulas

        if extrLeft is not None:
            pieces.append((extrLeft, x < self.xx[0]))

        for i in range(len(self.xx) - 2):
            pieces.append((self.poly(i), x < self.xx[i+1]))

        if extrRight is not None:
            pieces.append((self.poly(-1), x < self.xx[-1]))
            pieces.append((extrRight, sp.sympify(True)))
        else:
            pieces.append((self.poly(-1), sp.sympify(True)))

        with evaluate(False):
            if sbml:
                pieces = [
                    (
                        p.subs(amici_time_symbol, sbml_time_symbol),
                        c.subs(amici_time_symbol, sbml_time_symbol)
                    )
                    for (p, c) in pieces
                ]
            return sp.Piecewise(*pieces)

    @property
    def amiciAnnotation(self) -> str:
        """
        An SBML annotation describing the spline.
        """
        annotation = '<amici:spline xmlns:amici="http://github.com/AMICI-dev/AMICI"'
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
        dom = xml.dom.minidom.parseString(annotation)
        return dom.toprettyxml()

    def _annotation_attributes(self) -> Dict[str, Any]:
        attributes = {}

        attributes['spline_method'] = self.method

        if self.bc is not None:
            attributes['spline_bc'] = self.bc

        if self.extrapolate[0] == self.extrapolate[1]:
            extr = None if self.extrapolate is None else self.extrapolate[0]
        else:
            extr1, extr2 = self.extrapolate
            extr1 = 'no_extrapolation' if extr1 is None else extr1
            extr2 = 'no_extrapolation' if extr2 is None else extr2
            extr = f'({extr1}, {extr2})'
        if extr is not None:
            attributes['spline_extrapolate'] = extr

        return attributes

    def _annotation_children(self) -> Dict[str, Union[str, List[str]]]:
        children = {}

        with evaluate(False):
            x = self.x.subs(amici_time_symbol, sbml_time_symbol)
        children['spline_parameter'] = sbmlMathML(x)

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
            auto_add: bool = True,
            x_nominal: Sequence[float] = None,
            y_nominal: Optional[Union[Sequence[float], float]] = None,
            y_units: Optional[str] = None,
            x_units: Optional[str] = None
        ):
        """
        Function add the spline to an SBML model using an assignment rule
        with AMICI-specific annotations.

            :param sbmlImporter:
            A `SbmlImporter` object
            (used for testing whether formulas are adimissible).

            :param model:
            SBML model to which the spline is to be added.

            :param auto_add:
            Automatically missing parameters to the SBML model.
            Only used for expressions consisting in a single symbol.

            :param x_nominal:
            Values used when auto-adding parameters from `xx`.

            :param y_nominal:
            Values used when auto-adding parameters from `yy`.

            :param x_units:
            Units used when auto-adding parameters from `xx`.

            :param y_units:
            Units used when auto-adding parameters from `yy`.
        """
        # Convert time from AMICI to SBML naming
        with evaluate(False):
            x = self.x.subs(amici_time_symbol, sbml_time_symbol)

        # Try to autodetermine units
        if x_units is None:
            x_units = getUnits(model, x)
            for _x in self.xx:
                if x_units is not None:
                    break
                x_units = getUnits(model, _x)
        if y_units is None:
            for _y in self.yy:
                y_units = getUnits(model, _y)
                if y_units is not None:
                    break

        # Autoadd parameters
        if auto_add:
            if not hasParameter(model, self.sbmlId):
                addParameter(model, self.sbmlId, constant=False, units=y_units)

            if isinstance(x_nominal, collections.abc.Sequence):
                if len(x_nominal) != len(self.xx):
                    raise ValueError(
                        'if x_nominal is a list, then it must have '
                        'the same length as the spline grid'
                    )
            else:
                x_nominal = len(self.xx) * [x_nominal]
            for (_x, _val) in zip(self.xx, x_nominal):
                if _x.is_Symbol and not hasParameter(model, _x.name):
                    addParameter(model, _x.name, value=_val, units=x_units)

            if isinstance(y_nominal, collections.abc.Sequence):
                if len(y_nominal) != len(self.yy):
                    raise ValueError(
                        'if y_nominal is a list, then it must have '
                        'the same length as the spline values'
                    )
            else:
                y_nominal = len(self.yy) * [y_nominal]
            for (_y, _val) in zip(self.yy, y_nominal):
                if _y.is_Symbol and not hasParameter(model, _y.name):
                    addParameter(model, _y.name, value=_val, units=y_units)

        # Create assignment rule for spline
        rule = addAssignmentRule(model, self.sbmlId, self.sbmlFormula)

        # Add annotation specifying spline method
        rule.setAnnotation(self.amiciAnnotation)

        # Create additional assignment rule for periodic extrapolation
        # NB mod is not in the subset of MathML supported by SBML
        if self.extrapolate[0] == 'periodic' or self.extrapolate[1] == 'periodic':
            parameterId = self.sbmlId + '_x_in_fundamental_period'
            T = self.xx[-1] - self.xx[0]
            x0 = self.xx[0]
            s = 2 * sp.pi * ((x - x0) / T - sp.sympify(1)/4)
            k = sp.Piecewise((3, sp.cos(s) < 0), (1, True))
            formula = x0 + T * (sp.atan(sp.tan(s)) / (2 * sp.pi) + k/4)
            assert amici_time_symbol not in formula.free_symbols
            addParameter(model, parameterId, constant=False, units=x_units)
            addAssignmentRule(model, parameterId, formula)

    @classmethod
    def fromAnnotation(cls, annotation: libsbml.XMLNode):
        pass


class CubicHermiteSpline(AbstractSpline):
    def __init__(
            self,
            sbmlId: Union[str, sp.Symbol],
            x: Union[str, sp.Basic],
            xx: Sequence,
            yy: Sequence,
            dd: Sequence,
            *,
            periodic: Optional[bool] = None,
            extrapolate: Union[
                None, str,
                Tuple[Union[None, str], Union[None, str]]
            ] = None
        ):
        """
        Constructor for `CubicHermiteSpline` objects.

            :param sbmlId:
            The SBML ID of the parameter associated to the spline as a string or
            a SymPy symbol.

            :param x:
            The symbolic argument at which the spline is evaluated.
            It will be sympified.

            :param xx:
            The points at which the spline values are known.
            Currently they can only contain depend on constant parameters.
            It should be strictly increasing.
            It will be sympified.

            :param yy:
            The spline values at each of the points in `xx`.
            They may not depend on model species.
            It will be sympified.

            :param dd:
            The spline derivatives at each of the points in `xx`.
            They may not depend on model species.
            It will be sympified.

            :param periodic:
            Whether the spline satisfies periodic boundary conditions.
            If `None` (default) it is autodetected from `yy` and `dd`.

            Whether to extrapolate the spline outside the base interval
            defined by `(xx[0], xx[-1])`.
            It is a tuple of extrapolation methods, one for each side of the
            base interval.
            `None` or `'no_extrapolation'`
            means no extrapolation is to be performed
            (evaluating the spline outside the base interval is not guaranteed
            to work).
            If it is not a tuple, then the same extrapolation will be applied
            on both sides.
        """

        yy = np.asarray([sympify_noeval(y) for y in yy])
        dd = np.asarray([sympify_noeval(d) for d in dd])

        if len(xx) != len(dd):
            raise ValueError(
                'length of xx and dd must be the same '
                f'(instead len(xx) = {len(xx)} and len(dd) = {len(dd)})'
            )

        if periodic is None:
            periodic = yy[0] == yy[-1] and dd[0] == dd[-1]
        elif periodic:
            if yy[0] != yy[-1] or dd[0] != dd[-1]:
                raise ValueError(
                    'periodic=True but given yy and dd do not satisfy '
                    'periodic boundary conditions!'
                )

        super().__init__(
            sbmlId, x, xx, yy,
            periodic=periodic,
            extrapolate=extrapolate
        )

        self._dd = dd

    @property
    def dd(self) -> np.ndarray:
        "The spline derivatives at each of the points in `xx`."
        return self._dd

    @property
    def method(self) -> str:
        "Spline method (cubic Hermite spline)"
        return 'cubic_hermite'

    def check_if_valid(self, importer: SbmlImporter):
        """
        Check if the spline described by this object can be correctly
        be implemented by AMICI. E.g., check whether the formulas
        for spline grid points, values, ... contain species symbols.
        """
        # TODO this is very much a draft

        species = list(importer.symbols['species']['identifier'])
        for d in self.dd:
            if len(d.free_symbols.intersection(species)) != 0:
                raise ValueError('dd should not depend on model species')

        super().check_if_valid(importer)

    def _poly(self, x, i) -> sp.Basic:
        """
        Return the symbolic expression for polynomial interpolant
        on the `(xx[i], xx[i+1])` interval using `x` as the spline parameter.
        """
        assert 0 <= i < len(self.xx) - 1

        dx = self.xx[i+1] - self.xx[i]
        t = (x - self.xx[i]) / dx

        h00 = 2*t**3 - 3*t**2 + 1
        h10 = t**3 - 2*t**2 + t
        h01 = -2*t**3 + 3*t**2
        h11 = t**3 - t**2

        y0 = self.yy[i]
        y1 = self.yy[i+1]
        dy0 = self.dd[i]
        dy1 = self.dd[i+1]

        with evaluate(False):
            return h00 * y0 + h10 * dx * dy0 + h01 * y1 + h11 * dx * dy1

    def _annotation_children(self) -> Dict[str, Union[str, List[str]]]:
        children = super()._annotation_children()
        children['spline_derivatives'] = [sbmlMathML(d) for d in self.dd]
        return children


class CatmullRomSpline(CubicHermiteSpline):
    def __init__(
            self,
            sbmlId: Union[str, sp.Symbol],
            x: Union[str, sp.Basic],
            xx: UniformGrid,
            yy: Sequence,
            *,
            periodic: bool = False,
            extrapolate: Union[
                None, str,
                Tuple[Union[None, str], Union[None, str]]
            ] = None
        ):
        """
        Constructor for `CatmullRomSpline` objects.

            :param sbmlId:
            The SBML ID of the parameter associated to the spline as a string or
            a SymPy symbol.

            :param x:
            The symbolic argument at which the spline is evaluated.
            It will be sympified.

            :param xx:
            The points at which the spline values are known.
            They must be a `UniformGrid` object
            (uniformly spaced points).

            :param yy:
            The spline values at each of the points in `xx`.
            They may not depend on model species.
            It will be sympified.

            :param periodic:
            Whether the spline satisfies periodic boundary conditions.

            Whether to extrapolate the spline outside the base interval
            defined by `(xx[0], xx[-1])`.
            It is a tuple of extrapolation methods, one for each side of the
            base interval.
            `None` or `'no_extrapolation'`
            means no extrapolation is to be performed
            (evaluating the spline outside the base interval is not guaranteed
            to work).
            If it is not a tuple, then the same extrapolation will be applied
            on both sides.
        """

        if not isinstance(xx, UniformGrid):
            raise TypeError('Catmull-Rom splines are only defined for uniformly-space grid points')

        extrapolate = \
            AbstractSpline._normalize_extrapolate(extrapolate, periodic)

        yy = np.asarray([sympify_noeval(y) for y in yy])

        dd = []

        if periodic:
            dd.append(sp.Mul(
                1 / (2 * xx.step),
                yy[1] - yy[-2],
                evaluate=False
            ))
        elif extrapolate[0] == 'constant':
            dd.append(0)
        else:
            dd.append(sp.Mul(
                1 / xx.step,
                yy[1] - yy[0],
                evaluate=False
            ))

        for i in range(1, len(xx) - 1):
            dd.append(sp.Mul(
                1 / (2 * xx.step),
                yy[i+1] - yy[i-1],
                evaluate=False
            ))

        if periodic:
            dd.append(sp.Mul(
                1 / (2 * xx.step),
                yy[1] - yy[-2],
                evaluate=False
            ))
        elif extrapolate[1] == 'constant':
            dd.append(0)
        else:
            dd.append(sp.Mul(
                1 / xx.step,
                yy[-1] - yy[-2],
                evaluate=False
            ))

        super().__init__(
            sbmlId, x, xx, yy, dd,
            periodic=periodic, extrapolate=extrapolate
        )

    @property
    def method(self) -> str:
        "Spline method (Catmull-Rom spline)"
        return 'catmull_rom'

    def _annotation_children(self) -> Dict[str, Union[str, List[str]]]:
        children = super()._annotation_children()
        del children['spline_derivatives']  # not needed
        return children

    def _annotation_attributes(self) -> Dict[str, Any]:
        attributes = super()._annotation_attributes()
        attributes['spline_periodic'] = self.periodic
        return attributes
