"""
Test AMICI's Python spline implementation by
creating splines with different properties,
evaluating them and comparing them with
the true analytical values.
"""

import math

import amici
import sympy as sp
from amici.testing import skip_on_valgrind


@skip_on_valgrind
def test_SplineUniform():
    spline = amici.splines.CubicHermiteSpline(
        sbml_id="f",
        evaluate_at=amici.sbml_utils.amici_time_symbol,
        nodes=amici.splines.UniformGrid(0, 1, number_of_nodes=4),
        values_at_nodes=[0.0, 2.0, 0.5, 1.0],
    )
    assert math.isclose(float(spline.evaluate(0.0)), 0.0)
    assert math.isclose(float(spline.evaluate(0.25)), 1.74609375)
    assert math.isclose(float(spline.evaluate(1.0 / 3)), 2.0)
    assert math.isclose(float(spline.evaluate(0.50)), 1.3437499999999996)
    assert math.isclose(float(spline.evaluate(2.0 / 3)), 0.5)
    assert math.isclose(float(spline.evaluate(0.75)), 0.484375)
    assert math.isclose(float(spline.evaluate(1.00)), 1.0)


@skip_on_valgrind
def test_SplineNonUniform():
    spline = amici.splines.CubicHermiteSpline(
        sbml_id="f",
        evaluate_at=amici.sbml_utils.amici_time_symbol,
        nodes=[0.0, 0.1, 0.5, 1.0],
        values_at_nodes=[0.0, 2.0, 0.5, 1.0],
    )
    assert math.isclose(float(spline.evaluate(0.00)), 0.0)
    assert math.isclose(float(spline.evaluate(0.05)), 1.1484375)
    assert math.isclose(float(spline.evaluate(0.10)), 2.0)
    assert math.isclose(float(spline.evaluate(0.25)), 2.0498046875)
    assert math.isclose(float(spline.evaluate(0.50)), 0.5)
    assert math.isclose(float(spline.evaluate(0.75)), 0.6015625)
    assert math.isclose(float(spline.evaluate(1.00)), 1.0)


@skip_on_valgrind
def test_SplineExplicit():
    spline = amici.splines.CubicHermiteSpline(
        sbml_id="f",
        evaluate_at=amici.sbml_utils.amici_time_symbol,
        nodes=amici.splines.UniformGrid(0, 1, number_of_nodes=5),
        values_at_nodes=[0.0, 2.0, 0.5, 1.0, 0.75],
        derivatives_at_nodes=[1.0, 0.0, 0.1, -0.1, 0.0],
    )
    assert math.isclose(float(spline.evaluate(0.00)), 0.0)
    assert math.isclose(float(spline.evaluate(0.20)), 1.8000000000000003)
    assert math.isclose(float(spline.evaluate(0.25)), 2.0)
    assert math.isclose(float(spline.evaluate(0.40)), 1.02439999999999985)
    assert math.isclose(float(spline.evaluate(0.50)), 0.5)
    assert math.isclose(float(spline.evaluate(0.60)), 0.6819999999999999)
    assert math.isclose(float(spline.evaluate(0.75)), 1.0)
    assert math.isclose(float(spline.evaluate(0.80)), 0.9707999999999999)
    assert math.isclose(float(spline.evaluate(1.00)), 0.75)


@skip_on_valgrind
def test_SplineZeroBC():
    spline = amici.splines.CubicHermiteSpline(
        sbml_id="f",
        evaluate_at=amici.sbml_utils.amici_time_symbol,
        nodes=amici.splines.UniformGrid(0, 1, number_of_nodes=4),
        values_at_nodes=[0.0, 2.0, 0.5, 1.0],
        bc="zeroderivative",
    )
    assert math.isclose(float(spline.evaluate(0.00)), 0.0)
    assert math.isclose(float(spline.evaluate(0.25)), 1.65234375)
    assert math.isclose(float(spline.evaluate(0.50)), 1.3437499999999996)
    assert math.isclose(float(spline.evaluate(0.75)), 0.5078125)
    assert math.isclose(float(spline.evaluate(1.00)), 1.0)


@skip_on_valgrind
def test_SplineLogarithmic():
    spline = amici.splines.CubicHermiteSpline(
        sbml_id="f",
        evaluate_at=amici.sbml_utils.amici_time_symbol,
        nodes=amici.splines.UniformGrid(0, 1, number_of_nodes=5),
        values_at_nodes=[0.2, 2.0, 0.5, 1.0, 0.75],
        logarithmic_parametrization=True,
    )
    assert math.isclose(float(spline.evaluate(0.00)), 0.2)
    assert math.isclose(float(spline.evaluate(0.20)), 2.07939779651678)
    assert math.isclose(float(spline.evaluate(0.25)), 2.0)
    assert math.isclose(float(spline.evaluate(0.40)), 0.947459046694449)
    assert math.isclose(float(spline.evaluate(0.50)), 0.5)
    assert math.isclose(float(spline.evaluate(0.60)), 0.545987404053269)
    assert math.isclose(float(spline.evaluate(0.75)), 1.0)
    assert math.isclose(float(spline.evaluate(0.80)), 0.996753014029391)
    assert math.isclose(float(spline.evaluate(1.00)), 0.75)


@skip_on_valgrind
def test_SplineUniformConstantExtrapolation():
    spline = amici.splines.CubicHermiteSpline(
        sbml_id="f",
        evaluate_at=amici.sbml_utils.amici_time_symbol,
        nodes=amici.splines.UniformGrid(0, 1, number_of_nodes=4),
        values_at_nodes=[0.0, 2.0, 0.5, 1.0],
        extrapolate="constant",
    )
    assert math.isclose(float(spline.evaluate(-2.00)), 0.0)
    assert math.isclose(float(spline.evaluate(-1.00)), 0.0)
    assert math.isclose(float(spline.evaluate(0.00)), 0.0)
    assert math.isclose(float(spline.evaluate(0.25)), 1.65234375)
    assert math.isclose(float(spline.evaluate(1.0 / 3)), 2.0)
    assert math.isclose(float(spline.evaluate(0.50)), 1.3437499999999996)
    assert math.isclose(float(spline.evaluate(2.0 / 3)), 0.5)
    assert math.isclose(float(spline.evaluate(0.75)), 0.5078125)
    assert math.isclose(float(spline.evaluate(1.00)), 1.0)
    assert math.isclose(float(spline.evaluate(2.00)), 1.0)
    assert math.isclose(float(spline.evaluate(3.00)), 1.0)


@skip_on_valgrind
def test_SplineUniformLinearExtrapolation():
    spline = amici.splines.CubicHermiteSpline(
        sbml_id="f",
        evaluate_at=amici.sbml_utils.amici_time_symbol,
        nodes=amici.splines.UniformGrid(0, 1, number_of_nodes=4),
        values_at_nodes=[0.0, 2.0, 0.5, 1.0],
        extrapolate="linear",
    )
    assert math.isclose(float(spline.evaluate(-2.00)), -12.0)
    assert math.isclose(float(spline.evaluate(-1.00)), -6.0)
    assert math.isclose(float(spline.evaluate(0.00)), 0.0)
    assert math.isclose(float(spline.evaluate(0.25)), 1.74609375)
    assert math.isclose(float(spline.evaluate(1.0 / 3)), 2.0)
    assert math.isclose(float(spline.evaluate(0.50)), 1.3437499999999996)
    assert math.isclose(float(spline.evaluate(2.0 / 3)), 0.5)
    assert math.isclose(float(spline.evaluate(0.75)), 0.484375)
    assert math.isclose(float(spline.evaluate(1.00)), 1.0)
    assert math.isclose(float(spline.evaluate(2.00)), 2.5)
    assert math.isclose(float(spline.evaluate(3.00)), 4.0)


@skip_on_valgrind
def test_SplineUniformPolynomialExtrapolation():
    spline = amici.splines.CubicHermiteSpline(
        sbml_id="f",
        evaluate_at=amici.sbml_utils.amici_time_symbol,
        nodes=amici.splines.UniformGrid(0, 1, number_of_nodes=4),
        values_at_nodes=[0.0, 2.0, 0.5, 1.0],
        extrapolate="polynomial",
    )
    assert math.isclose(float(spline.evaluate(-2.00)), 429.0)
    assert math.isclose(float(spline.evaluate(-1.00)), 57.0)
    assert math.isclose(float(spline.evaluate(0.00)), 0.0)
    assert math.isclose(float(spline.evaluate(0.25)), 1.74609375)
    assert math.isclose(float(spline.evaluate(1.0 / 3)), 2.0)
    assert math.isclose(float(spline.evaluate(0.50)), 1.3437499999999996)
    assert math.isclose(float(spline.evaluate(2.0 / 3)), 0.5)
    assert math.isclose(float(spline.evaluate(0.75)), 0.484375)
    assert math.isclose(float(spline.evaluate(1.00)), 1.0)
    assert math.isclose(float(spline.evaluate(2.00)), -33.5)
    assert math.isclose(float(spline.evaluate(3.00)), -248.0)


@skip_on_valgrind
def test_SplineUniformPeriodicExtrapolation():
    spline = amici.splines.CubicHermiteSpline(
        sbml_id="f",
        evaluate_at=amici.sbml_utils.amici_time_symbol,
        nodes=amici.splines.UniformGrid(0, 1, number_of_nodes=4),
        values_at_nodes=[1.0, 2.0, 0.5, 1.0],
        extrapolate="periodic",
    )
    assert math.isclose(float(spline.evaluate(-4 / 3)), 0.5)
    assert math.isclose(float(spline.evaluate(-0.5)), 1.2812499999999996)
    assert math.isclose(float(spline.evaluate(0.00)), 1.0)
    assert math.isclose(float(spline.evaluate(0.25)), 1.9140625)
    assert math.isclose(float(spline.evaluate(1 / 3)), 2.0)
    assert math.isclose(float(spline.evaluate(0.50)), 1.2812499999999996)
    assert math.isclose(float(spline.evaluate(2 / 3)), 0.5)
    assert math.isclose(float(spline.evaluate(0.75)), 0.47265625)
    assert math.isclose(float(spline.evaluate(1.00)), 1.0)
    assert math.isclose(float(spline.evaluate(1.25)), 1.9140625)
    assert math.isclose(float(spline.evaluate(2.75)), 0.47265625)


@skip_on_valgrind
def test_SplineNonUniformPeriodicExtrapolation():
    spline = amici.splines.CubicHermiteSpline(
        sbml_id="f",
        evaluate_at=amici.sbml_utils.amici_time_symbol,
        nodes=[0.0, 0.1, 0.5, 1.0],
        values_at_nodes=[1.0, 2.0, 0.5, 1.0],
        extrapolate="periodic",
    )
    assert math.isclose(float(spline.evaluate(-1.90)), 2.0)
    assert math.isclose(float(spline.evaluate(-0.25)), 0.3203125)
    assert math.isclose(float(spline.evaluate(0.00)), 1.0)
    assert math.isclose(float(spline.evaluate(0.05)), 1.5296875)
    assert math.isclose(float(spline.evaluate(0.10)), 2.0)
    assert math.isclose(float(spline.evaluate(0.25)), 1.7568359375)
    assert math.isclose(float(spline.evaluate(0.50)), 0.5)
    assert math.isclose(float(spline.evaluate(0.75)), 0.3203125)
    assert math.isclose(float(spline.evaluate(1.00)), 1.0)
    assert math.isclose(float(spline.evaluate(1.50)), 0.5)
    assert math.isclose(float(spline.evaluate(2.05)), 1.5296875)


@skip_on_valgrind
def check_gradient(spline, t, params, params_values, expected, rel_tol=1e-9):
    value = spline.evaluate(t)
    subs = {
        pname: pvalue
        for (pname, pvalue) in zip(params, params_values, strict=True)
    }
    for p, exp in zip(params, expected, strict=True):
        assert math.isclose(
            float(value.diff(p).subs(subs)), exp, rel_tol=rel_tol
        )


@skip_on_valgrind
def test_SplineUniformSensitivity():
    params = (a, b, c) = sp.symbols("a b c")
    params_values = [0.5, 1.0, 2.5]
    spline = amici.splines.CubicHermiteSpline(
        sbml_id="f",
        evaluate_at=amici.sbml_utils.amici_time_symbol,
        nodes=amici.splines.UniformGrid(0, 1, number_of_nodes=4),
        values_at_nodes=[3 * a + b, c**2 - 3, 1, sp.log(b) + 3 * c - 6 * a],
    )
    check_gradient(spline, 0.00, params, params_values, [3.0, 1.0, 0.0])
    check_gradient(
        spline,
        0.25,
        params,
        params_values,
        [0.539062, 0.179688, 4.45312],
        rel_tol=1e-5,
    )
    check_gradient(spline, 1.0 / 3, params, params_values, [0.0, 0.0, 5.0])
    check_gradient(
        spline, 0.50, params, params_values, [0.1875, -0.125, 2.625]
    )
    check_gradient(spline, 2.0 / 3, params, params_values, [0.0, 0.0, 0.0])
    check_gradient(
        spline,
        0.75,
        params,
        params_values,
        [-1.07812, 0.179688, 0.1875],
        rel_tol=1e-5,
    )
    check_gradient(spline, 1.00, params, params_values, [-6.0, 1.0, 3.0])


@skip_on_valgrind
def test_SplineNonUniformSensitivity():
    params = (a, b, c) = sp.symbols("a b c")
    params_values = [0.5, 1.0, 2.5]
    spline = amici.splines.CubicHermiteSpline(
        sbml_id="f",
        evaluate_at=amici.sbml_utils.amici_time_symbol,
        nodes=[0.0, 0.1, 0.5, 1.0],
        values_at_nodes=[3 * a + b, c**2 - 3, 1, sp.log(b) + 3 * c - 6 * a],
    )
    check_gradient(spline, 0.00, params, params_values, [3.0, 1.0, 0.0])
    check_gradient(
        spline,
        0.05,
        params,
        params_values,
        [1.3125, 0.4375, 2.89062],
        rel_tol=1e-5,
    )
    check_gradient(spline, 0.10, params, params_values, [0.0, 0.0, 5.0])
    check_gradient(spline, 0.30, params, params_values, [-0.45, -0.3, 3.6])
    check_gradient(spline, 0.50, params, params_values, [0.0, 0.0, 0.0])
    check_gradient(
        spline, 0.75, params, params_values, [-2.625, 0.4375, 0.921875]
    )
    check_gradient(spline, 1.00, params, params_values, [-6.0, 1.0, 3.0])


@skip_on_valgrind
def test_SplineExplicitSensitivity():
    params = (a, b, c) = sp.symbols("a b c")
    params_values = [0.5, 1.0, 2.5]
    spline = amici.splines.CubicHermiteSpline(
        sbml_id="f",
        evaluate_at=amici.sbml_utils.amici_time_symbol,
        nodes=amici.splines.UniformGrid(0, 1, number_of_nodes=4),
        values_at_nodes=[3 * a + b, c**2 - 3, 1, sp.log(b) + 3 * c - 6 * a],
        derivatives_at_nodes=[
            c**3 - 2,
            sp.sqrt(b) * sp.log(b) + 3 * c,
            4 * a - sp.sin(b),
            1,
        ],
    )
    check_gradient(spline, 0.00, params, params_values, [3.0, 1.0, 0.0])
    check_gradient(
        spline,
        0.25,
        params,
        params_values,
        [0.46875, 0.109375, 4.37109],
        rel_tol=1e-6,
    )
    check_gradient(spline, 1.0 / 3, params, params_values, [0.0, 0.0, 5.0])
    check_gradient(
        spline,
        0.50,
        params,
        params_values,
        [-0.166667, 0.0641793, 2.625],
        rel_tol=1e-5,
    )
    check_gradient(spline, 2.0 / 3, params, params_values, [0.0, 0.0, 0.0])
    check_gradient(
        spline,
        0.75,
        params,
        params_values,
        [-0.75, 0.130923, 0.46875],
        rel_tol=1e-5,
    )
    check_gradient(spline, 1.00, params, params_values, [-6.0, 1.0, 3.0])


@skip_on_valgrind
def test_SplineLogarithmicSensitivity():
    params = (a, b, c) = sp.symbols("a b c")
    params_values = [0.5, 1.0, 2.5]
    spline = amici.splines.CubicHermiteSpline(
        sbml_id="f",
        evaluate_at=amici.sbml_utils.amici_time_symbol,
        nodes=amici.splines.UniformGrid(0, 1, number_of_nodes=4),
        values_at_nodes=[3 * a + b, c**2 - 3, 1, sp.log(b) + 3 * c - 6 * a],
        logarithmic_parametrization=True,
    )
    check_gradient(spline, 0.00, params, params_values, [3.0, 1.0, 0.0])
    check_gradient(
        spline,
        0.25,
        params,
        params_values,
        [0.585881, 0.195294, 4.38532],
        rel_tol=1e-5,
    )
    check_gradient(spline, 1.0 / 3, params, params_values, [0.0, 0.0, 5.0])
    check_gradient(
        spline,
        0.50,
        params,
        params_values,
        [0.514003, -0.132395, 1.52044],
        rel_tol=1e-5,
    )
    check_gradient(spline, 2.0 / 3, params, params_values, [0.0, 0.0, 0.0])
    check_gradient(
        spline,
        0.75,
        params,
        params_values,
        [-0.820743, 0.13679, -0.0577988],
        rel_tol=1e-5,
    )
    check_gradient(spline, 1.00, params, params_values, [-6.0, 1.0, 3.0])
