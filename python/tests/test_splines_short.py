"""
Test AMICI's C++ spline implementation by comparing
the results of simulations of simple SBML models
containing simple splines with a symbolically-computed
ground truth.
"""

import numpy as np
import sympy as sp
from amici.splines import CubicHermiteSpline, UniformGrid
from amici.testing import skip_on_valgrind
from splines_utils import check_splines_full, example_spline_1


@skip_on_valgrind
def test_spline_piecewise(**kwargs):
    """
    Test a SBML model containing a single spline.
    AMICI's behaviour in absence of spline annotations is also tested.
    """
    spline, params, tols = example_spline_1()
    check_splines_full(spline, params, tols, **kwargs)


@skip_on_valgrind
def test_two_splines(**kwargs):
    """
    Test a SBML model containing two splines.
    """
    spline0, params0, tols0 = example_spline_1(
        0, num_nodes=4, fixed_values=[0, 2], extrapolate="linear"
    )
    spline1, params1, tols1 = example_spline_1(
        1, num_nodes=5, scale=1.5, offset=5, extrapolate="linear"
    )

    splines = [spline0, spline1]

    params = dict(params0)
    params.update(params1)

    if isinstance(tols0, dict):
        tols0 = (tols0, tols0, tols0)
    if isinstance(tols1, dict):
        tols1 = (tols1, tols1, tols1)

    tols = []
    for t0, t1 in zip(tols0, tols1, strict=True):
        keys = set().union(t0.keys(), t1.keys())
        t = {
            key: max(
                t0.get(key, 0.0),
                t1.get(key, 0.0),
            )
            for key in keys
        }
        tols.append(t)

    tols[1]["x_rtol"] = max(1e-9, tols[1].get("x_rtol", -np.inf))
    tols[1]["x_atol"] = max(5e-9, tols[1].get("x_atol", -np.inf))
    tols[1]["sx_rtol"] = max(1e-5, tols[1].get("llh_rtol", -np.inf))
    tols[1]["sx_atol"] = max(5e-9, tols[1].get("sx_atol", -np.inf))
    tols[1]["llh_rtol"] = max(5e-14, tols[1].get("llh_rtol", -np.inf))
    tols[1]["sllh_atol"] = max(5e-5, tols[1].get("sllh_atol", -np.inf))

    tols[2]["x_rtol"] = max(5e-10, tols[2].get("x_rtol", -np.inf))
    tols[2]["x_atol"] = max(1e-8, tols[2].get("x_atol", -np.inf))
    tols[2]["llh_rtol"] = max(5e-14, tols[2].get("llh_rtol", -np.inf))
    tols[2]["sllh_atol"] = max(5e-5, tols[2].get("sllh_atol", -np.inf))

    check_splines_full(splines, params, tols, check_piecewise=False, **kwargs)


@skip_on_valgrind
def test_splines_plist():
    """
    Test if AMICI's spline implementation
    handles correctly a change in the parameter list.
    """
    # Dummy spline #1
    xx = UniformGrid(0, 5, number_of_nodes=3)
    yy = np.asarray([0.0, 1.0, 0.5])
    spline1 = CubicHermiteSpline(
        "y1",
        nodes=xx,
        values_at_nodes=yy,
        bc="auto",
        extrapolate=(None, "constant"),
    )
    # Dummy spline #2
    xx = UniformGrid(0, 5, number_of_nodes=4)
    yy = np.asarray([0.0, 0.5, -0.5, 0.5])
    spline2 = CubicHermiteSpline(
        "y2",
        nodes=xx,
        values_at_nodes=yy,
        bc="auto",
        extrapolate=(None, "constant"),
    )
    # Real spline #3
    xx = UniformGrid(0, 5, number_of_nodes=6)
    p1, p2, p3, p4, p5 = sp.symbols("p1 p2 p3 p4 p5")
    yy = np.asarray(
        [p1 + p2, p2 * p3, p4, sp.cos(p1 + p3), p4 * sp.log(p1), p3]
    )
    dd = np.asarray([-0.75, -0.875, p5, 0.125, 1.15057181, 0.0])
    params = {p1: 1.0, p2: 0.5, p3: 1.5, p4: -0.25, p5: -0.5}
    # print([y.subs(params).evalf() for y in yy])
    spline3 = CubicHermiteSpline(
        "y3",
        nodes=xx,
        values_at_nodes=yy,
        derivatives_at_nodes=dd,
        bc="auto",
        extrapolate=(None, "constant"),
    )
    # Dummy spline 4
    xx = UniformGrid(0, 5, number_of_nodes=3)
    yy = np.asarray([0.0, -0.5, 0.5])
    spline4 = CubicHermiteSpline(
        "y4",
        nodes=xx,
        values_at_nodes=yy,
        bc="auto",
        extrapolate=(None, "constant"),
    )
    tols = dict(
        x_rtol=1e-6,
        x_atol=1e-11,
        sx_rtol=1e-6,
        sx_atol=5e-11,
        llh_rtol=1e-14,
        sllh_atol=5e-9,
    )
    check_splines_full(
        [spline1, spline2, spline3, spline4],
        params,
        tols,
        check_piecewise=False,
        check_forward=False,
        check_adjoint=True,  # plist cannot be checked, but complex parameter dependence can
        parameter_lists=[[0, 1, 4], [2, 3]],
    )
