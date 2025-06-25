"""
Test AMICI's C++ spline implementation by comparing
the results of simulations of simple SBML models
containing simple splines with a symbolically-computed
ground truth. NB: The test in this file takes a long
time to complete.
"""

import os

import numpy as np
from amici.testing import skip_on_valgrind
from splines_utils import (
    check_splines_full,
    example_spline_1,
    example_spline_2,
    example_spline_3,
)


@skip_on_valgrind
def test_multiple_splines(**kwargs):
    """
    Test a SBML model containing multiple splines.
    """
    spline0, params0, tols0 = example_spline_1(
        0, num_nodes=9, fixed_values=[0, 2], extrapolate="linear"
    )
    spline1, params1, tols1 = example_spline_1(
        1, num_nodes=14, scale=1.5, offset=5, extrapolate="linear"
    )
    spline2, params2, tols2 = example_spline_1(
        2, num_nodes=5, scale=0.5, offset=-5, extrapolate="linear"
    )
    spline3, params3, tols3 = example_spline_1(
        3, fixed_values="all", extrapolate="linear"
    )
    spline4, params4, tols4 = example_spline_2(4)
    spline5, params5, tols5 = example_spline_3(5)

    splines = [spline0, spline1, spline2, spline3, spline4, spline5]

    params = dict(params0)
    params.update(params1)
    params.update(params2)
    params.update(params3)
    params.update(params4)
    params.update(params5)

    if isinstance(tols0, dict):
        tols0 = (tols0, tols0, tols0)
    if isinstance(tols1, dict):
        tols1 = (tols1, tols1, tols1)
    if isinstance(tols2, dict):
        tols2 = (tols2, tols2, tols2)
    if isinstance(tols3, dict):
        tols3 = (tols3, tols3, tols3)
    if isinstance(tols4, dict):
        tols4 = (tols4, tols4, tols4)
    if isinstance(tols5, dict):
        tols5 = (tols5, tols5, tols5)

    tols = []
    for t0, t1, t2, t3, t4, t5 in zip(
        tols0, tols1, tols2, tols3, tols4, tols5, strict=True
    ):
        keys = set().union(
            t0.keys(), t1.keys(), t2.keys(), t3.keys(), t4.keys(), t5.keys()
        )
        t = {
            key: max(
                t0.get(key, 0.0),
                t1.get(key, 0.0),
                t2.get(key, 0.0),
                t3.get(key, 0.0),
                t4.get(key, 0.0),
                t5.get(key, 0.0),
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

    if os.name == "nt":
        tols[2]["sllh_atol"] = max(5e-4, tols[2]["sllh_atol"])

    # Load precomputed results
    # They be computed again by
    #     groundtruth = test_multiple_splines(return_groundtruth=True)
    # They should be recomputed only if the splines used in the test change
    precomputed_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "test_splines_precomputed.npz",
    )
    with np.load(precomputed_path) as data:
        kwargs["groundtruth"] = dict(data)

    return check_splines_full(splines, params, tols, **kwargs)
