import numpy as np
import sympy as sp
from numpy.testing import assert_almost_equal

from amici.conserved_moieties2 import nullspace_by_rref, pivots, rref


def test_rref():
    """Create some random matrices and compare out ``rref`` to sympy"""
    rng = np.random.default_rng(12345)
    for rows, cols in rng.integers(0, 10, [200, 2]):
        M = np.random.rand(rows, cols)

        M_rref = rref(M)
        exp_rref, exp_pivots = sp.Matrix(M).rref()
        exp = np.asarray(exp_rref, dtype=float)
        print(M)
        print(exp)
        print(M_rref)
        assert list(exp_pivots) == pivots(M_rref)
        assert np.allclose(M_rref, exp)


def test_nullspace_by_rref():
    """Test ``nullspace_by_rref`` on a number of random matrices and compare
    to sympy results"""
    rng = np.random.default_rng(12345)
    for rows, cols in rng.integers(0, 50, [100, 2]):
        M = np.random.rand(rows, cols)

        actual = nullspace_by_rref(M)

        expected = sp.Matrix(M).nullspace()
        expected = np.hstack(expected).T if len(expected) else np.array([])

        if len(actual):
            assert np.allclose(M.dot(actual.T), 0)

        assert_almost_equal(actual, expected)
