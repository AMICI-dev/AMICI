import numpy as np
import sympy as sp
from numpy.testing import assert_almost_equal

from amici.conserved_moieties2 import nullspace_by_rref, pivots, rref


def test_rref():
    """Create some random matrices and compare out ``rref`` to sympy"""
    rng = np.random.default_rng(12345)
    for rows, cols in rng.integers(0, 10, [200, 2]):
        mat = np.random.rand(rows, cols)

        mat_rref = rref(mat)
        exp_rref, exp_pivots = sp.Matrix(mat).rref()
        exp = np.asarray(exp_rref, dtype=float)
        print(mat)
        print(exp)
        print(mat_rref)
        assert list(exp_pivots) == pivots(mat_rref)
        assert np.allclose(mat_rref, exp)


def test_nullspace_by_rref():
    """Test ``nullspace_by_rref`` on a number of random matrices and compare
    to sympy results"""
    rng = np.random.default_rng(12345)
    for rows, cols in rng.integers(0, 50, [50, 2]):
        mat = np.random.rand(rows, cols)

        actual = nullspace_by_rref(mat)

        if len(actual):
            assert np.allclose(mat.dot(actual.T), 0)

        expected = sp.Matrix(mat).nullspace()
        expected = np.hstack(expected).T if len(expected) else np.array([])

        assert_almost_equal(actual, expected)
