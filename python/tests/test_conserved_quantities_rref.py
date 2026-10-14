import numpy as np
import pytest
import sympy as sp
from amici.importers.sbml.conserved_quantities_rref import (
    nullspace_by_rref,
    pivots,
    rref,
)
from amici.testing import skip_on_valgrind


def random_matrix_generator(min_dim, max_dim, count):
    """Generate random 2D matrices"""
    rng = np.random.default_rng(12345)
    for rows, cols in rng.integers(min_dim, max_dim, [count, 2]):
        yield rng.random((rows, cols))


def random_stoichiometry_generator(max_dim, count):
    """Generate random sparse integer reaction x species matrices,
    i.e., transposed stoichiometric matrices"""
    rng = np.random.default_rng(12345)
    for n_species, n_reactions in rng.integers(2, max_dim, [count, 2]):
        stoich = np.zeros((n_species, n_reactions))
        for j in range(n_reactions):
            n_participants = min(rng.integers(1, 4), n_species)
            idx = rng.choice(n_species, size=n_participants, replace=False)
            stoich[idx, j] = rng.choice([-2, -1, 1, 2], size=n_participants)
        yield stoich.T


def assert_is_nullspace_basis(mat, basis, expected_dim):
    """Check that the rows of ``basis`` span a nullspace of ``mat`` of the
    expected dimension"""
    assert basis.shape[0] == expected_dim
    if not expected_dim:
        return
    assert basis.shape[1] == mat.shape[1]
    assert np.linalg.matrix_rank(basis) == expected_dim
    # the residual scales with the magnitude of the matrix entries and
    #  the basis vectors
    residual = np.abs(mat @ basis.T).max()
    scale = np.abs(mat).max() * np.abs(basis).max()
    assert residual <= 1e-10 * scale


@skip_on_valgrind
@pytest.mark.parametrize("mat", tuple(random_matrix_generator(0, 10, 200)))
def test_rref(mat):
    """Create some random matrices and compare output of ``rref`` and
    ``pivots`` to that of sympy"""
    actual_rref = rref(mat)
    expected_rref, expected_pivots = sp.Matrix(mat).rref()
    expected_rref = np.asarray(expected_rref, dtype=float)

    assert list(expected_pivots) == pivots(actual_rref)
    assert np.allclose(expected_rref, actual_rref)


@skip_on_valgrind
@pytest.mark.parametrize("mat", tuple(random_matrix_generator(0, 50, 50)))
def test_nullspace_by_rref(mat):
    """Test ``nullspace_by_rref`` on a number of random matrices and compare
    to sympy results"""
    actual = nullspace_by_rref(mat)

    expected_dim = mat.shape[1] - (
        np.linalg.matrix_rank(mat) if mat.size else 0
    )
    assert_is_nullspace_basis(mat, actual, expected_dim)

    expected = sp.Matrix(mat).nullspace()
    expected = (
        np.hstack(np.asarray(expected, dtype=float)).T
        if len(expected)
        else np.array([])
    )

    assert np.allclose(actual, expected, rtol=1e-8)


@skip_on_valgrind
@pytest.mark.parametrize("mat", tuple(random_stoichiometry_generator(25, 50)))
def test_nullspace_by_rref_stoichiometry(mat):
    """Test ``nullspace_by_rref`` on random stoichiometric matrices,
    which, unlike the dense random matrices above, are usually
    rank-deficient"""
    actual = nullspace_by_rref(mat)

    expected_dim = mat.shape[1] - sp.Matrix(mat.astype(int)).rank()
    assert_is_nullspace_basis(mat, actual, expected_dim)


@pytest.mark.xfail(
    reason="rref is not numerically robust, see "
    "https://github.com/AMICI-dev/AMICI/issues/3366",
    strict=True,
)
@pytest.mark.parametrize(
    "mat",
    [
        # misses one of three nullspace vectors
        [
            [0, 0, 0, 2, -1, 2, 0, 0, 0],
            [0, -2, -2, 0, 0, 0, 0, 0, 2],
            [0, -1, 2, 2, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, -2, 0, -1, 0, 0, 1, 0],
            [0, 0, 0, 0, -1, 0, 0, -2, 0],
            [0, 2, 0, 0, 0, 0, 0, 0, 0],
        ],
        # three vectors, but spanning only one dimension
        [
            [0, 0, 0, 0, 0, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0],
            [0, 0, -2, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, -2, 0, 0, 0, 0, 0, 0, 0, -2],
            [0, 0, -2, -1, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0],
            [0, 0, 0, 1, -2, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            [0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 1, 0, 0],
            [0, -2, 0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 1, 0, 0],
            [-2, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 2, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        ],
    ],
    ids=["wrong_dim", "wrong_basis"],
)
def test_nullspace_by_rref_stoichiometry_known_failures(mat):
    """Known failures of ``nullspace_by_rref`` on stoichiometric matrices"""
    mat = np.array(mat, dtype=float)
    actual = nullspace_by_rref(mat)

    expected_dim = mat.shape[1] - sp.Matrix(mat.astype(int)).rank()
    assert_is_nullspace_basis(mat, actual, expected_dim)
