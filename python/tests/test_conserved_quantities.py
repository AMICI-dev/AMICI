"""Tests for the identification of conservation laws from the stoichiometric
matrix"""

import math

import libsbml
import numpy as np
import pytest
import sympy as sp
from amici.importers.sbml import (
    SbmlImporter,
    _conservation_laws_from_stoichiometry,
)
from amici.testing import skip_on_valgrind


def random_stoichiometry_generator(max_dim, count):
    """Generate random sparse integer species x reactions matrices"""
    rng = np.random.default_rng(12345)
    for n_species, n_reactions in rng.integers(2, max_dim, [count, 2]):
        stoich = np.zeros((n_species, n_reactions), dtype=int)
        for j in range(n_reactions):
            n_participants = min(rng.integers(1, 4), n_species)
            idx = rng.choice(n_species, size=n_participants, replace=False)
            stoich[idx, j] = rng.choice([-2, -1, 1, 2], size=n_participants)
        yield stoich


def cls_to_matrix(raw_cls, n_species):
    """Assemble conservation laws into a matrix with one row per law"""
    mat = sp.zeros(len(raw_cls), n_species)
    for i_cl, (_, state_idxs, coefficients) in enumerate(raw_cls):
        for i, coeff in zip(state_idxs, coefficients, strict=True):
            mat[i_cl, i] = coeff
    return mat


def assert_exact_conservation_laws(stoichiometric_matrix, raw_cls):
    """Check that ``raw_cls`` is the reduced row echelon form of a basis of
    the left nullspace of ``stoichiometric_matrix``"""
    S = sp.Matrix(stoichiometric_matrix)
    n_species = S.shape[0]
    expected_dim = n_species - S.rank()
    assert len(raw_cls) == expected_dim

    cls = cls_to_matrix(raw_cls, n_species)
    assert all(isinstance(x, sp.Rational) for x in cls)
    assert (cls * S).is_zero_matrix
    assert cls.rank() == expected_dim

    pivots = [pivot for pivot, _, _ in raw_cls]
    assert pivots == sorted(set(pivots))
    for i_cl, pivot in enumerate(pivots):
        assert raw_cls[i_cl][1][0] == pivot
        assert list(cls[:, pivot]) == [
            int(i == i_cl) for i in range(len(pivots))
        ]


@skip_on_valgrind
@pytest.mark.parametrize(
    "stoichiometric_matrix", tuple(random_stoichiometry_generator(25, 50))
)
def test_conservation_laws_random(stoichiometric_matrix):
    """Conservation laws for random stoichiometric matrices are exact."""
    raw_cls = _conservation_laws_from_stoichiometry(
        sp.Matrix(stoichiometric_matrix)
    )
    assert_exact_conservation_laws(stoichiometric_matrix, raw_cls)


@pytest.mark.parametrize(
    "transposed_stoichiometric_matrix",
    [
        [
            [0, 0, 0, 2, -1, 2, 0, 0, 0],
            [0, -2, -2, 0, 0, 0, 0, 0, 2],
            [0, -1, 2, 2, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, -2, 0, -1, 0, 0, 1, 0],
            [0, 0, 0, 0, -1, 0, 0, -2, 0],
            [0, 2, 0, 0, 0, 0, 0, 0, 0],
        ],
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
)
def test_conservation_laws_previous_failures(transposed_stoichiometric_matrix):
    """Stoichiometric matrices for which the previous floating-point
    implementation missed conservation laws or returned wrong ones
    (https://github.com/AMICI-dev/AMICI/issues/3366)."""
    S = sp.Matrix(transposed_stoichiometric_matrix).T
    raw_cls = _conservation_laws_from_stoichiometry(S)
    assert_exact_conservation_laws(S, raw_cls)


def test_conservation_laws_float_stoichiometry():
    """Floating-point stoichiometries are interpreted as their decimal
    representation."""
    # A -> 0.1 B
    S = sp.Matrix([[-1], [sp.Float(0.1)]])
    assert _conservation_laws_from_stoichiometry(S) == [(0, [0, 1], [1, 10])]

    # A -> 0.3 B; 0.1 B -> C
    S = sp.Matrix([[-1, 0], [sp.Float(0.3), -sp.Float(0.1)], [0, 1]])
    raw_cls = _conservation_laws_from_stoichiometry(S)
    assert raw_cls == [
        (0, [0, 1, 2], [1, sp.Rational(10, 3), sp.Rational(1, 3)])
    ]
    assert_exact_conservation_laws(
        [[-1, 0], [sp.Rational(3, 10), -sp.Rational(1, 10)], [0, 1]],
        raw_cls,
    )


def test_conservation_laws_none():
    """No conservation laws for full-rank stoichiometric matrices."""
    assert (
        _conservation_laws_from_stoichiometry(sp.Matrix([[1, 0], [0, -1]]))
        == []
    )


@skip_on_valgrind
def test_conservation_laws_demartino2014(data_demartino2014):
    """Conservation laws for the E. coli network from De Martino et al.
    (2014)"""
    S_list, row_names = data_demartino2014
    n_species = len(row_names)
    S = np.asarray(S_list, dtype=int).reshape(-1, n_species).T
    # (constructing a SparseMatrix from a dense array is slow)
    S_sparse = sp.SparseMatrix(
        *S.shape,
        {(int(i), int(j)): int(S[i, j]) for i, j in zip(*np.nonzero(S))},
    )
    raw_cls = _conservation_laws_from_stoichiometry(S_sparse)
    assert len(raw_cls) == 38

    # rows scaled to integers, verified with integer arithmetic
    cls = np.zeros((len(raw_cls), n_species), dtype=np.int64)
    for i_cl, (pivot, state_idxs, coefficients) in enumerate(raw_cls):
        assert coefficients[0] == 1 and state_idxs[0] == pivot
        scale = math.lcm(*(c.q for c in coefficients))
        cls[i_cl, state_idxs] = [int(c * scale) for c in coefficients]
    assert not (cls @ S).any()
    # pivot columns are unit vectors (up to scaling), so the rows are
    #  linearly independent
    pivots = [pivot for pivot, _, _ in raw_cls]
    assert np.count_nonzero(cls[:, pivots]) == len(pivots)


def _sbml_variable_stoichiometry() -> str:
    """SBML model for ``A -> sr B``, with stoichiometry ``sr`` given by an
    assignment rule"""
    doc = libsbml.SBMLDocument(3, 2)
    model = doc.createModel()
    c = model.createCompartment()
    c.setId("c")
    c.setSize(1)
    c.setConstant(True)
    for species_id in "AB":
        s = model.createSpecies()
        s.setId(species_id)
        s.setCompartment("c")
        s.setInitialAmount(1)
        s.setHasOnlySubstanceUnits(True)
        s.setBoundaryCondition(False)
        s.setConstant(False)
    r = model.createReaction()
    r.setId("r")
    r.setReversible(False)
    reactant = r.createReactant()
    reactant.setSpecies("A")
    reactant.setStoichiometry(1)
    reactant.setConstant(True)
    product = r.createProduct()
    product.setSpecies("B")
    product.setId("sr")
    product.setStoichiometry(1)
    product.setConstant(False)
    r.createKineticLaw().setMath(libsbml.parseL3Formula("A"))
    rule = model.createAssignmentRule()
    rule.setVariable("sr")
    rule.setMath(libsbml.parseL3Formula("1 + time"))
    return libsbml.writeSBMLToString(doc)


def test_conservation_laws_variable_stoichiometry():
    """No conservation laws for non-numeric stoichiometries."""
    importer = SbmlImporter(_sbml_variable_stoichiometry(), from_file=False)
    importer._process_sbml(hardcode_symbols=[])
    assert importer._stoichiometric_matrix.free_symbols

    with pytest.warns(UserWarning, match="parameterized stoichiometric"):
        assert importer._get_conservation_laws_rref() == []
