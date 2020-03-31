"""Miscellaneous AMICI Python interface tests"""

import amici
import sympy as sp


def test_csc_matrix():
    """Test sparse CSC matrix creation"""
    matrix = sp.Matrix([[1, 0], [2, 3]])
    symbol_col_ptrs, symbol_row_vals, sparse_list, symbol_list, sparse_matrix \
        = amici.ode_export.csc_matrix(matrix, 'a')

    assert symbol_col_ptrs == [0, 2, 3]
    assert symbol_row_vals == [0, 1, 1]
    assert sparse_list == sp.Matrix([[1], [2], [3]])
    assert symbol_list == ['a0', 'a1', 'a2']
    assert str(sparse_matrix) == 'Matrix([[a0, 0], [a1, a2]])'


def test_csc_matrix_empty():
    """Test sparse CSC matrix creation for empty matrix"""
    matrix = sp.Matrix()
    symbol_col_ptrs, symbol_row_vals, sparse_list, symbol_list, sparse_matrix \
        = amici.ode_export.csc_matrix(matrix, 'a')

    assert symbol_col_ptrs == []
    assert symbol_row_vals == []
    assert sparse_list == sp.Matrix(0, 0, [])
    assert symbol_list == []
    assert str(sparse_matrix) == 'Matrix(0, 0, [])'


def test_csc_matrix_vector():
    """Test sparse CSC matrix creation from matrix slice"""

    matrix = sp.Matrix([[1, 0], [2, 3]])
    symbol_col_ptrs, symbol_row_vals, sparse_list, symbol_list, sparse_matrix \
        = amici.ode_export.csc_matrix(matrix[:, 0], 'a')

    assert symbol_col_ptrs == [0, 2]
    assert symbol_row_vals == [0, 1]
    assert sparse_list == sp.Matrix([[1], [2]])
    assert symbol_list == ['a0', 'a1']
    assert str(sparse_matrix) == 'Matrix([[a0], [a1]])'

    # Test continuation of numbering of symbols
    symbol_col_ptrs, symbol_row_vals, sparse_list, symbol_list, sparse_matrix \
        = amici.ode_export.csc_matrix(matrix[:, 1], 'a',
                                    base_index=len(symbol_list))

    assert symbol_col_ptrs == [0, 1]
    assert symbol_row_vals == [1]
    assert sparse_list == sp.Matrix([[3]])
    assert symbol_list == ['a2']
    assert str(sparse_matrix) == 'Matrix([[0], [a2]])'
