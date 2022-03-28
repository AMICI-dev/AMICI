"""Find conserved quantities deterministically"""

from typing import List, Literal, Optional, Union

import numpy as np


def rref(
        mat: np.array,
        round_ndigits: Optional[Union[Literal[False], int]] = None
) -> np.array:
    """
    Bring matrix ``mat`` to reduced row echelon form

    see https://en.wikipedia.org/wiki/Row_echelon_form

    :returns: ``mat`` in rref form.
    """
    # Rounding function
    if round_ndigits is False:
        # no-op
        def _round(mat):
            return mat
    else:
        if round_ndigits is None:
            # drop least significant digit (more or less)
            round_ndigits = - int(np.ceil(np.log10(np.spacing(1))))

        def _round(mat):
            mat = np.round(mat, round_ndigits)
            mat[np.abs(mat) <= 10**(-round_ndigits)] = 0
            return mat

    # create a copy  that will be modified
    mat = mat.copy()
    lead = 0
    n_rows, n_columns = mat.shape
    for r in range(n_rows):
        if n_columns <= lead:
            return mat

        i = r
        while mat[i, lead] == 0:
            i += 1
            if n_rows == i:
                i = r
                lead += 1
                if n_columns == lead:
                    return mat

        if i != r:
            # Swap rows
            mat[[i, r]] = mat[[r, i]]
        # Divide row
        mat[r] /= mat[r, lead]
        for i in range(n_rows):
            if i != r:
                # Subtract multiple
                mat[i] -= mat[i, lead] * mat[r]
                # mat[i, np.abs(mat[i]) < 1e-12] = 0
        mat = _round(mat)
        lead += 1
    return mat


def pivots(mat: np.array) -> List[int]:
    """Get indices of pivot columns in ``mat``, assumed to be in reduced row
    echelon form"""
    pivot_cols = []
    last_pivot_col = -1
    for i in range(mat.shape[0]):
        for j in range(last_pivot_col + 1, mat.shape[1]):
            if mat[i, j] != 0:
                pivot_cols.append(j)
                last_pivot_col = j
                break
    return pivot_cols


def nullspace_by_rref(mat: np.array):
    """Compute basis of the nullspace of ``mat`` based on the reduced row
    echelon form"""
    rref_mat = rref(mat)
    pivot_cols = pivots(rref_mat)
    rows, cols = mat.shape

    basis = []
    nonpivots = []
    for i in range(cols):
        if i in pivot_cols:
            continue
        nonpivots.append(i)
        vec = [1.0 if i == j else 0.0 for j in range(cols)]
        for ii, jj in enumerate(pivot_cols):
            vec[jj] -= rref_mat[ii][i]
        basis.append(vec)
    return np.array(basis)


