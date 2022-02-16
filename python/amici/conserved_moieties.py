import logging
import math
import random
import sys
from typing import List, MutableSequence, Sequence, Tuple

from .logging import get_logger

logger = get_logger(__name__, logging.ERROR)

# increase recursion limit for recursive quicksort
sys.setrecursionlimit(3000)


def compute_moiety_conservation_laws(
        stoichiometric_list: Sequence[float],
        num_species: int,
        num_reactions: int,
        max_num_monte_carlo: int = 20
) -> Tuple[List[List[int]], List[List[float]]]:
    """Compute moiety conservation laws.

    According to the algorithm proposed by De Martino et al. (2014)
    https://doi.org/10.1371/journal.pone.0100750

    :param stoichiometric_list:
        the stoichiometric matrix as a list (species x reactions,
        column-major ordering)
    :param num_species:
        total number of species in the reaction network
    :param num_reactions:
        total number of reactions in the reaction network
    :param max_num_monte_carlo:
        maximum number of MonteCarlo steps before changing to relaxation
    :returns:
        Integer MCLs as list of lists of indices of involved species and
        list of lists of corresponding coefficients.
    """
    # compute semi-positive conservation laws
    kernel_dim, engaged_species, int_kernel_dim, conserved_moieties, \
    cls_species_idxs, cls_coefficients = kernel(
        stoichiometric_list, num_species, num_reactions)
    done = (int_kernel_dim == kernel_dim)

    if not done:
        # construct interaction matrix
        J, J2, fields = fill(stoichiometric_list, engaged_species, num_species)

        timer = 0
        # maximum number of montecarlo search before starting relaxation
        while not done:
            yes, int_kernel_dim, conserved_moieties = monte_carlo(
                engaged_species, J, J2, fields, conserved_moieties,
                int_kernel_dim, cls_species_idxs, cls_coefficients,
                num_species, max_iter=max_num_monte_carlo
            )
            done = (int_kernel_dim == kernel_dim)
            timer = 0 if yes else timer + 1

            if timer == max_num_monte_carlo:
                done = relax(stoichiometric_list, conserved_moieties,
                             num_reactions, num_species)
                timer = 0
    reduce(int_kernel_dim, cls_species_idxs, cls_coefficients, num_species)

    return cls_species_idxs[:int_kernel_dim], cls_coefficients[:int_kernel_dim]


def _qsort(
        k: int,
        km: int,
        order: MutableSequence[int],
        pivots: Sequence[int]
) -> None:
    """Quicksort

    Recursive implementation of the quicksort algorithm

    :param k:
        number of elements to sort
    :param km:
        current center element
    :param order:
        ordering of the elements
    :param pivots:
        corresponding pivot elements from scaled partial pivoting strategy
    """
    # TODO: Rewrite into an iterative algorithm with pivoting strategy

    if k - km < 1:
        # nothing to sort
        return

    pivot = km + int((k - km) / 2)
    l = 0
    p = k - km - 1
    new_order = [0] * (k - km)
    for i in range(km, k):
        if i != pivot:
            if pivots[order[i]] < pivots[order[pivot]]:
                new_order[l] = order[i]
                l += 1
            else:
                new_order[p] = order[i]
                p -= 1
    new_order[p] = order[pivot]
    order[km:k] = new_order

    centre = p + km
    _qsort(k, centre + 1, order, pivots)
    _qsort(centre, km, order, pivots)


def kernel(
        stoichiometric_list: Sequence[float],
        num_species: int,
        num_reactions: int
) -> Tuple[float, List[int], int, List[int],
           List[List[int]], List[List[float]]]:
    """
    Kernel (left nullspace of :math:`S`) calculation by Gaussian elimination

    To compute the left nullspace of the stoichiometric matrix :math:`S`,
    a Gaussian elimination method with partial scaled pivoting is used to deal
    effectively with a possibly ill-conditioned stoichiometric matrix
    :math:`S`.

    Note that this is the Python reimplementation of the algorithm proposed
    by De Martino et al. (2014) https://doi.org/10.1371/journal.pone.0100750
    and thus a direct adaption of the original implementation in C/C++.

    :param stoichiometric_list:
        the stoichiometric matrix as a list (species x reactions,
        col-major ordering)
    :param num_species:
        total number of species in the reaction network
    :param num_reactions:
        total number of reactions in the reaction network
    :returns:
        kernel dimension, MCLs, integer kernel dimension, integer MCLs and
        indices to species and reactions in the preceding order as a tuple
    """
    MAX = 1e9
    MIN = 1e-9

    matrix: List[List[int]] = [[] for _ in range(num_species)]
    matrix2: List[List[float]] = [[] for _ in range(num_species)]

    i1 = 0
    j1 = 0
    for val in stoichiometric_list:
        if val != 0:
            matrix[j1].append(i1)
            matrix2[j1].append(val)
        j1 += 1
        if j1 == num_species:
            j1 = 0
            i1 += 1
    for i in range(num_species):
        matrix[i].append(num_reactions + i)
        matrix2[i].append(1)

    order: List[int] = list(range(num_species))
    pivots = [matrix[i][0] if len(matrix[i]) else MAX
              for i in range(num_species)]

    done = False
    while not done:
        _qsort(num_species, 0, order, pivots)
        for j in range(num_species - 1):
            if pivots[order[j + 1]] == pivots[order[j]] != MAX:
                min1 = 100000000
                if len(matrix[order[j]]) > 1:
                    for i in range(len(matrix[order[j]])):
                        min1 = min(min1, abs(matrix2[order[j]][0]
                                             / matrix2[order[j]][i]))

                min2 = 100000000
                if len(matrix[order[j + 1]]) > 1:
                    for i in range(len(matrix[order[j + 1]])):
                        min2 = min(min2, abs(matrix2[order[j + 1]][0]
                                             / matrix2[order[j + 1]][i]))

                if min2 > min1:
                    # swap
                    k2 = order[j + 1]
                    order[j + 1] = order[j]
                    order[j] = k2
        done = True

        for j in range(num_species - 1):
            if pivots[order[j + 1]] == pivots[order[j]] != MAX:
                k1 = order[j + 1]
                k2 = order[j]
                column = [0] * (num_species + num_reactions)
                g = matrix2[k2][0] / matrix2[k1][0]
                for i in range(1, len(matrix[k1])):
                    column[matrix[k1][i]] = matrix2[k1][i] * g

                for i in range(1, len(matrix[k2])):
                    column[matrix[k2][i]] -= matrix2[k2][i]

                matrix[k1] = []
                matrix2[k1] = []
                for col_idx, col_val in enumerate(column):
                    if abs(col_val) > MIN:
                        matrix[k1].append(col_idx)
                        matrix2[k1].append(col_val)

                done = False
                if len(matrix[order[j + 1]]):
                    pivots[order[j + 1]] = matrix[order[j + 1]][0]
                else:
                    pivots[order[j + 1]] = MAX

    RSolutions = [[] for _ in range(num_species)]
    RSolutions2 = [[] for _ in range(num_species)]
    kernel_dim = 0

    for i in range(num_species):
        done = all(matrix[i][j] >= num_reactions
                   for j in range(len(matrix[i])))
        if done and len(matrix[i]):
            for j in range(len(matrix[i])):
                RSolutions[kernel_dim].append(matrix[i][j] - num_reactions)
                RSolutions2[kernel_dim].append(matrix2[i][j])
            kernel_dim += 1
    del matrix, matrix2

    matched = []
    int_matched = []
    cls_species_idxs = [[] for _ in range(num_species)]
    cls_coefficients = [[] for _ in range(num_species)]

    i2 = 0
    for i in range(kernel_dim):
        ok2 = True
        for j in range(len(RSolutions[i])):
            if RSolutions2[i][j] * RSolutions2[i][0] < 0:
                ok2 = False
            if not matched or all(
                    cur_matched != RSolutions[i][j] for cur_matched in
                    matched
            ):
                matched.append(RSolutions[i][j])
        if ok2 and len(RSolutions[i]):
            min_value = MAX
            for j in range(len(RSolutions[i])):
                cls_species_idxs[i2].append(RSolutions[i][j])
                cls_coefficients[i2].append(abs(RSolutions2[i][j]))
                min_value = min(min_value, abs(RSolutions2[i][j]))
                if not int_matched or all(
                        cur_int_matched != cls_species_idxs[i2][j]
                        for cur_int_matched in int_matched
                ):
                    int_matched.append(cls_species_idxs[i2][j])
            for j in range(len(cls_species_idxs[i2])):
                cls_coefficients[i2][j] /= min_value
            i2 += 1
    int_kernel_dim = i2

    assert int_kernel_dim <= kernel_dim
    assert len(cls_species_idxs) == len(cls_coefficients), \
        "Inconsistent number of conserved quantities in coefficients and " \
        "species"
    return (kernel_dim, matched, int_kernel_dim, int_matched, cls_species_idxs,
            cls_coefficients)


def fill(
        stoichiometric_list: Sequence[float],
        matched: Sequence[int],
        num_species: int
) -> Tuple[List[List[int]], List[List[int]], List[int]]:
    """Construct interaction matrix

    Construct the interaction matrix out of the given stoichiometric matrix
    :math:`S`.

    :param stoichiometric_list:
        the stoichiometric matrix given as a flat list
    :param matched:
        found and independent moiety conservation laws (MCL)
    :param num_species:
        number of rows in :math:`S`
    :returns:
        interactions of metabolites and reactions, and matrix of interaction
    """
    dim = len(matched)
    MIN = 1e-9
    matrix = [[] for _ in range(dim)]
    matrix2 = [[] for _ in range(dim)]

    J = [[] for _ in range(num_species)]
    J2 = [[] for _ in range(num_species)]

    fields = [0] * num_species
    i1 = 0
    j1 = 0
    for val in stoichiometric_list:
        if val != 0:
            prendo = dim
            if dim > 0:
                for i in range(dim):
                    if j1 == matched[i]:
                        prendo = i
            if prendo < dim:
                matrix[prendo].append(i1)
                matrix2[prendo].append(val)
        j1 += 1
        if j1 == num_species:
            j1 = 0
            i1 += 1

    for i in range(dim):
        for j in range(i, dim):
            interactions = 0
            for po in range(len(matrix[i])):
                for pu in range(len(matrix[j])):
                    if matrix[i][po] == matrix[j][pu]:
                        interactions += matrix2[i][po] * matrix2[j][pu]
            if j == i:
                fields[i] = interactions
            elif abs(interactions) > MIN:
                J[i].append(j)
                J2[i].append(interactions)
                J[j].append(i)
                J2[j].append(interactions)
    return J, J2, fields


def _is_linearly_dependent(
        vector: Sequence[float],
        int_kernel_dim: int,
        cls_species_idxs: Sequence[Sequence[int]],
        cls_coefficients: Sequence[Sequence[float]],
        matched: Sequence[int],
        num_species: int
) -> bool:
    """Check for linear dependence between MCLs

    Check if the solutions found with Monte Carlo are linearly independent
    with respect to the previous found solution for all MCLs involved

    :param vector:
        found basis
    :param int_kernel_dim:
        number of integer conservative laws
    :param cls_species_idxs:
        NSolutions contains the species involved in the MCL
    :param cls_coefficients:
        NSolutions2 contains the corresponding coefficients in the MCL
    :param matched:
        actual found MCLs
    :param num_species:
        number of rows in :math:`S`
    :returns:
        boolean indicating linear dependence (true) or not (false)
    """
    K = int_kernel_dim + 1
    MIN = 1e-9
    MAX = 1e+9
    # TODO simplify - copy + []
    matrix: List[List[int]] = [[] for _ in range(K)]
    matrix2: List[List[float]] = [[] for _ in range(K)]
    for i in range(K - 1):
        for j in range(len(cls_species_idxs[i])):
            matrix[i].append(cls_species_idxs[i][j])
            matrix2[i].append(cls_coefficients[i][j])

    order2 = list(range(len(matched)))
    pivots2 = matched[:]
    _qsort(len(matched), 0, order2, pivots2)

    for i in range(len(matched)):
        if vector[order2[i]] > MIN:
            matrix[K - 1].append(matched[order2[i]])
            matrix2[K - 1].append(vector[order2[i]])

    order = list(range(K))
    pivots = [matrix[i][0] if len(matrix[i]) else MAX for i in range(K)]

    ok = False
    while not ok:
        _qsort(K, 0, order, pivots)
        for j in range(K - 1):
            if pivots[order[j + 1]] == pivots[order[j]] != MAX:
                min1 = MAX
                if len(matrix[order[j]]) > 1:
                    for i in range(len(matrix[order[j]])):
                        min1 = min(min1, abs(matrix2[order[j]][0]
                                             / matrix2[order[j]][i]))
                min2 = MAX
                if len(matrix[order[j + 1]]) > 1:
                    for i in range(len(matrix[order[j + 1]])):
                        min2 = min(min2, abs(matrix2[order[j + 1]][0]
                                             / matrix2[order[j + 1]][i]))
                if min2 > min1:
                    # swap
                    k2 = order[j + 1]
                    order[j + 1] = order[j]
                    order[j] = k2
        ok = True
        for j in range(K - 1):
            if pivots[order[j + 1]] == pivots[order[j]] != MAX:
                k1 = order[j + 1]
                k2 = order[j]
                column = [0] * num_species
                g = matrix2[k2][0] / matrix2[k1][0]
                for i in range(1, len(matrix[k1])):
                    column[matrix[k1][i]] = matrix2[k1][i] * g
                for i in range(1, len(matrix[k2])):
                    column[matrix[k2][i]] -= matrix2[k2][i]

                matrix[k1] = []
                matrix2[k1] = []
                for i in range(num_species):
                    if abs(column[i]) > MIN:
                        matrix[k1].append(i)
                        matrix2[k1].append(column[i])
                ok = False
                pivots[k1] = matrix[k1][0] if len(matrix[k1]) else MAX
    K1 = sum(len(matrix[i]) > 0 for i in range(K))
    return K == K1


def monte_carlo(
        matched: Sequence[int],
        J: Sequence[Sequence[int]],
        J2: Sequence[Sequence[float]],
        fields: Sequence[float],
        int_matched: MutableSequence[int],
        int_kernel_dim: int,
        cls_species_idxs: MutableSequence[MutableSequence[int]],
        cls_coefficients: MutableSequence[MutableSequence[float]],
        num_species: int,
        initial_temperature: float = 1,
        cool_rate: float = 1e-3,
        max_iter: int = 10
) -> Tuple[bool, int, Sequence[int]]:
    """MonteCarlo simulated annealing for finding integer MCLs

    Finding integer solutions for the MCLs by Monte Carlo, see step (b) in
    the De Martino (2014) paper and Eqs. 11-13 in the publication

    :param matched:
        matched
    :param J:
        J
    :param J2:
        J2
    :param fields:
        fields
    :param int_matched:
        actual matched MCLs
    :param int_kernel_dim:
        number of MCLs found in :math:`S`
    :param cls_species_idxs:
        Modified in-place.
    :param cls_coefficients:
        Modified in-place.
    :param initial_temperature:
        initial temperature
    :param cool_rate:
        cooling rate of simulated annealing
    :param max_iter:
        maximum number of MonteCarlo steps before changing to relaxation
    :returns:
        status of MC iteration, number of integer MCLs, number of MCLs,
        metabolites and reaction indices, MCLs and integer MCLs as a tuple
    """
    # TODO: doc: what does value of status indicate

    MIN = 1e-9
    dim = len(matched)
    num = [int(2 * random.uniform(0, 1)) if len(J[i]) else 0
           for i in range(dim)]
    numtot = sum(num)

    H = 0
    for i in range(dim):
        H += fields[i] * num[i] ** 2
        for j in range(len(J[i])):
            H += J2[i][j] * num[i] * num[J[i][j]]

    count = 0
    howmany = 0
    T1 = initial_temperature
    e = math.exp(-1 / T1)
    while True:
        en = int(random.uniform(0, 1) * dim)
        while not len(J[en]):
            en = int(random.uniform(0, 1) * dim)

        p = -1 if num[en] > 0 and random.uniform(0, 1) < 0.5 else 1
        delta = fields[en] * num[en]
        for i in range(len(J[en])):
            delta += J2[en][i] * num[J[en][i]]
        delta = 2 * p * delta + fields[en]

        if delta < 0 or random.uniform(0, 1) < math.pow(e, delta):
            num[en] += p
            numtot += p
            H += delta

        count += 1
        if count % dim == 0:
            T1 -= cool_rate
            if T1 <= 0:
                T1 = cool_rate
            e = math.exp(-1 / T1)

        if count == dim // cool_rate:
            count = 0
            T1 = initial_temperature
            e = math.exp(-1 / T1)
            en = int(random.uniform(0, 1) * dim)
            while not len(J[en]):
                en = int(random.uniform(0, 1) * dim)
            num = [0] * dim
            num[en] = 1
            numtot = 1

            H = 0
            for i in range(dim):
                H += fields[i] * num[i] ** 2
                for j in range(len(J[i])):
                    H += J2[i][j] * num[i] * num[J[i][j]]
            howmany += 1

        if (H < MIN and numtot > 0) or (howmany == 10 * max_iter):
            break

    if howmany < 10 * max_iter:
        # founds MCLS? need to check for linear independence
        if len(int_matched):
            yes = _is_linearly_dependent(
                num, int_kernel_dim, cls_species_idxs,
                cls_coefficients, matched, num_species)
        else:
            yes = True

        # reduce by MC procedure
        if yes:
            order2 = list(range(len(matched)))
            pivots2 = matched[:]
            _qsort(len(matched), 0, order2, pivots2)
            for i in range(len(matched)):
                if num[order2[i]] > 0:
                    cls_species_idxs[int_kernel_dim].append(
                        matched[order2[i]])
                    cls_coefficients[int_kernel_dim].append(num[order2[i]])
            int_kernel_dim += 1
            reduce(int_kernel_dim, cls_species_idxs, cls_coefficients,
                   num_species)
            min_value = 1000
            for i in range(len(cls_species_idxs[int_kernel_dim - 1])):
                if not len(int_matched) \
                        or all(cur_int_matched
                               != cls_species_idxs[int_kernel_dim - 1][i]
                               for cur_int_matched in int_matched):
                    int_matched.append(cls_species_idxs[int_kernel_dim - 1][i])

                min_value = min(min_value,
                                cls_coefficients[int_kernel_dim - 1][i])
            for i in range(len(cls_species_idxs[int_kernel_dim - 1])):
                cls_coefficients[int_kernel_dim - 1][i] /= min_value
            logger.debug(
                f"Found linearly independent moiety, now there are "
                f"{int_kernel_dim} engaging {len(int_matched)} species")
        else:
            logger.debug(
                "Found a moiety but it is linearly dependent... next.")
    else:
        yes = False
    return yes, int_kernel_dim, int_matched


def relax(
        stoichiometric_list: Sequence[float],
        int_matched: Sequence[int],
        num_reactions: int,
        num_species: int,
        relaxation_max: float = 1e6,
        relaxation_step: float = 1.9
) -> bool:
    """Relaxation scheme for Monte Carlo final solution

    Checking for completeness using Motzkin's theorem. See Step (c) in
    De Martino (2014) and the Eqs. 14-16 in the corresponding publication

    :param stoichiometric_list:
        stoichiometric matrix :math:`S` as a flat list (column-major ordering)
    :param int_matched:
        intmatched
    :param num_reactions:
        number of reactions in reaction network
    :param num_species:
        number of species in reaction network
    :param relaxation_max:
        maximum relaxation step
    :param relaxation_step:
        relaxation step width
    :returns:
        boolean indicating if relaxation has succeeded (``True` ) or not
        (``False``)
    """
    MIN = 1e-9
    MAX = 1e9

    i1 = 0
    j1 = 0
    K = len(int_matched)
    matrix: List[List[int]] = [[] for _ in range(K)]
    matrix2: List[List[float]] = [[] for _ in range(K)]

    for val in stoichiometric_list:
        if val != 0:
            prendo = K
            if K > 0:
                for i in range(K):
                    if j1 == int_matched[i]:
                        prendo = i
            if prendo < K:
                matrix[prendo].append(i1)
                matrix2[prendo].append(val)
        j1 += 1
        if j1 == num_species:
            j1 = 0
            i1 += 1

    # reducing the stoichiometric matrix of conserved moieties to row echelon
    # form by Gaussian elimination
    order = list(range(K))
    pivots = [matrix[i][0] if len(matrix[i]) else MAX for i in range(K)]
    done = False
    while not done:
        _qsort(K, 0, order, pivots)
        for j in range(K - 1):
            if pivots[order[j + 1]] == pivots[order[j]] != MAX:
                min1 = MAX
                if len(matrix[order[j]]) > 1:
                    for i in range(len(matrix[order[j]])):
                        min1 = min(min1, abs(matrix2[order[j]][0]
                                             / matrix2[order[j]][i]))
                min2 = MAX
                if len(matrix[order[j + 1]]) > 1:
                    for i in range(len(matrix[order[j + 1]])):
                        min2 = min(min2, abs(matrix2[order[j + 1]][0]
                                             / matrix2[order[j + 1]][i]))
                if min2 > min1:
                    # swap
                    k2 = order[j + 1]
                    order[j + 1] = order[j]
                    order[j] = k2
        done = True
        for j in range(K - 1):
            if pivots[order[j + 1]] == pivots[order[j]] != MAX:
                k1 = order[j + 1]
                k2 = order[j]
                column = [0] * num_reactions
                g = matrix2[k2][0] / matrix2[k1][0]
                for i in range(1, len(matrix[k1])):
                    column[matrix[k1][i]] = matrix2[k1][i] * g
                for i in range(1, len(matrix[k2])):
                    column[matrix[k2][i]] -= matrix2[k2][i]

                matrix[k1] = []
                matrix2[k1] = []
                for col_idx, col_val in enumerate(column):
                    if abs(col_val) > MIN:
                        matrix[k1].append(col_idx)
                        matrix2[k1].append(col_val)
                done = False
                if len(matrix[order[j + 1]]):
                    pivots[order[j + 1]] = matrix[order[j + 1]][0]
                else:
                    pivots[order[j + 1]] = MAX

    for i in range(K):
        if len(matrix[i]):
            norm = matrix2[i][0]
            for j in range(len(matrix[i])):
                matrix2[i][j] /= norm

    for k1 in reversed(range(K - 1)):
        k = order[k1]
        if len(matrix[k]) <= 1:
            continue

        i = 0
        while i < len(matrix[k]):
            for j1 in range(k1 + 1, K):
                j = order[j1]
                print("matrix[j]", matrix[j])
                print("matrix[k]", matrix[k], i)
                if not len(matrix[j]) or matrix[j][0] != matrix[k][i]:
                    continue

                row_k = [0] * num_reactions
                for a in range(len(matrix[k])):
                    row_k[matrix[k][a]] = matrix2[k][a]
                for a in range(len(matrix[j])):
                    row_k[matrix[j][a]] -= matrix2[j][a] * matrix2[k][i]
                matrix[k] = []
                matrix2[k] = []
                for row_idx, row_val in enumerate(row_k):
                    if row_val != 0:
                        matrix[k].append(row_idx)
                        matrix2[k].append(row_val)
            i += 1

    indip = [K + 1] * num_reactions
    for i in range(K):
        if len(matrix[i]):
            indip[matrix[i][0]] = i
    M1 = 0
    for i in range(num_reactions):
        if indip[i] == K + 1:
            indip[i] = K + M1
            M1 += 1

    matrixAus = [[] for _ in range(M1)]
    matrixAus2 = [[] for _ in range(M1)]
    i1 = 0
    for i in range(num_reactions):
        if indip[i] >= K:
            matrixAus[i1].append(i)
            matrixAus2[i1].append(1)
            i1 += 1
        else:
            t = indip[i]
            if len(matrix[t]) > 1:
                for k in range(1, len(matrix[t])):
                    quelo = indip[matrix[t][k]] - K
                    matrixAus[quelo].append(i)
                    matrixAus2[quelo].append(-matrix2[t][k])
    del matrix

    N1 = num_species - K
    matrix_aus = [[] for _ in range(N1)]
    matrix_aus2 = [[] for _ in range(N1)]
    k1 = 0
    i1 = 0
    j1 = 0
    for val in stoichiometric_list:
        prendo = 1
        for i in range(len(int_matched)):
            if j1 == int_matched[i]:
                prendo -= 1
        if val != 0 and prendo == 1:
            matrix_aus[k1].append(i1)
            matrix_aus2[k1].append(val)
        j1 += 1
        k1 += prendo
        if j1 == num_species:
            j1 = 0
            k1 = 0
            i1 += 1

    matrixb = [[] for _ in range(N1)]
    matrixb2 = [[] for _ in range(N1)]
    for i in range(M1):
        for j in range(N1):
            if len(matrix_aus[j]) * len(matrixAus[i]):
                prod = 0
                for ib in range(len(matrixAus[i])):
                    for jb in range(len(matrix_aus[j])):
                        if matrixAus[i][ib] == matrix_aus[j][jb]:
                            prod += matrixAus2[i][ib] * matrix_aus2[j][jb]
                if abs(prod) > MIN:
                    matrixb[j].append(i)
                    matrixb2[j].append(prod)
    del matrixAus, matrixAus2, matrix_aus, matrix_aus2

    var = [MIN] * M1
    time = 0
    cmin_idx = 0
    while True:
        cmin = 1000
        for j in range(N1):
            constr = 0
            if len(matrixb[j]):
                for i in range(len(matrixb[j])):
                    constr += matrixb2[j][i] * var[matrixb[j][i]]
                if constr < cmin:
                    cmin_idx = j
                    cmin = constr
        if cmin >= 0:
            # constraints satisfied
            break

        # Motzkin relaxation
        alpha = -relaxation_step * cmin
        fact = sum(val ** 2 for val in matrixb2[cmin_idx])
        alpha /= fact
        alpha = max(1e-9 * MIN, alpha)
        for j in range(len(matrixb[cmin_idx])):
            var[matrixb[cmin_idx][j]] += alpha * matrixb2[cmin_idx][j]

        time += 1
        if time >= relaxation_max:
            # timeout
            break

    return done


def reduce(
        int_kernel_dim: int,
        cls_species_idxs: MutableSequence[MutableSequence[int]],
        cls_coefficients: MutableSequence[MutableSequence[float]],
        num_species: int
) -> None:
    """Reducing the solution which has been found by the Monte Carlo process

    In case of superpositions of independent MCLs one can reduce by
    iteratively subtracting the other independent MCLs, taking care
    to maintain then non-negativity constraint, see Eq. 13 in De Martino (2014)

    :param int_kernel_dim:
        number of found MCLs
    :param cls_species_idxs:
        Species indices involved in each of the conservation laws.
        Modified in-place.
    :param cls_coefficients:
        Coefficients for each of the species involved in each of the
        conservation laws. Modified in-place.
    :param num_species:
        number of species / rows in :math:`S`
    """
    K = int_kernel_dim
    MIN = 1e-9
    order = list(range(K))
    pivots = [-len(cls_species_idxs[i]) for i in range(K)]

    done = False
    while not done:
        _qsort(K, 0, order, pivots)
        done = True
        for i in range(K - 1):
            for j in range(i + 1, K):
                k1 = order[i]
                k2 = order[j]
                column = [0] * num_species
                ok1 = True
                for species_idx, coefficient \
                        in zip(cls_species_idxs[k1], cls_coefficients[k1]):
                    column[species_idx] = coefficient
                for species_idx, coefficient \
                        in zip(cls_species_idxs[k2], cls_coefficients[k2]):
                    column[species_idx] -= coefficient
                    if column[species_idx] < -MIN:
                        ok1 = False
                if ok1:
                    done = False
                    cls_species_idxs[k1] = []
                    cls_coefficients[k1] = []
                    for col_idx, col_val in enumerate(column):
                        if abs(col_val) > MIN:
                            cls_species_idxs[k1].append(col_idx)
                            cls_coefficients[k1].append(col_val)
                    pivots[k1] = -len(cls_species_idxs[k1])
