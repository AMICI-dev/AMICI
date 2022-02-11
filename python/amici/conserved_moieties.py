import logging
import math
import random
import sys
from numbers import Number
from typing import List, MutableSequence, Sequence, Tuple

from .logging import get_logger

sys.setrecursionlimit(3000)
logger = get_logger(__name__, logging.ERROR)


def compute_moiety_conservation_laws(
        stoichiometric_list: Sequence[Number],
        num_species: int,
        num_reactions: int
) -> Tuple[List[List[int]], List[List[Number]]]:
    """Compute moiety conservation laws.

    According to the algorithm proposed by De Martino et al. (2014)
    https://doi.org/10.1371/journal.pone.0100750

    :param stoichiometric_list:
        the stoichiometric matrix as a list (species x reactions,
        row-major ordering)
    :param num_species:
        total number of species in the reaction network
    :param num_reactions:
        total number of reactions in the reaction network
    :returns:
        Integer MCLs as list of lists of indices of involved species and
        list of lists of corresponding coefficients.
    """
    # compute semi-positive conservation laws
    kernel_dim, engaged_species, int_kernel_dim, conserved_moieties, \
        cls_species_idxs, cls_coefficients = kernel(
        stoichiometric_list, num_species, num_reactions)

    # construct interaction matrix
    J, J2, fields = fill(stoichiometric_list, engaged_species, num_species)

    done = (int_kernel_dim == kernel_dim)
    timer = 0
    # maximum number of montecarlo search before starting relaxation
    max_num_monte_carlo = 10
    while not done:
        yes, int_kernel_dim, engaged_species, conserved_moieties = monte_carlo(
            engaged_species, J, J2, fields, conserved_moieties,
            int_kernel_dim, cls_species_idxs, cls_coefficients, num_species,
            max_iter=max_num_monte_carlo
        )
        done = (int_kernel_dim == kernel_dim)
        if yes:
            timer += 1
        else:
            timer = 0

        if timer == max_num_monte_carlo:
            done = relax(stoichiometric_list, conserved_moieties,
                         num_reactions, num_species)
            if not done:
                timer = 0
    reduce(int_kernel_dim, cls_species_idxs, cls_coefficients, num_species)
    return cls_species_idxs, cls_coefficients


def _qsort(
        k: int,
        km: int,
        orders: MutableSequence[int],
        pivots: Sequence[int]
) -> None:
    """Quicksort

    Recursive implementation of the quicksort algorithm
    TODO: Rewrite into an iterative algorithm with pivoting strategy

    :param k:
        number of elements to sort
    :param km:
        current center element
    :param orders:
        ordering of the elements
    :param pivots:
        corresponding pivot elements from scaled partial pivoting strategy
    """
    if k - km < 1:
        # nothing to do
        return

    pivot = km + int((k - km) / 2)
    l = 0
    p = k - km - 1
    neworders = [None] * (k - km)
    for i in range(km, k):
        if i != pivot:
            if pivots[orders[i]] < pivots[orders[pivot]]:
                neworders[l] = orders[i]
                l += 1
            else:
                neworders[p] = orders[i]
                p -= 1
    neworders[p] = orders[pivot]
    for i in range(km, k):
        orders[i] = neworders[i - km]

    centre = p + km
    _qsort(k, centre + 1, orders, pivots)
    _qsort(centre, km, orders, pivots)


def kernel(
        stoichiometric_list: Sequence[Number],
        num_species: int,
        num_reactions: int
) -> Tuple[Number, List[int], int, List[int],
           List[List[int]], List[List[Number]]]:
    """
    Kernel (left nullspace of :math:`S`) calculation by Gaussian elimination

    To compute the left nullspace of the stoichiometric matrix :math:`S`,
    a Gaussian elimination method with partial scaled pivoting is used to deal
    effectively with a possibly ill-conditioned stoichiometric matrix :math:`S`.

    Note that this is the Python reimplementation of the algorithm proposed
    by De Martino et al. (2014) https://doi.org/10.1371/journal.pone.0100750
    and thus a direct adaption of the original implementation in C/C++.

    :param stoichiometric_list:
        the stoichiometric matrix as a list (species x reactions,
        row-major ordering)
    :param num_species:
        total number of species in the reaction network
    :param num_reactions:
        total number of reactions in the reaction network
    :returns:
        kernel dimension, MCLs, integer kernel dimension, integer MCLs and
        indices to species and reactions in the preceding order as a tuple
    """
    il = 0
    jl = 0
    N = num_species
    M = num_reactions
    MAX = 1e9
    MIN = 1e-9

    matrix = [[] for _ in range(N)]
    matrix2 = [[] for _ in range(N)]
    matched = []
    int_matched = []
    cls_species_idxs = [[] for _ in range(N)]
    cls_coefficients = [[] for _ in range(N)]

    for _, val in enumerate(stoichiometric_list):
        if val != 0:
            matrix[jl].append(il)
            matrix2[jl].append(val)
        jl += 1
        if jl == N:
            jl = 0
            il += 1

    for i in range(N):
        matrix[i].append(M + i)
        matrix2[i].append(1)

    orders = list(range(N))
    pivots = [matrix[i][0] if len(matrix[i]) > 0 else MAX for i in range(N)]

    done = False
    while not done:
        _qsort(N, 0, orders, pivots)
        for j in range(N - 1):
            if pivots[orders[j + 1]] == pivots[orders[j]] \
                    and pivots[orders[j]] != MAX:
                min1 = 100000000
                if len(matrix[orders[j]]) > 1:
                    for i in range(len(matrix[orders[j]])):
                        if abs(matrix2[orders[j]][0]
                               / matrix2[orders[j]][i]) < min1:
                            min1 = abs(
                                matrix2[orders[j]][0] / matrix2[orders[j]][i])

                min2 = 100000000
                if len(matrix[orders[j + 1]]) > 1:
                    for i in range(len(matrix[orders[j + 1]])):
                        if abs(matrix2[orders[j + 1]][0] /
                               matrix2[orders[j + 1]][i]) < min2:
                            min2 = abs(matrix2[orders[j + 1]][0] /
                                       matrix2[orders[j + 1]][i])

                if min2 > min1:
                    k2 = orders[j + 1]
                    orders[j + 1] = orders[j]
                    orders[j] = k2
        done = 1

        for j in range(N - 1):
            if pivots[orders[j + 1]] == pivots[orders[j]] \
                    and pivots[orders[j]] != MAX:
                k1 = orders[j + 1]
                k2 = orders[j]
                column = [0 for _ in range(N + M)]
                g = matrix2[k2][0] / matrix2[k1][0]
                for i in range(1, len(matrix[k1])):
                    column[matrix[k1][i]] = matrix2[k1][i] * g

                for i in range(1, len(matrix[k2])):
                    column[matrix[k2][i]] -= matrix2[k2][i]

                matrix[k1] = []
                matrix2[k1] = []
                for i in range(N + M):
                    if abs(column[i]) > MIN:
                        matrix[k1].append(i)
                        matrix2[k1].append(column[i])

                done = 0
                if len(matrix[orders[j + 1]]) > 0:
                    pivots[orders[j + 1]] = matrix[orders[j + 1]][0]
                else:
                    pivots[orders[j + 1]] = MAX

    RSolutions = [[] for _ in range(N)]
    RSolutions2 = [[] for _ in range(N)]
    kernel_dim = 0

    for i in range(N):
        done = True
        if len(matrix[i]):
            for j in range(len(matrix[i])):
                if matrix[i][j] < M:
                    done = False
                    break
        if done and len(matrix[i]):
            for j in range(len(matrix[i])):
                RSolutions[kernel_dim].append(matrix[i][j] - M)
                RSolutions2[kernel_dim].append(matrix2[i][j])
            kernel_dim += 1

    i2 = 0
    for i in range(kernel_dim):
        ok2 = 1
        if (len(RSolutions[i])) > 0:
            for j in range(len(RSolutions[i])):
                if RSolutions2[i][j] * RSolutions2[i][0] < 0:
                    ok2 = 0
                if len(matched) == 0:
                    matched.append(RSolutions[i][j])
                elif all(cur_matched != RSolutions[i][j]
                         for cur_matched in matched):
                    matched.append(RSolutions[i][j])
        if ok2 == 1 and len(RSolutions[i]) > 0:
            min = MAX
            for j in range(len(RSolutions[i])):
                cls_species_idxs[i2].append(RSolutions[i][j])
                cls_coefficients[i2].append(abs(RSolutions2[i][j]))
                if min > abs(RSolutions2[i][j]):
                    min = abs(RSolutions2[i][j])
                if len(int_matched) == 0:
                    int_matched.append(cls_species_idxs[i2][j])
                elif all(cur_int_matched != cls_species_idxs[i2][j]
                         for cur_int_matched in int_matched):
                    int_matched.append(cls_species_idxs[i2][j])
            for j in range(len(cls_species_idxs[i2])):
                cls_coefficients[i2][j] /= min
            i2 += 1
    int_kernel_dim = i2

    assert int_kernel_dim <= kernel_dim
    assert len(cls_species_idxs) == len(cls_coefficients), \
        "Inconsistent number of conserved quantities in coefficients and " \
        "species"
    return (kernel_dim, matched, int_kernel_dim, int_matched, cls_species_idxs,
            cls_coefficients)


def fill(
        stoichiometric_list: Sequence[Number],
        matched: Sequence[int],
        num_rows: int
) -> Tuple[List[List[int]], List[List[int]], List[int]]:
    """Construct interaction matrix

    Construct the interaction matrix out of the given stoichiometric matrix
    :math:`S`.

    :param stoichiometric_list:
        the stoichiometric matrix given as a flat list
    :param matched:
        found and independent moiety conservation laws (MCL)
    :param num_rows:
        number of rows in :math:`S`
    :returns:
        interactions of metabolites and reactions, and matrix of interaction
    """
    dim = len(matched)
    MIN = 1e-9
    matrix = [[] for _ in range(dim)]
    matrix2 = [[] for _ in range(dim)]

    J = [[] for _ in range(num_rows)]
    J2 = [[] for _ in range(num_rows)]

    fields = [0] * num_rows
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
        if j1 == num_rows:
            j1 = 0
            i1 += 1

    for i in range(dim):
        for j in range(i, dim):
            if len(matrix[i]) > 0:
                for po in range(len(matrix[i])):
                    interactions = 0
                    if len(matrix[j]) > 0:
                        for pu in range(len(matrix[j])):
                            if matrix[i][po] == matrix[j][pu]:
                                interactions += (
                                        matrix2[i][po] * matrix2[j][pu])
                    if j == i:
                        fields[i] = interactions
                    elif abs(interactions) > MIN:
                        J[i].append(j)
                        J2[i].append(interactions)
                        J[j].append(i)
                        J2[j].append(interactions)
    return J, J2, fields


def _is_linearly_dependent(
        vectors: Sequence[Number],
        int_kernel_dim: int,
        cls_species_idxs: Sequence[Sequence[int]],
        cls_coefficients: Sequence[Sequence[Number]],
        matched: Sequence[int],
        num_rows: int
) -> bool:
    """Check for linear dependence between MCLs

    Check if the solutions found with Monte Carlo are linearly independent
    with respect to the previous found solution for all MCLs involved

    :param vectors:
        found basis
    :param int_kernel_dim:
        number of integer conservative laws
    :param cls_species_idxs:
        NSolutions contains the species involved in the MCL
    :param cls_coefficients:
        NSolutions2 contains the corresponding coefficients in the MCL
    :param matched:
        actual found MCLs
    :param num_rows:
        number of rows in :math:`S`
    :returns:
        boolean indicating linear dependence (true) or not (false)
    """
    K = int_kernel_dim + 1
    MIN = 1e-9
    MAX = 1e+9
    matrix = [[] for _ in range(K)]
    matrix2 = [[] for _ in range(K)]
    for i in range(K - 1):
        for j in range(len(cls_species_idxs[i])):
            matrix[i].append(cls_species_idxs[i][j])
            matrix2[i].append(cls_coefficients[i][j])

    orders2 = list(range(len(matched)))
    pivots2 = matched[:]

    _qsort(len(matched), 0, orders2, pivots2)
    for i in range(len(matched)):
        if vectors[orders2[i]] > MIN:
            matrix[K - 1].append(matched[orders2[i]])
            matrix2[K - 1].append(float(vectors[orders2[i]]))

    ok = 0
    orders = list(range(K))

    pivots = [matrix[i][0] if len(matrix[i]) else MAX for i in range(K)]

    while ok == 0:
        _qsort(K, 0, orders, pivots)
        for j in range(K - 1):
            if pivots[orders[j + 1]] == pivots[orders[j]] != MAX:
                min1 = MAX
                if len(matrix[orders[j]]) > 1:
                    for i in range(len(matrix[orders[j]])):
                        if (abs(matrix2[orders[j]][0]
                                / matrix2[orders[j]][i])) < min1:
                            min1 = abs(matrix2[orders[j]][0]
                                       / matrix2[orders[j]][i])
                min2 = MAX
                if len(matrix[orders[j + 1]]) > 1:
                    for i in range(len(matrix[orders[j + 1]])):
                        if (abs(matrix2[orders[j + 1]][0] /
                                matrix2[orders[j + 1]][i])) < min2:
                            min2 = abs(matrix2[orders[j + 1]][0] /
                                       matrix2[orders[j + 1]][i])
                if min2 > min1:
                    k2 = orders[j + 1]
                    orders[j + 1] = orders[j]
                    orders[j] = k2
        ok = 1
        for j in range(K - 2):
            if pivots[orders[j + 1]] == pivots[orders[j]] != MAX:
                k1 = orders[j + 1]
                k2 = orders[j]
                column = [0] * num_rows
                g = matrix2[k2][0] / matrix2[k1][0]
                for i in range(1, len(matrix[k1])):
                    column[matrix[k1][i]] = matrix2[k1][i] * g
                for i in range(1, len(matrix[k2])):
                    column[matrix[k2][i]] -= matrix2[k2][i]

                matrix[k1] = []
                matrix2[k1] = []
                for i in range(num_rows):
                    if abs(column[i]) > MIN:
                        matrix[k1].append(i)
                        matrix2[k1].append(column[i])
                ok = 0
                pivots[k1] = matrix[k1][0] if len(matrix[k1]) > 0 else MAX
    K1 = sum(len(matrix[i]) > 0 for i in range(K))
    return K == K1


def monte_carlo(
        matched,
        J,
        J2,
        fields,
        int_matched: Sequence[int],
        int_kernel_dim: int,
        cls_species_idxs: Sequence[Sequence[int]],
        cls_coefficients: Sequence[Sequence[Number]],
        num_rows: int,
        initial_temperature: float = 1,
        cool_rate: float = 1e-3,
        max_iter: int = 10
) -> Tuple[bool, int, Sequence[int], Sequence[int]]:
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
    num = [int(2 * random.uniform(0, 1)) if len(J[i]) > 0 else 0
           for i in range(dim)]
    numtot = sum(num)

    H = 0
    for i in range(dim):
        H += fields[i] * num[i] * num[i]
        if len(J[i]):
            for j in range(len([J[i]])):
                H += J2[i][j] * num[i] * num[J[i][j]]

    count = 0
    T1 = initial_temperature
    howmany = 0
    e = math.exp(-1 / T1)

    while True:
        en = int(random.uniform(0, 1) * dim)
        # Note: Bug in original c++ code (while loop without any side effect
        # changed to if statement to prevent a possibly infinite loop)
        if len(J[en]) == 0:
            en = int(random.uniform(0, 1) * dim)
        p = 1
        if num[en] > 0 and random.uniform(0, 1) < 0.5:
            p = -1
        delta = fields[en] * num[en]
        for i in range(len(J[en])):
            delta += J2[en][i] * num[J[en][i]]
        delta = 2 * p * delta + fields[en]

        if delta < 0 or random.uniform(0, 1) < math.pow(e, delta):
            num[en] += p
            numtot += p
            H += delta

        count += 1

        if count % int(dim) == 0:
            T1 -= cool_rate
            if T1 <= 0:
                T1 = cool_rate
                e = math.exp(-1 / T1)

        if count == int(float(dim) / cool_rate):
            T1 = initial_temperature
            e = math.exp(-1 / T1)
            count = 0
            for i in range(dim):
                num[i] = 0
            en = int(random.uniform(0, 1) * dim)
            # Note: Bug in original c++ code (while loop without any side
            #  effect changed to if statement to prevent a possibly infinite
            #  loop)
            if len(J[en]):
                en = int(random.uniform(0, 1) * dim)
            num[en] = 1
            numtot = 1
            H = 0
            for i in range(dim):
                H += fields[i] * num[i] * num[i]
                if len(J[i]) > 0:
                    for j in range(len(J[i])):
                        H += J2[i][j] * num[i] * num[J[i][j]]
            howmany += 1

        if (H < MIN and numtot > 0) or (howmany == (10 * max_iter)):
            break

    if howmany < 10 * max_iter:
        if len(int_matched) > 0:
            yes = _is_linearly_dependent(num, int_kernel_dim, cls_species_idxs,
                                         cls_coefficients, matched, num_rows)
            assert yes, "Not true!"
        else:
            yes = True
        if yes:
            orders2 = list(range(len(matched)))
            pivots2 = matched[:]
            _qsort(len(matched), 0, orders2, pivots2)
            for i in range(len(matched)):
                if num[orders2[i]] > 0:
                    cls_species_idxs[int_kernel_dim].append(matched[orders2[i]])
                    cls_coefficients[int_kernel_dim].append(num[orders2[i]])
            int_kernel_dim += 1
            # side-effects on num vector
            _is_linearly_dependent(num, int_kernel_dim, cls_species_idxs,
                                   cls_coefficients, matched, num_rows)
            reduce(int_kernel_dim, cls_species_idxs, cls_coefficients, num_rows)
            min = 1000
            for i in range(len(cls_species_idxs[int_kernel_dim - 1])):
                if len(int_matched) == 0:
                    int_matched.append(cls_species_idxs[int_kernel_dim - 1][i])
                elif all(cur_int_matched
                         != cls_species_idxs[int_kernel_dim - 1][i]
                         for cur_int_matched in int_matched):
                    int_matched.append(cls_species_idxs[int_kernel_dim - 1][i])

                if cls_coefficients[int_kernel_dim - 1][i] < min:
                    min = cls_coefficients[int_kernel_dim - 1][i]
            for i in range(len(cls_species_idxs[int_kernel_dim - 1])):
                cls_coefficients[int_kernel_dim - 1][i] /= min
            logger.debug(
                f"Found linearly independent moiety, now there are "
                f"{int_kernel_dim} engaging {len(int_matched)} species")
        else:
            logger.debug(
                "Found a moiety but it is linearly dependent... next.")
    else:
        yes = False
    return yes, int_kernel_dim, matched, int_matched


def relax(
        stoichiometric_list: Sequence[Number],
        int_matched: Sequence[int],
        M: int,
        N: int,
        relaxation_max: float = 1e6,
        relaxation_step: float = 1.9
) -> bool:
    """Relaxation scheme for Monte Carlo final solution

    Checking for completeness using Motzkin's theorem. See Step (c) in
    De Martino (2014) and the Eqs. 14-16 in the corresponding publication

    :param stoichiometric_list:
        stoichiometric matrix :math:`S` as a flat list
    :param int_matched:
        intmatched
    :param M:
        number of species in reaction network
    :param N:
        number of reactions in reaction network
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
    matrix = [[] for _ in range(N)]
    matrix2 = [[] for _ in range(N)]

    i1 = 0
    j1 = 0
    K = len(int_matched)
    for _, val in enumerate(stoichiometric_list):
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
        if j1 == N:
            j1 = 0
            i1 += 1

    orders = [i for i in range(N)]
    pivots = [matrix[i][0] if len(matrix[i]) > 0 else MAX for i in range(N)]

    done = False
    while not done:
        _qsort(K, 0, orders, pivots)
        for j in range(K):
            if pivots[orders[j + 1]] == pivots[orders[j]] \
                    and pivots[orders[j]] != MAX:
                min1 = MAX
                if len(matrix[orders[j]]) > 1:
                    for i in range(len(matrix[orders[j]])):
                        if abs(matrix2[orders[j]][0]
                               / matrix2[orders[j]][i]) < min1:
                            min1 = matrix2[orders[j]][0] \
                                   / matrix2[orders[j]][i]
                min2 = MAX
                if len(matrix[orders[j + 1]]) > 1:
                    for i in range(len(matrix[orders[j]])):
                        if abs(matrix2[orders[j + 1]][0] /
                               matrix2[orders[j + 1]][i]) < min2:
                            min2 = abs(matrix2[orders[j + 1]][0]) \
                                   / matrix2[orders[j + 1]][i]
                if min2 > min1:
                    k2 = orders[j + 1]
                    orders[j + 1] = orders[j]
                    orders[j] = k2
        for j in range(K):
            if pivots[orders[j + 1]] == pivots[orders[j]] \
                    and pivots[orders[j]] != MAX:
                k1 = orders[j + 1]
                k2 = orders[j]
                column = [0] * M
                g = matrix2[k2][0] / matrix2[k1][0]
                for i in range(1, len(matrix[k1])):
                    column[matrix[k1][i]] = matrix2[k1][i] * g
                for i in range(1, len(matrix[k2])):
                    column[matrix[k2][i]] -= matrix2[k2][i]

                matrix[k1] = []
                matrix2[k1] = []
                for i in range(M):
                    if abs(column[i]) > MIN:
                        matrix[k1].append(i)
                        matrix2[k1].append(column[i])
                if len(matrix[orders[j + 1]]) > 0:
                    pivots[orders[j + 1]] = matrix[orders[j + 1]][0]
                else:
                    pivots[orders[j + 1]] = MAX

        for i in range(K):
            if len(matrix[i]) > 0:
                norm = matrix2[i][0]
                for j in range(len(matrix[i])):
                    matrix2[i][j] /= norm

        for k1 in reversed(range(K - 1)):
            k = orders[k1]
            if len(matrix[k]) > 1:
                for i in range(1, len(matrix[k])):
                    for j1 in range(k1 + 1, K):
                        j = orders[j1]
                        if len(matrix[j]) > 0:
                            if matrix[j][0] == matrix[k][i]:
                                row_k = [0] * M
                                for a in range(len(matrix[k])):
                                    row_k[matrix[k]][a] = matrix2[k][a]
                                for a in range(len(matrix[j])):
                                    row_k[matrix[j]][a] -= matrix2[j][a] * \
                                                           matrix2[k][i]
                                matrix = []
                                matrix2 = []
                                for a in range(M):
                                    if row_k[a] != 0:
                                        matrix[k].append(a)
                                        matrix2[k].append(row_k[a])

        indip = [K + 1] * M

        for i in range(K):
            if len(matrix[i]) > 0:
                indip[matrix[i]][0] = i

        M1 = 0
        for i in range(M):
            if indip[i] == K + 1:
                indip[i] = K + M1
                M1 += 1

        matrixAus = [[] for _ in range(M1)]
        matrixAus2 = [[] for _ in range(M1)]
        i1 = 0
        for i in range(M):
            if indip[i] >= K:
                matrixAus[i1].append(i)
                matrixAus2[i1].append(1)
                i1 += 1
            else:
                t = indip[i]
                if len(matrix[t]) > 1:
                    for k in range(1, len(matrix[t])):
                        quelo = indip[matrix[t]][k] - K
                        matrixAus[quelo].append(i)
                        matrixAus2[quelo].append(-matrix2[t][k])

        for i in range(K):
            matrix[i] = []

        N1 = N - K
        matrix_aus = [[] for _ in range(N1)]
        matrix_aus2 = [[] for _ in range(N1)]

        k1 = 0
        i1 = 0
        j1 = 0

        for _, val in enumerate(stoichiometric_list):
            prendo = 1
            if len(int_matched) > 0:
                for i in range(len(int_matched)):
                    if j1 == int_matched[i]:
                        prendo -= 1
            if val != 0:
                if prendo == 1:
                    matrix_aus[k1].append(i1)
                    matrix_aus2[k2].append(val)
            j1 += 1
            k1 += prendo
            if j1 == N:
                j1 = 0
                k1 = 0

        matrixb = [[] for _ in range(N1)]
        matrixb2 = [[] for _ in range(N1)]
        for i in range(M1):
            for j in range(N1):
                prod = 0
                if len(matrix_aus[j]) * len(matrixAus[i]) > 0:
                    for ib in range(len(matrixAus[i])):
                        for jb in range(len(matrix_aus[j])):
                            if matrixAus[i][ib] == matrix_aus[j][jb]:
                                prod += matrixAus2[i][ib] * matrix_aus2[j][jb]
                    if abs(prod) > MIN:
                        matrixb[j].append(i)
                        matrixb2[j].append(prod)

        for i in range(M1):
            matrixAus[i] = []
            matrixAus2[i] = []

        for i in range(N1):
            matrix_aus[i] = []
            matrix_aus2[i] = []

        var = [MIN] * M1
        done = False
        time = 0
        while True:
            cmin = 1000
            for j in range(N1):
                constr = 0
                if len(matrixb[j]) > 0:
                    for i in range(len(matrixb[j])):
                        constr += matrixb2[j][i] * var[matrixb][j][i]
                    if constr < cmin:
                        min = j
                        cmin = constr
            time += 1
            if cmin >= 0:
                done = True
            else:
                alpha = -relaxation_step * cmin  # Motzkin relaxation
                fact = 0
                for j in range(len(matrixb[min])):
                    fact += matrixb2[min][j] * matrixb2[min][j]
                alpha /= fact
                if alpha < 1e-9 * MIN:
                    alpha = 1e-9 * MIN
                for j in range(len(matrixb[min])):
                    var[matrixb[min][j]] += alpha * matrixb2[min][j]

            if done or time >= relaxation_max:
                break

        if done:
            return True
    return False


def reduce(
        int_kernel_dim: int,
        cls_species_idxs: MutableSequence[MutableSequence[int]],
        cls_coefficients: MutableSequence[MutableSequence[Number]],
        num_rows: int
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
    :param num_rows:
        number of rows in :math:`S`
    """
    K = int_kernel_dim
    MIN = 1e-9
    orders = list(range(K))
    pivots = [-len(cls_species_idxs[i]) for i in range(K)]

    while True:
        _qsort(K, 0, orders, pivots)
        ok = True
        for i in range(K - 2):
            for j in range(i + 1, K):
                k1 = orders[i]
                k2 = orders[j]
                column = [0] * num_rows
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
                    ok = False
                    cls_species_idxs[k1] = []
                    cls_coefficients[k1] = []
                    for l in range(num_rows):
                        if abs(column[l]) > MIN:
                            cls_species_idxs[k1].append(l)
                            cls_coefficients[k1].append(column[l])
                    pivots[k1] = -len(cls_species_idxs[k1])
        if ok:
            break
