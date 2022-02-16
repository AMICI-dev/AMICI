"""Tests for conservation laws / conserved moieties"""
import os
from time import perf_counter

import sympy as sp
import numpy as np
import pytest

from amici.conserved_moieties import *
from amici.logging import get_logger, log_execution_time

logger = get_logger(__name__)


@pytest.fixture(scope="session")
def data_demartino2014():
    """Get tests from DeMartino2014 Suppl. Material"""
    import urllib.request
    import io
    import gzip

    # stoichiometric matrix
    response = urllib.request.urlopen(
        r'https://chimera.roma1.infn.it/SYSBIO/test-ecoli.dat.gz')
    data = gzip.GzipFile(fileobj=io.BytesIO(response.read()))
    S = [int(item) for sl in
         [entry.decode('ascii').strip().split('\t')
          for entry in data.readlines()] for item in sl]

    # metabolite / row names
    response = urllib.request.urlopen(
        r'https://chimera.roma1.infn.it/SYSBIO/test-ecoli-met.txt')
    row_names = [entry.decode('ascii').strip()
                 for entry in io.BytesIO(response.read())]

    return S, row_names


def output(
        int_kernel_dim, kernel_dim, int_matched, NSolutions, NSolutions2,
        row_names, verbose=False
):
    """
    Output solution

    :param int_kernel_dim:
        intKernelDim
    :param kernel_dim:
        kernelDim
    :param int_matched:
        intmatched
    :param NSolutions:
        NSolutions
    """
    print(
        f"There are {int_kernel_dim} linearly independent conserved moieties, "
        f"engaging {len(int_matched)} metabolites\n")
    if int_kernel_dim == kernel_dim:
        print("They generate all the conservation laws")
    else:
        print(f"They don't generate all the conservation laws, "
              f"{kernel_dim - int_kernel_dim} of them are not reducible to "
              "moieties")
    # print all conservation laws
    if verbose:
        for i, coefficients, engaged_species_idxs \
                in enumerate(zip(NSolutions2, NSolutions)):
            print(f"Moiety number {i + 1} engages {len(engaged_species_idxs)} "
                  "metabolites:")
            for species_idx, coefficient \
                    in zip(engaged_species_idxs, coefficients):
                print(f"\t{row_names[species_idx]}\t{coefficient}")


@log_execution_time("Detecting conservation laws", logger)
def test_detect_cl(data_demartino2014, quiet=False):
    """Invoke test case and benchmarking for De Martino's published results
    for E. coli network"""
    stoichiometric_list, row_names = data_demartino2014
    num_species = 1668
    num_reactions = 2381
    # Expected number of metabolites per conservation law found after kernel()
    expected_num_species = \
        [53] + [2] * 11 + [6] + [3] * 2 + [2] * 15 + [3] + [2] * 5

    assert len(stoichiometric_list) == num_species * num_reactions,\
        "Unexpected dimension of stoichiometric matrix"

    start = perf_counter()
    # TODO: remove redundancy with compute_moiety_conservation_laws()
    kernel_dim, engaged_species, int_kernel_dim, conserved_moieties, \
        cls_species_idxs, cls_coefficients = kernel(
        stoichiometric_list, num_species, num_reactions)

    if not quiet:
        output(int_kernel_dim, kernel_dim, engaged_species, cls_species_idxs,
               cls_coefficients, row_names)

    # There are 38 conservation laws, engaging 131 metabolites
    # 36 are integers (conserved moieties), engaging 128 metabolites (from C++)
    assert kernel_dim == 38, "Not all conservation laws found"
    assert int_kernel_dim == 36, "Not all conserved moiety laws found"
    assert len(engaged_species) == 131, \
        "Wrong number of engaged metabolites reported"
    assert len(conserved_moieties) == 128, \
        "Wrong number of conserved moieties reported"

    # Assert that each conserved moiety has the correct number of metabolites
    for i in range(int_kernel_dim - 2):
        assert (len(cls_species_idxs[i]) == expected_num_species[i]), \
            f"Moiety #{i + 1} failed for test case (De Martino et al.)"

    if int_kernel_dim != kernel_dim:
        J, J2, fields = fill(stoichiometric_list, engaged_species, num_species)

        timer = 0
        counter = 1
        finish = 0
        max_num_monte_carlo = 10
        while finish == 0:
            if not quiet:
                print(f"Monte Carlo call #{counter}")
            yes, int_kernel_dim, conserved_moieties = \
                monte_carlo(engaged_species, J, J2, fields, conserved_moieties,
                            int_kernel_dim, cls_species_idxs, cls_coefficients,
                            num_species=num_species, max_iter=max_num_monte_carlo)
            if not quiet:
                output(int_kernel_dim, kernel_dim, engaged_species,
                       cls_species_idxs, cls_coefficients, row_names)

            counter += 1
            if int_kernel_dim == kernel_dim:
                finish = 1
            if yes == 0:
                timer += 1
            else:
                timer = 0
            if timer == max_num_monte_carlo:
                if not quiet:
                    print("Relaxation...")
                finish = relax(stoichiometric_list, conserved_moieties, num_reactions, num_species)
                timer = 0
    reduce(int_kernel_dim, cls_species_idxs, cls_coefficients, num_species)

    if not quiet:
        output(int_kernel_dim, kernel_dim, engaged_species, cls_species_idxs,
               cls_coefficients, row_names)

    runtime = perf_counter() - start
    if not quiet:
        print(f"Execution time: {runtime} [s]")
    return runtime


@log_execution_time("Detecting moiety conservation laws", logger)
def test_cl_detect_execution_time(data_demartino2014):
    """Test execution time stays within a certain predefined bound.
    As the algorithm is non-deterministic, allow for some retries.
    Only one has to succeed."""
    max_tries = 3
    # <5s on modern hardware, but leave some slack
    max_time_seconds = 30 if "GITHUB_ACTIONS" in os.environ else 10

    runtime = np.Inf

    for _ in range(max_tries):
        runtime = test_detect_cl(data_demartino2014, quiet=True)
        if runtime < max_time_seconds:
            break
    assert runtime < max_time_seconds, "Took too long"


def test_compute_moiety_conservation_laws_simple():
    """Test a simple example, ensure the conservation laws are identified
     reliably. Requires the Monte Carlo to identify all."""
    stoichiometric_matrix = sp.Matrix([
        [-1.0, 1.0],
        [-1.0, 1.0],
        [1.0, -1.0],
        [1.0, -1.0]]
    )
    stoichiometric_list = [
        float(entry) for entry in stoichiometric_matrix.T.flat()
    ]

    num_tries = 1000
    found_all_n_times = 0
    for _ in range(num_tries):
        cls_state_idxs, cls_coefficients = compute_moiety_conservation_laws(
            stoichiometric_list, *stoichiometric_matrix.shape)

        assert cls_state_idxs in ([[0, 3], [1, 2], [1, 3]],
                                  [[0, 3], [1, 2], [0, 2]],
                                  # should happen rarely
                                  [[0, 3], [1, 2]])
        assert cls_coefficients in ([[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]],
                                    [[1.0, 1.0], [1.0, 1.0]])

        num_cls_found = len(cls_state_idxs)
        if num_cls_found == 3:
            found_all_n_times += 1
    # sometimes we don't find all conservation laws, but this should be rare
    assert found_all_n_times / num_tries > 0.999

