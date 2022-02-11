"""Tests for conservation laws / conserved moieties"""
import os
from time import perf_counter

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
        intKernelDim, kernelDim, intmatched, NSolutions, NSolutions2,
        row_names, verbose=False
):
    """
    Output solution

    :param intKernelDim:
        intKernelDim
    :param kernelDim:
        kernelDim
    :param intmatched:
        intmatched
    :param NSolutions:
        NSolutions
    """
    print(f"There are {intKernelDim} linearly independent conserved moieties, "
          f"engaging {len(intmatched)} metabolites\n")
    if intKernelDim == kernelDim:
        print("They generate all the conservation laws")
    else:
        print(f"They don't generate all the conservation laws, "
              f"{kernelDim - intKernelDim} of them are not reducible to "
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
    S, row_names = data_demartino2014
    N = 1668
    M = 2381
    # Expected number of metabolites per conservation law
    knownValuesFromDeMartino = \
        [53] + [2] * 11 + [6] + [3] * 2 + [2] * 15 + [3] + [2] * 5

    assert len(S) == N * M, "Unexpected dimension of stoichiometric matrix"

    start = perf_counter()
    kernelDim, engagedMetabolites, intKernelDim, conservedMoieties, \
    NSolutions, NSolutions2 = kernel(S, N, M)

    if not quiet:
        output(intKernelDim, kernelDim, engagedMetabolites, NSolutions,
               NSolutions2, row_names)

    # There are 38 conservation laws, engaging 131 metabolites
    # 36 are integers (conserved moieties), engaging 128 metabolites (from C++)
    assert kernelDim == 38, "Not all conservation laws found"
    assert intKernelDim == 36, "Not all conserved moiety laws found"
    assert len(engagedMetabolites) == 131, \
        "Wrong number of engaged metabolites reported"
    assert len(conservedMoieties) == 128, \
        "Wrong number of conserved moieties reported"

    J, J2, fields = fill(S, engagedMetabolites, N)

    if intKernelDim != kernelDim:
        timer = 0
        counter = 1
        maxIter = 10
        finish = 0
        while finish == 0:
            if not quiet:
                print(f"Monte Carlo call #{counter} (maxIter: {maxIter})")
            yes, intKernelDim, kernelDim, NSolutions, NSolutions2, \
            engagedMetabolites, conservedMoieties = monte_carlo(
                engagedMetabolites, J, J2, fields, conservedMoieties,
                intKernelDim, NSolutions, NSolutions2, kernelDim, N)
            if not quiet:
                output(intKernelDim, kernelDim, engagedMetabolites, NSolutions,
                       NSolutions2, row_names)

            counter += 1
            if intKernelDim == kernelDim:
                finish = 1
            if yes == 0:
                timer += 1
            if timer == max:
                if not quiet:
                    print("Relaxation...")
                finish = relax(S, conservedMoieties, M, N)
                if finish == 1:
                    timer = 0
        old = NSolutions
        old2 = NSolutions2
        intKernelDim, kernelDim, NSolutions, NSolutions2 = reduce(intKernelDim,
                                                                  kernelDim,
                                                                  NSolutions,
                                                                  NSolutions2,
                                                                  N)
        for i in range(len(old)):
            assert (set(old[i]) == set(NSolutions[i]))
            assert (set(old2[i]) == set(NSolutions2[i]))

    # Assert that each conserved moiety has the correct number of metabolites
    for i in range(intKernelDim - 2):
        assert (len(NSolutions[i]) == knownValuesFromDeMartino[i]), \
            f"Moiety #{i + 1} failed for test case (De Martino et al.)"

    if not quiet:
        output(intKernelDim, kernelDim, engagedMetabolites, NSolutions,
               NSolutions2, row_names)

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
    max_time_seconds = 25 if "GITHUB_ACTIONS" in os.environ else 10

    runtime = np.Inf

    for _ in range(max_tries):
        runtime = test_detect_cl(data_demartino2014, quiet=True)
        if runtime < max_time_seconds:
            break
    assert runtime < max_time_seconds, "Took too long"
