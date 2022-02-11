"""Tests for conservation laws / conserved moieties"""

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
    row_names = [entry.decode('ascii')
                 for entry in io.BytesIO(response.read())]

    return S, row_names


def output(
        intKernelDim, kernelDim, intmatched, NSolutions, NSolutions2, row_names
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

    for i in range(intKernelDim):
        print(f"Moiety number {i + 1} engages {len(NSolutions[i])} "
              "metabolites.")
        print('\t')
        for j in range(len(NSolutions[i])):
            print(f"{row_names[NSolutions[i][j]]} \t {NSolutions2[i][j]}")


@log_execution_time("Detecting conservation laws", logger)
def test_detect_cl(data_demartino2014):
    """Invoke test case and benchmarking for De Martino's published results
    for E. coli network"""
    S, row_names = data_demartino2014
    N = 1668
    M = 2381
    knownValuesFromDeMartino = \
        [53] + [2] * 11 + [6] + [3] * 2 + [2] * 15 + [3] + [2] * 5

    assert len(S) == N * M, "Unexpected dimension of stoichiometric matrix"

    start = perf_counter()
    logger.debug(f"Kernel calculation of S ({N} x {M})...\n")
    kernelDim, engagedMetabolites, intKernelDim, conservedMoieties, \
        NSolutions, NSolutions2 = kernel(S, N, M)
    logger.debug(f"There are {kernelDim} conservation laws, engaging, "
                 f"{len(engagedMetabolites)} metabolites, {intKernelDim} are "
                 f"integers (conserved moieties), engaging "
                 f"{len(conservedMoieties)} metabolites...""")
    logger.debug("Kernel calc")
    output(intKernelDim, kernelDim, engagedMetabolites, NSolutions,
           NSolutions2, row_names)
    logger.debug("Kernel calc")

    # There are 38 conservation laws, engaging 131 metabolites
    # 36 are integers (conserved moieties), engaging 128 metabolites (from C++)
    assert kernelDim == 38, "Not all conservation laws found"
    assert intKernelDim == 36, "Not all conserved moiety laws found"
    assert len(engagedMetabolites) == 131, \
        "Wrong number of engaged metabolites reported"
    assert len(conservedMoieties) == 128, \
        "Wrong number of conserved moieties reported"

    logger.debug("-" * 80)
    logger.debug("-" * 80)
    logger.debug("-" * 80)

    logger.debug("Filling interaction matrix...\n")
    J, J2, fields = fill(S, len(engagedMetabolites), engagedMetabolites, N)

    logger.debug("after fill")

    finish = 0
    if intKernelDim == kernelDim:
        finish = 1
    else:
        timer = 0
        counter = 1
        maxIter = 10
        while finish == 0:
            logger.debug(f"MonteCarlo call #{counter} (maxIter: {maxIter})")
            yes, intKernelDim, kernelDim, NSolutions, NSolutions2, \
                engagedMetabolites, conservedMoieties = MonteCarlo(
                    engagedMetabolites, J, J2, fields, conservedMoieties,
                    intKernelDim, NSolutions, NSolutions2, kernelDim, N)
            output(intKernelDim, kernelDim, engagedMetabolites, NSolutions,
                   NSolutions2, row_names)

            counter += 1
            if intKernelDim == kernelDim:
                finish = 1
            if yes == 0:
                timer += 1
            if timer == max:
                logger.debug("Relaxation...")
                finish = Relaxation(S, conservedMoieties, M, N)
                if finish == 1:
                    timer = 0
        old = NSolutions
        old2 = NSolutions2
        intKernelDim, kernelDim, NSolutions, NSolutions2 = Reduce(intKernelDim,
                                                                  kernelDim,
                                                                  NSolutions,
                                                                  NSolutions2,
                                                                  N)
        for i in range(len(old)):
            assert (set(old[i]) == set(NSolutions[i]))
            assert (set(old2[i]) == set(NSolutions2[i]))

    # Assert that each conserved moiety has the correct number of metabolites
    for i in range(intKernelDim - 2):
        assert (len(NSolutions[i]) == knownValuesFromDeMartino[i]),\
            f"Moiety #{i + 1} failed for test case (De Martino et al.)"

    logger.debug("*" * 80)
    logger.debug("Details about conserved moieties:")
    logger.debug("*" * 80)
    output(intKernelDim, kernelDim, engagedMetabolites, NSolutions,
           NSolutions2, row_names)
    logger.debug("-" * 80)
    logger.debug("-" * 80)
    logger.debug("-" * 80)

    end = perf_counter()
    logger.debug(f"Execution time: {end - start} [s]")

    return end - start


@log_execution_time("Detecting moiety conservation laws", logger)
def test_cl_detect_execution_time(data_demartino2014):
    """Test execution time stays within a certain predefined bound"""
    max_tries = 3
    # <5s on modern hardware, but leave some slack
    max_time_seconds = 10
    runtime = np.Inf

    for _ in range(max_tries):
        runtime = test_detect_cl(data_demartino2014)
        if runtime < max_time_seconds:
            break
    assert runtime < max_time_seconds, "Took too long"
