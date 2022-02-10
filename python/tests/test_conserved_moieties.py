"""Tests for conservation laws / conserved moieties"""

import time
from time import perf_counter

from amici.conserved_moieties import *
from amici.logging import get_logger, log_execution_time

logger = get_logger(__name__)


def GetRemoteInput():
    import urllib.request
    import io
    import gzip
    response = urllib.request.urlopen(
        r'http://chimera.roma1.infn.it/SYSBIO/test-ecoli.dat.gz')
    data = gzip.GzipFile(fileobj=io.BytesIO(response.read()))
    return [int(item) for sl in
            [entry.decode('ascii').strip().split('\t') for entry in
             data.readlines()] for item in sl]


def GetRemoteNames():
    import urllib.request
    import io
    response = urllib.request.urlopen(
        r'http://chimera.roma1.infn.it/SYSBIO/test-ecoli-met.txt')
    return [entry.decode('ascii') for entry in io.BytesIO(response.read())]


def Input():
    """ Read input from test stoichiometric matrix """
    with open('matrix.dat', 'r') as f:
        return [int(item) for sl in [entry.strip().split('\t') for entry in f]
                for item in sl]


def Names():
    """ Get names of metabolites from list of metabolites example"""
    with open('metabolites.txt', 'r') as f:
        return list(f)


def Output(
        intKernelDim, kernelDim, intmatched, NSolutions, NSolutions2,
        IsRemoteFile=False
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
              "moeities")

    names = IsRemoteFile and Names() or GetRemoteNames()
    for i in range(intKernelDim):
        print(f"Moiety number {i + 1} engages {len(NSolutions[i])} "
              "metabolites.")
        print('\t')
        for j in range(len(NSolutions[i])):
            print(f"{names[NSolutions[i][j]]} \t {NSolutions2[i][j]}")


@log_execution_time("Detecting conservation laws", logger)
def test_detect_cl():
    S = GetRemoteInput()
    N = 1668
    M = 2381
    knownValuesFromDeMartino = [53] + [2] * 11 + [6] + [3] * 2 + [2] * 15 + [
        3] + [2] * 5

    if len(S) != N * M:
        logger.debug("Stoichiometric matrix inconsistent")

    start = perf_counter()
    logger.debug(f"Kernel calculation of S ({N} x {M})...\n")
    kernelDim, engagedMetabolites, intKernelDim, conservedMoieties, \
    NSolutions, NSolutions2 = kernel(S, N, M)
    logger.debug(f"There are {kernelDim} conservation laws, engaging, "
                 f"{len(engagedMetabolites)} metabolites, {intKernelDim} are "
                 f"integers (conserved moieties), engaging "
                 f"{len(conservedMoieties)} metabolites...""")
    logger.debug("Kernel calc")
    Output(intKernelDim, kernelDim, engagedMetabolites, NSolutions,
           NSolutions2, IsRemoteFile=False)
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
            yes, intKernelDim, kernelDim, NSolutions, NSolutions2, engagedMetabolites, conservedMoieties = MonteCarlo(
                engagedMetabolites, J, J2, fields, conservedMoieties,
                intKernelDim, NSolutions, NSolutions2, kernelDim, N)
            Output(intKernelDim, kernelDim, engagedMetabolites, NSolutions,
                   NSolutions2)

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
    for i in range(0, intKernelDim - 2):
        assert (len(NSolutions[i]) == knownValuesFromDeMartino[
            i]), f"Moiety #{i + 1} failed for test case (De Martino et al.)"

    logger.debug("*" * 80)
    logger.debug("Details about conserved moieties:")
    logger.debug("*" * 80)
    Output(intKernelDim, kernelDim, engagedMetabolites, NSolutions,
           NSolutions2)
    logger.debug("-" * 80)
    logger.debug("-" * 80)
    logger.debug("-" * 80)

    end = perf_counter()
    logger.debug(f"Execution time: {end - start} [s]")

    return end - start


@log_execution_time("Detecting moeity conservation laws", logger)
def test_cl_detect_execution_time():
    """ Test execution time stays within a certain predefined bound """
    numIterations = 2
    thresholdForTimeout = 5
    timings = [test_detect_cl() for _ in range(numIterations)]

    for timing in timings:
        assert timing < thresholdForTimeout


def test_some_other_test():
    """Invoke test case and benchmarking for De Martino's published results
    for E. coli network"""
    start = time.time()
    print("*" * 80)
    print("Conserved moieties test cases")
    print("*" * 80)

    # Hard-coded stoichiometric matrix as test case containing _ONE_
    # conservative law (i.e. one conserved moiety)
    S = [-1, 0, 0, 0, 0, 0,
         1, -1, 1, 0, 0, 0,
         0, 1, -1, -1, -1, -1,
         0, 0, 0, 1, -1, 0,
         0, 0, 0, 0, 0, 1]
    N = 6  # number of metabolites (columns)
    M = 5  # number of reactions (rows)

    S = Input()
    N = 1668
    M = 2381
    knownValuesFromDeMartino = [53] + [2] * 11 + [6] + [3] * 2 + [2] * 15 + [
        3] + [2] * 5

    if len(S) != N * M:
        print("Stoichiometric matrix inconsistent")

    print(f"Kernel calculation of S ({N} x {M})...\n")
    kernelDim, engagedMetabolites, intKernelDim, conservedMoieties, \
    NSolutions, NSolutions2 = kernel(S, N, M)
    print(f"""There are {kernelDim} conservation laws, engaging, 
    	{len(engagedMetabolites)} metabolites, {intKernelDim} are integers (conserved 
    	moeities), engaging {len(conservedMoieties)} metabolites...""")
    print("Kernel calc")
    Output(intKernelDim, kernelDim, engagedMetabolites, NSolutions,
           NSolutions2)
    print("Kernel calc")

    # There are 38 conservation laws, engaging 131 metabolites
    # 36 are integers (conserved moieties), engaging 128 metabolites (from C++)
    assert (kernelDim == 38), "Not all conservation laws found"
    assert (intKernelDim == 36), "Not all conserved moiety laws found"
    assert (
            len(engagedMetabolites) == 131), "Wrong number of engaged metabolites reported"
    assert (
            len(conservedMoieties) == 128), "Wrong number of conserved moieties reported"

    print("-" * 80)
    print("-" * 80)
    print("-" * 80)

    print("Filling interaction matrix...\n")
    J, J2, fields = fill(S, len(engagedMetabolites), engagedMetabolites, N)

    print("after fill")

    finish = 0
    if intKernelDim == kernelDim:
        finish = 1
    else:
        timer = 0
        counter = 1
        maxIter = 10
        while finish == 0:
            # TODO: Sometimes MC runs into the wrong solution: Why?
            print(f"MonteCarlo call #{counter} (maxIter: {maxIter})")
            yes, intKernelDim, kernelDim, NSolutions, NSolutions2, \
            engagedMetabolites, conservedMoieties = MonteCarlo(
                engagedMetabolites, J, J2, fields, conservedMoieties,
                intKernelDim, NSolutions, NSolutions2, kernelDim, N)
            Output(intKernelDim, kernelDim, engagedMetabolites, NSolutions,
                   NSolutions2)

            counter += 1
            if intKernelDim == kernelDim:
                finish = 1
            if yes == 0:
                timer += 1
            if timer == max:
                print("Relaxation...")
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
    # (TODO: last two moeities fluctuate in DeMartino C++ implementation,
    #  likewise our implementation fluctuates, thus excluding
    #  -> Figure out how to avoid this...)
    for i in range(intKernelDim - 2):
        assert (len(NSolutions[i]) == knownValuesFromDeMartino[
            i]), f"Moiety #{i + 1} failed for test case (De Martino et al.)"

    print("*" * 80)
    print("Details about conserved moieties:")
    print("*" * 80)
    Output(intKernelDim, kernelDim, engagedMetabolites, NSolutions,
           NSolutions2)
    print("-" * 80)
    print("-" * 80)
    print("-" * 80)

    end = time.time()
    print(f"Execution time: {end - start} [s]")
