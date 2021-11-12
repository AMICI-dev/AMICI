""" Tests for conservative laws / conserved moities """

from amici.conserved_moieties import *
from time import perf_counter
from amici.logging import get_logger, log_execution_time

logger = get_logger(__name__)

@log_execution_time("Detecting moeity conservation laws", logger)
def test_detect_cl():
	S = GetRemoteInput()
	N = 1668
	M = 2381
	knownValuesFromDeMartino = [53] + [2] * 11 + [6] + [3] * 2 + [2] * 15 + [3] + [2] * 5 

	if len(S) != N*M:
		logger.debug("Stoichiometric matrix inconsistent")

	start = perf_counter()
	logger.debug(f"Kernel calculation of S ({N} x {M})...\n")
	kernelDim, engagedMetabolites, intKernelDim, conservedMoieties, NSolutions, NSolutions2 = kernel(S, N, M)
	logger.debug(f"""There are {kernelDim} conservation laws, engaging, 
	{len(engagedMetabolites)} metabolites, {intKernelDim} are integers (conserved 
	moeities), engaging {len(conservedMoieties)} metabolites...""")
	logger.debug("Kernel calc")
	Output(intKernelDim, kernelDim, engagedMetabolites, NSolutions, IsRemoteFile=False)
	logger.debug("Kernel calc")

	# There are 38 conservation laws, engaging 131 metabolites 
	# 36 are integers (conserved moieties), engaging 128 metabolites (from C++)
	assert (kernelDim == 38), "Not all conservation laws found"
	assert (intKernelDim == 36), "Not all conserved moiety laws found"
	assert (len(engagedMetabolites) == 131), "Wrong number of engaged metabolites reported"
	assert (len(conservedMoieties) == 128), "Wrong number of conserved moeities reported"

	logger.debug("".join(['-' for _ in range(0, 80)]))
	logger.debug("".join(['-' for _ in range(0, 80)]))
	logger.debug("".join(['-' for _ in range(0, 80)]))

	logger.debug("Filling interaction matrix...\n")
	J, J2, fields = fill(S, len(engagedMetabolites), engagedMetabolites)

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
			yes, intKernelDim, kernelDim, NSolutions, NSolutions2, engagedMetabolites, conservedMoieties = MonteCarlo(engagedMetabolites, J, J2, fields, conservedMoieties, intKernelDim, NSolutions, NSolutions2, kernelDim)
			Output(intKernelDim, kernelDim, engagedMetabolites, NSolutions)

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
		old=NSolutions
		old2=NSolutions2
		intKernelDim, kernelDim, NSolutions, NSolutions2 = Reduce(intKernelDim, kernelDim, NSolutions, NSolutions2)
		for i in range(0, len(old)):
			assert(set(old[i]) == set(NSolutions[i]))
			assert(set(old2[i]) == set(NSolutions2[i]))


	# Assert that each conserved moeity has the correct number of metabolites 
	for i in range(0, intKernelDim-2): 
		assert (len(NSolutions[i]) == knownValuesFromDeMartino[i]), f"Moeity #{i+1} failed for test case (De Martino et al.)"

	logger.debug("".join(['*' for _ in range(0, 80)]))
	logger.debug("Details about conserved moeities:")
	logger.debug("".join(['*' for _ in range(0, 80)]))
	Output(intKernelDim, kernelDim, engagedMetabolites, NSolutions)
	logger.debug("".join(['-' for _ in range(0, 80)]))
	logger.debug("".join(['-' for _ in range(0, 80)]))
	logger.debug("".join(['-' for _ in range(0, 80)]))
	
	end = perf_counter()
	logger.debug(f"Execution time: {end-start} [s]")
	return end-start

@log_execution_time("Detecting moeity conservation laws", logger)
def test_cl_detect_execution_time():
	""" Test execution time stays within a certain predefined bound """
	numIterations = 100
	thresholdForTimeout = 5
	timings = [test_detect_cl() for _ in range(0, numIterations)]

	for timing in timings:
		assert timing < thresholdForTimeout
