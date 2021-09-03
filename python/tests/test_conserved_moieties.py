""" Tests for conservative laws / conserved moities """

from amici.conserved_moeties import *

def test_detect_cl():
	S = GetRemoteInput()
	N = 1668
	M = 2381
	knownValuesFromDeMartino = [53] + [2] * 11 + [6] + [3] * 2 + [2] * 15 + [3] + [2] * 5 

	if len(S) != N*M:
		print("Stoichiometric matrix inconsistent")

	print(f"Kernel calculation of S ({N} x {M})...\n")
	kernelDim, engagedMetabolites, intKernelDim, conservedMoieties, NSolutions, NSolutions2 = kernel(S, N, M)
	print(f"""There are {kernelDim} conservation laws, engaging, 
	{len(engagedMetabolites)} metabolites, {intKernelDim} are integers (conserved 
	moeities), engaging {len(conservedMoieties)} metabolites...""")
	print("Kernel calc")
	Output(intKernelDim, kernelDim, engagedMetabolites, NSolutions, IsRemoteFile=False)
	print("Kernel calc")

	# There are 38 conservation laws, engaging 131 metabolites 
	# 36 are integers (conserved moieties), engaging 128 metabolites (from C++)
	assert (kernelDim == 38), "Not all conservation laws found"
	assert (intKernelDim == 36), "Not all conserved moiety laws found"
	assert (len(engagedMetabolites) == 131), "Wrong number of engaged metabolites reported"
	assert (len(conservedMoieties) == 128), "Wrong number of conserved moeities reported"

	print("".join(['-' for _ in range(0, 80)]))
	print("".join(['-' for _ in range(0, 80)]))
	print("".join(['-' for _ in range(0, 80)]))

	print("Filling interaction matrix...\n")
	J, J2, fields = fill(S, len(engagedMetabolites), engagedMetabolites)

	print("after fill")

	finish = 0
	if intKernelDim == kernelDim:
		finish = 1
	else:
		timer = 0
		counter = 1
		maxIter = 10
		while finish == 0:
			print(f"MonteCarlo call #{counter} (maxIter: {maxIter})")
			yes, intKernelDim, kernelDim, NSolutions, NSolutions2, engagedMetabolites, conservedMoieties = MonteCarlo(engagedMetabolites, J, J2, fields, conservedMoieties, intKernelDim, NSolutions, NSolutions2, kernelDim)
			Output(intKernelDim, kernelDim, engagedMetabolites, NSolutions)

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
		old=NSolutions
		old2=NSolutions2
		intKernelDim, kernelDim, NSolutions, NSolutions2 = Reduce(intKernelDim, kernelDim, NSolutions, NSolutions2)
		for i in range(0, len(old)):
			assert(set(old[i]) == set(NSolutions[i]))
			assert(set(old2[i]) == set(NSolutions2[i]))


	# Assert that each conserved moeity has the correct number of metabolites (TODO: last two moeities fluctuate in DeMartino C++ implementation, likewise our implementation fluctuates, thus excluding -> Figure out how to avoid this...)
	for i in range(0, intKernelDim-2): 
		assert (len(NSolutions[i]) == knownValuesFromDeMartino[i]), f"Moeity #{i+1} failed for test case (De Martino et al.)"

	print("".join(['*' for _ in range(0, 80)]))
	print("Details about conserved moeities:")
	print("".join(['*' for _ in range(0, 80)]))
	Output(intKernelDim, kernelDim, engagedMetabolites, NSolutions)
	print("".join(['-' for _ in range(0, 80)]))
	print("".join(['-' for _ in range(0, 80)]))
	print("".join(['-' for _ in range(0, 80)]))
	
	end = time.time()
	print(f"Execution time: {end-start} [s]")