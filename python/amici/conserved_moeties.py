import random
import math
import sys
import time
sys.setrecursionlimit(3000)

# TODO: Add proper Pysb unit tests for the example and all algorithmic parts

# TODO: Refactoring to more idiomatic Python code

# TODO: Write proper Pydoc and apply coding conventions of the group

# TODO: Rewrite quicksort algorithm (recursive) to more efficient iterative variant

# TODO: Refactor to a class

def qsort(k, km, orders, pivots):
	""" 
	Quicksort
	
	:param k:
		k
	:param km:
		km
	:param orders:
		orders
	:param pivots:
		pivots
	"""
	centre = 0
	if k-km >= 1:
		pivot = km+int((k-km)/2) # was // 2 without int
		l = 0
		p = k - km - 1
		neworders = [None] * (k-km)
		for i in range(km, k):
			if i != pivot:
				if (pivots[orders[i]] < pivots[orders[pivot]]):
					neworders[l] = orders[i]
					l += 1
				else:
					neworders[p] = orders[i]
					p -= 1
		neworders[p] = orders[pivot]
		for i in range(km, k):
			orders[i] = neworders[i-km]

		centre = p + km
		qsort(k, centre+1, orders, pivots)
		qsort(centre, km, orders, pivots)

def kernel(stoichiometricMatrixAsList, numberOfMetabolites, numberOfReactions):
	""" 
	Kernel calculation by Gaussian elimination
	
	:param stoichiometricMatrixAsList:
		the stoichiometric matrix as a list (reactions x metabolites)
	:param numberOfMetabolites
		total number of metabolites
	:param numberOfReactions
		total number of reactions
	"""
	il=0
	jl=0
	N = numberOfMetabolites
	M = numberOfReactions
	MAX = 1e9
	MIN = 1e-9
	
	matrix = [ [] for _ in range(N) ]
	matrix2 = [ [] for _ in range(N) ]
	matched = [ ]
	intmatched = [ ]
	NSolutions = [ [] for _ in range(N)]
	NSolutions2 = [ [] for _ in range(N)]

	for _, val in enumerate(stoichiometricMatrixAsList):
		if val != 0:
			matrix[jl].append(il)
			matrix2[jl].append(val)
		jl += 1
		if jl == N:
			jl=0
			il += 1

	for i in range(0, N):
		matrix[i].append(M+i)
		matrix2[i].append(1)
	
	ok = 0
	orders = [ i for i in range(0, N) ]
	pivots = [ matrix[i][0] if len(matrix[i]) > 0 else MAX for i in range(0, N)]

	while ok == 0:
		qsort(N, 0, orders, pivots)
		for j in range(0, N-1):
			if pivots[orders[j+1]]==pivots[orders[j]] and pivots[orders[j]] != MAX:
				min1 = 100000000
				if len(matrix[orders[j]]) > 1:
					for i in range(0, len(matrix[orders[j]])):
						if abs(matrix2[orders[j]][0]/matrix2[orders[j]][i]) < min1:
							min1 = abs(matrix2[orders[j]][0]/matrix2[orders[j]][i])

				min2 = 100000000
				if len(matrix[orders[j+1]]) > 1:
					for i in range(0, len(matrix[orders[j+1]])):
						if abs(matrix2[orders[j+1]][0]/matrix2[orders[j+1]][i]) < min2:
							min2 = abs(matrix2[orders[j+1]][0]/matrix2[orders[j+1]][i])
						
				if min2 > min1:
					k2 = orders[j+1]
					orders[j+1] = orders[j]
					orders[j] = k2
		ok = 1

		for j in range(0, N-1):
			if pivots[orders[j+1]] == pivots[orders[j]] and pivots[orders[j]] != MAX:
				k1 = orders[j+1]
				k2 = orders[j]
				colonna = [ 0 for _ in range(N+M) ]
				g  = matrix2[k2][0]/matrix2[k1][0]
				for i in range(1, len(matrix[k1])):
					colonna[matrix[k1][i]] = matrix2[k1][i]*g

				for i in range(1, len(matrix[k2])):
					colonna[matrix[k2][i]] -= matrix2[k2][i]

				matrix[k1] = []
				matrix2[k1] = []
				for i in range(0, N+M):
					if abs(colonna[i]) > MIN:
						matrix[k1].append(i)
						matrix2[k1].append(colonna[i])

				ok = 0
				if len(matrix[orders[j+1]]) > 0:
					pivots[orders[j+1]] = matrix[orders[j+1]][0]
				else:
					pivots[orders[j+1]] = MAX

	RSolutions = [ [] for _ in range(N)]
	RSolutions2 = [ [] for _ in range(N)]
	kernelDim = 0

	for i in range(0, N):
		ok = 1
		if len(matrix[i]) > 0:
			for j in range(0, len(matrix[i])):
				if matrix[i][j] < M:
					ok = 0
		if ok == 1 and len(matrix[i]) > 0:
			for j in range(0, len(matrix[i])):
				RSolutions[kernelDim].append(matrix[i][j]-M)
				RSolutions2[kernelDim].append(matrix2[i][j])
			kernelDim += 1

	for i in range(0, N):
		matrix[i] = []
		matrix2[i] = []

	i2 = 0
	for i in range(0, kernelDim):
		ok2 = 1
		if (len(RSolutions[i])) > 0:
			for j in range(0, len(RSolutions[i])):
				if RSolutions2[i][j]*RSolutions2[i][0] < 0:
					ok2 = 0
				if len(matched) == 0:
					matched.append(RSolutions[i][j])
				else:
					ok3 = 1
					for k in range(0, len(matched)):
						if matched[k] == RSolutions[i][j]:
							ok3 = 0
					if ok3 == 1:
						matched.append(RSolutions[i][j])
		if ok2 == 1 and len(RSolutions[i]) > 0:
			min = MAX
			for j in range(0, len(RSolutions[i])):
				NSolutions[i2].append(RSolutions[i][j])
				NSolutions2[i2].append(abs(RSolutions2[i][j]))
				if min > abs(RSolutions2[i][j]):
					min = abs(RSolutions2[i][j])
				if len(intmatched) == 0:
					intmatched.append(NSolutions[i2][j])
				else:
					ok3 = 1
					for k in range(0, len(intmatched)):
						if intmatched[k] == NSolutions[i2][j]:
							ok3 = 0
					if ok3 == 1:
						intmatched.append(NSolutions[i2][j])
			for j in range(0, len(NSolutions[i2])):
				NSolutions2[i2][j] /= min
			i2 += 1
	intKernelDim=i2

	return kernelDim, matched, intKernelDim, intmatched, NSolutions, NSolutions2

def fill(stoichiometricMatrixAsList, matched_size, matched):
	""" 
	Interaction matrix construction 
	
	:param stoichiometricMatrixAsList:
		the stoichiometric matrix as a list
	:param matched_size:
		found MCLs in the matrix S
	:param matched
		actual found MCLs
	"""
	dim = matched_size
	MIN = 1e-9
	matrix = [ [] for _ in range(dim) ]
	matrix2 = [ [] for _ in range(dim) ]

	J = [ [] for _ in range(N) ]
	J2 = [ [] for _ in range(N) ]

	fields = [0]*N
	i1 = 0
	j1 = 0
	for _, val in enumerate(stoichiometricMatrixAsList):
		if val != 0:
			prendo=dim
			if dim > 0:
				for i in range(0, dim):
					if j1 == matched[i]:
						prendo = i
			if prendo < dim:
				matrix[prendo].append(i1)
				matrix2[prendo].append(val)		
		j1 += 1
		if j1 == N:
			j1 = 0 
			i1 += 1

	for i in range(0, dim):
		for j in range(i, dim):
			interactions = 0
			if len(matrix[i]) > 0:
				for po in range(0, len(matrix[i])):
					if len(matrix[j]) > 0:
						for pu in range(0, len(matrix[j])):
							if matrix[i][po] == matrix[j][pu]:
								interactions += (matrix2[i][po]*matrix2[j][pu])
					if j == i:
						fields[i] = interactions
					else:
						if abs(interactions) > MIN:
							J[i].append(j)
							J2[i].append(interactions)
							J[j].append(i)
							J2[j].append(interactions)
	return J, J2, fields

def LinearDependence(vectors, intkerneldim, NSolutions, NSolutions2, matched):
	""" 
	Check if the solution found with MonterCarlo is linearly independent 
	with respect to the previous found solution

	:param vectors:
		vectors
	:param intkerneldim:
		number of integer conservative laws
	:param NSolutions:
		NSolutions
	:param NSolutions2:
		NSolutions2
	:param matched:
		actual found MCLs
	"""
	K = intkerneldim+1
	MIN = 1e-9
	MAX = 1e+9
	matrix = [ [] for _ in range(K) ]
	matrix2 = [ [] for _ in range(K) ]
	for i in range(0, K-1):
		for j in range(0, len(NSolutions[i])):
			matrix[i].append(NSolutions[i][j])
			matrix2[i].append(NSolutions2[i][j])
	orders2 = [None] * len(matched)
	pivots2 = [None] * len(matched)

	for i in range(0, len(matched)):
		orders2[i] = i
		pivots2[i] = matched[i]

	qsort(len(matched), 0, orders2, pivots2)
	for i in range(0, len(matched)):
		if vectors[orders2[i]] > MIN:
			matrix[K-1].append(matched[orders2[i]])
			matrix2[K-1].append(float(vectors[orders2[i]]))

	ok = 0
	orders = [None] * K
	for i in range(0, K):
		orders[i] = i

	pivots = [None] * K
	for i in range(0, K):
		if len(matrix[i]) > 0:
			pivots[i] = matrix[i][0]
		else:
			pivots[i] = MAX

	while ok == 0:
		qsort(K, 0, orders, pivots)
		for j in range(0, K-1):
			if (pivots[orders[j+1]] == pivots[orders[j]]) and (pivots[orders[j]] != MAX):
				min1=MAX
				if len(matrix[orders[j]]) > 1:
					for i in range(0, len(matrix[orders[j]])):
						if (abs(matrix2[orders[j]][0] / matrix2[orders[j]][i])) < min1:
							min1 = abs(matrix2[orders[j]][0]/matrix2[orders[j]][i])
				min2=MAX
				if len(matrix[orders[j+1]]) > 1:
					for i in range(0, len(matrix[orders[j+1]])):
						if (abs(matrix2[orders[j+1]][0]/matrix2[orders[j+1]][i])) < min2:
							min2 = abs(matrix2[orders[j+1]][0] / matrix2[orders[j+1]][i])
				if min2 > min1:
					k2 = orders[j+1]
					orders[j+1] = orders[j]
					orders[j] = k2
		ok = 1
		for j in range(0, K-2):
			if (pivots[orders[j+1]] == pivots[orders[j]]) and (pivots[orders[j]] != MAX):
				k1 = orders[j+1]
				k2 = orders[j]
				colonna = [None] * N
				for i in range(0, N):
					colonna[i] = 0
				g = matrix2[k2][0]/matrix2[k1][0]
				for i in range(1, len(matrix[k1])):
					colonna[matrix[k1][i]] = matrix2[k1][i]*g
				for i in range(1, len(matrix[k2])):
					colonna[matrix[k2][i]] -= matrix2[k2][i]

				matrix[k1] = []
				matrix2[k1] = []
				for i in range(0, N):
					if abs(colonna[i]) > MIN:
						matrix[k1].append(i)
						matrix2[k1].append(colonna[i])
				ok = 0
				if len(matrix[k1]) > 0:
					pivots[k1] = matrix[k1][0]
				else:
					pivots[k1] = MAX
		
		K1 = 0
		for i in range(0, K):
			if len(matrix[i]) > 0:
				K1 += 1
		yes = 0
		if (K == K1):
			yes = 1
		return yes

def MonteCarlo(matched, J, J2, fields, intmatched, intkerneldim, NSolutions, NSolutions2, kerneldim, initT=1, coolrate=1e-3, maxIter=10):
	"""
	MonteCarlo simulated annealing 
	
	:param matched:
		matched
	:param J:
		J
	:param J2:
		J2
	:param fields:
		fields
	:param intmatched:
		actual matched MCLs 
	:param intkerneldim:
		number of MCLs found in S
	:param NSolutions:
		NSolutions
	:param NSolutions2:
		NSolutions2
	:param kerneldim:
		kerneldim
	:param initT:
		initial temperature
	:param coolrate:
		cooling rate of simulated annealing
	:param maxIter:
		maximum number of MonteCarlo steps before changing to relaxation

	"""
	MIN = 1e-9
	dim = len(matched)
	num = [0] * dim
	numtot = 0
	for i in range(0, dim):
		if (len(J[i])) > 0:
			num[i] = int(2*random.uniform(0, 1))
		else: 
			num[i] = 0
		numtot += num[i]
	
	H = 0
	for i in range(0, dim):
		H += fields[i]*num[i]*num[i]
		if len(J[i]) > 0:
			for j in range(0, len([J[i]])):
				H += J2[i][j] * num[i] * num[J[i][j]]

	stop = 0
	count = 0
	T1 = initT
	howmany = 0
	e = math.exp(-1/T1)

	while True:
		en = int(random.uniform(0, 1) * dim)
		if len(J[en]) == 0: # TODO: bug in original c++ code (while loop without any side effect changed to if statement to prevent infinite loop)
			en = int(random.uniform(0, 1)*dim)
		p=1
		if num[en] > 0 and random.uniform(0, 1) < 0.5:
			p = -1
		delta = fields[en] * num[en]
		for i in range(0, len(J[en])):
			delta += J2[en][i] * num[J[en][i]]
		delta = 2*p*delta + fields[en]

		if delta < 0 or random.uniform(0, 1) < math.pow(e, delta):
			num[en] += p
			numtot += p
			H += delta

		count += 1

		if count % int(dim) == 0:
			T1 -= coolrate
			if (T1 <= 0):
				T1 = coolrate
				e = math.exp(-1/T1)
		
		if count == int(float(dim)/coolrate):
			T1 = initT
			e = math.exp(-1/T1)
			count = 0
			numtot = 0
			for i in range(0, dim):
				num[i] = 0
			en = int(random.uniform(0, 1) * dim)
			if len(J[en]) > 0: # TODO: bug in original c++ code (while loop without any side effect changed to if statement to prevent infinite loop)
				en = int(random.uniform(0, 1) * dim)
			num[en] = 1
			numtot = 1
			H = 0
			for i in range(0, dim):
				H += fields[i] * num[i] * num[i]
				if len(J[i]) > 0:
					for j in range(0, len(J[i])):
						H += J2[i][j] * num[i] * num[J[i][j]]
			howmany += 1

		if (H < MIN and numtot > 0) or (howmany == (10*maxIter)):
			stop = 1

		if stop != 0:
			break

	yes = 0
	if howmany < 10*maxIter:
		if len(intmatched) > 0:
			yes = LinearDependence(num, intkerneldim, NSolutions, NSolutions2, matched)
			assert (yes == True), "Not true!"
		else:
			yes = 1
		if yes == 1:
			orders2 = [None] * len(matched)
			pivots2 = [None] * len(matched)
			for i in range(0, len(matched)):
				orders2[i] = i
				pivots2[i] = matched[i]
			print("sort...")
			qsort(len(matched), 0, orders2, pivots2)
			print("end")
			for i in range(0, len(matched)):
				if num[orders2[i]] > 0:
					NSolutions[intkerneldim].append(matched[orders2[i]])
					NSolutions2[intkerneldim].append(num[orders2[i]])
			intkerneldim += 1
			yes2  = LinearDependence(num, intkerneldim, NSolutions, NSolutions2, matched) # TODO: yes2 never used, does LinearDependence have side effects?
			intkerneldim, kerneldim, NSolutions, NSolutions2 = Reduce(intkerneldim, kerneldim, NSolutions, NSolutions2)
			min = 1000
			for i in range(0, len(NSolutions[intkerneldim-1])):
				if len(intmatched) == 0:
					intmatched.append(NSolutions[intkerneldim-1][i])
				else:
					ok3 = 1
					for k in range(0, len(intmatched)):
						if (intmatched[k] == NSolutions[intkerneldim-1][i]):
							ok3 = 0
					if ok3 == 1:
						intmatched.append(NSolutions[intkerneldim-1][i])
				if NSolutions2[intkerneldim-1][i] < min:
					min = NSolutions2[intkerneldim-1][i]
			for i in range(0, len(NSolutions[intkerneldim-1])):
				NSolutions2[intkerneldim-1][i] /= min
			print(f"Found linearly independent moiety, now there are {intkerneldim} engaging {len(intmatched)} metabolites")
		else:
			print("Found a moiety but it is linearly dependent... next.")
	else:
		yes = 0
	return yes, intkerneldim, kernelDim, NSolutions, NSolutions2, matched, intmatched
					


def Relaxation(stoichiometricMatrixAsList, intmatched, M, N, relaxationmax=1e6, relaxation_step=1.9):
	""" 
	Relaxation scheme for MonteCarlo final solution

	:param stoichiometricMatrixAsList:
		stoichiometric matrix as a list
	:param intmatched:
		intmatched
	:param M:
		number of metabolites
	:param N:
		number of reactions
	:param relaxationmax:
		maximum relaxation
	:param relaxation_step:
		relaxation step width
	"""
	MIN = 1e-9
	MAX = 1e9
	matrix = [ [] for _ in range(N) ]
	matrix2 = [ [] for _ in range(N) ]

	i1 = 0
	j1 = 0
	K = len(intmatched)
	for _, val in enumerate(stoichiometricMatrixAsList):
		if val != 0:
			prendo=K
			if K > 0:
				for i in range(0, K):
					if j1 == intmatched[i]:
						prendo = i
			if prendo < K:
				matrix[prendo].append(i1)
				matrix2[prendo].append(val)		
		j1 += 1
		if j1 == N:
			j1=0
			i1 += 1

	ok = 0

	orders = [ i for i in range(0, N) ]
	pivots = [ matrix[i][0] if len(matrix[i]) > 0 else MAX for i in range(0, N)]

	while ok == 0:
		qsort(K, 0, orders, pivots)
		for j in range(0, K):
			if pivots[orders[j+1]] == pivots[orders[j]] and pivots[orders[j]] != MAX:
				min1 = MAX
				if len(matrix[orders[j]]) > 1:
					for i in range(0, len(matrix[orders[j]])):
						if abs(matrix2[orders[j]][0]/matrix2[orders[j]][i]) < min1:
							min1 = matrix2[orders[j]][0]/matrix2[orders[j]][i]
				min2 = MAX
				if len(matrix[orders[j+1]]) > 1:
					for i in range(0, len(matrix[orders[j]])):
						if abs(matrix2[orders[j+1]][0]/matrix2[orders[j+1]][i]) < min2:
							min2 = abs(matrix2[orders[j+1]][0]) / matrix2[orders[j+1]][i]
				if min2 > min1:
					k2 = orders[j+1]
					orders[j+1] = orders[j]
					orders[j] = k2
		ok = 1
		j = 0
		for j in range(0, K):
			if pivots[orders[j+1]] == pivots[orders[j]] and pivots[orders[j]] != MAX:
				k1 = orders[j+1]
				k2 = orders[j]
				colonna = [None] * M 
				for i in range(0, M):
					colonna[i] = 0
				g = matrix2[k2][0]/matrix2[k1][0]
				for i in range(1, len(matrix[k1])):
					colonna[matrix[k1][i]] = matrix2[k1][i]*g
				for i in range(1, len(matrix[k2])):
					colonna[matrix[k2][i]] -= matrix2[k2][i]

				matrix[k1] = []
				matrix2[k1] = []
				for i in range(0, M):
					if abs(colonna[i]) > MIN:
						matrix[k1].append(i)
						matrix2[k1].append(colonna[i])
				ok = 0
				if len(matrix[orders[j+1]]) > 0:
					pivots[orders[j+1]] = matrix[orders[j+1]][0]
				else:
					pivots[orders[j+1]] = MAX

		for i in range(0, K):
			if len(matrix[i]) > 0:
				norm = matrix2[i][0]
				for j in range(0, len(matrix[i])):
					matrix2[i][j] /= norm

		for k1 in reversed(range(0, K-1)):
			k = orders[k1]
			if len(matrix[k]) > 1:
				for i in range(1, len(matrix[k])):
					for j1 in range(k1+1, K):
						j = orders[j1]
						if len(matrix[j]) > 0:
							if matrix[j][0] == matrix[k][i]:
								rigak = [None] * M
								for a in range(0, M):
									rigak[a] = 0
								for a in range(0, len(matrix[k])):
									rigak[matrix[k]][a] = matrix2[k][a]
								for a in range(0, len(matrix[j])):
									rigak[matrix[j]][a] -= matrix2[j][a] * matrix2[k][i]
								matrix = []
								matrix2 = []
								for a in range(0, M):
									if rigak[a] != 0:
										matrix[k].append(a)
										matrix2[k].append(rigak[a])

		indip = [None] * M
		for i in range(0, M):
			indip[i] = K+1
		
		for i in range(0, K):
			if len(matrix[i]) > 0:
				indip[matrix[i]][0] = i
		
		M1 = 0
		for i in range(0, M):
			if indip[i] == K+1:
				indip[i] = K + M1
				M1 += 1
		
		matrixAus = [ [] for _ in range(M1)]
		matrixAus2 = [ [] for _ in range(M1)]
		i1 = 0
		for i in range(0, M):
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
		
		for i in range(0, K):
			matrix[i] = []

		N1 = N-K
		matrix_aus = [ [] for _ in range(N1)]
		matrix_aus2 = [ [] for _ in range(N1)]

		k1 = 0
		i1 = 0
		j1 = 0
		
		for _, val in enumerate(stoichiometricMatrixAsList):
			prendo = 1
			if len(intmatched) > 0:
				for i in range(0, len(intmatched)):
					if j1 == intmatched[i]:
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

		matrixb = [ [] for _ in range(N1)]
		matrixb2 = [ [] for _ in range(N1)]
		for i in range(0, M1):
			for j in range(0, N1):
				prod = 0
				if len(matrix_aus[j])*len(matrixAus[i]) > 0:
					for ib in range(0, len(matrixAus[i])):
						for jb in range(0, len(matrix_aus[j])):
							if matrixAus[i][ib] == matrix_aus[j][jb]:
								prod += matrixAus2[i][ib] * matrix_aus2[j][jb]
					if abs(prod) > MIN:
						matrixb[j].append(i)
						matrixb2[j].append(prod)

		for i in range(0, M1):
			matrixAus[i] = []
			matrixAus2[i] = []

		for i in range(0, N1):
			matrix_aus[i] = []
			matrix_aus2[i] = []

		var = [None] * M1
		for i in range(0, M1):
			var[i] = MIN
		ok = 0
		time = 0
		while True:
			cmin = 1000
			for j in range(0, N1):
				constr = 0
				if len(matrixb[j]) > 0:
					for i in range(0, len(matrixb[j])):
						constr += matrixb2[j][i] * var[matrixb][j][i]
					if constr < cmin:
						min = j
						cmin = constr
			time += 1
			if cmin >= 0:
				ok = 1
			else:
				alpha = -relaxation_step*cmin # Motzkin relaxation
				fact = 0
				for j in range(0, len(matrixb[min])):
					fact += matrixb2[min][j]*matrixb2[min][j]
				alpha /= fact
				if alpha < 1e-9 * MIN:
					alpha = 1e-9 * MIN
				for j in range(0, len(matrixb[min])):
					var[matrixb[min][j]] += alpha * matrixb2[min][j]

			if not (ok == 0 and time < relaxationmax):
				break

		yes = 0
		if ok == 1:
			yes = 1
		return yes


def Reduce(intKernelDim, kernelDim, NSolutions, NSolutions2):
	""" 
	Reducing the solution found by MonteCarlo
	
	:param intKernelDim:
		number of found MCLs
	:param kernelDim:
		number of found conservative laws
	:param NSolutions:
		NSolutions
	:Param NSolutions2:
		NSolutions2
	"""
	K = intKernelDim
	MIN = 1e-9
	ok = 0
	orders = [None] * K
	for i in range(0, K):
		orders[i] = i
	pivots = [None] * K
	for i in range(0, K):
		pivots[i] = -len(NSolutions[i])

	while True:
		qsort(K, 0, orders, pivots)
		ok = 1
		for i in range(0, K-2):
			for j in range(i+1, K):
				k1 = orders[i]
				k2 = orders[j]
				colonna = [None] * N
				for l in range(0, N):
					colonna[l] = 0
				ok1 = 1
				for l in range(0, len(NSolutions[k1])):
					colonna[NSolutions[k1][l]] = NSolutions2[k1][l]
				for l in range(0, len(NSolutions[k2])):
					colonna[NSolutions[k2][l]] -= NSolutions2[k2][l]
					if colonna[NSolutions[k2][l]] < -MIN:
						ok1 = 0
				if ok1 == 1:
					ok = 0
					NSolutions[k1] = []
					NSolutions2[k1] = []
					for l in range(0, N):
						if abs(colonna[l]) > MIN:
							NSolutions[k1].append(l)
							NSolutions2[k1].append(colonna[l])
					pivots[k1] = -len(NSolutions[k1])
		if ok != 0:
			break
	return intKernelDim, kernelDim, NSolutions, NSolutions2

def Input():
	""" Read input from test stoichiometric matrix """
	with open('matrix.dat', 'r') as f:
		return [int(item) for sl in [entry.strip().split('\t') for entry in f] for item in sl]

def Names():
	""" Get names of metabolites from list of metabolites example"""
	with open('metabolites.txt', 'r') as f:
		return [entry for entry in f]
			
def Output(intKernelDim, kernelDim, intmatched, NSolutions):
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
	print(f"There are {intKernelDim} linearly independent conserved moieties, engaging {len(intmatched)} metabolites\n")
	if intKernelDim == kernelDim:
		print("They generate all the conservation laws")
	else:
		print(f"They don't generate all the conservation laws, {kernelDim-intKernelDim} of the mare not reducible to moeities")

	names = Names()
	for i in range(0, intKernelDim):
		print(f"Moeity number {i+1} engages {len(NSolutions[i])} metabolites.")
		print(f"\t")
		for j in range(0, len(NSolutions[i])):
			print(f"{names[NSolutions[i][j]]} \t {NSolutions2[i][j]}")

if __name__ == "__main__":
	""" Invoke test case and benchmarking for De Martino's publishes results for E-Coli network """
	start = time.time()
	print("".join(['*' for _ in range(0, 80)]))
	print("Conserved moeties test cases")
	print("".join(['*' for _ in range(0, 80)]))

	# Hard-coded stoichiometric matrix as test case containing _ONE_ conservative law (i.e. one conserved moiety)
	S  = [ -1, 0, 0, 0, 0, 0,
			1, -1, 1, 0, 0, 0,
			0, 1, -1, -1, -1, -1,
			0, 0, 0, 1, -1, 0,
			0, 0, 0, 0, 0, 1]
	N = 6 # number of metabolites (columns)
	M = 5 # number of reactions (rows)

	S = Input()
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
			# TODO: Sometimes MC runs into the wrong solution: Why?
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
