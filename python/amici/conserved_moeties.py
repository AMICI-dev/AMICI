def qsort(k, km, orders, pivots):
	centre = 0
	if k-km >= 1:
		pivot = km+(k-km)//2
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
	pivots = [  matrix[i][0] if len(matrix[i]) > 0 else MAX for i in range(0, N)]

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
	for i in range(0, N-1):
		ok = 1
		if len(matrix[i]) > 0:
			for j in range(0, len(matrix[i])-1):
				if matrix[i][j] < M:
					ok = 0
		if ok == 1 and len(matrix[i]) > 0:
			for j in range(0, len(matrix[i])-1):
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
			for j in range(0, len(RSolutions[i])-1):
				if RSolutions2[i][j]*RSolutions2[i][0] < 0:
					ok2 = 0
				if len(matched) == 0:
					matched.append(RSolutions[i][j])
				else:
					ok3 = 1
					for k in range(0, len(matched)-1):
						if matched[k] == RSolutions[i][j]:
							ok3 = 0
						if ok3 == 1:
							matched.append(RSolutions[i][j])
		if ok2 == 1 and len(RSolutions[i]) > 0:
			min = MAX
			for j in range(0, len(RSolutions[i])-1):
				NSolutions[i2].append(RSolutions[i][j])
				NSolutions2[i2].append(abs(RSolutions2[i][j]))
				if min > abs(RSolutions2[i][j]):
					min = abs(RSolutions2[i][j])
				if len(intmatched) == 0:
					intmatched.append(NSolutions[i2][j])
				else:
					ok3 = 1
					for k in range(0, len(intmatched)-1):
						if intmatched[k] == NSolutions[i2][j]:
							ok3 = 0
			for j in range(0, len(NSolutions[i2])-1):
				NSolutions2[i2][j] /= min
			i2 += 1
	intKernelDim=i2

	return kernelDim, matched, intKernelDim, intmatched

def fill(stoichiometricMatrixAsList, matched_size, matched):
	MIN = 1e-9
	matrix = [ [] for _ in range(N) ]
	matrix2 = [ [] for _ in range(N) ]

	J = [ [] for _ in range(N) ]
	J2 = [ [] for _ in range(N) ]

	fields = [None]*N
	i1 = 0
	j1 = 0
	dim = matched_size
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
			j1=0
			i1 += 1

	for i in range(0, dim):
		for j in range(i, dim):
			interactions = 0
			if len(matrix[i]) > 0:
				for po in range(0, len(matrix[i])):
					if len(matrix[j]) > 0:
						for pu in range(0, len(matrix[j])):
							if matrix[i][po] == matrix[j][pu]:
								interactions += matrix2[i][po]*matrix2[j][pu]

				if j == 1:
					fields[i] = interactions
			else:
				if abs(interactions) > MIN:
					J[i].append(j)
					J2[i].append(interactions)
					J[j].append(i)
					J2[j].append(interactions)

def MonteCarlo():
	""" TODO: Implement """
	pass

def Relaxation():
	""" TODO: Implement """
	pass

def Reduce():
	""" TODO: Implement """
	pass

def Output():
	""" TODO: Implement """
	pass

if __name__ == "__main__":
	print("Conserved moeties test case...")

	# Hard-coded stoichiometric matrix as test case
	S = [1, 0, 3, 4,
		 1, 5, 1, 5, 
		-10, 3, 1, 6]
	N = 4 # number of metabolites (columns)
	M = 3 # number of reactions (rows)

	if len(S) != N*M:
		print("Stoichiometric matrix inconsistent")
	else:
		print("Kernel...")
		kernelDim, engagedMetabolites, intKernelDim, conservedMoieties = kernel(S, N, M)
		print(f"""There are {kernelDim} conservation laws, engaging, 
		{len(engagedMetabolites)} metabolites, {intKernelDim} are integers (conserved 
		moeities), engaging {len(conservedMoieties)} metabolites...""")

		print("Filling...")
		fill(S, len(engagedMetabolites), engagedMetabolites)

		print("MonteCarlo...")
		finish = 0
		if intKernelDim == kernelDim:
			finish = 1
		timer = 0
		while finish == 0:
			yes, intKernelDim, kernelDim = MonteCarlo()
			if intKernelDim == kernelDim:
				finish = 1
			if yes == 0:
				timer += 1
			if timer == max:
				print("Relaxation...")
				finish = Relaxation()
				if finish == 1:
					timer = 0

		print("Reduce...")
		Reduce()

		print("Pretty print output...")
		Output()
