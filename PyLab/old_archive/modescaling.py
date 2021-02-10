"Aleksi Bossart, 28.3.2019"
"Plot the number of modes as a function of size for a given 2x2 primitive cell"

import numpy as np
from scipy import linalg as algebra
from matplotlib import pyplot as plt


"load functions defined in mechadef.py and autopenta.py"

from mechadef import drawlattice
from mechadef import compatibility
from mechadef import kernel
from mechadef import posadj


"initialise the mode number matrix"

maxdim = 8

modenumber = np.zeros(shape = (maxdim-1, maxdim-1))
prediction = np.zeros(shape = (maxdim-1, maxdim-1))


"specify the primitive cell"

tiling = 'p4g'

primcell = np.transpose(np.genfromtxt('lattices/aperiodicity/penta/primcell/' + tiling, delimiter=','))[:, ::-1]

def lowerbound(n,m):
	return 4;



"load the position file for the primitive cell"

cellpos = np.genfromtxt('lattices/positions/primitivecell/lieb', delimiter=',')


"prediction of the number of modes from vertex model"

for n in range(1, maxdim):

	for m in range(1, maxdim):

		prediction[m-1, n-1] = lowerbound(2*n, 2*m)


"loop over n and m, produce appropriate aperiodicity and generate positions and adjacency from there"

for n in range(1, maxdim):

	for m in range(1, maxdim):

		aperio = np.tile(primcell, (n,m))
		positions, adjacency = posadj(aperio, cellpos)

		"compute the compatibility matrix"

		cmat = compatibility(positions, adjacency)


		"compute the number of floppy modes and save it"

		modenumber[m-1, n-1] = kernel(cmat).shape[0] - 3


		if m == maxdim-1 and n == 3:
			drawlattice(positions, adjacency)
			plt.show()


"plot something"

syssize = 2 * np.arange(1, maxdim)
	

plt.title('Mode scaling for the ' + tiling + ' tiling')
plt.xlabel = ('system size')
plt.xlabel = ('number of modes')
plt.plot(syssize, prediction[2, :], '-', color = 'k', label='vertex prediction')
plt.plot(syssize, prediction[:, 2], '--', color='k', label='vertex prediction')
plt.plot(syssize, modenumber[2, :], marker='o', lw=0, color='r', label='size (6,n)')
plt.plot(syssize, modenumber[:, 2], marker='^', lw=0, color='b', label='size (m,6)')
plt.legend()
plt.savefig('results/figures/sanitycheck/' + tiling)
plt.show()


"save the data as a csv"

np.savetxt('results/modescaling/' + tiling, modenumber.astype(int), delimiter = ',', fmt="%d")


