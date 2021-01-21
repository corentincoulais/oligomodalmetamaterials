"Aleksi Bossart, 6.8.2019, last modified 21.1.21"
"Plot the number of modes as a function of size for tilings specified in folder fullaperio, or averaging over random ones"

import numpy as np
from scipy import linalg as algebra
from matplotlib import pyplot as plt


"load functions to compute adacency, compatibility matrix and draw lattices"

from mechadef import drawlattice
from mechadef import compatibility
from mechadef import posadj
from mechadef import twist

"initialise the mode number matrix: give the minimal and maximal tiling dimensions for which the number of zero modes is to be computed"

mindim = 1
maxdim = 16

"number of trials over which to average if investigating random tilings"

trials = 1

modenumber = np.zeros(maxdim)


"specify the primitive cell, with numbers from 0 to 3 corresponding to the four possible orientations of the unit cell."

tiling = '1203'


"load the position file for the vertices primitive cell"

cellpos = np.genfromtxt('primitivecell', delimiter=',')


"load orientations of the maximal lattice, comment if investigating random lattices"

fullaperio = np.transpose(np.genfromtxt('fullaperio/'+tiling, delimiter=','))[:, ::-1]


"loop over n and m, produce appropriate aperiodicity and generate positions and adjacency from there"

progress = ''
remaining = maxdim * ' \u2235 '


print('{' + progress + remaining + '}')

for n in range(mindim, maxdim+1, 1):

	for k in range(trials):

		"comment for random lattices"

		aperio = fullaperio[:n, (32-n):]
		
		"uncomment for random lattices"
		
		#aperio = np.random.randint(4, size=(n,n))

		positions, adjacency = posadj(aperio, cellpos)
		
		"uncomment to pre-twist the lattice by a certain amount"
		
		#positions = twist(positions, 0.24)

		"compute the compatibility matrix"

		cmat = compatibility(positions, adjacency)


		"compute the number of floppy modes and save it"

		modenumber[n-1] += np.shape(cmat)[1] - np.linalg.matrix_rank(cmat) - 3
		

		"draw the lattice if needed, set if n == 5, for example"
		
		if False:
			drawlattice(positions, adjacency)
			plt.show()


	progress += ' \u22C8 '
	remaining = remaining[:-3]
	print('{' + progress + remaining + '}')


modenumber /= trials


"plot the number of zero modes against tiling side length"

np.savetxt('data/'+tiling+'_test.out', modenumber, delimiter=',')

syssize = np.arange(1, maxdim+1)
	

plt.title('Mode scaling')
plt.xlabel('system size')
plt.ylabel('number of modes')

plt.plot(syssize.astype(int), modenumber, '.', color='darkblue')

plt.savefig('data/'+tiling+'_test')
plt.show()


