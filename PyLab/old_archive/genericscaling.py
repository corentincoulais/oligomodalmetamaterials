"Aleksi Bossart, 9.4.2019"
"Plot the average number of modes as a function of size for random configurations"

import numpy as np
from scipy import linalg as algebra
from matplotlib import pyplot as plt


"load functions defined in mechadef.py and autopenta.py"

from mechadef import drawlattice
from mechadef import compatibility
from mechadef import kernel
from mechadef import posadj


"initialise the mode number matrix"

maxdim = 12
trials = 300
modenumber = np.zeros(maxdim-1)


"specify the primitive cell"

tiling = 'fullspectrum'


"load the position file for the primitive cell"

cellpos = np.genfromtxt('lattices/positions/primitivecell/lieb', delimiter=',')



"loop over n and m, produce appropriate aperiodicity and generate positions and adjacency from there"

progress = ''
remaining = maxdim * ' \u2235 '

print('{' + progress + remaining + '}')

for n in range(1, maxdim):

	for k in range(trials):

		aperio = np.random.randint(4, size = (n, n))
		positions, adjacency = posadj(aperio, cellpos)


		"compute the compatibility matrix"

		cmat = compatibility(positions, adjacency)


		"compute the number of floppy modes and save it"

		modenumber[n-1] += kernel(cmat).shape[0] - 3

		#drawlattice(positions, adjacency)
		#plt.show()

	progress += ' \u22C8 '
	remaining = remaining[:-3]
	print('{' + progress + remaining + '}')


modenumber /= trials


"plot something"

syssize = np.arange(1, maxdim)
	

plt.title('Mode scaling for various tilings')
plt.xlabel = ('system size')
plt.xlabel = ('number of modes')

plt.plot(syssize[1::2].astype(int), 2 * syssize[1::2] / 1 + 0, color='saddlebrown', label='p4m', marker = 'o', lw = 0)
plt.plot(syssize[1::2].astype(int), 3 * syssize[1::2] / 2 + 0, color='sienna', label='pm', marker = 'o', lw = 0)
plt.plot(syssize[1::2].astype(int), 1 * syssize[1::2] / 1 + 1, color='peru', label='cm:2,cm:3, pmg:2', marker = 'o', lw = 0)
plt.plot(syssize[1::2].astype(int), 1 * syssize[1::2] / 2 + 2, color='burlywood', label='p2', marker = 'o', lw = 0)
plt.plot(syssize[1::2].astype(int), 1 * syssize[1::2] / 2 + 1, color='navajowhite', label='p1', marker = 'o', lw = 0)

plt.plot(syssize[1::2].astype(int), 0 * syssize[1::2] + 4, color='lightskyblue', label='p4g', marker = 'o', lw = 0)
plt.plot(syssize[1::2].astype(int), 0 * syssize[1::2] + 2, color='deepskyblue', label='cm:4', marker = 'o', lw = 0)
plt.plot(syssize[1::2].astype(int), 0 * syssize[1::2] + 1, color='dodgerblue', label='pmg', marker = 'o', lw = 0)

plt.plot(syssize.astype(int), modenumber, color='darkblue', label='generic', marker = 'o', lw = 0)
plt.legend()
plt.savefig('results/figures/modescaling/' + tiling)
plt.show()


