"Aleksi Bossart, 15.4.2019"
"Find interesting 3x3 configurations"

import numpy as np
from scipy import linalg as algebra
from matplotlib import pyplot as plt


"load functions defined in mechadef.py and autopenta.py"

from mechadef import drawlattice
from mechadef import compatibility
from mechadef import kernel
from mechadef import posadj


"initialise the mode number matrix"
trials = 30
modenumber = 1


"load the position file for the primitive cell"

cellpos = np.genfromtxt('lattices/positions/primitivecell/lieb', delimiter=',')



"Produce random 3x3 aperiodicity and generate positions and adjacency from there"

progress = ''
remaining = trials * ' \u2235 '

print('{' + progress + remaining + '}')

for k in range(trials):

	aperio = np.tile(np.random.randint(4, size = (3, 3)), (4,4))
	positions, adjacency = posadj(aperio, cellpos)


	"compute the compatibility matrix"

	cmat = compatibility(positions, adjacency)


	"compute the number of floppy modes and plot if interesting"

	if kernel(cmat).shape[0] == 4:
		drawlattice(positions, adjacency)
		plt.show()

	progress += ' \u22C8 '
	remaining = remaining[:-3]
	print('{' + progress + remaining + '}')




