"Aleksi Bossart, 23.4.2019"
"Generate the next level of L-substituted lattice"


import numpy as np
from matplotlib import pyplot as plt

from mechadef import drawlattice
from mechadef import posadj


"load the aperiodicity matrix"

nw = np.transpose(np.genfromtxt('lattices/aperiodicity/penta/toy', delimiter=','))[:, ::-1]
#nw = np.genfromtxt('lattices/aperiodicity/penta/L', delimiter=',')

"generate the 3 other block matrices"

se = np.copy(nw)

ne = np.flip(np.copy(nw), axis = 0)
ne[ne == 0] = 1.5
ne[ne == 1] = 0.5
ne[ne == 2] = 3.5
ne[ne == 3] = 2.5
ne -= 0.5

sw = np.flip(np.copy(nw), axis = 1)
sw[sw == 0] = 3.5
sw[sw == 1] = 2.5
sw[sw == 2] = 1.5
sw[sw == 3] = 0.5
sw -= 0.5

"concatenate and save"

aperio = np.block([[sw, nw], [se, ne]])
np.savetxt('lattices/aperiodicity/penta/L', aperio, delimiter = ',')


"load the position file for the primitive cell"

cellpos = np.genfromtxt('lattices/positions/primitivecell/lieb', delimiter=',')


"compute posadj"

positions, adjacency = posadj(aperio, cellpos)


"draw the lattice"

plt.rcParams["figure.figsize"] = [18,10]
drawlattice(positions, adjacency, rayon = 0.08)

plt.show()


