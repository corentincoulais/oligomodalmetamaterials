"Aleksi Bossart, 19.2.2019"
"this program plots all the floppy modes of an arbitrary mechanical lattice"

import numpy as np
from scipy import linalg as algebra
from matplotlib import pyplot as plt

"load functions defined in mechadef.py"

from mechadef import drawlattice
from mechadef import compatibility
from mechadef import kernel
from mechadef import gramaway
from mechadef import hourglass
from mechadef import colorcode
from mechadef import rotatebasis
from mechadef import drawstress

positions = np.genfromtxt('lattices/positions/gen/pentagone', delimiter=',')
adjacency = np.genfromtxt('lattices/adjacency/gen/pentagone', delimiter=',')


"compute the kernel of the transposed compatibility matrix to find states of self-stress"

qmat = np.transpose(compatibility(positions, adjacency))
selfstress = kernel(qmat)

n1, n2 = np.shape(selfstress)


"save a state of self-stress"

np.savetxt('lattices/selfstress/all1', selfstress[0], delimiter = ',')


"print/draw the state of self-stress"

print("There are %d states of self-stress" % n1)

plt.rcParams["figure.figsize"] = [17,9]

for i in range(0,n1,1):

	drawstress(positions, adjacency, selfstress[i])
	plt.show()
	print((np.abs(selfstress[i])>1e-15)*selfstress[i])



