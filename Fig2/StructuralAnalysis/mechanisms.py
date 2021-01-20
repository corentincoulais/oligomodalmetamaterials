"Aleksi Bossart, 19.2.2019"
"this program plots all the floppy modes of an arbitrary mechanical lattice"

import numpy as np
from scipy import linalg as algebra
from matplotlib import pyplot as plt

import sys
from sys import platform as _platform

path_sys={"darwin":"/Users/coulais/science/Shared/Metacombinatorial/",
          "win32":"TBD",
          "linux":"TBD"}
#sys.path.append(path_sys[_platform])
sys.path.append(path_sys[_platform]+'PyLab')

workingpath=path_sys[_platform]+"StructuralAnalysis/"


"load functions defined in mechadef.py"

from mechadef import drawlattice
from mechadef import compatibility
from mechadef import kernel
from mechadef import gramaway
from mechadef import hourglass
from mechadef import colorcode
from mechadef import rotatebasis



positions = np.genfromtxt(workingpath+'lattices/positions/gen/pentagone', delimiter=',')
adjacency = np.genfromtxt(workingpath+'lattices/adjacency/gen/pentagone', delimiter=',')


"draw the original lattice"

#plt.rcParams["figure.figsize"] = [16,9]
#drawlattice(positions, adjacency)
#plt.close()


"compute the kernel of the compatibility matrix with some fixed points, to avoid rigid body motions"

fix = 4

cmat = compatibility(positions, adjacency)
nullvec = kernel(cmat[:, fix:])

n1, n2 = np.shape(nullvec)

zeromodes = np.zeros(n1*(n2 + fix))
zeromodes = zeromodes.reshape(n1, n2 + fix)
zeromodes[:, fix:] += nullvec

n2 = n2 + fix

"save cmat"

#np.savetxt('results/compatibility/pentagone', cmat.astype(int), delimiter = ',', fmt="%d")


"remove known floppy modes, typically counter-rotations, then re-orthogonalize the result"

#counter = np.genfromtxt('lattices/zeromodes/gen/bigcorner1', delimiter=',')
#zeromodes = gramaway(zeromodes, counter)
#zeromodes = np.transpose(algebra.orth(np.transpose(zeromodes)))
#n1, n2 = np.shape(zeromodes)


"enforce periodic boundary conditions by orthogonalizing against incoherent moves"
"input: boundary vectors"
#TODO: protect against wrong inputs (true for all of the code, actually)

#top = np.genfromtxt('lattices/boundaries/gen/penta/top', delimiter=',')
#bottom = np.genfromtxt('lattices/boundaries/gen/penta/bottom', delimiter=',')

#bottop = np.size(top)

#left = np.genfromtxt('lattices/boundaries/gen/penta/left', delimiter=',')
#right = np.genfromtxt('lattices/boundaries/gen/penta/right', delimiter=',')

#lefrite = np.size(left)

#divmat = np.zeros(n2 * 2 * (bottop + lefrite))
#divmat = np.reshape(divmat, newshape = (2 * (bottop + lefrite), n2))

#for i in range(bottop):

	#divmat[i, 2 * int(top[i])] = 1
	#divmat[i, 2 * int(bottom[i])] = -1

	#divmat[i, 2 * int(top[i]) + 1] = 1
	#divmat[i, 2 * int(bottom[i]) + 1] = -1

#for i in range(lefrite):

	#divmat[2 * bottop + i, 2 * int(left[i])] = 1
	#divmat[2 * bottop + i, 2 * int(right[i])] = -1

	#divmat[2 * bottop + i, 2 * int(left[i]) + 1] = 1
	#divmat[2 * bottop + i, 2 * int(right[i]) + 1] = -1

#print(divmat)

#for divec in divmat[0:3]:

	#zeromodes = gramaway(zeromodes, divec)
	
#zeromodes = np.transpose(algebra.orth(np.transpose(zeromodes)))
#n1, n2 = np.shape(zeromodes)


"save a floppy mode for later use"

np.savetxt(workingpath+'lattices/zeromodes/gen/coins', zeromodes, delimiter = ',')

"draw the corresponding zero modes with some small displacement delta"

floppy = zeromodes.reshape((n1, int(n2/2), 2))
delta = 5

print("There are %d floppy modes" % n1)

plt.rcParams["figure.figsize"] = [17,9]

for i in range(0,n1,1):

	modpos = positions - delta*floppy[i]

	drawlattice(positions, adjacency, colour = 'orange', name = "axes%d" % i)
	drawlattice(modpos, adjacency, name = "axes%d" % i)
	
	"add color code for hourglass deformations"

	#poi = np.arange(75,144,6)
	#for k in np.arange(11):
		#poi = np.concatenate(( poi , np.arange(k * 144 + 150, k * 144 + 211,6), np.arange(k * 144 + 219, k * 144 + 286,6) ))
	#sabliers = np.array([modpos[i] for i in poi])

	#deformations = hourglass(modpos, sabliers + np.array([1,1]))
	#colorcode(sabliers  + np.array([1,1]), deformations, rayon = 0.8, name = "axes%d" % i)

	#print(deformations)

	#plt.show()
	#plt.close()
	
	#diagonale = np.array([210, 348, 486, 624, 762, 900, 1038, 1176, 1314, 1452, 1590])
	#diagonale = np.array([75, 225, 375, 525, 675, 825, 975, 1125, 1275, 1425, 1575, 1725])

	#deformations = hourglass(modpos, np.array([modpos[i] for i in diagonale]))
	#colorcode(np.array([modpos[i] for i in diagonale]), deformations, rayon = 0.8, name = "axes%d" % i)

	plt.show()
	#  plt.close()

	#horloges = np.array([modpos[i] for i in diagonale])
	#decaysab = hourglass(modpos, horloges)
	#xpos = np.arange(12)

	#hourslope, hourigin = np.polyfit(np.log(xpos[4:]), np.log(decaysab[4:]), 1)
	#print('The polynomial exponent for the hourglass parameter is %f' % hourslope)

	#print('The sum-hourglass parameter is %f' % np.sum(decaysab))

	#plt.plot(xpos, decaysab, 'o')
	#plt.plot(np.log(xpos), hourslope * np.log(xpos) + hourigin, 'r')
	#plt.show()


