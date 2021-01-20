"Aleksi Bossart, 26.2.2019"
"semi-automated generation of arbitrary aperiodic lattices"


import numpy as np
from matplotlib import pyplot as plt
import sys
from sys import platform as _platform

path_sys={"darwin":"/Users/coulais/science/Shared/Metacombinatorial/",
          "win32":"TBD",
          "linux":"TBD"}
#sys.path.append(path_sys[_platform])
sys.path.append(path_sys[_platform]+'PyLab')
#import mechadef
from mechadef import drawlattice
from mechadef import posadj

workingpath=path_sys[_platform]+"StructuralAnalysis/"


"load the aperiodicity matrix"

aperio = np.transpose(np.genfromtxt(workingpath+'lattices/aperiodicity/penta/toy', delimiter=','))[:, ::-1]


"load the position file for the primitive cell"

cellpos = np.genfromtxt(workingpath+'lattices/positions/primitivecell/lieb', delimiter=',')


"compute posadj"

positions, adjacency = posadj(aperio, cellpos,fileout_gen=workingpath+'lattices/adjacency/gen/pentagone',fileout_pos=workingpath+'lattices/positions/gen/pentagone')


"draw the lattice"

plt.rcParams["figure.figsize"] = [18,10]
drawlattice(positions, adjacency, rayon = 0.08)

plt.show()


