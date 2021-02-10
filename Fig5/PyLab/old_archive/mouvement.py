
"""
Matplotlib Animation Example

author: Jake Vanderplas
email: vanderplas@astro.washington.edu
website: http://jakevdp.github.com
license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!
"""

"Modified by Aleksi Bossart, 16.4.2019"
"Interactive multimodal mechanism"


from numpy import sin, cos
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation


"load functions defined in mechadef.py"

from mechadef import drawlattice
from mechadef import compatibility
from mechadef import kernel

positions = np.genfromtxt('lattices/positions/gen/pentagone', delimiter=',')
posx, posy = positions.reshape((2, int(np.size(positions)/2) ))[0], positions.reshape((2, int(np.size(positions)/2) ))[1]
adjacency = np.genfromtxt('lattices/adjacency/gen/pentagone', delimiter=',')

#------------------------------------------------------------
# set up initial state and global variables

dt = 1./30 # 30 fps

#------------------------------------------------------------
# set up figure and animation
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(-2, 2), ylim=(-2, 2))
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)

def init():
    """initialize animation"""
    line.set_data([], [])
    return line

def animate(i):
    """perform animation step"""
    line.set_data(posx, posy)
    return line


# choose the interval based on dt and the time to animate one step

from time import time

t0 = time()
animate(0)
t1 = time()
interval = 1000 * dt - (t1 - t0)

ani = animation.FuncAnimation(fig, animate, frames=300, interval=interval, blit=True, init_func=init)

plt.show()
