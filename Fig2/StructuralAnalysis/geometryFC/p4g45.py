#Modified by Aleksi on 3/4/19 to produce geometries for multimodal tilings

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import matplotlib as mpl

def p4gtile(ax=None, a=8., b=8., Nx=5, Ny=7, theta0=np.pi/4, t=0.3, blow=0.0):
    
    
    ########### Quadrilaterials ###########
    
    a1 = [0.0 ,-a-blow]
    a2 = [a+blow , 0.0]
    a3 = [0.0 , a+blow]
    a4 = [-a-blow, 0.0]
    AA = np.array([a1,a2,a3,a4])
    
    
    ########### Triangles ###########
    
    b1 = [0.0-blow, 0.0-blow]
    b2 = [b+blow, 0.0-blow]
    b3 = [0.0-blow, b+blow]
    BB = np.array([b1,b2,b3])
    
    
    ########### rotation matrices #############
    
    r01 = [np.cos(theta0), -np.sin(theta0)]
    r02 = [np.sin(theta0), np.cos(theta0) ]
    R0 = np.array([r01, r02])
        
    r11 = [np.cos(np.pi/2.+theta0), -np.sin(np.pi/2.+theta0)]
    r12 = [np.sin(np.pi/2.+theta0), np.cos(np.pi/2.+theta0) ]
    R1 = np.array([r11, r12])
        
    r21 = [np.cos(np.pi+theta0), -np.sin(np.pi+theta0)]
    r22 = [np.sin(np.pi+theta0), np.cos(np.pi+theta0) ]
    R2 = np.array([r21, r22])
        
    r31 = [np.cos(3*np.pi/2.+theta0), -np.sin(3*np.pi/2.+theta0)]
    r32 = [np.sin(3*np.pi/2.+theta0), np.cos(3*np.pi/2.+theta0) ]
    R3 = np.array([r31, r32])
        
    r41 = [np.cos(np.pi/4.+theta0), -np.sin(np.pi/4.+theta0)]
    r42 = [np.sin(np.pi/4.+theta0), np.cos(np.pi/4.+theta0) ]
    R4 = np.array([r41, r42])
    
    r51 = [np.cos(-np.pi/4.+theta0), -np.sin(-np.pi/4.+theta0)]
    r52 = [np.sin(-np.pi/4.+theta0), np.cos(-np.pi/4.+theta0) ]
    R5 = np.array([r51, r52])
    
    
    ########## Connections ##########
    
    t1 = [-a/4.,-t]
    t2 = [a/4.,-t]
    t3 = [a/4., t]
    t4 = [-a/4., t]
    T = np.array([t1,t2,t3,t4])
    
    
    # tile the plane with squares and triangles
    
    step = 4*a*np.cos(theta0)
    
    TT = [np.dot(AA, R0) + [step*i, step*j] for j in range(0, Ny) for i in range(0, Nx)]

    TTr = [np.dot(BB, R0) + [step*(i - 1/2), step*(j - 1/2)] for j in range(0, Ny + 1) for i in range(0, Nx) if (i + j) % 2 and (i, j) != (0, Ny) and (i,j) != (Nx, Ny)]
    TTr += [np.dot(BB, R2) + [step*(i - 1/2), step*(j - 1/2)] for j in range(0, Ny + 1) for i in range(1, Nx + 1) if (i + j) % 2 and (i, j) != (0, Ny) and (i, j) != (Nx, 0) and (i,j) != (Nx, Ny)]
    
    TTr += [np.dot(BB, R1) + [step*(i - 1/2), step*(j - 1/2)] for j in range(1, Ny + 1) for i in range(0, Nx + 1) if (i + j + 1) % 2 and (i, j) != (0, 0) and (i,j) != (Nx, Ny) and (i, j) != (0, Ny) and (i, j) != (Nx, 0)]
    TTr += [np.dot(BB, R3) + [step*(i - 1/2), step*(j - 1/2)] for j in range(0, Ny) for i in range(0, Nx + 1) if (i + j + 1) % 2 and (i, j) != (0, 0) and (i,j) != (Nx, Ny) and (i, j) != (0, Ny) and (i, j) != (Nx, 0)]

    
    # add horizontal struts
    
    TT2 = [np.dot(np.dot(T, R0.transpose()) + [0, 0], R0) + [step*(i - 1/2), step*(j - 1/2)] for j in range(0, Ny + 1) for i in range(0, Nx + 1) if (i + j) % 2 and (i, j) != (0, Ny) and (i, j) != (Nx, 0) and (i, j) != (Nx, Ny)]
    
    # add vertical struts  
    
    TT3 = [np.dot(np.dot(T, R1.transpose()) + [0, 0], R0) + [step*(i - 1/2), step*(j - 1/2)] for j in range(0, Ny + 1) for i in range(0, Nx + 1) if (i + j + 1) % 2 and (i, j) != (0, 0) and (i, j) != (Nx, 0) and (i, j) != (0, Ny) and (i, j) != (Nx, Ny)]

    
    # add diagonal struts  
    
    TTd = [np.dot(np.dot(T, R5.transpose()) + [0, 0], R0) + [step*i / 2 - step / 4, step*j / 2 - step / 4] for j in range(0, 2 * Ny) for i in range(0, 2 * Nx) if (i + j) % 2]
    
    
    # add antidiagonal struts  
    
    TTd += [np.dot(np.dot(T, R4.transpose()) + [0, 0], R0) + [step*i / 2 - step / 4, step*j / 2 - step / 4] for j in range(0, 2 * Ny) for i in range(0, 2 * Nx) if (i + j + 1) % 2]
    
    
    if ax is not None:
        
        #coll = PolyCollection(TT, cmap=mpl.cm.jet, edgecolors='none', color="b")
        #ax.add_collection(coll)

        #coll = PolyCollection(TTr, cmap=mpl.cm.jet, edgecolors='none', color="b")
        #ax.add_collection(coll)
        
        #coll = PolyCollection(TT2[0:1][:][:], cmap=mpl.cm.jet, edgecolors='none', color="r")
        #ax.add_collection(coll)
        
        coll = PolyCollection(TT3, cmap=mpl.cm.jet, edgecolors='none', color="g")
        ax.add_collection(coll)
        
        coll = PolyCollection(TTd, cmap=mpl.cm.jet, edgecolors='none', color="y")
        ax.add_collection(coll)
        
        ax.set_aspect("equal")
        ax.set_xlim(-30, 150)
        ax.set_ylim(-30, 160)
    
    
    return TT, TTr, TT2, TT3, TTd, step



# plot the result
fig = plt.figure(1, figsize = (10, 10))
ax = fig.add_subplot(111)
TT, TTr, TT2, TT3, TTd, step = p4gtile(ax)
n = 5
vecx = TT2[n:n+1][:][:]
vecx = vecx[0][:,0]
vecy = TT2[n:n+1][:][:]
vecy = vecy[0][:,0]
plt.plot(vecx, vecy)
#plt.plot(TT2[0][0][0],(TT2[0][0][1]+TT2[0][1][1])/2.)
plt.show()