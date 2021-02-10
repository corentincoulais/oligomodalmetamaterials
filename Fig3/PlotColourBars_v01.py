# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 15:15:30 2015

@author: coulais, dykstra

# On newer versions of MySQL, after starting server,
# type "SET GLOBAL local_infile = true;" in MySQL Command

"""
###############################################################################
############################## import packages ################################
###############################################################################
import sys
from sys import platform as _platform
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.patches import Ellipse,Rectangle
from matplotlib.collections import PatchCollection
import matplotlib as mpl
import time 
import cv2
import copy
path_sys={"darwin":"/Users/coulais/science/David_Corentin/Multimode",
          "win32":r"C:/Users/David/Documents/0_PhD/15_Viscoelastic_Metamaterial/Multimode2_v02/",    #@@@ CHANGE MADE   26-07-2019
          "linux":"/home/aleksi/Desktop/Metacombinatorial"}

sys.path.append(path_sys[_platform])
sys.path.append(path_sys[_platform]+'PyLab')
import Multimode2


def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)


plt.close('all')


mpl.rcParams['axes.labelsize'] = 24
mpl.rcParams['xtick.labelsize'] =  24
mpl.rcParams['ytick.labelsize'] = 24
mpl.rcParams['figure.subplot.left'] = 0.25
mpl.rcParams['figure.subplot.right'] = 0.95
mpl.rcParams['figure.subplot.bottom'] = 0.15
mpl.rcParams['figure.subplot.top'] = 0.95
#mpl.rcParams['figure.figsize']=[ 8,8]
mpl.rcParams['font.family']='Calibri'
mpl.rcParams['lines.linewidth']=10
mpl.rcParams['axes.linewidth']= 3
mpl.rcParams['xtick.major.size']= 8
mpl.rcParams['xtick.major.width']=3
mpl.rcParams['xtick.major.pad']=3
mpl.rcParams['ytick.major.pad']=3
mpl.rcParams['ytick.major.size']= 8
mpl.rcParams['ytick.major.width']=3
mpl.rcParams['xtick.minor.size']= 0#7.5
mpl.rcParams['xtick.minor.width']=0
mpl.rcParams['ytick.minor.size']= 0#7.5
mpl.rcParams['ytick.minor.width']=0      
    






#%% Plot for Tests of 20190913


plt.close('all')

Omlims=(-1,1)
deltaOmega = Omlims[1]-Omlims[0]
N = 1001
Omega = np.linspace(Omlims[0],Omlims[1],N)
dOmega = Omega[1]-Omega[0]
icol = (Omega+1)/2
icol=1-icol
dicol = icol[1]-icol[0]
AR_cb = 0.1

ialpha=np.abs(Omega)*5
ialpha[np.where(ialpha>1)]=1

w_cb = AR_cb*deltaOmega


fig1=plt.figure(1,figsize=(3.4*0.25,3.4))
ax1=fig1.add_axes([0.3,0.1,0.8,0.8])
Rect_back = Rectangle((0,Omlims[0]),w_cb,deltaOmega)
Rect_back.set_facecolor('k')#plt.cm.PuOr(icol))
Rect_back.set_edgecolor(None)
ax1.add_artist(Rect_back)



for ind in range(N-1):
    Rect = Rectangle((0,Omlims[0]+ind*dOmega),w_cb,dOmega)
    Rect.set_facecolor(plt.cm.bwr(icol[ind]))
    Rect.set_edgecolor(None)
    Rect.set_alpha(ialpha[ind])
    
    ax1.add_artist(Rect)

ax1.axis('scaled')

ax1.set_xlim([0,w_cb])
ax1.set_ylim(Omlims)
ax1.set_xticks([])
ax1.set_yticks([Omlims[0],0,Omlims[1]])


ax1.set_xticks([])

ax1.set_xlabel('$\Omega$')

fig1.savefig('Colorbar_20190913.pdf')
fig1.savefig('Colorbar_20190913.png',dpi=600,transparent=True)


#%% Plot for Tests of 20190620




Omlims=(-1,1)
deltaOmega = Omlims[1]-Omlims[0]
N = 1001
Omega = np.linspace(Omlims[0],Omlims[1],N)
dOmega = Omega[1]-Omega[0]
icol = (Omega+1)/2
icol=1-icol
dicol = icol[1]-icol[0]
AR_cb = 0.1

ialpha=np.abs(Omega)*2
ialpha[np.where(ialpha>1)]=1

w_cb = AR_cb*deltaOmega


fig2=plt.figure(2,figsize=(45/25.4*2*0.25,45/25.4*2))
ax2=fig2.add_axes([0.3,0.1,0.8,0.8])
Rect_back = Rectangle((0,Omlims[0]),w_cb,deltaOmega)
Rect_back.set_facecolor('k')#plt.cm.PuOr(icol))
Rect_back.set_edgecolor(None)
ax2.add_artist(Rect_back)



for ind in range(N-1):
    Rect = Rectangle((0,Omlims[0]+ind*dOmega),w_cb,dOmega)
    Rect.set_facecolor(plt.cm.bwr(icol[ind]))
    Rect.set_edgecolor(None)
    Rect.set_alpha(ialpha[ind])
    
    ax2.add_artist(Rect)

ax2.axis('scaled')

ax2.set_xlim([0,w_cb])
ax2.set_ylim(Omlims)
ax2.set_xticks([])
ax2.set_yticks([Omlims[0],0,Omlims[1]])


ax2.set_xticks([])

ax2.set_xlabel('$\Omega$')

fig2.savefig('Colorbar_20190620.pdf')
fig2.savefig('Colorbar_20190620.png',dpi=600,transparent=True)
