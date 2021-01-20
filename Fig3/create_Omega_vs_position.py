# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 15:15:30 2015

@author: coulais
"""
###############################################################################
############################## import packages ################################
###############################################################################
import sys
from sys import platform as _platform
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.patches import Ellipse
from matplotlib.collections import PatchCollection
import cv2
import matplotlib as mpl
import pandas as pd


mpl.rcParams['axes.labelsize'] = 24
mpl.rcParams['xtick.labelsize'] =  24
mpl.rcParams['ytick.labelsize'] = 24
mpl.rcParams['figure.subplot.left'] = 0.25
mpl.rcParams['figure.subplot.right'] = 0.95
mpl.rcParams['figure.subplot.bottom'] = 0.15
mpl.rcParams['figure.subplot.top'] = 0.95
#mpl.rcParams['figure.figsize']=[ 8,8]
#mpl.rcParams['font.family']='Times'
mpl.rcParams['lines.linewidth']=5
mpl.rcParams['axes.linewidth']= 2
mpl.rcParams['xtick.major.size']= 10
mpl.rcParams['xtick.major.width']=2
mpl.rcParams['xtick.major.pad']=8
mpl.rcParams['ytick.major.pad']=8
mpl.rcParams['ytick.major.size']= 10
mpl.rcParams['ytick.major.width']=2
mpl.rcParams['xtick.minor.size']= 0#7.5
mpl.rcParams['xtick.minor.width']=0
mpl.rcParams['ytick.minor.size']= 0#7.5
mpl.rcParams['ytick.minor.width']=0
#mpl.rcParams['text.usetex'] = True 
#mpl.rcParams['font.sans-serif'] = 'helvet'


#mpl.rcParams['text.latex.preamble'] = [
#       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
#       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
#       r'\usepackage{sfmath}',
#       r'\usepackage{helvet}',    # set the normal font here
#       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
#       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
#]  


mpl.rcParams['legend.fontsize']=22;

path_sys={"darwin":"/Users/coulais/science/Shared/Metacombinatorial/",
          "win32":"TBD",
          "linux":"TBD"}
sys.path.append(path_sys[_platform])
sys.path.append(path_sys[_platform]+'PyLab')
import Multimode


###############################################################################
################## extract the position of the tip ############################
###############################################################################

fig1=plt.figure(1,figsize=(10,10))
ax1=fig1.add_axes([0.15,0.15,0.8,0.8])

for sub in [6]:#,4,5,6,9]:
    F=Multimode.Focus('20190620',sub)
    F.load()
    #F.plot_ellipses_orderparam(ax1,1,0)
    #ts=range(1,F.infos.tmax,1)
    #t=50
    Omega1=np.zeros(16)
    Omega2=np.zeros(16)
    for nx,ny in zip(range(1,17),range(1,17)):
        F.extract_ellipses_orderparam(F.infos.tmax-1,nx,ny)
        Omega1[nx-1]=F.Time.data.Omega*(-1)**(nx+ny)
    for nx,ny in zip(range(1,17),range(16,0,-1)):
        F.extract_ellipses_orderparam(F.infos.tmax-1,nx,ny)
        Omega2[nx-1]=F.Time.data.Omega*(-1)**(nx+ny)

ax1.plot(range(1,17),-Omega1,color="b",ls="--")
ax1.plot(range(1,17),-Omega2,color="b")

dic_out={}
dic_out['Omega1, x'] = np.array(range(1,17))
dic_out['Omega1, y'] = -Omega1

dic_out['Omega2, x'] = np.array(range(1,17))
dic_out['Omega2, y'] = -Omega2

DF_out = pd.DataFrame(dic_out)
DF_out.to_csv('PositionOmega_20190620_Test6.csv')

for sub in [9]:#,4,5,6,9]:
    F=Multimode.Focus('20190620',sub)
    F.load()
    #F.plot_ellipses_orderparam(ax1,1,0)
    #ts=range(1,F.infos.tmax,1)
    #t=50
    Omega1=np.zeros(16)
    Omega2=np.zeros(16)
    for nx,ny in zip(range(1,17),range(1,17)):
        F.extract_ellipses_orderparam(F.infos.tmax-1,nx,ny)
        Omega1[nx-1]=F.Time.data.Omega*(-1)**(nx+ny)
    for nx,ny in zip(range(1,17),range(16,0,-1)):
        F.extract_ellipses_orderparam(F.infos.tmax-1,nx,ny)
        Omega2[nx-1]=F.Time.data.Omega*(-1)**(nx+ny)

ax1.plot(range(1,17),-Omega1,color="r",ls="--")
ax1.plot(range(1,17),-Omega2,color="r")

dic_out={}
dic_out['Omega1, x'] = np.array(range(1,17))
dic_out['Omega1, y'] = -Omega1

dic_out['Omega2, x'] = np.array(range(1,17))
dic_out['Omega2, y'] = -Omega2

DF_out = pd.DataFrame(dic_out)
DF_out.to_csv('PositionOmega_20190620_Test9.csv')



ax1.set_xlim(0,17)
ax1.set_xticks(range(1,17))
ax1.set_ylim(-1,1)
ax1.set_yticks([-0.5,0,0.5,1])

ax1.set_xlabel(r"x position")
ax1.set_ylabel(r"polarisation $P$")
ax1.yaxis.set_label_coords(-0.1,0.5)

#fig1.savefig("/Users/coulais/science/Shared/Metacombinatorial/Paperdraft/draft_v2/Figs/Fig3/Omega_vs_position_2.pdf")