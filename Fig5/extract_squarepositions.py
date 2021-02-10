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
import matplotlib.pyplot as plt

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
mpl.rcParams['text.usetex'] = True 
mpl.rcParams['font.sans-serif'] = 'helvet'


mpl.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{sfmath}',
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]  


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
typeD="ellipses"

for sub,c,lab in zip([6,9],["b","r"],["counterrotations","lineardecay"]):    
    F=Multimode.Focus('20190620',sub)
    F.load()
    
    fid=open("/Users/coulais/science/Shared/Metacombinatorial/Paperdraft/draft_v2/Figs/Fig3/corners_"+lab+"_"+typeD+".txt","w")
    fid.write("frame;Xbottomleft;Ybottomleft;Xbottomright;Ybottomright;Xtopright;Ytopright;Xtopleft;Ytopleft;\n")
    for t in range(F.infos.tmax-1):
        fid.write("%d;" % t)        
        for nx,ny in zip([2,15,15,2],[2,2,15,15]):
            if typeD=="squares":
                F.extract_strain(t,nx,ny)
            elif typeD=="ellipses":
                F.extract_ellipses_orderparam(t,nx,ny)
            fid.write("%f;%f;" % (F.Time.data.X,F.Time.data.Y))
        #fid.write("%d\t" % t)
        fid.write("\n")
    fid.close()
        
