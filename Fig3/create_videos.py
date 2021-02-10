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
path_sys={"darwin":"/Users/coulais/science/Shared/Metacombinatorial/",
          "win32":"TBD",
          "linux":"TBD"}
sys.path.append(path_sys[_platform])
sys.path.append(path_sys[_platform]+'PyLab')
import Multimode


###############################################################################
################## extract the position of the tip ############################
###############################################################################

for sub in [3,4,5,6,9]:
    F=Multimode.Focus('20190620',sub)
    F.load()
    #F.plot_ellipses_orderparam(ax1,1,0)
    ts=range(1,F.infos.tmax,1)
    for it in range(len(ts)):
        fig1=plt.figure(1,figsize=(8,8))
        ax1=fig1.add_axes([0.,0.,1,1])
        ax1.set_axis_off()
        F.plot_ellipses_orderparam(ax1,ts[it],it)
        fig1.savefig(F.movie.file+"frame_%05d.png" % (ts[it]))
        plt.close()
    
    os.system("/Applications/ffmpeg -framerate 30 -i "+ F.movie.file+"frame_%05d.png -s '800x778' -c:v libx264 -pix_fmt yuv420p "+F.movie.file+"../video" +str(sub)+".mp4")

