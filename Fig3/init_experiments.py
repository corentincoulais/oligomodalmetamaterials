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

'''
for sub in [4,5,6,9]:#,2,3,4,5,6]:#range(10):#,1,3,4]:#range(4,5):#[15]:#,25,26]:
    F=Multimode.Focus('20190620',sub)
    F.detect_islands(display=False,save=True,displaypic=1)
    F.track_islands()


'''
#F.time(0)
#F.Time.fields(["X","Y"])
#plt.scatter(F.Time.data.X,F.Time.data.Y)
#fndflk

'''
for sub in [4,5,6,9]:#range(10):#,1,3,4]:#range(4,5):#[15]:#,25,26]:
    S=Multimode.CreateDB('20190620',sub)
    S.load(dat="infos")
    S.load(dat="imagedata")
    S.makeDB()
#    #F.forcedisp()
 #   #plt.figure(3)

#Assign labels to particles


fig1=plt.figure(1,figsize=(8,8))
ax1=fig1.add_axes([0.1,0.1,0.8,0.8])


sub=9
D=Multimode.DBReshape('20190620',sub)
D.load(dat="infos")
D.time(0)
D.Time.fields(['t','Y','X','p'])
ax1.scatter(D.Time.data.X,D.Time.data.Y)
#D.assign_nxny_auto(ax1,offx=0.4,offy=0.4)
#D.assign_nxny_auto(ax1,offx=0.3,offy=0.3)
D.assign_nxny_auto(ax1,offx=-0.3,offy=-0.5)
'''

fig1=plt.figure(1,figsize=(8,8))
ax1=fig1.add_axes([0.1,0.1,0.8,0.8])
#Define squares based on particles and initial positions
sub=9
D=Multimode.CreateDB('20190620',sub)
D.load(dat="infos")
#D.Time.fields(['t','Y','X','p'])
#ax1.scatter(D.Time.data.X,D.Time.data.Y,c="b")
D.writetxtSQUARE(ax1=ax1)
D.writetxtELLIPSE(ax1=ax1)
