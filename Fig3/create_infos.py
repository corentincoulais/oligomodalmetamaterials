# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 18:56:45 2015

@author: coulais
"""

#import matplotlib as plt
#import cv2
import csv
import pickle
from scipy import ndimage
import os
import matplotlib as mpl
import numpy as np
import sys
#reload(sys)
#sys.setdefaultencoding("UTF8")
import os
from scipy.optimize import curve_fit
import shutil
import numpy as np
from sys import platform as _platform
#Region of interes for calculating angles
path_sys={"darwin":"/Users/coulais/science/Shared/Metacombinatorial",
          "win32":"TBD",
          "linux":"/home/aleksi/Desktop/Metacombinatorial"}
          
path_rawdata={"darwin":"/Volumes/LaCie/Aleksi",
          "win32":"TBD",
          "linux":"/home/aleksi/Desktop/Metacombinatorial"}

sys.path.append(path_sys[_platform])
sys.path.append(path_sys[_platform]+'PyLab')
import Multimode

mpl.rcParams['savefig.dpi'] = 50

suffs=[]
paths=[]
paths_rawdata=[]
pathOs=[]
instron_file=[]
network=[]
BBox=[]
ths=[]
minArea=[]
maxArea=[]
instronpush="TBD"
good=[]

squaresize=12 #./np.sqrt(2)#measured from the the picture at rest.

suffs+=["1"];instron_file+=["20190620.is_ccyclic_Exports/20190620_1_1.csv"];network+=[0];BBox+=[[580,2748,0,2748]];ths+=[150];minArea+=[80];maxArea+=[300];good+=["no"]
suffs+=["2"];instron_file+=["20190620.is_ccyclic_Exports/20190620_2_1.csv"];network+=[0];BBox+=[[570,2748,0,2748]];ths+=[150];minArea+=[80];maxArea+=[300];good+=["no"]
suffs+=["3"];instron_file+=["20190620.is_ccyclic_Exports/20190620_3_1.csv"];network+=[0];BBox+=[[580,2748,0,2748]];ths+=[150];minArea+=[80];maxArea+=[300];good+=["yes"]
suffs+=["4"];instron_file+=["20190620.is_ccyclic_Exports/20190620_4_1.csv"];network+=[0];BBox+=[[590,2748,0,2748]];ths+=[150];minArea+=[80];maxArea+=[300];good+=["yes"]
suffs+=["5"];instron_file+=["20190620.is_ccyclic_Exports/20190620_5_1.csv"];network+=[0];BBox+=[[610,2748,0,2748]];ths+=[150];minArea+=[80];maxArea+=[300];good+=["yes"]
suffs+=["6"];instron_file+=["20190620.is_ccyclic_Exports/20190620_6_1.csv"];network+=[0];BBox+=[[524,2748,0,2748]];ths+=[150];minArea+=[80];maxArea+=[300];good+=["yes"]
suffs+=["7"];instron_file+=["20190620.is_ccyclic_Exports/20190620_7_1.csv"];network+=[0];BBox+=[[570,2748,0,2748]];ths+=[150];minArea+=[80];maxArea+=[300];good+=["no"]
suffs+=["8"];instron_file+=["20190620.is_ccyclic_Exports/20190620_8_1.csv"];network+=[0];BBox+=[[570,2748,0,2748]];ths+=[150];minArea+=[80];maxArea+=[300];good+=["no"]
suffs+=["9"];instron_file+=["20190620.is_ccyclic_Exports/20190620_9_1.csv"];network+=[0];BBox+=[[520,2748,0,2748]];ths+=[150];minArea+=[80];maxArea+=[300];good+=["yes"]

N=len(suffs)
paths+=[path_sys[_platform]+'/Experiments/20190620/']*N
paths_rawdata+=[path_rawdata[_platform]+'/20190620/']*N
if os.path.isdir("Files") is False:
    os.mkdir("Files")
for i in range(N):
    #os.chdir("//storage01//data//AMOLF//groups//hecke-group//Corentin//rawdata//")
    F=Multimode.Focus('20190620',i+1)
    print(F.path)
    F.infos.path_rawdata=paths_rawdata[i]
    F.infos.suffix=suffs[i]
    F.infos.instron_file=instron_file[i]
    F.infos.network=network[i]
    F.infos.path=paths[i]
    F.infos.BBox=BBox[i]
    F.infos.ths=ths[i]
    F.infos.square=squaresize
    F.infos.tracking_maxdisp=20
    F.infos.pix2mm=16.*12./2634.#
    F.infos.perim_lim=20
    F.infos.minArea=minArea[i]
    F.infos.maxArea=maxArea[i]
    F.infos.tmax=360
    #F.infos.maxaspect=3
    F.infos.N=17
    F.infos.size=squaresize*F.infos.N
    #F.infos.ROI=ROI()
    F.save(dat='infos')
    
    
