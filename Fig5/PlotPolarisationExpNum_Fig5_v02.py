# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 15:15:30 2015

@author: dykstra


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
import pandas as pd


def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)

titles=['SlowRS','FastBD']
CSVFiles_Exp = [r"C:\\Users\david\Documents\0_PhD\15_Viscoelastic_Metamaterial\Python\AveragePolarisationExperimental\EpsOmega4.csv",
                r"C:\\Users\david\Documents\0_PhD\15_Viscoelastic_Metamaterial\Python\AveragePolarisationExperimental\EpsOmega5.csv"]


CSVFiles_Num=[r"G:\\15_Viscoelastic_Metamaterial\AbaqusNew\Model_v013\batch_v13_42_1deg_Visco-Agilus_improved_norelax19mm\OutputOmega\Omega-minmaxavg_MAT-18.csv",
              r"G:\\15_Viscoelastic_Metamaterial\AbaqusNew\Model_v013\batch_v13_42_1deg_Visco-Agilus_improved_norelax19mm\OutputOmega\Omega-minmaxavg_MAT-1.csv"]


CSVFiles_Num=[r"G:\\15_Viscoelastic_Metamaterial\AbaqusNew\Model_v013\batch_v13_42_1deg_Visco-Agilus_improved_norelax19mm\OutputOmega\Omega-minmaxavg_MAT-34.csv",
              r"G:\\15_Viscoelastic_Metamaterial\AbaqusNew\Model_v013\batch_v13_42_1deg_Visco-Agilus_improved_norelax19mm\OutputOmega\Omega-minmaxavg_MAT-33.csv"]
N=2

Data_Exp={}
Data_Num={}

for indr in range(N):
    Data_Exp[indr]=pd.read_csv(CSVFiles_Exp[indr],sep=',')
    Data_Num[indr]=pd.read_csv(CSVFiles_Num[indr],sep=',')




#%%

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
mpl.rcParams['xtick.minor.size']= 5#7.5
mpl.rcParams['xtick.minor.width']=3
mpl.rcParams['ytick.minor.size']= 4#7.5
mpl.rcParams['ytick.minor.width']=3      
    






#%% Plot color for Tests of 20190913

plt.close('all')


figsz=(4,4)
#axshape = [0.25,0.3,0.65,0.65]
axshape = [0.25,0.25,0.7,0.7]

lw=4
colors=['dodgerblue','darkorange']

for indr in range(N):
    

    fig1=plt.figure(100+indr,figsize=figsz)
    ax1=fig1.add_axes(axshape)  
    
    ax1.fill_between(-Data_Exp[indr]['Epsilon'],Data_Exp[indr]['Omega_min'],Data_Exp[indr]['Omega_max'],color=colors[0],alpha=0.3)    
    ax1.fill_between(Data_Num[indr]['Epsilon'],Data_Num[indr]['Omega_min'],Data_Num[indr]['Omega_max'],color=colors[1],alpha=0.3)
    
    ax1.plot(-Data_Exp[indr]['Epsilon'],Data_Exp[indr]['Omega_avg'],color=colors[0],Linewidth=lw)   
    ax1.plot(Data_Num[indr]['Epsilon'],Data_Num[indr]['Omega_avg'],color=colors[1],Linewidth=lw)   
    
    
    ax1.set_xlabel(r'$\varepsilon$')
    ax1.set_ylabel(r'$\Omega$')#,rotation=0)
    ax1.set_xlim([0,0.05])
    ax1.set_ylim([-0.5,0.5])
    ax1.set_xticks([0.,0.025,0.05])
    ax1.set_yticks([-0.5,0.,0.5])
    ax1.set_xticklabels(['0','0.025','0.5'])
    ax1.set_yticklabels(['-0.5','0','0.5'])
    
    #ax1.xaxis.set_label_coords(0.5, -0.2)
    #ax1.yaxis.set_label_coords(-0.23, 0.5)
    ax1.xaxis.set_label_coords(0.5, -0.17)
    ax1.yaxis.set_label_coords(-0.13, 0.5)
     
    #ax1.xaxis.set_label_coords(0.5, -0.16)
    #ax1.yaxis.set_label_coords(0., 1.04)

    
    FigName = 'Fig5_EpsilonOmegaExpNum_'+titles[indr]+'_v02'
    fig1.savefig(FigName+'.png',dpi=600,transparent=True)
    fig1.savefig(FigName+'.pdf')

    
