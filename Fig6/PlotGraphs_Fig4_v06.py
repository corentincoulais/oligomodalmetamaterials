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
import pandas as pd


def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)


CurrentFolder = os.getcwd()
MainFolder =  os.path.dirname(os.path.dirname(CurrentFolder))
MainFolder = r"G:\15_Viscoelastic_Metamaterial"
ImageAnalysisMainFolder = MainFolder+'\\Multimode2_v02\\Experiments\\20190913'

Test_OmegaFolder=2*[None]
Test_OmegaFolder[0] = ImageAnalysisMainFolder+'\\OmegaDiag_average\\4\\CSVfiles'
Test_OmegaFolder[1] = ImageAnalysisMainFolder+'\\OmegaDiag_average\\5\\CSVfiles'

TestOmegaCSV =2*[None] 
TestOmegaCSV[0] = Test_OmegaFolder[0]+'\\frame_00301.csv'
TestOmegaCSV[1] = Test_OmegaFolder[1]+'\\frame_00049.csv'

Test_PoissonFolder=2*[None]
Test_PoissonFolder[0] = ImageAnalysisMainFolder+'\\Poisson\\4\\CSVfiles'
Test_PoissonFolder[1] = ImageAnalysisMainFolder+'\\Poisson\\5\\CSVfiles'

TestPoissonCSV =2*[None] 
TestPoissonCSV[0] = Test_PoissonFolder[0]+'\\RelativePoisson_avgside.csv'
TestPoissonCSV[1] = Test_PoissonFolder[1]+'\\RelativePoisson_avgside.csv'

Final_Data_Folder = MainFolder


#%%

Numerics_PoissonFolder = r"G:\15_Viscoelastic_Metamaterial\AbaqusNew\Model_v013\batch_v13_42_1deg_Visco-Agilus_improved_norelax19mm\OutputOmega"
NumericsPoissonCSV =2*[None] 
NumericsPoissonCSV[0] = Numerics_PoissonFolder+'\\EpsilonPoissonAverage_MAT-34.csv'
NumericsPoissonCSV[1] = Numerics_PoissonFolder+'\\EpsilonPoissonAverage_MAT-33.csv'


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

figsz=(4.,4.)
#axshape = [0.29,0.3,0.65,0.65]
axshape = [0.25,0.25,0.7,0.7]


lw=4
colors=['dodgerblue','darkorange']

Poisson_DF = {}
for ind in range(len(TestPoissonCSV)):
    Poisson_DF[ind] = pd.read_csv(TestPoissonCSV[ind])    


PoissonNumerics_DF = {}
for ind in range(len(NumericsPoissonCSV)):
    PoissonNumerics_DF[ind] = pd.read_csv(NumericsPoissonCSV[ind])    



fig1=plt.figure(1,figsize=figsz)
ax1=fig1.add_axes(axshape)  

ax1.plot(Poisson_DF[0]['- epsilon_y RelativePoisson_avgside'],Poisson_DF[0]['nu RelativePoisson_avgside'],'--',color=colors[0],marker='o',markersize=8,markeredgewidth=2.,Linewidth=lw, markerfacecolor='white')   
ax1.plot(Poisson_DF[1]['- epsilon_y RelativePoisson_avgside'],Poisson_DF[1]['nu RelativePoisson_avgside'],'--',color=colors[1],marker='s',markersize=8,markeredgewidth=2.,Linewidth=lw, markerfacecolor='white')


ax1.set_xlabel(r'$\varepsilon$')
ax1.set_ylabel(r'$\nu$')
ax1.set_xlim([0,0.05])
ax1.set_ylim([-0.5,0.5])
ax1.set_xticks([0.,0.025,0.05])
ax1.set_yticks([-0.5,0,0.5])
ax1.set_xticklabels(['0','0.025','0.5'])
ax1.set_yticklabels(['-0.5','0','0.5'])

#ax1.xaxis.set_label_coords(0.5, -0.2)
#ax1.yaxis.set_label_coords(-0.27, 0.5)
ax1.xaxis.set_label_coords(0.5, -0.17)
ax1.yaxis.set_label_coords(-0.13, 0.5)


FigName = '//Fig4_AveragePoisson'
fig1.savefig(CurrentFolder+FigName+'.png',dpi=600,transparent=True)
fig1.savefig(CurrentFolder+FigName+'.pdf')



#%%

fig2=plt.figure(2,figsize=figsz)
ax2=fig2.add_axes(axshape)  

Omega_DF = {}
for ind in range(len(TestOmegaCSV)):
    Omega_DF[ind] = pd.read_csv(TestOmegaCSV[ind])
    
    
ax2.plot(Omega_DF[0]['Position x']+1,Omega_DF[0]['Omega x, centre'],'--',color=colors[0],Linewidth=lw)#,marker='o',markersize=20, markerfacecolor='none')        
ax2.plot(Omega_DF[0]['Position y'],Omega_DF[0]['Omega y, centre'],'-',color=colors[0],Linewidth=lw)#,marker='o',markersize=20, markerfacecolor='none')        


ax2.plot(Omega_DF[1]['Position x']+1,Omega_DF[1]['Omega x, centre'],'--',color=colors[1],Linewidth=lw)#,marker='s',markersize=20, markerfacecolor='none')        
ax2.plot(Omega_DF[1]['Position y'],Omega_DF[1]['Omega y, centre'],'-',color=colors[1],Linewidth=lw)#,marker='s',markersize=20, markerfacecolor='none')        

    

ax2.set_xlabel('position')
ax2.set_ylabel('$\Omega$')

ax2.set_xlim([0,10])
ax2.set_ylim([-0.6,0.6])
ax2.set_yticks([-0.6,-0.3,0,0.3,0.6])

ax2.xaxis.set_label_coords(0.5, -0.2)
ax2.yaxis.set_label_coords(-0.27, 0.5)

#ax2.minorticks_on()
#ax2.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())


FigName = '//Fig4_Omega'
fig2.savefig(CurrentFolder+FigName+'.png',dpi=600,transparent=True)
fig2.savefig(CurrentFolder+FigName+'.pdf')


#%% Seperate Poisson Ratio plots


#axshape = [0.25,0.3,0.65,0.65]
#axshape = [0.20,0.19,0.7,0.7]

lw=4
colors=['dodgerblue','darkorange']
colors=['dodgerblue','dodgerblue']
colorsnum='darkorange'

Poisson_DF = {}
for ind in range(len(TestPoissonCSV)):
    Poisson_DF[ind] = pd.read_csv(TestPoissonCSV[ind])    



fig1=plt.figure(10,figsize=figsz)
ax1=fig1.add_axes(axshape)  

ax1.fill_between(PoissonNumerics_DF[0]['Epsilon'],PoissonNumerics_DF[0]['PoissonMin'],PoissonNumerics_DF[0]['PoissonMax'],color=colorsnum,alpha=0.3)    
ax1.plot(PoissonNumerics_DF[0]['Epsilon'],PoissonNumerics_DF[0]['PoissonAverage'],'-',color=colorsnum,Linewidth=lw)   

ax1.errorbar(Poisson_DF[0]['- epsilon_y RelativePoisson_avgside'],Poisson_DF[0]['nu RelativePoisson_avgside'],yerr=0.1,elinewidth=2,capsize=6,fmt='--',color=colors[0],marker='o',markersize=8,markeredgewidth=2.,Linewidth=lw, markerfacecolor='white')   
#ax1.plot(Poisson_DF[1]['- epsilon_y RelativePoisson_avgside'],Poisson_DF[1]['nu RelativePoisson_avgside'],'--',color=colors[1],marker='s',markersize=8,markeredgewidth=2.,Linewidth=lw, markerfacecolor='white')



ax1.set_xlabel(r'$\varepsilon$')
ax1.set_ylabel(r'$\nu$')#,rotation=0)
ax1.set_xlim([0,0.05])
ax1.set_ylim([-0.8,0.8])
ax1.set_xticks([0.,0.025,0.05])
ax1.set_yticks([-0.8,0,0.8])
ax1.set_xticklabels(['0','0.025','0.5'])
ax1.set_yticklabels(['-0.8','0','0.8'])

#ax1.xaxis.set_label_coords(0.5, -0.2)
#ax1.yaxis.set_label_coords(-0.27, 0.5)
ax1.xaxis.set_label_coords(0.5, -0.17)
ax1.yaxis.set_label_coords(-0.13, 0.5)

FigName = '//Fig4_AveragePoisson_CR_v06'
fig1.savefig(CurrentFolder+FigName+'.png',dpi=600,transparent=True)
fig1.savefig(CurrentFolder+FigName+'.pdf')

#%%# =============================================================================
# ax1.set_xlim([0,0.1])
# ax1.set_ylim([-2,2])
# ax1.set_xticks([0.,0.05,0.1])
# ax1.set_yticks([-2,-1,0,1,2])
# 
# 
# FigName = '//Fig4_AveragePoisson_CR_largedomain'
# fig1.savefig(CurrentFolder+FigName+'.png',dpi=600,transparent=True)
# fig1.savefig(CurrentFolder+FigName+'.pdf')
# =============================================================================




lw=4

Poisson_DF = {}
for ind in range(len(TestPoissonCSV)):
    Poisson_DF[ind] = pd.read_csv(TestPoissonCSV[ind])    



fig1=plt.figure(11,figsize=figsz)
ax1=fig1.add_axes(axshape)  


ax1.fill_between(PoissonNumerics_DF[1]['Epsilon'],PoissonNumerics_DF[1]['PoissonMin'],PoissonNumerics_DF[1]['PoissonMax'],color=colorsnum,alpha=0.3)    
ax1.plot(PoissonNumerics_DF[1]['Epsilon'],PoissonNumerics_DF[1]['PoissonAverage'],'-',color=colorsnum,Linewidth=lw)   

#ax1.plot(Poisson_DF[0]['- epsilon_y RelativePoisson_avgside'],Poisson_DF[0]['nu RelativePoisson_avgside'],'--',color=colors[0],marker='o',markersize=8,markeredgewidth=2.,Linewidth=lw, markerfacecolor='white')   
#ax1.plot(Poisson_DF[1]['- epsilon_y RelativePoisson_avgside'],Poisson_DF[1]['nu RelativePoisson_avgside'],'--',color=colors[1],marker='o',markersize=8,markeredgewidth=2.,Linewidth=lw, markerfacecolor='white')
ax1.errorbar(Poisson_DF[1]['- epsilon_y RelativePoisson_avgside'],Poisson_DF[1]['nu RelativePoisson_avgside'],yerr=0.1,elinewidth=2,capsize=6,fmt='--',color=colors[0],marker='o',markersize=8,markeredgewidth=2.,Linewidth=lw, markerfacecolor='white')   


ax1.set_xlabel(r'$\varepsilon$')
ax1.set_ylabel(r'$\nu$')#,rotation=0)
ax1.set_xlim([0,0.05])
ax1.set_ylim([-0.8,0.8])
ax1.set_xticks([0.,0.025,0.05])
ax1.set_yticks([-0.8,0,0.8])
ax1.set_xticklabels(['0','0.025','0.5'])
ax1.set_yticklabels(['-0.8','0','0.8'])

#ax1.xaxis.set_label_coords(0.5, -0.2)
#ax1.yaxis.set_label_coords(-0.27, 0.5)
ax1.xaxis.set_label_coords(0.5, -0.17)
ax1.yaxis.set_label_coords(-0.13, 0.5)


FigName = '//Fig4_AveragePoisson_LD_v06'
fig1.savefig(CurrentFolder+FigName+'.png',dpi=600,transparent=True)
fig1.savefig(CurrentFolder+FigName+'.pdf')

# =============================================================================
# ax1.set_xlim([0,0.1])
# ax1.set_ylim([-2,2])
# ax1.set_xticks([0.,0.05,0.1])
# ax1.set_yticks([-2,-1,0,1,2])
# 
# FigName = '//Fig4_AveragePoisson_LD_largedomain'
# fig1.savefig(CurrentFolder+FigName+'.png',dpi=600,transparent=True)
# fig1.savefig(CurrentFolder+FigName+'.pdf')
# =============================================================================

        