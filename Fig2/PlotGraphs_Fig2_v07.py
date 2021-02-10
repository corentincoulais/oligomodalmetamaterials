# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 15:15:30 2015

@author: coulais, dykstra


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
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)




def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)


CurrentFolder = os.getcwd()
MainFolder =  os.path.dirname(os.path.dirname(CurrentFolder))

#DataAleksiFolder = MainFolder+'\\Python\\PlotGraphFig2\\DataAleksi\\fig2data\\'
DataAleksiFolder = MainFolder+'\\Python\\PlotGraphFig2\\DataAleksi\\fig2data_v02\\scalingdata\\'

#ScalingFileNames = ['0110.out', '1111.out', '3113.out','quasi.out']
#ScalingFileNames = ['3113.out',         '0110.out', '0310.out','1010.out','1112.out','2310.out', '3210.out'                 ,'1311.out','1111.out']

ScalingFileNames = ['pmg.out', 'random.out','pm.out','p2.out','p1.out', 'pmg2.out'                 ,'cm1.out','cm2.out' ,'p4m.out',         'cm3.out', 'p4g.out']
markers = ['x' ,'2','v','<','>','x','s','o','+','^','D' ]#['x','^','D','v','<','>','x','s','o','+','o']


ScalingCorrespondence = ['Linear','Constant','1','quasi']

N = 11#%len(ScalingFileNames)
ScalingFiles=N*[None]
for ind in range(N):
    ScalingFiles[ind]=DataAleksiFolder+ScalingFileNames[ind]

Final_Data_Folder = MainFolder


 

#%%

plt.close('all')


mpl.rcParams['axes.labelsize'] = 24
mpl.rcParams['xtick.labelsize'] =  24
mpl.rcParams['ytick.labelsize'] = 24
mpl.rcParams['figure.subplot.left'] = 0.28
mpl.rcParams['figure.subplot.right'] = 0.94
mpl.rcParams['figure.subplot.bottom'] = 0.18
mpl.rcParams['figure.subplot.top'] = 0.94
#mpl.rcParams['figure.figsize']=[ 8,8]
mpl.rcParams['font.family']='Calibri'
mpl.rcParams['lines.linewidth']=5
mpl.rcParams['axes.linewidth']= 2
mpl.rcParams['xtick.major.size']= 10
mpl.rcParams['xtick.major.width']=2
mpl.rcParams['xtick.major.pad']=8
mpl.rcParams['ytick.major.pad']=8
mpl.rcParams['ytick.major.size']= 10
mpl.rcParams['ytick.major.width']=2
mpl.rcParams['xtick.minor.size']= 5#7.5
mpl.rcParams['xtick.minor.width']=2
mpl.rcParams['ytick.minor.size']= 5#7.5
mpl.rcParams['ytick.minor.width']=2
#mpl.rcParams['text.usetex'] = True 
#mpl.rcParams['font.sans-serif'] = 'helvet'






#%% Plot Fig2 Main



plt.close('all')



col_mono = (1,0.5,0.5)
col_oligo = (0.5,1,0.5)
col_pluri = (0.5,0.5,1)


colorvec=2*[col_mono]+7*[col_pluri]+2*[col_oligo]#[col_mono,    'green','green','green','green','green','green',    'blue','blue']
lw_vec = [1.5,1.5, 1.5, 1.5,0., 1.5,0.,0., 1.5, 5., 5.]  #2*[1.5]+7*[1.5]+2*[4.]




plot_vec=1*[1]+1*[2]+9*[1]#[1,        1,0,1,1,1,0,        1,1]
mksize_vec = 2*[7.]+7*[7.]+2*[7.]
#colorvec=['green','red','orange','blue','cyan']
#markers = ['v','s','^','o','+']


figsz=(4,4)
axshape = [0.22,0.21,0.73,0.73]


Scaling_DF = {}
N_mat=N*[None]
N_mat2=N*[None]

inds1_mat = N*[None]
for ind in range(N):
    Scaling_DF[ind] = pd.read_csv(ScalingFiles[ind],header=None)    

    N_mat[ind] = range(2,33,2)
    N_mat[ind]=np.array(N_mat[ind])
    inds1_mat[ind] = N_mat[ind]-1
    N_mat2[ind] = np.array(range(len(N_mat[ind])))
    

    
markers = ['x' ,'2','v','<','>','*','s','o','+','^','D' ]#['x','^','D','v','<','>','x','s','o','+','o']

N_mat2[0] = [0,1,2,3,5,7,9,11,13,15]
N_mat2[1] = [0,1,2,4,6,8,10,12,14]

N_mat2[3] = [0,1,3,5,7,9,11,13,15]
N_mat2[4] = [0,2,4,6,8,10,12,14]

N_mat2[5] = [0,3,6,9,12,15,12,14]
N_mat2[6] = [0,1,4,7,10,13]
N_mat2[7] = [0,2,5,8,11,14]


fig1=plt.figure(1,figsize=figsz)
ax1=fig1.add_axes(axshape)  

for ind in range(N):
    if plot_vec[ind]==1:
        xdata = N_mat[ind][N_mat2[ind]]
        ydata = np.array(Scaling_DF[ind][0][inds1_mat[ind]])[N_mat2[ind]]
        ax1.plot(xdata,ydata,color=colorvec[ind],linewidth=lw_vec[ind],linestyle='-',marker=markers[ind],markersize=mksize_vec[ind],markeredgewidth=1.5, markerfacecolor='white')#colorvec[ind])        

    elif plot_vec[ind]==2:
        xdata = N_mat[ind][N_mat2[ind]]
        ydata = np.array(Scaling_DF[ind][0][inds1_mat[ind]])[N_mat2[ind]]
        ydata[np.isnan(ydata)]=1.
        ax1.plot(xdata,ydata,color=colorvec[ind],linewidth=lw_vec[ind],linestyle='-',marker=markers[ind],markersize=mksize_vec[ind],markeredgewidth=1.5, markerfacecolor='white')#colorvec[ind])        


ax1.set_xlabel('$n$')
ax1.set_ylabel('$N_M$')
ax1.xaxis.set_label_coords(0.5, -0.18)
ax1.yaxis.set_label_coords(-0.18, 0.5)


ax1.tick_params(which='both',direction='out')
ax1.set_xlim([0,32])
ax1.set_ylim([0,16])
ax1.set_xticks([0,8,16,24,32])
ax1.set_yticks([0,4,8,12,16])
ax1.minorticks_on()
#ax1.tick_params(bottom=True,top=True,left=True,right=True)

#ax1.grid(which='both')



FigName = '//Fig2_MainGraph'
fig1.savefig(CurrentFolder+FigName+'.png',dpi=600,transparent=True)
fig1.savefig(CurrentFolder+FigName+'.pdf')






#%% Plot Fig2 Appendix



plt.close('all')



col_mono = (1,0.5,0.5)
col_oligo = (0.5,1,0.5)
col_pluri = (0.5,0.5,1)


colorvec=2*[col_mono]+7*[col_pluri]+2*[col_oligo]#[col_mono,    'green','green','green','green','green','green',    'blue','blue']
lw_vec = [1.,1., 1., 1.,0., 1.,0.,0., 1., 1., 1.]  #2*[1.5]+7*[1.5]+2*[4.]




plot_vec=1*[1]+1*[2]+9*[1]#[1,        1,0,1,1,1,0,        1,1]
mksize_vec = 2*[3.]+7*[3.]+2*[3.]
#colorvec=['green','red','orange','blue','cyan']
#markers = ['v','s','^','o','+']


figsz=(4,4)
axshape = [0.22,0.21,0.73,0.73]


Scaling_DF = {}
N_mat=N*[None]
N_mat2=N*[None]

inds1_mat = N*[None]
for ind in range(N):
    Scaling_DF[ind] = pd.read_csv(ScalingFiles[ind],header=None)    

    N_mat[ind] = range(2,33,2)
    N_mat[ind]=np.array(N_mat[ind])
    inds1_mat[ind] = N_mat[ind]-1
    N_mat2[ind] = np.array(range(len(N_mat[ind])))
    

    
markers = ['x' ,'2','v','<','>','*','s','o','+','^','D' ]#['x','^','D','v','<','>','x','s','o','+','o']

N_mat2[0] = [0,1,2,3,5,7,9,11,13,15]
N_mat2[1] = [0,1,2,4,6,8,10,12,14]

N_mat2[3] = [0,1,3,5,7,9,11,13,15]
N_mat2[4] = [0,2,4,6,8,10,12,14]

N_mat2[5] = [0,3,6,9,12,15,12,14]
N_mat2[6] = [0,1,4,7,10,13]
N_mat2[7] = [0,2,5,8,11,14]


fig1=plt.figure(2,figsize=figsz)
ax1=fig1.add_axes(axshape)  

# =============================================================================
# for ind in range(N):
#     if colorvec[ind]==col_pluri:
#         pass
#     else:
#             
#         if plot_vec[ind]==1:
#             xdata = N_mat[ind][N_mat2[ind]]
#             ydata = np.array(Scaling_DF[ind][0][inds1_mat[ind]])[N_mat2[ind]]
#             ax1.plot(xdata,ydata,color=colorvec[ind],linewidth=lw_vec[ind],linestyle='-',marker=markers[ind],markersize=mksize_vec[ind],markeredgewidth=1., markerfacecolor='white')#colorvec[ind])        
#     
#         elif plot_vec[ind]==2:
#             xdata = N_mat[ind][N_mat2[ind]]
#             ydata = np.array(Scaling_DF[ind][0][inds1_mat[ind]])[N_mat2[ind]]
#             ydata[np.isnan(ydata)]=1.
#             ax1.plot(xdata,ydata,color=colorvec[ind],linewidth=lw_vec[ind],linestyle='-',marker=markers[ind],markersize=mksize_vec[ind],markeredgewidth=1., markerfacecolor='white')#colorvec[ind])        
# 
# 
# =============================================================================

#%% quasi-crystal data

xdata=[1,2,4,8,16,32]
ydata=[2,3,4,5,6,7]
col = col_oligo
lw=5.
mk = 'p'
mksize=12.

ax1.plot(xdata,ydata,color=col,linewidth=lw,linestyle='--',marker=mk,markersize=mksize,markeredgewidth=1.5, markerfacecolor='white')#colorvec[ind])        


ax1.set_xlabel('$n$')
ax1.set_ylabel('$N_M$')
ax1.xaxis.set_label_coords(0.5, -0.18)
ax1.yaxis.set_label_coords(-0.18, 0.5)


ax1.tick_params(which='both',direction='out')
ax1.set_xlim([0,32])
ax1.set_ylim([0,16])
ax1.set_xticks([0,8,16,24,32])
ax1.set_yticks([0,4,8,12,16])
ax1.minorticks_on()
#ax1.tick_params(bottom=True,top=True,left=True,right=True)

#ax1.grid(which='both')



FigName = '//Fig2_AppendixQuasi'
fig1.savefig(CurrentFolder+FigName+'.png',dpi=600,transparent=True)
fig1.savefig(CurrentFolder+FigName+'.pdf')






#%% Plot Fig2 Inlet


lw=5.
mk = 'p'
mksize=10.


figsz=(6/3,5/3)
axshape = [0.26,0.32,0.65,0.6]





fig2=plt.figure(3,figsize=figsz)
ax2=fig2.add_axes(axshape)  


ax2.plot(xdata,ydata,color=col,linewidth=lw,linestyle='-',marker=mk,markersize=mksize,markeredgewidth=1.5, markerfacecolor='white')#colorvec[ind])        

#for ind in [3]:#range(N):
#    ax2.plot(N_mat[ind],Scaling_DF[ind][0][inds1_mat[ind]],color=colorvec[ind],linestyle='--',marker=markers[ind],markersize=10,markeredgewidth=2., markerfacecolor='white')#colorvec[ind])        



#ax2.set_xlabel('$N$')
#ax2.set_ylabel('$M$')
ax2.xaxis.set_label_coords(0.5, -0.1)
ax2.yaxis.set_label_coords(-0.1, 0.5)

ax2.set_xlim([1,32])
ax2.set_ylim([0,8])
ax2.set_yticks([0,4,8])

ax2.set_xscale('log',basex=2)
ax2.set_xticks((1,4,16))


FigName = '//Fig2_AppendixQuasiInlet'
fig2.savefig(CurrentFolder+FigName+'.png',dpi=600,transparent=True)
fig2.savefig(CurrentFolder+FigName+'.pdf')




        