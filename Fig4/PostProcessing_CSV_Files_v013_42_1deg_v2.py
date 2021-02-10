# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 09:53:49 2020

@author: David
"""


import numpy as np
import pandas as pd
import os
import scipy as sc
import matplotlib as mpl
from matplotlib import pyplot as plt
from copy import deepcopy as copy

plt.close('all')

savefigures=True


cmapOmega = 'seismic'
cmapPoisson = 'PRGn'


Drive = 'G'
# =============================================================================
# figsz=(4.,4.)
# axshape = [0.28,0.28,0.62,0.62]
# axshapeCB = [0.28,0.28,0.4,0.62]
# 
# =============================================================================
figsz=(7.5,4.)
axshape = [0.18,0.18,0.75,0.8]
axshapeCB = [0.15,0.15,0.75,0.8]
xtickvals = [0,0.02,0.04,0.06,0.08,0.1]
ytickvals = [1e-6,1e-4,1e-2]

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
    



#%% Define functions


def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)
        
def SmoothStepAbaqus(t,t0=0,t1=1,A0=0,A1=1): #Normalised, assuming always between 0 and 1
    if len(t)==1:
        if t<t0:    Amp=A0
        elif t>t1:    Amp=A1
        else:
            ksi = (t-t0)/(t1-t0)
            Amp = A0+(A1-A0)*ksi**3*(10-15*ksi+6*ksi**2)
    else:
        Amp = np.zeros(len(t))        
        for ind in range(len(t)):
            if t[ind]<t0:    Amp[ind]=A0
            elif t[ind]>t1:    Amp[ind]=A1
            else:
                ksi = (t[ind]-t0)/(t1-t0)
                Amp[ind] = A0+(A1-A0)*ksi**3*(10-15*ksi+6*ksi**2)
    return Amp

def is_odd(num):
   return num % 2 != 0
#%% Read output folder


CSVFolderName =Drive+r':\15_Viscoelastic_Metamaterial\AbaqusNew\Model_v013\batch_v13_42_1deg_Visco-Agilus_improved_norelax19mm\OutputOmega'

OutputFolderName = copy(CSVFolderName)

#%% Basic data
Loadingdist = 18.75
NormFacLoading = 180. # WRT point of load application

t_loading_MAT = 3*[1.5e-1,1.5e0,1.5e1,1.5e2,1.5e3,1.5e4]		+3*[1.5e5]		+3*[1.5*2.15e2,1.5*4.64e2,1.5*2.15e3,1.5*4.64e3]  	+[1.5e0*0.69/0.11,1.5e5*(9.3e-6)/(6.9e-7)]
E_TriH_vec = 6*[1.]+6*[2.]+6*[0.5]		+[1.,2.,0.5]		+4*[1.]+4*[2.]+4*[0.5]		+2*[1.]

eps_loading = Loadingdist/NormFacLoading
epsdot_loading_vec = np.array(t_loading_MAT)**(-1)*eps_loading
MaterialIDS = [1,2,3,  21,22,     4,    23,24,  5,   18]
MaterialIDS = [1,2,3,  21,22,     4,    23,24,     18]

N_IDs = len(MaterialIDS)



#%% Omega plots

CSV_DFs ={} 
#Index(['Unnamed: 0', 'Time', 'Epsilon', 'Omega_min', 'Omega_max', 'Omega_avg','domain_y_loc'],      dtype='object')


maxintcomp = 250
for ind3 in MaterialIDS:
    DFName = CSVFolderName+'\\Omega-minmaxavg_MAT-'+str(ind3)+'.csv'    
    CSV_DFs[ind3]=pd.read_csv(DFName)

maxintcomp = min(np.shape(CSV_DFs[1])[0],maxintcomp)
#CSVInput_Mat = np.zeros(((np.shape(CSV_DFs[1])[0],np.shape(CSV_DFs[1])[1],N_IDs)))
CSVInput_Mat = np.zeros(((maxintcomp,np.shape(CSV_DFs[1])[1],N_IDs)))


ind4=-1
for ind3 in MaterialIDS:
    ind4+=1
    CSVInput_Mat[:,:,ind4] = np.array(CSV_DFs[ind3])[:maxintcomp,:]
    
    
eps_mat = CSVInput_Mat[:,1,:]
epsdot_mat,eps_mat2 = np.meshgrid(epsdot_loading_vec[MaterialIDS],CSVInput_Mat[:,1,0])

Omega_absavgmat = np.absolute(CSVInput_Mat[:,2,:])
Omega_absmaxmat = np.max(   np.absolute(CSVInput_Mat[:,3:5,:]),axis=1)




stepscmap=101


nb=100
fig=plt.figure(nb,figsize=figsz )
ax=fig.add_axes(axshape)
ax = plt.gca()
ax.tick_params(axis = 'both', which = 'major')  
ax.set_xlabel('$\epsilon$')
ax.set_ylabel('$\dot{\epsilon}$ [s$^{-1}$]')
ax.set_yscale('log')
ax.set_xlim([0,0.1])    
ax.set_ylim([np.min(epsdot_loading_vec[MaterialIDS]),np.max(epsdot_loading_vec[MaterialIDS])])

ax.set_xticks(xtickvals)
ax.set_yticks(ytickvals)

ax.xaxis.set_label_coords(0.5, -0.12)
ax.yaxis.set_label_coords(-0.14, 0.5)
ax.contourf(eps_mat,epsdot_mat,Omega_absavgmat,levels=np.linspace(-1,1, stepscmap),cmap=cmapOmega,extend='both')

if savefigures:
    FigName = OutputFolderName+'\\Epsilon-Epsilondot-absOmega_avg-Contour'
    plt.savefig(FigName+'.pdf')
    plt.savefig(FigName+'.png',dpi=300,Transparent=True)


nb=101
fig=plt.figure(nb,figsize=figsz )
ax=fig.add_axes(axshape)
ax = plt.gca()
ax.tick_params(axis = 'both', which = 'major')  
ax.set_xlabel('$\epsilon$')
ax.set_ylabel('$\dot{\epsilon}$ [s$^{-1}$]')
ax.set_yscale('log')
ax.set_xlim([0,0.1])    
ax.set_ylim([np.min(epsdot_loading_vec[MaterialIDS]),np.max(epsdot_loading_vec[MaterialIDS])])

ax.set_xticks(xtickvals)
ax.set_yticks(ytickvals)

ax.xaxis.set_label_coords(0.5, -0.12)
ax.yaxis.set_label_coords(-0.14, 0.5)
ax.contourf(eps_mat,epsdot_mat,Omega_absmaxmat,levels=np.linspace(-1,1, stepscmap),cmap=cmapOmega,extend='both')

if savefigures:
    FigName = OutputFolderName+'\\Epsilon-Epsilondot-absOmega_max-Contour'
    plt.savefig(FigName+'.pdf')
    plt.savefig(FigName+'.png',dpi=300,Transparent=True)

nb=102
fig=plt.figure(nb,figsize=figsz )
ax=fig.add_axes(axshapeCB)
ax = plt.gca()
ax.tick_params(axis = 'both', which = 'major')  
ax.set_xlabel('$\epsilon$')
ax.set_ylabel('$\dot{\epsilon}$ [s$^{-1}$]')
ax.set_yscale('log')
ax.set_xlim([0,0.1])    
ax.set_ylim([np.min(epsdot_loading_vec[MaterialIDS]),np.max(epsdot_loading_vec[MaterialIDS])])

ax.set_xticks(xtickvals)
ax.set_yticks(ytickvals)

ax.xaxis.set_label_coords(0.5, -0.12)
ax.yaxis.set_label_coords(-0.14, 0.5)
im=ax.contourf(eps_mat,epsdot_mat,Omega_absmaxmat,levels=np.linspace(-1,1, stepscmap),cmap=cmapOmega,extend='neither')
fig.colorbar(im,ax=ax,ticks=[-1,-0.5,0,0.5,1])



if savefigures:
    FigName = OutputFolderName+'\\Epsilon-Epsilondot-absOmega_max-ContourCOLORBAR'
    plt.savefig(FigName+'.pdf')
    plt.savefig(FigName+'.png',dpi=300,Transparent=True)


#%%


axshapetest = [-0.7,0.1,1.,0.8]
nb=103
fig=plt.figure(nb,figsize=(figsz[1]/3,figsz[1]) )
ax2=fig.add_axes(axshapetest)
ax2.get_xaxis().set_visible(False)
ax2.get_yaxis().set_visible(False)
fig.colorbar(im,ax=ax2,ticks=[-1,-0.5,0,0.5,1])


if savefigures:
    FigName = OutputFolderName+'\\COLORBAR-OMEGA-Contour'
    plt.savefig(FigName+'.pdf')
    plt.savefig(FigName+'.png',dpi=300,Transparent=True)

#%% Poisson's ratio



stepscmap=1001

Poissonlims =2.

PoissonDFs = {}
# =============================================================================
# Index(['Unnamed: 0', 'Time', 'epsilon_loading', '-epsilon_y_holesavg',
#        'epsilon_x_holesavg', 'epsilon_y_holesmid', 'epsilon_x_holesmid',
#        'Poisson_avg_yloading', 'Poisson_avg_yholesavg'],
#       dtype='object')
# =============================================================================
for ind3 in MaterialIDS:
    PoissonDFName =     CSVFolderName+'\\Poisson-epsilon_MAT-'+str(ind3)+'.csv'
    PoissonDFs[ind3]=pd.read_csv(PoissonDFName)


PoissonInput_Mat = np.zeros(((maxintcomp,np.shape(PoissonDFs[ind3])[1],N_IDs)))


ind4=-1
for ind3 in MaterialIDS:
    ind4+=1
    PoissonInput_Mat[:,:,ind4] = np.array(PoissonDFs[ind3])[:maxintcomp,:]
    
    
    
# =============================================================================
#     
# eps_loading_mat = PoissonInput_Mat[:,2,:]
# eps_yholesavg_mat = PoissonInput_Mat[:,3,:]
# Poisson_avg_yloading = PoissonInput_Mat[:,7,:]
# Poisson_avg_yholesavg = PoissonInput_Mat[:,8,:]
# 
# epsdot_mat2,eps_mat2 = np.meshgrid(epsdot_loading_vec[MaterialIDS],PoissonInput_Mat[:,2,0])
# 
# =============================================================================

eps_loading_mat = PoissonInput_Mat[:,5,:]
eps_yholesavg_mat = PoissonInput_Mat[:,1,:]
Poisson_avg_yloading = PoissonInput_Mat[:,3,:]
Poisson_avg_yholesavg = PoissonInput_Mat[:,2,:]

epsdot_mat2,eps_mat2 = np.meshgrid(epsdot_loading_vec[MaterialIDS],PoissonInput_Mat[:,5,0])




nb=200
fig=plt.figure(nb,figsize=figsz )
ax=fig.add_axes(axshape)
ax = plt.gca()
ax.tick_params(axis = 'both', which = 'major')  
ax.set_xlabel('$\epsilon$')
ax.set_ylabel('$\dot{\epsilon}$ [s$^{-1}$]')
ax.set_yscale('log')
ax.set_xlim([0,0.1])    
ax.set_ylim([np.min(epsdot_loading_vec[MaterialIDS]),np.max(epsdot_loading_vec[MaterialIDS])])

ax.set_xticks(xtickvals)
ax.set_yticks(ytickvals)

ax.xaxis.set_label_coords(0.5, -0.12)
ax.yaxis.set_label_coords(-0.14, 0.5)
ax.contourf(eps_loading_mat,epsdot_mat2,Poisson_avg_yloading,levels=np.linspace(-Poissonlims,Poissonlims, stepscmap),cmap=cmapPoisson,extend='both')


if savefigures:
    FigName = OutputFolderName+'\\Epsilon-Epsilondot-Poisson_yloading-Contour'
    plt.savefig(FigName+'.pdf')
    plt.savefig(FigName+'.png',dpi=300,Transparent=True)


nb=201
fig=plt.figure(nb,figsize=figsz )
ax=fig.add_axes(axshape)
ax = plt.gca()
ax.tick_params(axis = 'both', which = 'major')  
ax.set_xlabel('$\epsilon$')
ax.set_ylabel('$\dot{\epsilon}$ [s$^{-1}$]')
ax.set_yscale('log')
ax.set_xlim([0,0.1])    
ax.set_ylim([np.min(epsdot_loading_vec[MaterialIDS]),np.max(epsdot_loading_vec[MaterialIDS])])

ax.set_xticks(xtickvals)
ax.set_yticks(ytickvals)

ax.xaxis.set_label_coords(0.5, -0.12)
ax.yaxis.set_label_coords(-0.14, 0.5)
ax.contourf(eps_yholesavg_mat,epsdot_mat2,Poisson_avg_yholesavg,cmap=cmapPoisson,levels=np.linspace(-Poissonlims,Poissonlims, stepscmap),extend='both')


if savefigures:
    FigName = OutputFolderName+'\\Epsilon-Epsilondot-Poisson_yholes-Contour'
    plt.savefig(FigName+'.pdf')
    plt.savefig(FigName+'.png',dpi=300,Transparent=True)



nb=202
fig=plt.figure(nb,figsize=figsz )
ax=fig.add_axes(axshapeCB)
ax = plt.gca()
ax.tick_params(axis = 'both', which = 'major')  
ax.set_xlabel('$\epsilon$')
ax.set_ylabel('$\dot{\epsilon}$ [s$^{-1}$]')
ax.set_yscale('log')
ax.set_xlim([0,0.1])    
ax.set_ylim([np.min(epsdot_loading_vec[MaterialIDS]),np.max(epsdot_loading_vec[MaterialIDS])])

ax.set_xticks(xtickvals)
ax.set_yticks(ytickvals)

ax.xaxis.set_label_coords(0.5, -0.12)
ax.yaxis.set_label_coords(-0.14, 0.5)
im=ax.contourf(eps_yholesavg_mat,epsdot_mat2,Poisson_avg_yholesavg,cmap=cmapPoisson,levels=np.linspace(-Poissonlims,Poissonlims, stepscmap),extend='both')
fig.colorbar(im,ax=ax,ticks=[-2.,-1.,0.,1.,2.])


if savefigures:
    FigName = OutputFolderName+'\\Epsilon-Epsilondot-Poisson_yholes-ContourCOLORBAR'
    plt.savefig(FigName+'.pdf')
    plt.savefig(FigName+'.png',dpi=300,Transparent=True)



nb=203
fig=plt.figure(nb,figsize=(figsz[1]/3,figsz[1]) )
ax2=fig.add_axes(axshapetest)
ax2.get_xaxis().set_visible(False)
ax2.get_yaxis().set_visible(False)
fig.colorbar(im,ax=ax2,ticks=[-2,-1,0,1,2])


if savefigures:
    FigName = OutputFolderName+'\\COLORBAR-POISSON-Contour'
    plt.savefig(FigName+'.pdf')
    plt.savefig(FigName+'.png',dpi=300,Transparent=True)










