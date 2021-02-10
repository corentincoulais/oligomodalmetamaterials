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

Drive='G'

cmapOmega = 'seismic'
cmapPoisson = 'PRGn'

figsz=(5.,5.)
axshape = [0.23,0.23,0.7,0.7]


OutputFolderName = os.getcwd()


mpl.rcParams['axes.labelsize'] = 24
mpl.rcParams['xtick.labelsize'] =  24
mpl.rcParams['ytick.labelsize'] = 24
mpl.rcParams['figure.subplot.left'] = 0.25
mpl.rcParams['figure.subplot.right'] = 0.95
mpl.rcParams['figure.subplot.bottom'] = 0.15
mpl.rcParams['figure.subplot.top'] = 0.95
#mpl.rcParams['figure.figsize']=[ 8,8]
mpl.rcParams['font.family']='Calibri'
mpl.rcParams['lines.linewidth']=5
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
        


#%%######################################################################################################
#  Polarisation, varying angle (different loading speeds)
#########################################################################################################


FoldersExperiments = r'G:\15_Viscoelastic_Metamaterial\Multimode2_v02\Experiments\20190913\tmpDB\\'
runs=[4,5]
Nruns=len(runs)

Nx = 9
Ny = 11

Ellipsetitles = ['ie','t', 'elXcm', 'elYcm', 'a_el', 'b_el', 'phi']
Ellipses={}
for indr in runs:
    Ellipses[indr] = pd.read_csv(FoldersExperiments+str(indr)+'\\Ellipse.txt',sep=';',names=Ellipsetitles)
    Ellipses[indr]['Omega'] = (1-Ellipses[indr]['a_el']/Ellipses[indr]['b_el'])*np.cos(2*(Ellipses[indr]['phi']-np.pi/4))



Ellipsesmat={}
Ellipsesarr={}
for indr in runs:
    nie = np.max(Ellipses[indr]['ie'])+1
    nt = np.max(Ellipses[indr]['t'])+1
    Niet = np.shape(Ellipses[indr])[0]
    
    Ellipsesarr[indr] = np.array(Ellipses[indr])

    
    Ellipsesmat[indr] = np.zeros(((nie,nt,len(Ellipses[indr].keys())-2)))
    for iniet in range(Niet):
        Ellipsesmat[indr][Ellipses[indr]['ie'][iniet],Ellipses[indr]['t'][iniet],:]=Ellipsesarr[indr][iniet,2:]

    
    
# Correspondence
Nxyvec = {}

indsbottomtop={}
sidesleftright={}
for indr in runs:
    nie = np.max(Ellipses[indr]['ie'])+1
    sidesleftright[indr]=np.zeros((Ny,2),dtype=int)
    
    Nxyvec[indr] = np.zeros((nie,4),dtype=int)
    
    indsbottomtop[indr] = 2*[None]
        
    vec = copy(Ellipsesmat[4][:,0,0])
    elXnorm0 = (vec-min(vec))/(max(vec)-min(vec))
    vec = copy(Ellipsesmat[4][:,0,1])
    elYnorm0 = (vec-min(vec))/(max(vec)-min(vec))
    

    Nxyvec[indr][:,0] = range(nie)
    Nxyvec[indr][:,1] = np.round(elXnorm0*(Nx-1))
    Nxyvec[indr][:,2] = np.round(elYnorm0*(Ny-1))
    Nxyvec[indr][:,3] = (-1)**np.round(elYnorm0*(Ny-1))

    for indie in range(nie):
        Ellipsesmat[indr][indie,:,5] = Ellipsesmat[indr][indie,:,5]*Nxyvec[indr][indie,3]
        
# =============================================================================
#     indsbottomtop[indr][0] = np.where(Nxyvec[indr][:,2]==0)[0]
#     indsbottomtop[indr][1] = np.where(Nxyvec[indr][:,2]==Ny-1)[0]
# =============================================================================
        
    indsbottomtop[indr][0] = np.where(Nxyvec[indr][:,2]==1)[0]
    indsbottomtop[indr][1] = np.where(Nxyvec[indr][:,2]==Ny-2)[0]
    
    for inds in range(0,Ny,2):#[0,2,4,6,8,10]:
        sidesleftright[indr][inds,0]=np.where((Nxyvec[indr][:,2]==inds) & (Nxyvec[indr][:,1]==1))[0]
        sidesleftright[indr][inds,1]=np.where((Nxyvec[indr][:,2]==inds) & (Nxyvec[indr][:,1]==Nx-2))[0]
        
    for inds in range(1,Ny,2):
        sidesleftright[indr][inds,0]=np.where((Nxyvec[indr][:,2]==inds) & (Nxyvec[indr][:,1]==0))[0]
        sidesleftright[indr][inds,1]=np.where((Nxyvec[indr][:,2]==inds) & (Nxyvec[indr][:,1]==Nx-1))[0]
        
        

Omega_avg={}
Omega_min={}
Omega_max={}
elYcm_min={}
elYcm_max={}
Ylen={}
eps={}
Ylenminmax={}
epsminmax={}

Xlen_local = {}
Xlen = {}
epsx_local={}
epsx_avg={}

Neps = 6
epslocs = np.linspace(0.00,0.05,Neps) 
epslocs_inds={}

epsylocs={}
epsx_locallocs={}
epsx_avglocs={}
epsx_cornerslocs={}
Poisson_avg={}
Poisson_corners={}
Poisson_local={}


for indr in runs:
    nt = np.max(Ellipses[indr]['t'])+1

    Omega_avg[indr] = np.average(Ellipsesmat[indr][:,:,5],axis=0)
    Omega_min[indr] = np.min(Ellipsesmat[indr][:,:,5],axis=0)
    Omega_max[indr] = np.max(Ellipsesmat[indr][:,:,5],axis=0)
    elYcm_min[indr] = np.min(Ellipsesmat[indr][:,:,1],axis=0)
    elYcm_max[indr] = np.max(Ellipsesmat[indr][:,:,1],axis=0)
    Ylenminmax[indr] = elYcm_max[indr]-elYcm_min[indr]
    epsminmax[indr] = (Ylenminmax[indr]-Ylenminmax[indr][0])/Ylenminmax[indr][0]
       
    Ylen[indr] = np.average(Ellipsesmat[indr][indsbottomtop[indr][1][range(0,5,1)],:,1],axis=0)-np.average(Ellipsesmat[indr][indsbottomtop[indr][0][range(0,5,1)],:,1],axis=0)
    eps[indr] = (Ylen[indr]-Ylen[indr][0])/Ylen[indr][0]
    
    Xlen_local[indr]=np.zeros((nt,Ny))
    epsx_local[indr]=np.zeros((nt,Ny))
    for  inds in range(Ny):
        Xlen_local[indr][:,inds] = Ellipsesmat[indr][sidesleftright[indr][inds,1],:,0]-Ellipsesmat[indr][sidesleftright[indr][inds,0],:,0]
        epsx_local[indr][:,inds]= (Xlen_local[indr][:,inds]-Xlen_local[indr][0,inds])/Xlen_local[indr][0,inds]
        #Poisson_local[indr][:,inds] = 
        
    plt.figure()
    plt.xlabel('$\epsilon_y$')
    plt.ylabel('$\epsilon_x$')
    for  inds in range(Ny):
        plt.plot(-eps[indr],-epsx_local[indr][:,inds])
        
    epslocs_inds[indr]=np.zeros(len(epslocs),dtype=int)
    for indeps in range(Neps):
        epslocs_inds[indr][indeps] = np.argmin(np.abs(-eps[indr]-epslocs[indeps]))

    epsylocs[indr] = eps[indr][epslocs_inds[indr]]
    epsx_locallocs[indr] = epsx_local[indr][epslocs_inds[indr],:]
    epsx_avglocs[indr] = np.average(epsx_locallocs[indr][:,range(0,11,2)],axis=1)
    
    depsylocs = epsylocs[indr][1:]-epsylocs[indr][:-1]
    depsx_locallocs = epsx_locallocs[indr][1:,:]-epsx_locallocs[indr][:-1,:]
    depsx_avglocs = epsx_avglocs[indr][1:]-epsx_avglocs[indr][:-1]


    depsylocsmat = np.repeat(depsylocs[:,np.newaxis],Ny,axis=1)
    Poisson_avg[indr]=-depsx_avglocs/depsylocs
    Poisson_local[indr]=    -depsx_locallocs/depsylocsmat
    
    
    plt.figure()
    plt.xlabel(r'$\epsilon_y$')
    plt.ylabel(r'$ \nu $')
    plt.ylim([-0.6,0.6])
    #plt.plot(epslocs[1:],Poisson_avg[indr],'--o')
    plt.errorbar(epslocs[1:],Poisson_avg[indr],yerr=0.1,fmt='o',ecolor='k',capthick=2,elinewidth=5,capsize=10)
    plt.plot(epslocs[1:],Poisson_avg[indr],'--o',markersize=10)
    



plt.figure()
plt.plot(-eps[4],np.abs(Omega_avg[4]))    
plt.plot(-eps[5][:50],np.abs(Omega_avg[5][:50]))    
plt.xlabel('$\epsilon$')
plt.ylabel('$\Omega$')
plt.savefig('EpsilonOmegaExperimental.pdf')
plt.savefig('EpsilonOmegaExperimental.png',dpi=300,transparent=True)




EpsOmega4={}
EpsOmega4['Epsilon'] = eps[4]
EpsOmega4['Omega_avg'] = Omega_avg[4]
EpsOmega4['Omega_min'] = Omega_min[4]
EpsOmega4['Omega_max'] = Omega_max[4]
EpsOmega4_DF = pd.DataFrame(EpsOmega4)
EpsOmega4_DF.to_csv('EpsOmega4.csv',sep = ',')

EpsOmega5={}
EpsOmega5['Epsilon'] = eps[5]
EpsOmega5['Omega_avg'] = Omega_avg[5]
EpsOmega5['Omega_min'] = Omega_min[5]
EpsOmega5['Omega_max'] = Omega_max[5]
EpsOmega5_DF = pd.DataFrame(EpsOmega5)
EpsOmega5_DF.to_csv('EpsOmega5.csv',sep = ',')


# =============================================================================
# #%% Plot Omega average
#     
# xlims   =   [0,0.1]
# ylims   =   [0.,0.6]
# xticksplot = [0,0.05,0.1]
# yticksplot = [0,0.3,0.6]
#     
# nb=100
# fig=plt.figure(nb,figsize=figsz )
# ax=fig.add_axes(axshape)
# ax = plt.gca()
# 
# for ind0 in range(N_angles):
#     ax.plot(Data_AngleOmega_fast[ind0]['Epsilon'],np.absolute(Data_AngleOmega_fast[ind0]['Omega_avg']),color = colors[ind0])
#     
# ax.tick_params(axis = 'both', which = 'major')  
# ax.set_xlabel(r'$\epsilon$')
# ax.set_ylabel(r'$\Omega$')
# ax.set_xlim(xlims)    
# ax.set_ylim(ylims)
# ax.set_xticks(xticksplot)
# ax.set_yticks(yticksplot)
# 
# ax.xaxis.set_label_coords(0.5, -0.2)
# ax.yaxis.set_label_coords(-0.2, 0.5)
# 
# 
# if savefigures:
#     FigName = OutputFolderName+'\\Epsilon-OmegaAvg_VarAngle_Fast'
#     plt.savefig(FigName+'.pdf')
#     plt.savefig(FigName+'.png',dpi=300,Transparent=True)
# 
# 
# 
# =============================================================================







