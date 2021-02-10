# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 12:02:26 2020

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


def ComputePolarisationMat(holesinds,Xcorneravg,Ycorneravg,plusminfac):    #Note: Square matrices only
    # plusminfac multiplies polarisation with +1 or -1 depending on requirement. Can be implemented either as integer or 2D array of same size as output
    Nholes = len(holesinds)
    a_mat=np.zeros((Nholes,Nholes));    b_mat=np.zeros((Nholes,Nholes));    angle_mat=np.zeros((Nholes,Nholes));    Omega_mat=np.zeros((Nholes,Nholes));    xcent_mat = np.zeros((Nholes,Nholes));    ycent_mat = np.zeros((Nholes,Nholes)); 
    
    for ind0 in range(Nholes):
        for ind1 in range(Nholes):
            holeind0 = holes1inds[ind0];       holeind1 = holes1inds[ind1]
            
            elldiag1 =np.sqrt( (Xcorneravg[holeind0+1,holeind1]-Xcorneravg[holeind0,holeind1+1])**2 + (Ycorneravg[holeind0+1,holeind1]-Ycorneravg[holeind0,holeind1+1])**2  )
            elldiag2 =np.sqrt( (Xcorneravg[holeind0+1,holeind1+1]-Xcorneravg[holeind0,holeind1])**2 + (Ycorneravg[holeind0+1,holeind1+1]-Ycorneravg[holeind0,holeind1])**2  )
            if elldiag1 >=elldiag2:
                a_el = 1.*elldiag1;          b_el = 1.*elldiag2
                angle = np.pi/4 + np.arctan((Ycorneravg[holeind0+1,holeind1]-Ycorneravg[holeind0,holeind1+1])      /        (Xcorneravg[holeind0+1,holeind1]-Xcorneravg[holeind0,holeind1+1]))
            elif elldiag2 >=elldiag1:
                a_el = 1.*elldiag2;          b_el = 1.*elldiag1
                angle = np.pi/4 + np.arctan((Ycorneravg[holeind0+1,holeind1+1]-Ycorneravg[holeind0,holeind1])      /        (Xcorneravg[holeind0+1,holeind1+1]-Xcorneravg[holeind0,holeind1]))
            Omega = (1-b_el/a_el)*np.cos(2*(angle))
            
            xcent_mat[ind0,ind1] = np.average(Xcorneravg[holeind0:holeind0+2,holeind1:holeind1+2])
            ycent_mat[ind0,ind1] = np.average(Ycorneravg[holeind0:holeind0+2,holeind1:holeind1+2])
            
            
            a_mat[ind0,ind1] = 1.*a_el;    b_mat[ind0,ind1] = 1.*b_el;        angle_mat[ind0,ind1] = 1.*angle; Omega_mat[ind0,ind1] = 1.*Omega; 
            Omega_mat = plusminfac *Omega_mat 
            
            #xcent_mat[]
    return(a_mat,b_mat,angle_mat,Omega_mat,xcent_mat,ycent_mat)



def Corneravg(NodesCornerDistIndex,Xnormnodes,Ynormnodes,U1_Data_Norm,U2_Data_Norm):

    Xcorneravg=-np.ones(((Regimes+1,Regimes+1,TimeSteps)))
    Ycorneravg=-np.ones(((Regimes+1,Regimes+1,TimeSteps)))
    
    for ind0 in range(Regimes+1):
        for ind1 in range(Regimes+1):
            indices= NodesCornerDistIndex[ind0][ind1]             
            for ind2 in range(TimeSteps):            
                if len(indices)>0.1:
                    Xarr = Xnormnodes+U1_Data_Norm[ind2,1:]
                    Yarr = Ynormnodes+U2_Data_Norm[ind2,1:]
                    Xcorneravg[ind0,ind1,ind2] = np.average(Xarr[indices])    
                    Ycorneravg[ind0,ind1,ind2] = np.average(Yarr[indices])    
    return Xcorneravg, Ycorneravg


def HoleSpecAnalysis(holes1inds,holes2inds,Xcorneravg,Ycorneravg,U1_Data_Norm,U2_Data_Norm):
    Nholes1 = len(holes1inds);    Nholes2 = len(holes2inds)
    TimeSteps = np.shape(U1_Data_Norm)[0]
    
    zeroarr = np.zeros(((Nholes1,Nholes1,TimeSteps))) 
    a_holes1=1.*zeroarr;    b_holes1=1.*zeroarr;    angle_holes1=1.*zeroarr;    Omega_holes1=1.*zeroarr;        xcent1=1.*zeroarr;   ycent1=1.*zeroarr    
    
    zeroarr = np.zeros(((Nholes2,Nholes2,TimeSteps))) 
    a_holes2=1.*zeroarr;    b_holes2=1.*zeroarr;     angle_holes2=1.*zeroarr;     Omega_holes2=1.*zeroarr;      xcent2=1.*zeroarr;   ycent2=1.*zeroarr    
    
    for ind2 in range(TimeSteps):
        a_holes1[:,:,ind2], b_holes1[:,:,ind2],angle_holes1[:,:,ind2], Omega_holes1[:,:,ind2], xcent1[:,:,ind2], ycent1[:,:,ind2] = ComputePolarisationMat(holes1inds,Xcorneravg[:,:,ind2],Ycorneravg[:,:,ind2],1)
        a_holes2[:,:,ind2], b_holes2[:,:,ind2],angle_holes2[:,:,ind2], Omega_holes2[:,:,ind2], xcent2[:,:,ind2], ycent2[:,:,ind2] = ComputePolarisationMat(holes2inds,Xcorneravg[:,:,ind2],Ycorneravg[:,:,ind2],1)


    return a_holes1, b_holes1, angle_holes1, Omega_holes1,      a_holes2, b_holes2, angle_holes2, Omega_holes2,     xcent1,ycent1,xcent2,ycent2


def OmegaDist(Omega_holes1,Omega_holes2):
    TimeSteps = np.shape(Omega_holes1)[2]
    zeroarr= 1.*np.zeros(TimeSteps)
    Omega_avg=1.*zeroarr;    Omega_min=1.*zeroarr;    Omega_max=1.*zeroarr           
    for ind2 in range(TimeSteps):            
        Omega_avg[ind2] = (np.average(Omega_holes1[:,:,ind2])*np.size(Omega_holes1[:,:,ind2])+np.average(Omega_holes2[:,:,ind2])*np.size(Omega_holes2[:,:,ind2])) / (np.size(Omega_holes1[:,:,ind2])+np.size(Omega_holes2[:,:,ind2])) 
        Omega_min[ind2] = np.min([np.min(Omega_holes1[:,:,ind2]),np.min(Omega_holes2[:,:,ind2])])
        Omega_max[ind2] = np.max([np.max(Omega_holes1[:,:,ind2]),np.max(Omega_holes2[:,:,ind2])])
        

    return Omega_avg,Omega_min,Omega_max
#%% Create output folder
    
OutputFolderName ='OutputOmega' 
createFolder(OutputFolderName)

#%% Load data
Loadingdist = 18.75
t_loading_MAT = 3*[1.5e-1,1.5e0,1.5e1,1.5e2,1.5e3,1.5e4]        +3*[1.5e5]        
E_TriH_vec = 6*[1.]+6*[2.]+6*[0.5]        +[1.,2.,0.5]        

MaterialIDS = list(range(0,18))


N_IDs = len(MaterialIDS)



for ind3 in MaterialIDS:
    
    t_loading = t_loading_MAT[ind3]
    E_TriH = E_TriH_vec[ind3]
    
    RPT_U_Name = 'Output_U_v013_041-Mat-'+str(ind3)+'.rpt'
    U_Data_RAW = pd.read_csv(RPT_U_Name,sep='\s{2,}',skiprows=2)
    csv_XY_Name = 'Output_XYstart_v013_041.csv'
    XY_Data_RAW = pd.read_csv(csv_XY_Name,sep=',')
    
    Regimes = 13
    
    s_Regime=1./Regimes
    U_Data_array = np.array(U_Data_RAW,dtype=float)
    TimeSteps = np.shape(U_Data_array)[0]
    Time = U_Data_array[:,0]
    IndsCompression = np.where(Time<t_loading)[0] 
     
    d_loading_vec = SmoothStepAbaqus(Time,0,t_loading,0,Loadingdist)
    UtitlesRAW = U_Data_RAW.keys()
    numnodes = int((len(UtitlesRAW)-1)/2)
    
    U1titlesRAW = UtitlesRAW[1:numnodes+1]
    U2titlesRAW = UtitlesRAW[1+numnodes:2*numnodes+1]
    
    Nodes = numnodes*[0]
    for ind0 in range(numnodes):
        Nodes[ind0] = int(U1titlesRAW[ind0][12:])
    Nodes = np.array(Nodes,dtype=int)
    
    Titles_tNodes = ['t']+list(Nodes)
    U1_Data_array=U_Data_array[:,range(0,numnodes+1)]
    U1_DF = pd.DataFrame(U1_Data_array,columns = Titles_tNodes)
        
    U2_Data_array=U_Data_array[:,[0]+list(range(1+numnodes,2*numnodes+1))]
    U2_DF = pd.DataFrame(U2_Data_array,columns = Titles_tNodes)
        
    
    
    
    
    #%% Process base XY data
    indiceslocations = np.where(XY_Data_RAW['Part Instance Name']=='ASSEMBLYPART-1')[0]
    all_nodes_titles = XY_Data_RAW['    Node Label'][indiceslocations]
    
    Xnodesfull = np.array(XY_Data_RAW[['X']])[indiceslocations][:,0]
    Ynodesfull = np.array(XY_Data_RAW[['Y']])[indiceslocations][:,0]
    Xnodes = Xnodesfull[Nodes]
    Ynodes = Ynodesfull[Nodes]
    
    XYbaseDFfull = pd.DataFrame(np.array(XY_Data_RAW[['X', 'Y']])[indiceslocations,:].transpose(), index = ['X','Y'],columns = all_nodes_titles)
    XYbaseDF = XYbaseDFfull[Nodes]
    
    XYbaserange = [np.min(Xnodesfull),np.max(Xnodesfull),np.min(Ynodesfull),np.max(Ynodesfull)]
    
    NormFac = (XYbaserange[3]-XYbaserange[2]) # WRT sample Size
    NormFacLoading = 180 # WRT point of load application
    
    Xnormnodes = (Xnodes-XYbaserange[0])/NormFac #Normalise by Y (load application)
    Ynormnodes = (Ynodes-XYbaserange[2])/NormFac #Normalise by Y (load application)
    
    U1_Data_Norm = U1_Data_array/NormFac
    U2_Data_Norm = U2_Data_array/NormFac
    
    SQindexX = np.floor(Xnormnodes*Regimes)
    SQindexY = np.floor(Ynormnodes*Regimes)
    
    eps_loading_vec = d_loading_vec/NormFacLoading
        
    
    
    #%% #########################################
    #Find nodes around corners of elipses
    
    distcorner = 0.01
    
    NodesCornerDist = (Regimes+1)*[None]
    for ind0 in range(Regimes+1):
        NodesCornerDist[ind0] = (Regimes+1)*[None]    
    NodesCornerDistIndex = copy(NodesCornerDist)
    
    Xcorneravg=-np.ones((Regimes+1,Regimes+1))
    Ycorneravg=-np.ones((Regimes+1,Regimes+1))
    
    
    
    for ind0 in range(Regimes+1):
        for ind1 in range(Regimes+1):
            indices = np.where(      (np.abs(Xnormnodes-s_Regime*ind0)<distcorner)       &       (np.abs(Ynormnodes-s_Regime*ind1)<distcorner)       )[0]
            NodesCornerDistIndex[ind0][ind1] = indices
    
    
    #%% Calculate ellipse properties across time
    Xcorneravg, Ycorneravg =  Corneravg(NodesCornerDistIndex,Xnormnodes,Ynormnodes,U1_Data_Norm,U2_Data_Norm)
    
    
    #%% Compute Polarisation of all holes
    
    holes1inds = [1,3,5,7,9,11]
    holes2inds = [2,4,6,8,10]
    
    [a_holes1, b_holes1, angle_holes1, Omega_holes1,      a_holes2, b_holes2, angle_holes2, Omega_holes2, xcent1,ycent1,xcent2,ycent2  ] = HoleSpecAnalysis(holes1inds,holes2inds,Xcorneravg,Ycorneravg,U1_Data_Norm,U2_Data_Norm)
    
    
    #%% Compute key Polarisation stats
    Omega_avg,Omega_min,Omega_max = OmegaDist(Omega_holes1,Omega_holes2)
    
    #%% Poisson ratio
    
    Widthsholes1 = xcent1[-1,:,:]-xcent1[0,:,:]  
    Widthsholes2 = xcent2[-1,:,:]-xcent2[0,:,:]  
    Widthsholes1con = np.repeat([Widthsholes1[:,0]],TimeSteps,axis=0).T
    Widthsholes2con = np.repeat([Widthsholes2[:,0]],TimeSteps,axis=0).T
    eps_x1 =  (Widthsholes1-Widthsholes1con)/Widthsholes1con
    eps_x2 =  (Widthsholes2-Widthsholes2con)/Widthsholes2con

    eps_x = np.zeros((len(holes1inds)+len(holes2inds),TimeSteps))
    eps_x[np.array(holes1inds)-1,:]=1.*eps_x1
    eps_x[np.array(holes2inds)-1,:]=1.*eps_x2
        
    eps_x1_avg = np.average(eps_x1,axis=0)
    eps_x2_avg = np.average(eps_x2,axis=0)
    eps_x_avg = np.average(eps_x,axis=0)
    eps_x_middle = eps_x[int((np.shape(eps_x)[0]-1)/2),:]
    
    Heightsholes1 = ycent1[:,-1,:]-ycent1[:,0,:]  
    Heightsholes2 = ycent2[:,-1,:]-ycent2[:,0,:]  
    Heightsholes1con = np.repeat([Heightsholes1[:,0]],TimeSteps,axis=0).T
    Heightsholes2con = np.repeat([Heightsholes2[:,0]],TimeSteps,axis=0).T
    eps_y1 =  (Heightsholes1-Heightsholes1con)/Heightsholes1con
    eps_y2 =  (Heightsholes2-Heightsholes2con)/Heightsholes2con
    
    eps_y = np.zeros((len(holes1inds)+len(holes2inds),TimeSteps))
    eps_y[np.array(holes1inds)-1,:]=1.*eps_y1
    eps_y[np.array(holes2inds)-1,:]=1.*eps_y2
    
    eps_y1_avg = np.average(eps_y1,axis=0)
    eps_y2_avg = np.average(eps_y2,axis=0)
    eps_y_avg = np.average(eps_y,axis=0)
    eps_y_middle = eps_y[int((np.shape(eps_y)[0]-1)/2),:]


    # based on applied strain
    eps_loading_conMAT = np.repeat([eps_loading_vec],len(holes1inds)+len(holes2inds),axis=0)
    deps_x_comp = eps_x[:,IndsCompression[1:]]-eps_x[:,IndsCompression[:-1]]
    deps_loading_conMAT_comp = eps_loading_conMAT[:,IndsCompression[1:]]-eps_loading_conMAT[:,IndsCompression[:-1]]
    Poisson_mat_comp = deps_x_comp/deps_loading_conMAT_comp
    Poisson_avg = np.average(deps_x_comp,axis=0)/np.average(deps_loading_conMAT_comp,axis=0)
       
    
    # based on average measure y-strin of the holes 
    eps_loading_conMAT2 = np.repeat([-eps_y_avg],len(holes1inds)+len(holes2inds),axis=0)
    deps_loading_conMAT2_comp = eps_loading_conMAT2[:,IndsCompression[1:]]-eps_loading_conMAT2[:,IndsCompression[:-1]]
    Poisson_mat_comp2 = deps_x_comp/deps_loading_conMAT2_comp
    Poisson_avg2 = np.average(deps_x_comp,axis=0)/np.average(deps_loading_conMAT2_comp,axis=0)
    
    
    
    #np.average(Poisson_mat_comp,axis=0)
    
    plt.close('all')
    plt.figure()
    plt.plot(Time,eps_loading_vec)
    plt.plot(Time,-eps_y_avg)
    plt.plot(Time,-eps_y_middle)
    #plt.plot(Time,-eps_y1_avg)
    #plt.plot(Time,-eps_y2_avg)
    plt.xlim([0,t_loading])
    
    
    plt.figure()
    plt.plot(Time,-eps_x_avg)
    plt.plot(Time,-eps_x_middle)
    #plt.plot(Time,-eps_x1_avg)
    #plt.plot(Time,-eps_x2_avg)
    plt.xlim([0,t_loading])
    
    
    plt.figure()
    plt.plot(eps_loading_vec[IndsCompression[1:]],Poisson_avg)
    plt.plot(-eps_y_avg[IndsCompression[1:]],Poisson_avg2)
    #plt.plot(Time,-eps_x1_avg)
    #plt.plot(Time,-eps_x2_avg)
    
    
    
    #%% Location domain wall
    Nholes = np.shape(eps_x)[0]
    
    Omega_avg_y =np.zeros((len(holes1inds)+len(holes2inds),TimeSteps))
    hole_y_locations = np.linspace(-5*s_Regime*NormFac,5*s_Regime*NormFac,11)
    hole_y_locationsNORM = hole_y_locations/NormFacLoading
    
    hole_y_locationsFINE = np.linspace(-5*s_Regime*NormFac,5*s_Regime*NormFac,101)
    hole_y_locationsNORMFINE = hole_y_locationsFINE/NormFacLoading
    
    polyfits_Omegay = TimeSteps*[None];    polyfits_Omegay_coeffs = TimeSteps*[None];    y_domainfit = TimeSteps*[None]
    
    domain_y_loc = np.zeros(TimeSteps)
    for ind2 in range(TimeSteps):
        Omega_avg_y[np.array(holes1inds)-1,ind2] = np.average(Omega_holes1[:,:,ind2],axis=0)
        Omega_avg_y[np.array(holes2inds)-1,ind2] = np.average(Omega_holes2[:,:,ind2],axis=0)
        
        polyfits_Omegay_coeffs[ind2] = np.polyfit(hole_y_locationsNORM, Omega_avg_y[:,ind2], 3)
        polyfits_Omegay[ind2] = np.poly1d(polyfits_Omegay_coeffs[ind2])
        y_domainfit[ind2] = np.roots(polyfits_Omegay[ind2])
        
        
        
        domain_y_loc[ind2] = 1.*y_domainfit[ind2][2]
    
        
        if domain_y_loc[ind2]> 1:
            domain_y_loc[ind2]=1
        elif domain_y_loc[ind2]< -1:
            domain_y_loc[ind2]=-1
            
        
    #%%
    """
    
    plt.figure()
    
    for ind2 in [50,100,200,250,298]:
        plt.plot(hole_y_locationsNORM,Omega_avg_y[:,ind2],'b')
        plt.plot(hole_y_locationsNORMFINE,polyfits_Omegay[ind2](hole_y_locationsNORMFINE),'r')
    
    """
    
    
    
    
    
    
    
    
    
    
    
    
    
    #%%###################################################################
    # #Create and save Omega epsilon plots
    ####################################################################
    plt.close('all')
    
    
    
    figsz=(4.,4.)
    axshape = [0.28,0.28,0.62,0.62]

    
    
    eps_mat, hole_y_locationsNORM_mat = np.meshgrid(eps_loading_vec,hole_y_locationsNORM)
    
    nb=100+ind3
    fig=plt.figure(nb,figsize=figsz )
    title='gamma_squareaxis'
    ax=fig.add_axes(axshape)
    ax = plt.gca()
    ax.tick_params(axis = 'both', which = 'major')  
    ax.set_xlabel(r'$\epsilon$')
    ax.set_ylabel(r'$y_{n}$')
    ax.set_xlim([0,0.1])    
    ax.set_ylim([-0.4,0.4])

    ax.set_xticks([0,0.05,0.1])
    ax.set_yticks([-0.4,-0.2,0,0.2,0.4])

    ax.xaxis.set_label_coords(0.5, -0.2)
    ax.yaxis.set_label_coords(-0.25, 0.5)
    ax.contourf(eps_mat,hole_y_locationsNORM_mat,Omega_avg_y,levels=np.linspace(-1,1, 101),cmap=cmapOmega,extend='both')

# =============================================================================
#     
#     FigName = OutputFolderName+'\\EpsilonOmegaContour_MAT-'+str(ind3)
#     if savefigures:
#         plt.savefig(FigName+'.pdf')
#         plt.savefig(FigName+'.png',dpi=300,Transparent=True)
# =============================================================================



    # Do the same with wide lines    
    figsz=(4.,4.)
    axshape = [0.28,0.28,0.62,0.62]
    
    nb=1100+ind3
    fig=plt.figure(nb,figsize=figsz )
    ax=fig.add_axes(axshape)
    ax = plt.gca()
    ax.tick_params(axis = 'both', which = 'major')  
    ax.set_xlabel('$\epsilon$')
    ax.set_ylabel(r'$y_{n}$')
    ax.set_xlim([0,0.07])    
    ax.set_ylim([-0.4,0.4])

    ax.set_xticks([0,0.035,0.07])
    ax.set_yticks([-0.5,0,0.5])

    ax.xaxis.set_label_coords(0.5, -0.2)
    ax.yaxis.set_label_coords(-0.25, 0.5)
    
    
    for ind in range(Nholes):
        x_mat = np.repeat([eps_mat[ind,:]],2,axis=0)
        y_mat = np.repeat([hole_y_locationsNORM_mat[ind,:]],2,axis=0)
        y_mat[0,:]+=-s_Regime/2*NormFac/NormFacLoading
        y_mat[1,:]+=+s_Regime/2*NormFac/NormFacLoading
        z_mat = np.repeat([Omega_avg_y[ind,:]],2,axis=0)
        ax.contourf(x_mat,y_mat,z_mat,levels=np.linspace(-0.8,0.8, 101),cmap=cmapOmega,extend='both')
  
    FigName = OutputFolderName+'\\EpsilonOmegaContour_MAT-'+str(ind3)
    if savefigures:
        plt.savefig(FigName+'.pdf')
        plt.savefig(FigName+'.png',dpi=300,Transparent=True)



    #%%
    
    
    Time_mat = np.repeat([Time-t_loading],len(holes1inds)+len(holes2inds),axis=0)
    
    nb=150+ind3
    fig=plt.figure(nb,figsize=figsz )
    title='gamma_squareaxis'
    ax=fig.add_axes(axshape)
    ax = plt.gca()
    ax.tick_params(axis = 'both', which = 'major')  
    ax.set_xlabel('$t$ [s]')
    ax.set_ylabel(r'$y_{n}$')
    ax.set_xlim([0,1000])    
    ax.set_ylim([-0.4,0.4])

    ax.set_xticks([0,500,1000])
    ax.set_yticks([-0.4,-0.2,0,0.2,0.4])

    ax.xaxis.set_label_coords(0.5, -0.2)
    ax.yaxis.set_label_coords(-0.25, 0.5)
    ax.contourf(Time_mat,hole_y_locationsNORM_mat,Omega_avg_y,levels=np.linspace(-1,1, 101),cmap=cmapOmega,extend='both')



    
# =============================================================================
# 
#     FigName = OutputFolderName+'\\TimeRelaxOmegaContour_MAT-'+str(ind3)
#     if savefigures:
#         plt.savefig(FigName+'.pdf')
#         plt.savefig(FigName+'.png',dpi=300,Transparent=True)
# =============================================================================


    OmegaDict ={}
    for indy in range(np.shape(eps_mat)[0]):
        OmegaDict['Time-'+str(indy)] = Time_mat[indy,:]
    for indy in range(np.shape(eps_mat)[0]):
        OmegaDict['epsilon-'+str(indy)] = eps_mat[indy,:]
    for indy in range(np.shape(eps_mat)[0]):
        OmegaDict['hole_y_locationsNORM-'+str(indy)] = hole_y_locationsNORM_mat[indy,:]
    for indy in range(np.shape(eps_mat)[0]):
        OmegaDict['Omega_avg-'+str(indy)] = Omega_avg_y[indy,:]

        
    OmegaDF = pd.DataFrame(OmegaDict)
    DFName = OutputFolderName+'\\Omega-avg_Time_Epsilon_MAT-'+str(ind3)+'.csv'
    OmegaDF.to_csv(DFName)
    
    
        
    

    nb=1150+ind3
    fig=plt.figure(nb,figsize=figsz )
    ax=fig.add_axes(axshape)
    ax = plt.gca()
    ax.tick_params(axis = 'both', which = 'major')  
    ax.set_xlabel('$t$ [s]')
    ax.set_ylabel(r'$y_{n}$')
    ax.set_xlim([0,1000])    
    ax.set_ylim([-0.4,0.4])

    ax.set_xticks([0,500,1000])
    ax.set_yticks([-0.5,0,0.5])

    ax.xaxis.set_label_coords(0.5, -0.2)
    ax.yaxis.set_label_coords(-0.25, 0.5)
    
    
    for ind in range(Nholes):
        x_mat = np.repeat([Time_mat[ind,:]],2,axis=0)
        y_mat = np.repeat([hole_y_locationsNORM_mat[ind,:]],2,axis=0)
        y_mat[0,:]+=-s_Regime/2*NormFac/NormFacLoading
        y_mat[1,:]+=+s_Regime/2*NormFac/NormFacLoading
        z_mat = np.repeat([Omega_avg_y[ind,:]],2,axis=0)
        ax.contourf(x_mat,y_mat,z_mat,levels=np.linspace(-1.,1., 101),cmap=cmapOmega,extend='both')
  

    FigName = OutputFolderName+'\\TimeRelaxOmegaContour_MAT-'+str(ind3)
    if savefigures:
        plt.savefig(FigName+'.pdf')
        plt.savefig(FigName+'.png',dpi=300,Transparent=True)    
    


    
    #%% #Create and save Poisson epsilon plots
    
    eps_mat_comp1, hole_y_locationsNORM_mat_comp1 = np.meshgrid(eps_loading_vec[IndsCompression[1:]],hole_y_locationsNORM)
    
    nb=200+ind3
    fig=plt.figure(nb,figsize=figsz )
    ax=fig.add_axes(axshape)
    ax = plt.gca()
    ax.tick_params(axis = 'both', which = 'major')  
    ax.set_xlabel('$\epsilon$')
    ax.set_ylabel(r'$\nu$')
    ax.set_xlim([0,0.07])    
    ax.set_ylim([-0.5,0.5])

    ax.set_xticks([0,0.035,0.07])
    ax.set_yticks([-0.5,0,0.5])

    ax.xaxis.set_label_coords(0.5, -0.2)
    ax.yaxis.set_label_coords(-0.25, 0.5)
    
    
    for ind in range(Nholes):
        x_mat = np.repeat([eps_mat_comp1[ind,:]],2,axis=0)
        y_mat = np.repeat([hole_y_locationsNORM_mat_comp1[ind,:]],2,axis=0)
        y_mat[0,:]+=-s_Regime/2*NormFac/NormFacLoading
        y_mat[1,:]+=+s_Regime/2*NormFac/NormFacLoading
        z_mat = np.repeat([Poisson_mat_comp[ind,:]],2,axis=0)
        ax.contourf(x_mat,y_mat,z_mat,levels=np.linspace(-2,2.,101),cmap=cmapPoisson,extend='both')
        
    #ax.contourf(eps_mat_comp1,hole_y_locationsNORM_mat_comp1,Poisson_mat_comp,levels=np.linspace(-2.,2, 1001),cmap='seismic',extend='both')

    #plt.title('$E_{tri}/E_{sq} = $'+str(E_TriH)[:3]+', '+'$t_{load} = $'+str(t_loading)[:5] )
    #plt.xlim([0,0.08])
    #plt.ylim([-1,1])
    
    FigName = OutputFolderName+'\\EpsilonPoissonContour_MAT-'+str(ind3)
    if savefigures:
        plt.savefig(FigName+'.pdf')
        plt.savefig(FigName+'.png',dpi=300,Transparent=True)  
        
    #%% Poisson's ratio average plots
    
    nb=300+ind3
    fig=plt.figure(nb,figsize=figsz )
    title='gamma_squareaxis'
    ax=fig.add_axes(axshape)
    ax = plt.gca()
    ax.tick_params(axis = 'both', which = 'major')  
    ax.set_xlabel('$\epsilon$')
    ax.set_ylabel(r'$y_{n}$')
    ax.set_xlim([0,0.1])    
    ax.set_ylim([-2,0.5])

    ax.set_xticks([0,0.05,0.1])
    #ax.set_yticks([-0.4,-0.2,0,0.2,0.4])

    ax.xaxis.set_label_coords(0.5, -0.2)
    ax.yaxis.set_label_coords(-0.25, 0.5)
    ax.plot(eps_loading_vec[IndsCompression[1:]][:-30],Poisson_avg[:-30])
    
    FigName = OutputFolderName+'\\EpsilonPoissonAverage_MAT-'+str(ind3)
    if savefigures:
        plt.savefig(FigName+'.pdf')
        plt.savefig(FigName+'.png',dpi=300,Transparent=True)       
    
    
    
    
        
    
    
    #%%
    
    
    plt.figure(figsize=(7,7))
    plt.plot(eps_loading_vec,Omega_avg)
    plt.plot(eps_loading_vec,Omega_min)
    plt.plot(eps_loading_vec,Omega_max)
    plt.xlabel('$\epsilon$')
    plt.ylabel('$\Omega$')
    plt.title('$E_{tri}/E_{sq} = $'+str(E_TriH)[:3]+', '+'$t_{load} = $'+str(t_loading)[:5] )
    plt.xlim([0,0.08])
    plt.ylim([-1,1])
    
    
    plt.figure(figsize=(7,7))
    plt.plot(eps_loading_vec,Omega_avg)
    plt.plot(eps_loading_vec,Omega_min)
    plt.plot(eps_loading_vec,Omega_max)
    plt.xlabel('$\epsilon$')
    plt.ylabel('$\Omega$')
    plt.title('$E_{tri}/E_{sq} = $'+str(E_TriH)[:3]+', '+'$t_{load} = $'+str(t_loading)[:5] )
    plt.xlim([0,0.08])
    plt.ylim([-1,1])
    
    FigName = OutputFolderName+'\\EpsilonOmega_MAT-'+str(ind3)
    if savefigures:
        plt.savefig(FigName+'.pdf')
        plt.savefig(FigName+'.png',dpi=300,Transparent=True)
    
    
    plt.figure(figsize=(7,7))
    plt.plot(eps_loading_vec,domain_y_loc)
    plt.xlabel('$\epsilon$')
    plt.ylabel('Domain wall y')
    plt.xlim([0,0.08])
    plt.ylim([-0.5,0.5])    
    FigName = OutputFolderName+'\\EpsilonDomain_MAT-'+str(ind3)    
    if savefigures:
        plt.savefig(FigName+'.pdf')
        plt.savefig(FigName+'.png',dpi=300,Transparent=True)    

    plt.figure(figsize=(7,7))
    plt.plot(Time,domain_y_loc)
    plt.xlabel('t [s]')
    plt.ylabel('Domain wall y')
    plt.ylim([-0.5,0.5])    
    FigName = OutputFolderName+'\\TimeDomain_MAT-'+str(ind3)
    if savefigures:
        plt.savefig(FigName+'.pdf')
        plt.savefig(FigName+'.png',dpi=300,Transparent=True)    
    
    
    

    OmegaminmaxavgDict ={}
    OmegaminmaxavgDict['Time'] = Time
    OmegaminmaxavgDict['Epsilon'] = eps_loading_vec
    OmegaminmaxavgDict['Omega_min'] = Omega_min
    OmegaminmaxavgDict['Omega_max'] = Omega_max
    OmegaminmaxavgDict['Omega_avg'] = Omega_avg
    OmegaminmaxavgDict['domain_y_loc'] = domain_y_loc
    
    OmegaminmaxavgDict['eps_x_avg'] = eps_x_avg
    OmegaminmaxavgDict['eps_x_middle'] = eps_x_middle
    OmegaminmaxavgDict['eps_y_avg'] = eps_y_avg
    OmegaminmaxavgDict['eps_y_middle'] = eps_y_middle
        
    OmegaminmaxavgDF = pd.DataFrame(OmegaminmaxavgDict)
    DFName = OutputFolderName+'\\Omega-minmaxavg_MAT-'+str(ind3)+'.csv'
    OmegaminmaxavgDF.to_csv(DFName)
    
    
    
    
    PoissonDict = {}
    PoissonDict['Time'] = Time[IndsCompression[1:]]
    PoissonDict['epsilon_loading'] = eps_loading_vec[IndsCompression[1:]]
    PoissonDict['-epsilon_y_holesavg'] = -eps_y_avg[IndsCompression[1:]]
    PoissonDict['epsilon_x_holesavg'] = eps_x_avg[IndsCompression[1:]]
    PoissonDict['epsilon_y_holesmid'] = eps_y_middle[IndsCompression[1:]]
    PoissonDict['epsilon_x_holesmid'] = eps_x_middle[IndsCompression[1:]]    
    PoissonDict['Poisson_avg_yloading'] =Poisson_avg
    PoissonDict['Poisson_avg_yholesavg'] =Poisson_avg2
    PoissonDF = pd.DataFrame(PoissonDict)
    DFName = OutputFolderName+'\\Poisson-epsilon_MAT-'+str(ind3)+'.csv'
    PoissonDF.to_csv(DFName)
    
    
    
    
    
    #%%
    """
    #%%
    
    
    domainXlinesX = [0,1]
    domainXlinesY = (Regimes+1)*[0]#[0,1]
    for ind0 in range(Regimes+1):
        domainXlinesY[ind0] = [1./Regimes*ind0,1./Regimes*ind0]
        
    domainYlinesX = copy(domainXlinesY)
    domainYlinesY = copy(domainXlinesX)
    
    
    plt.close('all')
    
    
    
    plt.figure()
    plt.plot(Xnormnodes,Ynormnodes,'ro',markersize=10)
    for ind0 in range(Regimes+1):
        plt.plot(domainXlinesX,domainXlinesY[ind0],'k')
        plt.plot(domainYlinesX[ind0],domainYlinesY,'k')
    
    for ind0 in range(Regimes):
        for ind1 in range(Regimes):     
            if is_odd(ind1-ind0):
                plt.fill_between(   [(ind0)*s_Regime,(ind0+1)*s_Regime],   [(ind1)*s_Regime,(ind1)*s_Regime]   ,[(ind1+1)*s_Regime,(ind1+1)*s_Regime]    ,color='b')
                
            
    
    plt.figure()
    plt.plot(Xnormnodes,SQindexX,'ro')
    plt.figure()
    plt.plot(Ynormnodes,SQindexY,'ro')
    """
    



