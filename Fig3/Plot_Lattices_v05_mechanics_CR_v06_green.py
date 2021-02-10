# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 12:50:12 2019

@author: David
"""

import numpy as np
from matplotlib import pyplot as plt
from copy import deepcopy as copy
import pandas as pd

plt.close('all')


ori=0
struc=0
nx=0
ny=0
lw=1
lw2=1


def plotlattice(ori,nx,ny,struc=1,lw=1,lw2=1,coltri='k'):
    
    
    #ori = orientation type (0,1,2,3)
    #struc: 0 = alll structure with vertex, 1 = real only, 2 = vertex only
    
    zo=1
    
    colmesh='k'
    colmesh_tri = False #    False if black grid, True if colour of triangles 
    
    coltritrue = True
    #coltri='k'
    
    
    #type 1
    if ori==0:    
        if coltritrue:
            coltri=(185/255,255/255,185/255)#'dodgerblue'

        
        if colmesh_tri:
            colmesh = coltri
        
        
        if struc==0 or struc==1:
            trianglex = [ 0+nx, 0.5+nx, 0+nx, 0+nx ] 
            triangley = [ 0+ny, 0+ny,0.5+ny,0+ny]    
            for i in range(3):
                plt.plot(trianglex, triangley, color=coltri,linewidth=lw,zorder=zo)
            plt.fill(trianglex, triangley,color=coltri)
            
            
            trianglex = [ 0+nx, 0.5+nx, 0+nx, 0+nx ] 
            triangley = [ 1+ny, 1+ny,0.5+ny,1+ny]    
            for i in range(3):
                plt.plot(trianglex, triangley, color=coltri,linewidth=lw,zorder=zo)
            plt.fill(trianglex, triangley,color=coltri)
            
            
            trianglex = [ 1+nx, 0.5+nx, 1+nx, 1+nx ] 
            triangley = [ 0+ny, 0+ny,0.5+ny,0+ny]    
            for i in range(3):
                plt.plot(trianglex, triangley, color=coltri,linewidth=lw,zorder=zo)
            plt.fill(trianglex, triangley,color=coltri)
            
            plt.plot([0.5+nx,1+nx],[1+ny,1+ny], color=coltri,linewidth=lw,zorder=zo)
            plt.plot([1+nx,1+nx],[0.5+ny,1+ny], color=coltri,linewidth=lw,zorder=zo)
        
        if struc==0 or struc==2:
    
            plt.plot([0.5+nx,0.5+nx],[0+ny,1+ny], colmesh,linewidth=lw2)
            plt.plot([0+nx,1+nx],[0.5+ny,0.5+ny], colmesh,linewidth=lw2)
            plt.plot([0.5+nx,1+nx],[0.5+ny,1+ny], colmesh,linewidth=lw2)
            
    #type 2
    if ori==1:    
        if coltritrue:
            coltri=(0,140/255,0)#'gold'
        if colmesh_tri:
            colmesh = coltri
            
        if struc==0 or struc==1:
    
            trianglex = [ 0+nx, 0.5+nx, 0+nx, 0+nx ] 
            triangley = [ 0+ny, 0+ny,0.5+ny,0+ny]    
            for i in range(3):
                plt.plot(trianglex, triangley, color=coltri,linewidth=lw,zorder=zo)
            plt.fill(trianglex, triangley,color=coltri)
            
            
            trianglex = [ 0+nx, 0.5+nx, 0+nx, 0+nx ] 
            triangley = [ 1+ny, 1+ny,0.5+ny,1+ny]    
            for i in range(3):
                plt.plot(trianglex, triangley, color=coltri,linewidth=lw,zorder=zo)
            plt.fill(trianglex, triangley,color=coltri)
            
            
            trianglex = [ 1+nx, 0.5+nx, 1+nx, 1+nx ] 
            triangley = [ 1+ny, 1+ny,0.5+ny,1+ny]    
            for i in range(3):
                plt.plot(trianglex, triangley, color=coltri,linewidth=lw,zorder=zo)
            plt.fill(trianglex, triangley,color=coltri)
            
            plt.plot([0.5+nx,1+nx],[0+ny,0+ny], color=coltri,linewidth=lw,zorder=zo)
            plt.plot([1+nx,1+nx],[0.5+ny,0+ny], color=coltri,linewidth=lw,zorder=zo)
            
        if struc==0 or struc==2:
            plt.plot([0.5+nx,0.5+nx],[0+ny,1+ny], colmesh,linewidth=lw2)
            plt.plot([0+nx,1+nx],[0.5+ny,0.5+ny], colmesh,linewidth=lw2)
            plt.plot([0.5+nx,1+nx],[0.5+ny,0+ny], colmesh,linewidth=lw2)
            
        
    
    
    
    
    
    
    #type 3
    if ori==2:
        if coltritrue:
            coltri=(0,240/255,0)#'darkorange'
        if colmesh_tri:
            colmesh = coltri
            
        if struc==0 or struc==1:
            
            
            trianglex = [ 1+nx, 0.5+nx, 1+nx, 1+nx ] 
            triangley = [ 0+ny, 0+ny,0.5+ny,0+ny]    
            for i in range(3):
                plt.plot(trianglex, triangley, color=coltri,linewidth=lw,zorder=zo)
            plt.fill(trianglex, triangley,color=coltri)
            
            
            trianglex = [ 0+nx, 0.5+nx, 0+nx, 0+nx ] 
            triangley = [ 1+ny, 1+ny,0.5+ny,1+ny]    
            for i in range(3):
                plt.plot(trianglex, triangley, color=coltri,linewidth=lw,zorder=zo)
            plt.fill(trianglex, triangley,color=coltri)
            
            
            trianglex = [ 1+nx, 0.5+nx, 1+nx, 1+nx ] 
            triangley = [ 1+ny, 1+ny,0.5+ny,1+ny]    
            for i in range(3):
                plt.plot(trianglex, triangley, color=coltri,linewidth=lw,zorder=zo)
            plt.fill(trianglex, triangley,color=coltri)
            
            plt.plot([0.5+nx,0+nx],[0+ny,0+ny], color=coltri,linewidth=lw,zorder=zo)
            plt.plot([0+nx,0+nx],[0+ny,0.5+ny], color=coltri,linewidth=lw,zorder=zo)
            
        if struc==0 or struc==2:
            plt.plot([0.5+nx,0.5+nx],[0+ny,1+ny], colmesh,linewidth=lw2)
            plt.plot([0+nx,1+nx],[0.5+ny,0.5+ny], colmesh,linewidth=lw2)
            plt.plot([0.5+nx,0+nx],[0.5+ny,0+ny], colmesh,linewidth=lw2)
        
        
    
    
    
    
    
    #type 4
    if ori==3:
        if coltritrue:
            coltri=(125/255,255/255,125/255)#'turquoise'
        if colmesh_tri:
            colmesh = coltri
            
            
        if struc==0 or struc==1:
            trianglex = [ 1+nx, 0.5+nx, 1+nx, 1+nx ] 
            triangley = [ 0+ny, 0+ny,0.5+ny,0+ny]    
            for i in range(3):
                plt.plot(trianglex, triangley, color=coltri,linewidth=lw,zorder=zo)
            plt.fill(trianglex, triangley,color=coltri)
            
            
            trianglex = [ 0+nx, 0.5+nx, 0+nx, 0+nx ] 
            triangley = [ 0+ny, 0+ny,0.5+ny,0+ny]    
            for i in range(3):
                plt.plot(trianglex, triangley, color=coltri,linewidth=lw,zorder=zo)
            plt.fill(trianglex, triangley,color=coltri)
            
            
            trianglex = [ 1+nx, 0.5+nx, 1+nx, 1+nx ] 
            triangley = [ 1+ny, 1+ny,0.5+ny,1+ny]    
            for i in range(3):
                plt.plot(trianglex, triangley, color=coltri,linewidth=lw,zorder=zo)
            plt.fill(trianglex, triangley,color=coltri)
            
            plt.plot([0.5+nx,0+nx],[1+ny,1+ny], color=coltri,linewidth=lw,zorder=zo)
            plt.plot([0+nx,0+nx],[1+ny,0.5+ny], color=coltri,linewidth=lw,zorder=zo)
            
        if struc==0 or struc==2:
            plt.plot([0.5+nx,0.5+nx],[0+ny,1+ny], colmesh,linewidth=lw2)
            plt.plot([0+nx,1+nx],[0.5+ny,0.5+ny], colmesh,linewidth=lw2)
            plt.plot([0.5+nx,0+nx],[0.5+ny,1+ny], colmesh,linewidth=lw2)
        
def plotarrows(x,y,angle_ar,arrowcount,w_ar,lw=1.,coltri='k',arrowdepth=0.5,relplacement = 0.5,countlist=[]):
    
    
    countlist+=[arrowcount]
    
    #plots multiple arrows at given location. Plots arrows in positive x-direction by default, with centre of all arrows at given location by default.
    if arrowcount != 0:
        if arrowcount<-0.1:
            angle_ar = angle_ar+180.
            arrowcount=-arrowcount
            
        xy_loc = [x,y]
        
        wtnorm = arrowcount*arrowdepth
        
        
        
        xline=[0.,0.]
        yline=[0.,0.]
        
        for ind in range(arrowcount):
            
            xnorm=np.array([0.,-arrowdepth,-arrowdepth])+relplacement*wtnorm-ind*arrowdepth
            ynorm=np.array([0.,0.5,-0.5])
            xynorm = np.array([xnorm,ynorm])
            #w_ar= 0.2
            
            xyref = w_ar*xynorm
            #xref = xyref[0]
            #yref = xyref[1]
            
        
            
            Rmat = np.array([[np.cos(angle_ar*np.pi/180),-np.sin(angle_ar*np.pi/180)],[np.sin(angle_ar*np.pi/180),np.cos(angle_ar*np.pi/180)]])
            
            xyrot = np.dot(Rmat,xyref)
            xrot = xyrot[0]
            yrot = xyrot[1]
            
            xfin = xrot  + xy_loc[0]
            yfin = yrot  + xy_loc[1]
            #xyfin = np.array([xfin,yfin])
            
            if ind == 0:
                xline[0] = 1.*xfin[0] 
                yline[0] = 1.*yfin[0] 
            if ind == range(arrowcount)[-1]:
                xline[1] = (xfin[1]+xfin[2])/2 
                yline[1] = (yfin[1]+yfin[2])/2 
                
                
                
            plt.fill(xfin, yfin,coltri,zorder=2)
        #plt.plot(xline,yline,coltri,linewidth=lw,zorder=2)
    return countlist
        
        
        
    
def plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=0.5,relplacement = 0.5,countlist=[]):
    angle_vec=np.array([90.,0.,-90.,180.])
    dxnormvec = np.array([0.5,1.,0.5,0.])
    dynormvec = np.array([1.,0.5,0.,0.5])
    
    #lw=0.
    for ind in range(len(arrowcountvec)):
        x = a*(nx+dxnormvec[ind])
        y = a*(ny+dynormvec[ind])
        countlist = plotarrows(x,y,angle_vec[ind],arrowcountvec[ind],w_ar,lw,coltri,arrowdepth,relplacement,countlist)
        plt.show()
    return countlist
        
        
                





#%%##########################################
#############################################
#############################################
# Perform all calculations
#############################################
#############################################
#############################################
w_ar=0.25
ardepth=3**0.5/2





plt.close('all')

plotlattices = True
if plotlattices:
    
    #%%###########################################
    # Define lattice 2x2
    ############################################
    size=2
    struc2 = np.array([[3,3],
                      [1,3]])
                      
    startname = '2x2Lattice_CR_struc'
    
    for struc in range(3):    
        plt.figure(figsize=(1.5*size,1.5*size))
        plt.hold(True)
        plt.axis('equal')
        
        for nx in range(size):
            for ny in range(size):
                plotlattice(struc2[ny,nx],nx,ny,struc,lw=lw,lw2=lw2)
        plt.axis('off' )
        plt.tight_layout()
    
        #%%###########################################
        # Add arrows
        ############################################    
    
        
        if  struc ==0 or struc ==2:
            
            if struc ==0:
                plt.savefig(startname+str(struc)+'_noarrows.pdf')

                
            countlist=[]
            
            w_ar=0.25
            
            
            
            ###########################
            # Add arrows to squares
            
            # order of locations on square: north, east, south, west
            
            #Orange squares
            
            
            
            for nx in range(0,size,2):
                for ny  in range(1,size,2):
                    nyinv = size-1-ny
                    arrowcountvec = np.array([1,-1,1,-1])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    
            # Other squares
            for nx in range(1,size,2):
                for ny  in range(0,size,2):
                    nyinv = size-1-ny
                    arrowcountvec = np.array([1,-1,1,-1])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
            
            # Edges
            for nx in [0]:
                ind=0
                for ny  in range(0,size,2):
                    nyinv = size-1-ny
                    arrowcountvec = np.array([0,0,0,1])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1
            for nx in [size-1]:
                ind=0
                for ny  in range(1,size,2):                    
                    nyinv = size-1-ny
                    arrowcountvec = np.array([0,1,0,0])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1
                    
            for ny  in [0]:
                ind=0
                for nx in range(0,size,2):
                    nxinv = size-1-nx
                    arrowcountvec = np.array([0,0,-1,0])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1
                    
            for ny  in [size-1]:
                ind=0
                for nx in range(1,size,2):                    
                    nxinv = size-1-nx
                    arrowcountvec = np.array([-1,0,0,0])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1    

        plt.savefig(startname+str(struc)+'.pdf')

    #%%###########################################
    # Define lattice 4x4
    ############################################
    size=4
    struc4 = np.array([[3,3,3,3],
                      [1,3,1,3],
                      [3,3,3,3],
                      [1,3,1,3]])
                      
    startname = '4x4Lattice_CR_struc'
    
    for struc in range(3):    
        plt.figure(figsize=(1.5*size,1.5*size))
        plt.hold(True)
        plt.axis('equal')
        
        for nx in range(size):
            for ny in range(size):
                plotlattice(struc4[ny,nx],nx,ny,struc,lw=lw,lw2=lw2)
        plt.axis('off' )
        plt.tight_layout()
    
        #%%###########################################
        # Add arrows
        ############################################    
    
        
        if  struc ==0 or struc ==2:
            
            if struc ==0:
                plt.savefig(startname+str(struc)+'_noarrows.pdf')

                
            countlist=[]
        
            w_ar=0.25
            
            
            
            ###########################
            # Add arrows to squares
            
            # order of locations on square: north, east, south, west
            
            #Orange squares
            
            
            
            for nx in range(0,size,2):
                for ny  in range(1,size,2):
                    nyinv = size-1-ny
                    arrowcountvec = np.array([1,-1,1,-1])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    
            # Other squares
            for nx in range(1,size,2):
                for ny  in range(0,size,2):
                    nyinv = size-1-ny
                    arrowcountvec = np.array([1,-1,1,-1])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
            
            # Edges
            for nx in [0]:
                ind=0
                for ny  in range(0,size,2):
                    nyinv = size-1-ny
                    arrowcountvec = np.array([0,0,0,1])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1
            for nx in [size-1]:
                ind=0
                for ny  in range(1,size,2):                    
                    nyinv = size-1-ny
                    arrowcountvec = np.array([0,1,0,0])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1
                    
            for ny  in [0]:
                ind=0
                for nx in range(0,size,2):
                    nxinv = size-1-nx
                    arrowcountvec = np.array([0,0,-1,0])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1
                    
            for ny  in [size-1]:
                ind=0
                for nx in range(1,size,2):                    
                    nxinv = size-1-nx
                    arrowcountvec = np.array([-1,0,0,0])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1    

        plt.savefig(startname+str(struc)+'.pdf')
        
    
        if struc ==0:
                
            #%%###########################################
            # Process energy mechanics
            ############################################
            
            countlist4=np.array(countlist)
            
            Efac4 = 1/2*np.sum(countlist4**2)   # *kappa*phi**2
            
            Efac4norm = Efac4    / size**2 
            # Energy per unit cell width = Efac4norm *size *         1/2*kappa*phi**2
            
            dx_arrow_t_min = 1.#np.min(dx_arrow_t_lin)
            
            
            dx_min_lin4 = dx_arrow_t_min / np.sqrt(2)  #Linear in a and phi a
            
            
            
            k_min4 = 2*Efac4/dx_min_lin4**2
            
            k_min4_norm=k_min4/size**2
            
            
    
    #%%###########################################
    # Define lattice 8x8
    ############################################
    
    size=8
    
    struc8 = np.zeros((size,size),dtype=int)
    
    startname = '8x8Lattice_CR_struc'

    
    for ix in range(2):
        for iy in range(2):
            struc8[0+4*iy:4+4*iy,0+4*ix:4+4*ix] = copy(struc4)
    
    
        
    for struc in range(3):    
        plt.figure(figsize=(1.5*size,1.5*size))
        plt.hold(True)
        plt.axis('equal')
        
        for nx in range(size):
            for ny in range(size):
                plotlattice(struc8[ny,nx],nx,ny,struc,lw=1,lw2=1)
        plt.axis('off' )
        plt.tight_layout()
        
    
        #%%###########################################
        # Add arrows
        ############################################    
        if  struc ==0 or struc ==2:
            
            if struc ==0:
                plt.savefig(startname+str(struc)+'_noarrows.pdf')

            countlist=[]
            
            
            
            
            ###########################
            # Add arrows to squares
            
            # order of locations on square: north, east, south, west
            
            #Orange squares
            
            
            
            for nx in range(0,size,2):
                for ny  in range(1,size,2):
                    nyinv = size-1-ny
                    arrowcountvec = -np.array([1,-1,1,-1])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    
            # Other squares
            for nx in range(1,size,2):
                for ny  in range(0,size,2):
                    nyinv = size-1-ny
                    arrowcountvec = -np.array([1,-1,1,-1])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
            
            # Edges
            for nx in [0]:
                ind=0
                for ny  in range(0,size,2):
                    nyinv = size-1-ny
                    arrowcountvec = -np.array([0,0,0,1])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1
            for nx in [size-1]:
                ind=0
                for ny  in range(1,size,2):                    
                    nyinv = size-1-ny
                    arrowcountvec = -np.array([0,1,0,0])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1
                    
            for ny  in [0]:
                ind=0
                for nx in range(0,size,2):
                    nxinv = size-1-nx
                    arrowcountvec = -np.array([0,0,-1,0])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1
                    
            for ny  in [size-1]:
                ind=0
                for nx in range(1,size,2):                    
                    nxinv = size-1-nx
                    arrowcountvec = -np.array([-1,0,0,0])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1    









        plt.savefig(startname+str(struc)+'.pdf')
    
    
        if struc ==0:

            #%%###########################################
            # Process energy mechanics
            ############################################
            
            countlist8=np.array(countlist)
            
            Efac8 = 1/2*np.sum(countlist8**2)   # *kappa*phi**2
            
            Efac8norm = Efac8    / size**2 
            # Energy per unit cell width = Efac8norm *size *         1/2*kappa*phi**2
            
            
            dx_arrow_t_min = 1.#np.min(dx_arrow_t_lin)
            
            
            dx_min_lin8 = dx_arrow_t_min / np.sqrt(2)  #Linear in a and phi
            
            
            
            k_min8 = 2*Efac8/dx_min_lin8**2
            k_min8_norm=k_min8/size**2
        
        
        
        








#%%
plot16 = True
if plot16:
    
    #%%###########################################
    # Define lattice 16x16
    ############################################
    
    size=16
    
    struc16 = np.zeros((size,size),dtype=int)
    
    startname = '16x16Lattice_CR_struc'

    
    for ix in range(4):
        for iy in range(4):
            struc16[0+4*iy:4+4*iy,0+4*ix:4+4*ix] = copy(struc4)
    
    
    for struc in range(3):    
            
        plt.figure(figsize=(1.5*size,1.5*size))
        plt.hold(True)
        plt.axis('equal')
        
        for nx in range(size):
            for ny in range(size):
                plotlattice(struc16[ny,nx],nx,ny,struc,lw=lw,lw2=lw2)
        plt.axis('off' )
        plt.tight_layout()

        #%%###########################################
        # Add arrows
        ############################################    
        # Add arrows to triangles
        
        if  struc ==0 or struc ==2:
            
            if struc ==0:
                plt.savefig(startname+str(struc)+'_noarrows.pdf')            
            countlist=[]
            
            
            ###########################
            # Add arrows to squares
            
            # order of locations on square: north, east, south, west
            
            #Orange squares
            
            
            
            for nx in range(0,size,2):
                for ny  in range(1,size,2):
                    nyinv = size-1-ny
                    arrowcountvec = np.array([1,-1,1,-1])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    
            # Other squares
            for nx in range(1,size,2):
                for ny  in range(0,size,2):
                    nyinv = size-1-ny
                    arrowcountvec = np.array([1,-1,1,-1])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
            
            # Edges
            for nx in [0]:
                ind=0
                for ny  in range(0,size,2):
                    nyinv = size-1-ny
                    arrowcountvec = np.array([0,0,0,1])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1
            for nx in [size-1]:
                ind=0
                for ny  in range(1,size,2):                    
                    nyinv = size-1-ny
                    arrowcountvec = np.array([0,1,0,0])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1
                    
            for ny  in [0]:
                ind=0
                for nx in range(0,size,2):
                    nxinv = size-1-nx
                    arrowcountvec = np.array([0,0,-1,0])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1
                    
            for ny  in [size-1]:
                ind=0
                for nx in range(1,size,2):                    
                    nxinv = size-1-nx
                    arrowcountvec = np.array([-1,0,0,0])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1    









    


        plt.savefig(startname+str(struc)+'.pdf')

        if struc ==0:

            #%%###########################################
            # Process energy mechanics
            ############################################
            
            countlist16=np.array(countlist)
            
            Efac16 = 1/2*np.sum(countlist16**2)   # *kappa*phi**2
            
            Efac16norm = Efac16    / size**2 
            # Energy per unit cell width = Efac16norm *size *         1/2*kappa*phi**2
            
            
            dx_arrow_t_min = 1.#np.min(dx_arrow_t_lin)
            
            dx_min_lin16 = dx_arrow_t_min / np.sqrt(2)  #Linear in a and phi
                
            
            k_min16 = 2*Efac16/dx_min_lin16**2
            k_min16_norm=k_min16/size**2










    #%% Save k norm
    
    k_dic = {}
    k_dic['k_min4'] = [k_min4] 
    k_dic['k_min8'] = [k_min8] 
    k_dic['k_min16'] = [k_min16] 
    k_dic['k_min4_norm'] = [k_min4_norm] 
    k_dic['k_min8_norm'] = [k_min8_norm ]
    k_dic['k_min16_norm'] = [k_min16_norm] 

    k_DF = pd.DataFrame(k_dic)
    
    k_DF.to_csv('k_norm_CR.csv')







    