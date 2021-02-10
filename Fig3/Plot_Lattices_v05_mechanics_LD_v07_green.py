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
    


plt.close('all')

plotlattices = True
if plotlattices:
    
    #%%###########################################
    # Define lattice 4x4
    ############################################
        
    w_ar=0.2
    ardepth=1/2

    size=4
    struc4 = np.array([[3,3,3,3],
                      [1,3,1,3],
                      [3,3,3,3],
                      [1,3,1,3]])
                      
    startname = '4x4Lattice_LD_struc'

    
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
        # Add arrows to triangles
    
        
        if  struc ==0 or struc ==2:
            
            if struc ==0:
                plt.savefig(startname+str(struc)+'_noarrows.pdf')

            countlist=[]
            trilocs=[]
            angle_ar = -45.
            arrowcount=-4
            xlocs = [1.,3.]
            ylocs = copy(xlocs)
            for x in xlocs:
                for y in ylocs:
                    trilocs+=[x,y]
                    countlist = plotarrows(x,y,angle_ar,arrowcount,w_ar,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
            
                    
                    
            
                
            
            ###########################
            # Add arrows to squares
            
            # order of locations on square: north, east, south, west
            
            #Orange squares
            
            for nx in [0,2]:
                for ny  in [1,3]:
                    nyinv = 3-ny
                    arrowcountvec = np.array([1,-3,-3,1])+int(0.5*(nx-nyinv))*np.array([4, -4, 4, -4])
                    arrowcountvec = -arrowcountvec
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=2.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    
            # Other squares
            for nx in [1,3]:
                for ny  in [0,2]:
                    nyinv = 3-ny
                    arrowcountvec = np.array([3,-1,-1,3])+int(0.5*(nx-nyinv))*np.array([4, -4, 4, -4])
                    arrowcountvec = -arrowcountvec
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=2.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
            
            # Edges
            for nx in [0]:
                ind=0
                for ny  in [2,0]:
                    nyinv = 3-ny
                    arrowcountvec = np.array([0,0,0,-3])+ind*np.array([0,0,0,-4])
                    arrowcountvec = -arrowcountvec
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=2.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1
            for nx in [3]:
                ind=0
                for ny  in [1,3]:                    
                    nyinv = 3-ny
                    arrowcountvec = np.array([0,3,0,0])+ind*np.array([0,4,0,0])
                    arrowcountvec = -arrowcountvec
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=2.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1
                    
            for ny  in [0]:
                ind=0
                for nx in [2,0]:
                    nxinv = 3-nx
                    arrowcountvec = np.array([0,0,3,0])+ind*np.array([0,0,4,0])
                    arrowcountvec = -arrowcountvec
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=2.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1
                    
            for ny  in [3]:
                ind=0
                for nx in [1,3]:                    
                    nxinv = 3-nx
                    arrowcountvec = np.array([-3,0,0,0])+ind*np.array([-4,0,0,0])
                    arrowcountvec = -arrowcountvec
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=2.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
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
            
            
            
            # Length change: taylor series 2*a*(sin(phi/2+pi/4)-sin(pi/4)) = a/sqrt(2)* phi - a/(4*sqrt2) * phi**2 +...

            dx_arrow = size*[None]
            dy_arrow = size*[None]

            arrow_vec = np.array(range(1,100,2))
            indsvec = np.array(range(1,len(arrow_vec),2))
            arrow_vec[indsvec]=-1*arrow_vec[indsvec]
            
            arrow_vec = np.append(arrow_vec[::-1],arrow_vec)
            larr = len(arrow_vec)
            
            

            prefac=1
            ind1=0
            for ind in range(size):
                dx_arrow[ind] = prefac*arrow_vec[int(larr/2-size+2*ind1):int(larr/2+2*ind1)]        
                if ind%2==0:# !=int(size/2-1):
                    ind1=ind1+1
                prefac=prefac*-1
            
            
            
            """
            dx_arrow[0] = [-7,5,-3,1]
            dx_arrow[1] = [3,-1,-1,3]
            dx_arrow[2] = [-3,1,1,-3]
            dx_arrow[3] = [-1,3,-5,7]
            """
            
            dx_arrow_t_lin = np.zeros(size)
            dx_arrow_t_quad = np.zeros(size)
            for ind in range(size):
                dx_arrow_t_lin[ind] = np.sum(dx_arrow[ind])
                dx_arrow_t_quad[ind] = np.sum(np.array(dx_arrow[ind])**2)
            dx_arrow_avg_lin = np.average(dx_arrow_t_lin)
            dx_arrow_avg_quad = np.average(dx_arrow_t_quad)
            
            dx_avg_lin4 = dx_arrow_avg_lin  / np.sqrt(2) #Linear in a and phi
            dx_avg_quad4 =-1/(4* np.sqrt(2))  *dx_arrow_avg_quad   #Linear in a, quadratic in phi
            
            # dxavg = dx_norm_avg_quad4 * a * size * phi**2            
            
            dx_arrow_t_min = np.min(dx_arrow_t_lin)
            
            
            dx_min_lin4 = dx_arrow_t_min / np.sqrt(2)  #Linear in a and phi
            
            
            
            k_min4 = 2*Efac4/dx_min_lin4**2
            
            k_min4_norm=k_min4/size**2
            
            
    #%%###########################################
    # Define lattice 8x8
    ############################################
        
    w_ar=0.15
    ardepth=1/2
    

    size=8
    
    struc8 = np.zeros((size,size),dtype=int)
    
    startname = '8x8Lattice_LD_struc'

    for ix in range(2):
        for iy in range(2):
            struc8[0+4*iy:4+4*iy,0+4*ix:4+4*ix] = copy(struc4)
    
    
        
    for struc in range(3):    
        plt.figure(figsize=(1.5*size,1.5*size))
        plt.hold(True)
        plt.axis('equal')
        
        for nx in range(size):
            for ny in range(size):
                plotlattice(struc8[ny,nx],nx,ny,struc,lw=lw,lw2=lw2)
        plt.axis('off' )
        plt.tight_layout()

    
        #%%###########################################
        # Add arrows
        ############################################    
        # Add arrows to triangles
        
        if  struc ==0 or struc ==2:
            
            if struc ==0:
                plt.savefig(startname+str(struc)+'_noarrows.pdf')            
            
            trilocs=[]
            angle_ar = -45.
            arrowcount=-4
            xlocs = [1.,3.,5.,7.]
            ylocs = copy(xlocs)
            for x in xlocs:
                for y in ylocs:
                    trilocs+=[x,y]
                    countlist = plotarrows(x,y,angle_ar,arrowcount,w_ar,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    plt.show()
            
                    
                    
            
                
            
            ###########################
            # Add arrows to squares
            
            # order of locations on square: north, east, south, west
            
            #Orange squares
            
            for nx in [0,2,4,6]:
                for ny  in [1,3,5,7]:
                    nyinv = 7-ny
                    arrowcountvec = np.array([1,-3,-3,1])+int(0.5*(nx-nyinv))*np.array([4, -4, 4, -4])
                    arrowcountvec = -arrowcountvec
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=2.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    
            # Other squares
            for nx in [1,3,5,7]:
                for ny  in [0,2,4,6]:
                    nyinv = 7-ny
                    arrowcountvec = np.array([3,-1,-1,3])+int(0.5*(nx-nyinv))*np.array([4, -4, 4, -4])
                    arrowcountvec = -arrowcountvec
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=2.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
            
            # Edges
            for nx in [0]:
                ind=0
                for ny  in [6,4,2,0]:
                    
                    nyinv = 7-ny
                    arrowcountvec = np.array([0,0,0,-3])+ind*np.array([0,0,0,-4])
                    arrowcountvec = -arrowcountvec
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=2.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1
            for nx in [7]:
                ind=0
                for ny  in [1,3,5,7]:                    
                    nyinv = 7-ny
                    arrowcountvec = np.array([0,3,0,0])+ind*np.array([0,4,0,0])
                    arrowcountvec = -arrowcountvec
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=2.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1
                    
            for ny  in [0]:
                ind=0
                for nx in [6,4,2,0]:
                    
                    nxinv = 7-nx
                    arrowcountvec = np.array([0,0,3,0])+ind*np.array([0,0,4,0])
                    arrowcountvec = -arrowcountvec
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=2.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1
                    
            for ny  in [7]:
                ind=0
                for nx in [1,3,5,7]:
                    
                    nxinv = 7-nx
                    arrowcountvec = np.array([-3,0,0,0])+ind*np.array([-4,0,0,0])
                    arrowcountvec = -arrowcountvec
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=2.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
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
            
            
            
            # Length change: taylor series 2*a*(sin(phi/2+pi/4)-sin(pi/4)) = a/sqrt(2)* phi - a/(4*sqrt2) * phi**2 +...

            dx_arrow = size*[None]
            dy_arrow = size*[None]

            arrow_vec = np.array(range(1,100,2))
            indsvec = np.array(range(1,len(arrow_vec),2))
            arrow_vec[indsvec]=-1*arrow_vec[indsvec]
            
            arrow_vec = np.append(arrow_vec[::-1],arrow_vec)
            larr = len(arrow_vec)
            
            
            prefac=1
            ind1=0
            for ind in range(size):
                dx_arrow[ind] = prefac*arrow_vec[int(larr/2-size+2*ind1):int(larr/2+2*ind1)]        
                if ind%2==0:# !=int(size/2-1):
                    ind1=ind1+1
                prefac=prefac*-1
            
            
            
            """
            dx_arrow[0] = [-7,5,-3,1]
            dx_arrow[1] = [3,-1,-1,3]
            dx_arrow[2] = [-3,1,1,-3]
            dx_arrow[3] = [-1,3,-5,7]
            """
            
            dx_arrow_t_lin = np.zeros(size)
            dx_arrow_t_quad = np.zeros(size)
            for ind in range(size):
                dx_arrow_t_lin[ind] = np.sum(dx_arrow[ind])
                dx_arrow_t_quad[ind] = np.sum(np.array(dx_arrow[ind])**2)
            dx_arrow_avg_lin = np.average(dx_arrow_t_lin)
            dx_arrow_avg_quad = np.average(dx_arrow_t_quad)
            
            dx_avg_lin8 = dx_arrow_avg_lin  / np.sqrt(2) #Linear in a and phi
            dx_avg_quad8 =-1/(4* np.sqrt(2))  *dx_arrow_avg_quad   #Linear in a, quadratic in phi
            
            # dxavg = dx_norm_avg_quad4 * a * size * phi**2            
            
            dx_arrow_t_min = np.min(dx_arrow_t_lin)
            
            
            dx_min_lin8 = dx_arrow_t_min / np.sqrt(2)  #Linear in a and phi
            
            
            
            k_min8 = 2*Efac8/dx_min_lin8**2
            k_min8_norm=k_min8/size**2


#%% 
plot16 = True
if plot16:

    #%%###########################################
    # Define lattice 16x16
    ############################################
            
    w_ar=0.15
    ardepth=1/2

    size=16
    
    struc16 = np.zeros((size,size),dtype=int)
    
    startname = '16x16Lattice_LD_struc'
    
    
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
            trilocs=[]
            angle_ar = -45.
            arrowcount=4
            xlocs = [1.,3.,5.,7.,9.,11,13.,15]
            ylocs = copy(xlocs)
            for x in xlocs:
                for y in ylocs:
                    trilocs+=[x,y]
                    countlist = plotarrows(x,y,angle_ar,arrowcount,w_ar,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    plt.show()
            
                    
        #%%
            
                
            
            ###########################
            # Add arrows to squares
            
            # order of locations on square: north, east, south, west
            
            #Orange squares
            
            for nx in [0,2,4,6,8,10,12,14]:
                for ny  in [1,3,5,7,9,11,13,15]:
                    nyinv = 15-ny
                    arrowcountvec = np.array([1,-3,-3,1])+int(0.5*(nx-nyinv))*np.array([4, -4, 4, -4])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=lw,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    
            # Other squares
            for nx in [1,3,5,7,9,11,13,15]:
                for ny  in [0,2,4,6,8,10,12,14]:
                    nyinv = 15-ny
                    arrowcountvec = np.array([3,-1,-1,3])+int(0.5*(nx-nyinv))*np.array([4, -4, 4, -4])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=lw,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
            
            # Edges
            for nx in [0]:
                ind=0
                for ny  in [14,12,10,8,   6,4,2,0]:
                    
                    nyinv = 15-ny
                    arrowcountvec = np.array([0,0,0,-3])+ind*np.array([0,0,0,-4])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1
            for nx in [15]:
                ind=0
                for ny  in [1,3,5,7,9,11,13,15]:                    
                    nyinv = 15-ny
                    arrowcountvec = np.array([0,3,0,0])+ind*np.array([0,4,0,0])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1
                    
            for ny  in [0]:
                ind=0
                for nx in [14,12,10,8,6,4,2,0]:                    
                    nxinv = 15-nx
                    arrowcountvec = np.array([0,0,3,0])+ind*np.array([0,0,4,0])
                    countlist = plotarrowsquare(arrowcountvec,nx,ny,w_ar,a=1.,lw=1.,coltri='k',arrowdepth=ardepth,relplacement = 0.5,countlist=countlist)
                    ind+=1
                    
            for ny  in [15]:
                ind=0
                for nx in [1,3,5,7,9,11,13,15]:
                    nxinv = 15-nx
                    arrowcountvec = np.array([-3,0,0,0])+ind*np.array([-4,0,0,0])
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
            
            
            
            # Length change: taylor series 2*a*(sin(phi/2+pi/4)-sin(pi/4)) = a/sqrt(2)* phi - a/(4*sqrt2) * phi**2 +...

            dx_arrow = size*[None]
            dy_arrow = size*[None]

            arrow_vec = np.array(range(1,100,2))
            indsvec = np.array(range(1,len(arrow_vec),2))
            arrow_vec[indsvec]=-1*arrow_vec[indsvec]
            
            arrow_vec = np.append(arrow_vec[::-1],arrow_vec)
            larr = len(arrow_vec)
            
            
            prefac=1
            ind1=0
            for ind in range(size):
                dx_arrow[ind] = prefac*arrow_vec[int(larr/2-size+2*ind1):int(larr/2+2*ind1)]        
                if ind%2==0:# !=int(size/2-1):
                    ind1=ind1+1
                prefac=prefac*-1
            
            
            
            """
            dx_arrow[0] = [-7,5,-3,1]
            dx_arrow[1] = [3,-1,-1,3]
            dx_arrow[2] = [-3,1,1,-3]
            dx_arrow[3] = [-1,3,-5,7]
            """
            
            dx_arrow_t_lin = np.zeros(size)
            dx_arrow_t_quad = np.zeros(size)
            for ind in range(size):
                dx_arrow_t_lin[ind] = np.sum(dx_arrow[ind])
                dx_arrow_t_quad[ind] = np.sum(np.array(dx_arrow[ind])**2)
            dx_arrow_avg_lin = np.average(dx_arrow_t_lin)
            dx_arrow_avg_quad = np.average(dx_arrow_t_quad)
            
            dx_avg_lin16 = dx_arrow_avg_lin  / np.sqrt(2) #Linear in a and phi
            dx_avg_quad16 =-1/(4* np.sqrt(2))  *dx_arrow_avg_quad   #Linear in a, quadratic in phi
            
            # dxavg = dx_norm_avg_quad4 * a * size * phi**2            
            
            dx_arrow_t_min = np.min(dx_arrow_t_lin)
            
            
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
    
    k_DF.to_csv('k_norm_LD.csv')







    